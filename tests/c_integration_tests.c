// Copyright (C) 2025- Jonas Greiner
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

/* Pure-C integration tests for the public C interface declared in opentrustregion.h
 *
 * The Fortran-side c_interface_unit_tests cover the bind(C) wrappers but never compile
 * against the C header itself. These tests do, and test_settings_layout and
 * test_state_layout check that the settings and state structs have the size of the
 * bind(C) types and that every field is read back under its own name, so any drift
 * between the bind(C) types (Fortran) and the structs (C) is caught here, except a
 * field left in a struct that fits into the padding of the bind(C) type, which changes
 * neither the size nor the offset of a named field. */

#include <math.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>

#include "opentrustregion.h"

#include "c_integration_tests.h"

/* ------------------------------------------------------------------
 * Hartmann 6D function
 * ------------------------------------------------------------------ */

#define N_PARAM 6
#define N_TERM 4

extern const c_int hartmann6d_n_param;
extern const c_int hartmann6d_n_terms;
extern const c_real hartmann6d_alpha[N_TERM];
static const c_real *alpha = hartmann6d_alpha;
extern const c_real hartmann6d_A[N_TERM * N_PARAM];
#define A(i, j) (hartmann6d_A[(i) + (j) * N_TERM])
extern const c_real hartmann6d_P[N_TERM * N_PARAM];
#define P(i, j) (hartmann6d_P[(i) + (j) * N_TERM])
extern const c_real hartmann6d_minimum1[N_PARAM];
static const c_real *minimum1 = hartmann6d_minimum1;
extern const c_real hartmann6d_saddle_point[N_PARAM];
static const c_real *saddle_point = hartmann6d_saddle_point;
extern const c_real hartmann6d_near_minimum[N_PARAM];

/* Host data handed to the callbacks through the context of the settings: the current
 * point and its Hessian, flags recording which callbacks were reached, and the state
 * records of the call with what the callbacks read from them. */
typedef struct {
  c_real curr_vars[N_PARAM];
  c_real hess[N_PARAM][N_PARAM];
  const solver_state_type *solver_state;
  const stability_state_type *stability_state;
  int n_update_orbs_calls;
  int state_not_current;
  c_int n_hess_x_seen;
  c_int nested_n_hess_x_seen;
  int precond_called;
  int project_called;
  int logger_called;
  int stability_precond_called;
  int stability_project_called;
  int stability_logger_called;
  int nested_check_ran;
  int nested_check_error;
  int in_nested;
  int outer_used_inner_hess_x;
} hartmann_context;

static void exp_terms(const c_real x[N_PARAM], c_real out[N_TERM]) {
  for (int i = 0; i < N_TERM; i++) {
    c_real s = 0.0;
    for (int j = 0; j < N_PARAM; j++) {
      c_real d = x[j] - P(i, j);
      s += A(i, j) * d * d;
    }
    out[i] = exp(-s);
  }
}

static c_real hartmann_func(const c_real x[N_PARAM]) {
  c_real e[N_TERM];
  exp_terms(x, e);
  c_real f = 0.0;
  for (int i = 0; i < N_TERM; i++)
    f -= alpha[i] * e[i];
  return f;
}

static void hartmann_grad(const c_real x[N_PARAM], c_real grad[N_PARAM]) {
  c_real e[N_TERM];
  exp_terms(x, e);
  for (int j = 0; j < N_PARAM; j++) {
    c_real g = 0.0;
    for (int i = 0; i < N_TERM; i++)
      g += 2.0 * alpha[i] * A(i, j) * (x[j] - P(i, j)) * e[i];
    grad[j] = g;
  }
}

static void hartmann_hess(const c_real x[N_PARAM], c_real hess[N_PARAM][N_PARAM]) {
  c_real e[N_TERM];
  exp_terms(x, e);
  for (int i = 0; i < N_PARAM; i++) {
    c_real h_ii = 0.0;
    for (int k = 0; k < N_TERM; k++) {
      c_real d = x[i] - P(k, i);
      h_ii += alpha[k] * A(k, i) * e[k] * (1.0 - 2.0 * A(k, i) * d * d);
    }
    hess[i][i] = 2.0 * h_ii;
    for (int j = 0; j < i; j++) {
      c_real h_ij = 0.0;
      for (int k = 0; k < N_TERM; k++) {
        h_ij +=
            alpha[k] * A(k, i) * A(k, j) * (x[i] - P(k, i)) * (x[j] - P(k, j)) * e[k];
      }
      hess[i][j] = -4.0 * h_ij;
      hess[j][i] = hess[i][j];
    }
  }
}

/* ------------------------------------------------------------------
 * Callbacks exposed to the Fortran solver via the C ABI
 * ------------------------------------------------------------------ */

static c_int hess_x_fun(const c_real *x, c_real *hx, void *context) {
  hartmann_context *ctx = context;
  if (!ctx)
    return 1;
  if (ctx->stability_state)
    ctx->n_hess_x_seen = ctx->stability_state->n_hess_x;
  /* the internal stability check of a solve calls this function too, the last time
   * after all its other Hessian linear transformations */
  if (ctx->solver_state)
    ctx->nested_n_hess_x_seen = ctx->solver_state->stability_state.n_hess_x;
  for (int i = 0; i < N_PARAM; i++) {
    c_real s = 0.0;
    for (int j = 0; j < N_PARAM; j++)
      s += ctx->hess[i][j] * x[j];
    hx[i] = s;
  }
  return 0;
}

static c_int update_orbs(const c_real *delta_vars, c_real *func, c_real *grad,
                         c_real *h_diag, hess_x_fp *hess_x_ptr, void *context) {
  hartmann_context *ctx = context;
  if (!ctx)
    return 1;
  /* the state is copied before every callback, so it counts the updates before this
   * one */
  if (ctx->solver_state && ctx->solver_state->n_update_orbs != ctx->n_update_orbs_calls)
    ctx->state_not_current = 1;
  ctx->n_update_orbs_calls++;
  for (int i = 0; i < N_PARAM; i++)
    ctx->curr_vars[i] += delta_vars[i];
  *func = hartmann_func(ctx->curr_vars);
  hartmann_grad(ctx->curr_vars, grad);
  hartmann_hess(ctx->curr_vars, ctx->hess);
  for (int i = 0; i < N_PARAM; i++)
    h_diag[i] = ctx->hess[i][i];
  *hess_x_ptr = hess_x_fun;
  return 0;
}

static c_int obj_func(const c_real *delta_vars, c_real *func, void *context) {
  hartmann_context *ctx = context;
  if (!ctx)
    return 1;
  c_real x[N_PARAM];
  for (int i = 0; i < N_PARAM; i++)
    x[i] = ctx->curr_vars[i] + delta_vars[i];
  *func = hartmann_func(x);
  return 0;
}

/* Identity preconditioners exercise the precond callback slots without risking a zero
 * vector when mu=0 (which would trip the Gram-Schmidt zero-vector guard). The
 * stability variant is set only on the nested stability settings, so the test can tell
 * whether the internal stability check reached its own callback slots rather than the
 * solver's. */
static c_int precond(const c_real *residual, const c_real *mu, c_real *precond_residual,
                     void *context) {
  (void)mu;
  hartmann_context *ctx = context;
  if (!ctx)
    return 1;
  ctx->precond_called = 1;
  for (int i = 0; i < N_PARAM; i++)
    precond_residual[i] = residual[i];
  return 0;
}

static c_int stability_precond(const c_real *residual, const c_real *mu,
                               c_real *precond_residual, void *context) {
  (void)mu;
  hartmann_context *ctx = context;
  if (!ctx)
    return 1;
  ctx->stability_precond_called = 1;
  for (int i = 0; i < N_PARAM; i++)
    precond_residual[i] = residual[i];
  return 0;
}

/* Identity projections, the stability variant again only for the nested settings. */
static c_int project(c_real *vector, void *context) {
  (void)vector;
  hartmann_context *ctx = context;
  if (!ctx)
    return 1;
  ctx->project_called = 1;
  return 0;
}

static c_int stability_project(c_real *vector, void *context) {
  (void)vector;
  hartmann_context *ctx = context;
  if (!ctx)
    return 1;
  ctx->stability_project_called = 1;
  return 0;
}

/* The nested stability check is given a Hessian-vector product of its own, distinct
 * from the one the outer solve is using. If anything were still shared between the two
 * calls, the outer solve would resume against this one, which it records. */
static c_int inner_hess_x(const c_real *x, c_real *hx, void *context) {
  hartmann_context *ctx = context;
  if (!ctx)
    return 1;
  if (!ctx->in_nested)
    ctx->outer_used_inner_hess_x = 1;
  return hess_x_fun(x, hx, context);
}

/* Convergence check that runs a whole stability check from inside the running solve,
 * with its own settings object; only nest once, and report convergence by the gradient
 * criterion with the default threshold it replaces, so the outer solve is
 * unaffected. */
static c_int conv_check_nested(c_bool *converged, void *context) {
  hartmann_context *ctx = context;
  if (!ctx)
    return 1;
  c_real grad[N_PARAM], norm_squared = 0.0;
  hartmann_grad(ctx->curr_vars, grad);
  for (int i = 0; i < N_PARAM; i++)
    norm_squared += grad[i] * grad[i];
  *converged = sqrt(norm_squared / N_PARAM) < solver_settings_init().conv_tol;
  if (!ctx->nested_check_ran) {
    c_real h_diag[N_PARAM];
    c_bool stable = false;
    stability_settings_type inner = stability_settings_init();
    inner.context = ctx;
    ctx->nested_check_ran = 1;
    for (int i = 0; i < N_PARAM; i++)
      h_diag[i] = ctx->hess[i][i];
    ctx->in_nested = 1;
    ctx->nested_check_error =
        stability_check(h_diag, inner_hess_x, N_PARAM, &stable, &inner, NULL, NULL);
    ctx->in_nested = 0;
  }
  return 0;
}

/* Loggers, the stability variant again only for the nested settings. A logger cannot
 * return an error, so a call without the host context is recorded here instead, since
 * there is no context to record it in. */
static int logger_without_context = 0;

static void logger(const char *message, void *context) {
  (void)message;
  hartmann_context *ctx = context;
  if (!ctx) {
    logger_without_context = 1;
    return;
  }
  ctx->logger_called = 1;
}

static void stability_logger(const char *message, void *context) {
  (void)message;
  hartmann_context *ctx = context;
  if (!ctx) {
    logger_without_context = 1;
    return;
  }
  ctx->stability_logger_called = 1;
}

/* Settings that are not given or were not initialized carry no host context, so the
 * calls without them use callbacks which keep their point in a static context. */
static hartmann_context settingless_ctx;

static c_int settingless_hess_x(const c_real *x, c_real *hx, void *context) {
  (void)context;
  return hess_x_fun(x, hx, &settingless_ctx);
}

static c_int settingless_update_orbs(const c_real *delta_vars, c_real *func,
                                     c_real *grad, c_real *h_diag,
                                     hess_x_fp *hess_x_ptr, void *context) {
  (void)context;
  c_int error =
      update_orbs(delta_vars, func, grad, h_diag, hess_x_ptr, &settingless_ctx);
  *hess_x_ptr = settingless_hess_x;
  return error;
}

static c_int settingless_obj_func(const c_real *delta_vars, c_real *func,
                                  void *context) {
  (void)context;
  return obj_func(delta_vars, func, &settingless_ctx);
}

/* ------------------------------------------------------------------
 * Reference settings provided by test_reference.f90, so that these tests need no
 * values of their own: the reference value of a field by its name, prefixed by
 * "stability_settings." for the nested settings, C settings filled with these values
 * field by field with only the named logical set if one is named instead of NULL, the
 * sizes of the bind(C) settings types, and comparisons with the default settings
 * ------------------------------------------------------------------ */

void reference_field(const char *name, c_real *value, char *keyword);
void get_reference_solver_values(solver_settings_type *settings,
                                 const char *true_logical);
size_t solver_settings_size(void);
size_t stability_settings_size(void);
bool is_default_solver_settings(const solver_settings_type *settings);
bool is_default_stability_settings(const stability_settings_type *settings);

/* ------------------------------------------------------------------
 * Reference states provided by test_reference.f90 in the same way: the reference
 * value of a field by its name, prefixed by "stability_state." for the nested state, a
 * C state filled with these values field by field, and the sizes of the bind(C) state
 * types
 * ------------------------------------------------------------------ */

c_real reference_state_field(const char *name);
void get_reference_solver_state_values(solver_state_type *state);
size_t solver_state_size(void);
size_t stability_state_size(void);

/* the flush of the standard error of the Fortran runtime, which buffers what the
 * library writes there while it is a file, provided by test_reference.f90 */
void flush_fortran_stderr(void);

/* ------------------------------------------------------------------
 * Helpers
 * ------------------------------------------------------------------ */

static bool hartmann_dimensions_match(const char *test_name) {
  if (hartmann6d_n_param == N_PARAM && hartmann6d_n_terms == N_TERM)
    return true;
  fprintf(stderr, "%s failed: Hartmann 6D dimensions differ from the library's.\n",
          test_name);
  return false;
}

/* the capture of the standard error declared in c_integration_tests.h */
static FILE *stderr_capture;
static int saved_stderr;

void start_stderr_capture(void) {
  fflush(stderr);
  flush_fortran_stderr();
  stderr_capture = tmpfile();
  saved_stderr = dup(STDERR_FILENO);
  dup2(fileno(stderr_capture), STDERR_FILENO);
}

bool stop_stderr_capture(const char *warning) {
  fflush(stderr);
  flush_fortran_stderr();
  dup2(saved_stderr, STDERR_FILENO);
  close(saved_stderr);
  char captured[1024] = {0};
  rewind(stderr_capture);
  size_t length = fread(captured, 1, sizeof(captured) - 1, stderr_capture);
  fclose(stderr_capture);
  return warning ? strstr(captured, warning) != NULL : length == 0;
}

/* check that a bool the library wrote holds 0 or 1, the only values of a C bool,
 * through its byte, since reading any other value as a bool is undefined */
static bool is_c_bool(const c_bool *value) {
  unsigned char byte;
  _Static_assert(sizeof(c_bool) == 1, "is_c_bool: c_bool must be a single byte");
  memcpy(&byte, value, sizeof byte);
  return byte <= 1;
}

static bool check_field(c_real value, const char *name) {
  c_real ref_value;
  char ref_keyword[OTR_KW_LEN + 1];
  reference_field(name, &ref_value, ref_keyword);
  if (value != ref_value) {
    fprintf(stderr, "test_settings_layout failed: Field %s misplaced.\n", name);
    return false;
  }
  return true;
}

static bool check_state_field(c_real value, const char *name) {
  if (value != reference_state_field(name)) {
    fprintf(stderr, "test_state_layout failed: Field %s misplaced.\n", name);
    return false;
  }
  return true;
}

static bool check_keyword(const char *keyword, const char *name) {
  c_real ref_value;
  char ref_keyword[OTR_KW_LEN + 1];
  reference_field(name, &ref_value, ref_keyword);
  if (strcmp(keyword, ref_keyword) != 0) {
    fprintf(stderr, "test_settings_layout failed: Field %s misplaced.\n", name);
    return false;
  }
  return true;
}

/* ------------------------------------------------------------------
 * Tests
 * ------------------------------------------------------------------ */

bool test_settings_layout(void) {
  bool ok = true;

  /* check that the settings structs have the size of the bind(C) settings types, so
   * that a field missing from either is detected, which the field checks below would
   * not read, and stop otherwise, since filling smaller settings would write past
   * their end */
  if (sizeof(solver_settings_type) != solver_settings_size()) {
    fprintf(stderr, "test_settings_layout failed: Size of solver settings differs "
                    "from the library's.\n");
    ok = false;
  }
  if (sizeof(stability_settings_type) != stability_settings_size()) {
    fprintf(stderr, "test_settings_layout failed: Size of stability check settings "
                    "differs from the library's.\n");
    ok = false;
  }
  if (!ok)
    return false;

  /* check that every logical is stored as a C bool and read back under its own name,
   * only one is set at a time so that swapped logicals can be told apart */
  static const char *logicals[] = {"stability", "line_search", "initialized",
                                   "stability_settings.initialized"};
  const int n_logicals = sizeof logicals / sizeof *logicals;
  for (int i = 0; i < n_logicals; i++) {
    solver_settings_type s = {0};
    get_reference_solver_values(&s, logicals[i]);
    const c_bool *read[] = {&s.stability, &s.line_search, &s.initialized,
                            &s.stability_settings.initialized};
    _Static_assert(sizeof read / sizeof *read == sizeof logicals / sizeof *logicals,
                   "test_settings_layout: every logical needs a name and a field");
    for (int j = 0; j < n_logicals; j++) {
      if (!is_c_bool(read[j])) {
        fprintf(stderr,
                "test_settings_layout failed: Field %s not stored as a "
                "C bool.\n",
                logicals[j]);
        ok = false;
      } else if (*read[j] != (i == j)) {
        fprintf(stderr, "test_settings_layout failed: Field %s misplaced.\n",
                logicals[j]);
        ok = false;
      }
    }
  }

  /* check that every other field is read back under its own name */
  solver_settings_type s = {0};
  get_reference_solver_values(&s, NULL);
  ok &= check_field((uintptr_t)s.precond, "precond");
  ok &= check_field((uintptr_t)s.project, "project");
  ok &= check_field((uintptr_t)s.conv_check, "conv_check");
  ok &= check_field((uintptr_t)s.logger, "logger");
  ok &= check_field((uintptr_t)s.context, "context");
  ok &= check_field(s.conv_tol, "conv_tol");
  ok &= check_field(s.start_trust_radius, "start_trust_radius");
  ok &= check_field(s.global_red_factor, "global_red_factor");
  ok &= check_field(s.local_red_factor, "local_red_factor");
  ok &= check_field(s.n_random_trial_vectors, "n_random_trial_vectors");
  ok &= check_field(s.n_macro, "n_macro");
  ok &= check_field(s.n_micro, "n_micro");
  ok &= check_field(s.jacobi_davidson_start, "jacobi_davidson_start");
  ok &= check_field(s.seed, "seed");
  ok &= check_field(s.verbose, "verbose");
  ok &= check_keyword(s.subsystem_solver, "subsystem_solver");

  /* check the nested stability settings, which share their type with the settings
   * of a standalone stability check */
  stability_settings_type ss = s.stability_settings;
  ok &= check_field((uintptr_t)ss.precond, "stability_settings.precond");
  ok &= check_field((uintptr_t)ss.project, "stability_settings.project");
  ok &= check_field((uintptr_t)ss.logger, "stability_settings.logger");
  ok &= check_field((uintptr_t)ss.context, "stability_settings.context");
  ok &= check_field(ss.conv_tol, "stability_settings.conv_tol");
  ok &= check_field(ss.n_random_trial_vectors,
                    "stability_settings.n_random_trial_vectors");
  ok &= check_field(ss.n_iter, "stability_settings.n_iter");
  ok &=
      check_field(ss.jacobi_davidson_start, "stability_settings.jacobi_davidson_start");
  ok &= check_field(ss.seed, "stability_settings.seed");
  ok &= check_field(ss.verbose, "stability_settings.verbose");
  ok &= check_keyword(ss.diag_solver, "stability_settings.diag_solver");

  return ok;
}

bool test_state_layout(void) {
  bool ok = true;

  /* check that the state structs have the size of the bind(C) state types, so that a
   * field missing from either is detected, which the field checks below would not
   * read, and stop otherwise, since filling a smaller state would write past its end */
  if (sizeof(solver_state_type) != solver_state_size()) {
    fprintf(stderr, "test_state_layout failed: Size of solver state differs from the "
                    "library's.\n");
    ok = false;
  }
  if (sizeof(stability_state_type) != stability_state_size()) {
    fprintf(stderr, "test_state_layout failed: Size of stability check state differs "
                    "from the library's.\n");
    ok = false;
  }
  if (!ok)
    return false;

  /* check that every field is read back under its own name */
  solver_state_type s = {0};
  get_reference_solver_state_values(&s);
  if (!is_c_bool(&s.max_precision_reached)) {
    fprintf(stderr, "test_state_layout failed: Field max_precision_reached not stored "
                    "as a C bool.\n");
    ok = false;
  } else {
    ok &= check_state_field(s.max_precision_reached, "max_precision_reached");
  }
  ok &= check_state_field(s.func, "func");
  ok &= check_state_field(s.grad_rms, "grad_rms");
  ok &= check_state_field(s.trust_radius, "trust_radius");
  ok &= check_state_field(s.relative_residual, "relative_residual");
  ok &= check_state_field(s.level_shift, "level_shift");
  ok &= check_state_field(s.step_size, "step_size");
  ok &= check_state_field(s.macro_iteration, "macro_iteration");
  ok &= check_state_field(s.micro_iteration, "micro_iteration");
  ok &= check_state_field(s.jacobi_davidson_iteration, "jacobi_davidson_iteration");
  ok &= check_state_field(s.n_update_orbs, "n_update_orbs");
  ok &= check_state_field(s.n_obj_func, "n_obj_func");
  ok &= check_state_field(s.n_hess_x, "n_hess_x");

  /* the nested state, which also describes a standalone stability check state */
  stability_state_type ss = s.stability_state;
  ok &= check_state_field(ss.min_eigval, "stability_state.min_eigval");
  ok &= check_state_field(ss.residual_rms, "stability_state.residual_rms");
  ok &= check_state_field(ss.iteration, "stability_state.iteration");
  ok &= check_state_field(ss.n_hess_x, "stability_state.n_hess_x");

  return ok;
}

bool test_solver_settings_init(void) {
  bool ok = true;

  /* call function and compare with the default settings without callback functions
   * and host contexts */
  solver_settings_type s = solver_settings_init();
  if (!is_default_solver_settings(&s)) {
    fprintf(stderr,
            "test_solver_settings_init failed: Settings not initialized to "
            "the default values without callback functions and host contexts.\n");
    ok = false;
  }

  return ok;
}

bool test_stability_settings_init(void) {
  bool ok = true;

  /* call function and compare with the default settings without callback functions
   * and host contexts */
  stability_settings_type s = stability_settings_init();
  if (!is_default_stability_settings(&s)) {
    fprintf(stderr,
            "test_stability_settings_init failed: Settings not initialized "
            "to the default values without callback functions and host context.\n");
    ok = false;
  }

  return ok;
}

bool test_solver_c(void) {
  bool ok = true;

  if (!hartmann_dimensions_match("test_solver_c"))
    return false;

  /* start in the quadratic region near the first minimum */
  hartmann_context ctx = {0};
  memcpy(ctx.curr_vars, hartmann6d_near_minimum, sizeof(ctx.curr_vars));

  solver_settings_type settings = solver_settings_init();
  settings.context = &ctx;
  settings.precond = precond;
  settings.project = project;
  settings.logger = logger;
  settings.conv_check = conv_check_nested;
  settings.stability = true;
  settings.stability_settings.precond = stability_precond;
  settings.stability_settings.project = stability_project;
  settings.stability_settings.logger = stability_logger;
  settings.verbose = 3; /* ensure the logger callbacks are exercised */
  logger_without_context = 0;

  /* the maximum precision flag starts opposite to the one the converging solve
   * returns, so that its write-back is detected */
  solver_state_type state = {0};
  state.max_precision_reached = true;
  ctx.solver_state = &state;

  c_int error = solver(update_orbs, obj_func, N_PARAM, &settings, &state);
  if (error != 0) {
    fprintf(stderr, "test_solver_c failed: Produced error.\n");
    ok = false;
  }
  if (state.max_precision_reached) {
    fprintf(stderr, "test_solver_c failed: Maximum precision reached flag was not "
                    "returned.\n");
    ok = false;
  }
  if (!ctx.precond_called) {
    fprintf(stderr, "test_solver_c failed: Preconditioner was not called.\n");
    ok = false;
  }
  if (!ctx.project_called) {
    fprintf(stderr, "test_solver_c failed: Projection was not called.\n");
    ok = false;
  }
  if (!ctx.logger_called) {
    fprintf(stderr, "test_solver_c failed: Logger was not called.\n");
    ok = false;
  }
  if (logger_without_context) {
    fprintf(stderr, "test_solver_c failed: Logger was called without the host "
                    "context.\n");
    ok = false;
  }
  if (state.n_update_orbs <= 0) {
    fprintf(stderr, "test_solver_c failed: Orbital update counter was not "
                    "populated.\n");
    ok = false;
  }
  if (state.n_obj_func <= 0) {
    fprintf(stderr, "test_solver_c failed: Objective function counter was not "
                    "populated.\n");
    ok = false;
  }
  if (state.n_hess_x <= 0) {
    fprintf(stderr, "test_solver_c failed: Hessian linear transformation counter was "
                    "not populated.\n");
    ok = false;
  }
  if (state.stability_state.n_hess_x <= 0) {
    fprintf(stderr, "test_solver_c failed: Hessian linear transformation counter of "
                    "the internal stability check was not populated.\n");
    ok = false;
  }
  /* every orbital update read the updates before it, and the last Hessian linear
   * transformation of the internal stability check the ones before it */
  if (ctx.state_not_current || ctx.n_update_orbs_calls != state.n_update_orbs ||
      ctx.nested_n_hess_x_seen != state.stability_state.n_hess_x - 1) {
    fprintf(stderr, "test_solver_c failed: State was not current in the callbacks.\n");
    ok = false;
  }
  if (!ctx.stability_precond_called) {
    fprintf(stderr, "test_solver_c failed: Preconditioner set on the nested stability "
                    "settings was not called by the internal stability check.\n");
    ok = false;
  }
  if (!ctx.stability_project_called) {
    fprintf(stderr, "test_solver_c failed: Projection set on the nested stability "
                    "settings was not called by the internal stability check.\n");
    ok = false;
  }
  if (!ctx.stability_logger_called) {
    fprintf(stderr, "test_solver_c failed: Logger set on the nested stability "
                    "settings was not called by the internal stability check.\n");
    ok = false;
  }
  if (!ctx.nested_check_ran) {
    fprintf(stderr, "test_solver_c failed: A stability check nested inside the "
                    "running solve did not run.\n");
    ok = false;
  }
  if (ctx.nested_check_error != 0) {
    fprintf(stderr, "test_solver_c failed: A stability check nested inside the "
                    "running solve produced an error.\n");
    ok = false;
  }
  if (ctx.outer_used_inner_hess_x) {
    fprintf(stderr, "test_solver_c failed: After the nested stability check the outer "
                    "solve resumed against the nested call's callbacks.\n");
    ok = false;
  }

  /* settings or nested stability check settings that were not initialized are
   * replaced by the default settings with a warning on the standard error, since the
   * former provide no logger and the latter none here, while no settings stand for the
   * default settings without a warning */
  solver_settings_type uninitialized_settings = {0};
  solver_settings_type nested_uninitialized_settings = solver_settings_init();
  nested_uninitialized_settings.stability_settings = (stability_settings_type){0};
  const solver_settings_type *case_settings[3] = {NULL, &uninitialized_settings,
                                                  &nested_uninitialized_settings};
  const char *case_warnings[3] = {NULL, test_settings_uninitialized_warning,
                                  test_stability_settings_uninitialized_warning};
  const char *case_names[3] = {"without settings",
                               "with settings that were not initialized",
                               "with nested settings that were not initialized"};
  for (int i_case = 0; i_case < 3; i_case++) {
    memset(&settingless_ctx, 0, sizeof(settingless_ctx));
    memcpy(settingless_ctx.curr_vars, hartmann6d_near_minimum,
           sizeof(settingless_ctx.curr_vars));
    start_stderr_capture();
    error = solver(settingless_update_orbs, settingless_obj_func, N_PARAM,
                   case_settings[i_case], NULL);
    if (!stop_stderr_capture(case_warnings[i_case])) {
      fprintf(stderr, "test_solver_c failed: Warning on the standard error wrong %s.\n",
              case_names[i_case]);
      ok = false;
    }
    if (error != 0) {
      fprintf(stderr, "test_solver_c failed: Produced error %s.\n", case_names[i_case]);
      ok = false;
    }
  }

  return ok;
}

bool test_stability_check_c(void) {
  bool ok = true;

  if (!hartmann_dimensions_match("test_stability_check_c"))
    return false;

  hartmann_context ctx = {0};

  stability_settings_type settings = stability_settings_init();
  settings.context = &ctx;
  settings.precond = precond;
  settings.project = project;
  settings.logger = logger;
  settings.verbose = 3; /* ensure the logger callback is exercised */
  logger_without_context = 0;

  /* at a minimum, expect stable */
  hartmann_hess(minimum1, ctx.hess);
  c_real h_diag[N_PARAM];
  for (int i = 0; i < N_PARAM; i++)
    h_diag[i] = ctx.hess[i][i];
  /* the direction starts nonzero so that the zeros written at a stable point are
   * detected */
  c_real direction[N_PARAM];
  for (int i = 0; i < N_PARAM; i++)
    direction[i] = 1.0;
  c_bool stable = false;
  stability_state_type state = {0};
  ctx.stability_state = &state;
  c_int error = stability_check(h_diag, hess_x_fun, N_PARAM, &stable, &settings, &state,
                                direction);
  if (error != 0) {
    fprintf(stderr, "test_stability_check_c failed: Produced error.\n");
    ok = false;
  }
  if (!is_c_bool(&stable)) {
    fprintf(stderr, "test_stability_check_c failed: Stability flag not returned as a C "
                    "bool.\n");
    ok = false;
  } else if (!stable) {
    fprintf(stderr, "test_stability_check_c failed: Stability incorrectly classifies "
                    "stability of minimum.\n");
    ok = false;
  }
  if (!ctx.precond_called) {
    fprintf(stderr, "test_stability_check_c failed: Preconditioner was not called.\n");
    ok = false;
  }
  if (!ctx.project_called) {
    fprintf(stderr, "test_stability_check_c failed: Projection was not called.\n");
    ok = false;
  }
  if (!ctx.logger_called) {
    fprintf(stderr, "test_stability_check_c failed: Logger was not called.\n");
    ok = false;
  }
  if (logger_without_context) {
    fprintf(stderr, "test_stability_check_c failed: Logger was called without the "
                    "host context.\n");
    ok = false;
  }
  if (state.n_hess_x <= 0) {
    fprintf(stderr, "test_stability_check_c failed: Hessian linear transformation "
                    "counter was not populated.\n");
    ok = false;
  }
  if (state.iteration <= 0) {
    fprintf(stderr, "test_stability_check_c failed: State was not returned.\n");
    ok = false;
  }
  /* the last Hessian linear transformation read the ones before it */
  if (ctx.n_hess_x_seen != state.n_hess_x - 1) {
    fprintf(stderr, "test_stability_check_c failed: State was not current in the "
                    "callbacks.\n");
    ok = false;
  }
  ctx.stability_state = NULL;
  for (int i = 0; i < N_PARAM; i++) {
    if (direction[i] != 0.0) {
      fprintf(stderr, "test_stability_check_c failed: Stability check does not return "
                      "a zero direction for minimum.\n");
      ok = false;
      break;
    }
  }

  /* at a saddle, expect unstable */
  hartmann_hess(saddle_point, ctx.hess);
  for (int i = 0; i < N_PARAM; i++)
    h_diag[i] = ctx.hess[i][i];

  stable = true;
  error =
      stability_check(h_diag, hess_x_fun, N_PARAM, &stable, &settings, NULL, direction);
  if (error != 0) {
    fprintf(stderr, "test_stability_check_c failed: Produced error near saddle.\n");
    ok = false;
  }
  if (stable) {
    fprintf(stderr, "test_stability_check_c failed: Stability incorrectly classifies "
                    "stability of saddle point.\n");
    ok = false;
  }

  /* the returned direction replaces the zero direction returned at the minimum and has
   * to be a normalized direction of negative curvature */
  c_real norm_squared = 0.0, curvature = 0.0;
  for (int i = 0; i < N_PARAM; i++) {
    norm_squared += direction[i] * direction[i];
    for (int j = 0; j < N_PARAM; j++)
      curvature += direction[i] * ctx.hess[i][j] * direction[j];
  }
  if (fabs(sqrt(norm_squared) - 1.0) > 1e-6) {
    fprintf(stderr, "test_stability_check_c failed: Stability check does not return a "
                    "normalized direction for saddle point.\n");
    ok = false;
  }
  if (curvature >= 0.0) {
    fprintf(stderr, "test_stability_check_c failed: Stability check does not return a "
                    "direction of negative curvature for saddle point.\n");
    ok = false;
  }

  /* also exercise the no-direction path */
  stable = true;
  error = stability_check(h_diag, hess_x_fun, N_PARAM, &stable, &settings, NULL, NULL);
  if (error != 0) {
    fprintf(stderr, "test_stability_check_c failed: Produced error when not passing "
                    "direction.\n");
    ok = false;
  }
  if (stable) {
    fprintf(stderr, "test_stability_check_c failed: Stability incorrectly classifies "
                    "stability of saddle point when not passing direction.\n");
    ok = false;
  }

  /* settings that were not initialized are replaced by the default settings with a
   * warning on the standard error, since they provide no logger, while no settings
   * stand for the default settings without a warning */
  stability_settings_type uninitialized_settings = {0};
  const stability_settings_type *case_settings[2] = {NULL, &uninitialized_settings};
  const char *case_warnings[2] = {NULL, test_settings_uninitialized_warning};
  const char *case_names[2] = {"without settings",
                               "with settings that were not initialized"};
  memcpy(settingless_ctx.hess, ctx.hess, sizeof(settingless_ctx.hess));
  for (int i_case = 0; i_case < 2; i_case++) {
    start_stderr_capture();
    error = stability_check(h_diag, settingless_hess_x, N_PARAM, &stable,
                            case_settings[i_case], NULL, NULL);
    if (!stop_stderr_capture(case_warnings[i_case])) {
      fprintf(
          stderr,
          "test_stability_check_c failed: Warning on the standard error wrong %s.\n",
          case_names[i_case]);
      ok = false;
    }
    if (error != 0) {
      fprintf(stderr, "test_stability_check_c failed: Produced error %s.\n",
              case_names[i_case]);
      ok = false;
    }
  }

  return ok;
}
