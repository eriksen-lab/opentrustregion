// Copyright (C) 2025- Jonas Greiner
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENTRUSTREGION_H
#define OPENTRUSTREGION_H

#include <stdbool.h>
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

/* ------------------------------------------------------------------
 * Type aliases matching Fortran kinds
 * ------------------------------------------------------------------ */
/* defined by the CMake target of a library with 64-bit integers, define it when
 * compiling against such a library without CMake */
#ifdef OTR_ILP64
typedef int64_t c_int; /* corresponds to integer(c_ip) */
#else
typedef int32_t c_int; /* corresponds to integer(c_ip) */
#endif
typedef double c_real; /* corresponds to real(c_rp) */
typedef bool c_bool;   /* corresponds to logical(c_bool) */

/* ------------------------------------------------------------------
 * Maximum keyword length
 * ------------------------------------------------------------------ */
#define OTR_KW_LEN 64

/* ------------------------------------------------------------------
 * Declarations for functions and function pointer types
 * ------------------------------------------------------------------ */

/* Hessian-vector product callback */
typedef c_int hess_x_fn(const c_real *x_c, c_real *hess_x, void *context);
typedef hess_x_fn *hess_x_fp;

/* Orbital update callback */
typedef c_int update_orbs_fn(const c_real *kappa, c_real *func, c_real *grad,
                             c_real *h_diag, hess_x_fp *hess_x_ptr, void *context);
typedef update_orbs_fn *update_orbs_fp;

/* Objective function callback */
typedef c_int obj_func_fn(const c_real *kappa, c_real *func, void *context);
typedef obj_func_fn *obj_func_fp;

/* Preconditioner callback */
typedef c_int precond_fn(const c_real *residual, const c_real *mu,
                         c_real *precond_residual, void *context);
typedef precond_fn *precond_fp;

/* Projection callback */
typedef c_int project_fn(c_real *vector, void *context);
typedef project_fn *project_fp;

/* Convergence check callback */
typedef c_int conv_check_fn(c_bool *converged, void *context);
typedef conv_check_fn *conv_check_fp;

/* Extra trial vector callback */
typedef c_int get_extra_trial_vectors_fn(c_real *trial_vectors,
                                         c_int n_extra_trial_vectors, void *context);
typedef get_extra_trial_vectors_fn *get_extra_trial_vectors_fp;

/* Logger callback */
typedef void logger_fn(const char *message, void *context);
typedef logger_fn *logger_fp;

/* ------------------------------------------------------------------
 * Structs corresponding to Fortran settings
 * ------------------------------------------------------------------ */

// Struct corresponding to Fortran type(stability_settings_type_c)
typedef struct {
  precond_fp precond;
  project_fp project;
  get_extra_trial_vectors_fp get_extra_trial_vectors;
  conv_check_fp conv_check;
  logger_fp logger;

  c_bool hess_symm;
  c_bool stop_on_instability;
  c_bool initialized;

  c_real conv_tol;
  c_real stability_thresh;

  c_int n_random_trial_vectors;
  c_int n_extra_trial_vectors;
  c_int n_iter;
  c_int jacobi_davidson_start;
  c_int seed;
  c_int verbose;

  char diag_solver[OTR_KW_LEN + 1];

  /* opaque host data handed back to every callback, must outlive the call */
  void *context;
} stability_settings_type;

// Fortran-callable init subroutine for stability check settings
void init_stability_settings(stability_settings_type *settings);

// Struct corresponding to Fortran type(solver_settings_type_c)
typedef struct {
  precond_fp precond;
  project_fp project;
  get_extra_trial_vectors_fp get_extra_trial_vectors;
  conv_check_fp conv_check;
  logger_fp logger;

  c_bool stability;
  c_bool line_search;
  c_bool grad_noise;
  c_bool hess_symm;
  c_bool initialized;

  c_real conv_tol;
  c_real start_trust_radius;
  c_real global_red_factor;
  c_real local_red_factor;

  c_int n_random_trial_vectors;
  c_int n_extra_trial_vectors;
  c_int n_macro;
  c_int n_micro;
  c_int jacobi_davidson_start;
  c_int seed;
  c_int verbose;

  char subsystem_solver[OTR_KW_LEN + 1];

  stability_settings_type stability_settings;

  /* opaque host data handed back to every callback, must outlive the call */
  void *context;
} solver_settings_type;

// Fortran-callable init subroutine for solver settings
void init_solver_settings(solver_settings_type *settings);

/* ------------------------------------------------------------------
 * Structs corresponding to Fortran states
 * ------------------------------------------------------------------ */

// Struct corresponding to Fortran type(stability_state_type_c)
typedef struct {
  c_real min_eigval;
  c_real residual_rms;

  c_int iteration;
  c_int n_hess_x;
} stability_state_type;

// Struct corresponding to Fortran type(solver_state_type_c)
typedef struct {
  c_bool max_precision_reached;

  c_real func;
  c_real grad_rms;
  c_real trust_radius;
  c_real relative_residual;
  c_real level_shift;
  c_real step_size;

  c_int macro_iteration;
  c_int micro_iteration;
  c_int jacobi_davidson_iteration;
  c_int n_update_orbs;
  c_int n_obj_func;
  c_int n_hess_x;

  stability_state_type stability_state;
} solver_state_type;

/* ------------------------------------------------------------------
 * Fortran wrappers
 * ------------------------------------------------------------------ */

/**
 * C interface to the solver
 *
 * @param update_orbs_ptr   Pointer to update_orbs callback
 * @param obj_func_ptr      Pointer to objective function callback
 * @param n_param           Number of parameters
 * @param settings          Pointer to struct of solver settings, may be NULL to use
 *                          the default settings
 * @param state             Pointer to struct that receives the state of the
 *                          optimization, may be NULL to skip returning it
 * @return                  Integer error code from Fortran
 */
c_int solver(update_orbs_fp update_orbs_ptr, obj_func_fp obj_func_ptr, c_int n_param,
             const solver_settings_type *settings, solver_state_type *state);

/**
 * C interface to the stability check
 *
 * @param h_diag_ptr        Pointer to diagonal Hessian array of n_param elements
 * @param hess_x_ptr        Pointer to Hessian–vector product callback
 * @param n_param           Number of parameters
 * @param stable            Pointer to bool that receives stability result
 * @param settings          Pointer to struct of stability solver settings, may be
 *                          NULL to use the default settings
 * @param state             Pointer to struct that receives the state of the
 *                          stability check, may be NULL to skip returning it
 * @param kappa_ptr         Pointer to array of n_param elements that receives the
 *                          eigenvector of the lowest Hessian eigenvalue, the direction
 *                          of instability at an unstable point, may be NULL to skip
 *                          returning it
 * @return                  Integer error code from Fortran
 */
c_int stability_check(const c_real *h_diag_ptr, hess_x_fp hess_x_ptr, c_int n_param,
                      c_bool *stable, const stability_settings_type *settings,
                      stability_state_type *state, c_real *kappa_ptr);

#ifdef __cplusplus
}
#endif

/* ------------------------------------------------------------------
 * Small C helper functions returning the default settings, which a declared Fortran
 * settings object holds
 * ------------------------------------------------------------------ */
static inline solver_settings_type solver_settings_init(void) {
  solver_settings_type s;
  init_solver_settings(&s);
  return s;
}

static inline stability_settings_type stability_settings_init(void) {
  stability_settings_type s;
  init_stability_settings(&s);
  return s;
}

#endif
