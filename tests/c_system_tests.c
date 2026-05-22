// Copyright (C) 2025- Jonas Greiner
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
//
// Pure-C system test for the public C interface declared in opentrustregion.h.
//
// The Fortran-side c_interface_unit_tests cover the bind(C) wrappers but never compile 
// against the C header itself. This test does, so any drift between
// solver_settings_type_c (Fortran) and solver_settings_type (C) is caught here.
//

#include <math.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "opentrustregion.h"

// ---------------------------------------------------------------------------
// Compile-time layout checks for the C structs.
//
// These only verify that the C header is self-consistent: each field sits where the 
// field order claims it does, with no surprise padding before the pointer block. 
// Cross-language drift (Fortran vs. C) is caught at runtime by the default-value test 
// below.
// ---------------------------------------------------------------------------

_Static_assert(offsetof(solver_settings_type, precond) == 0,
               "solver_settings_type: precond must be the first field");
_Static_assert(offsetof(solver_settings_type, project) == 1 * sizeof(void*),
               "solver_settings_type: project must follow precond");
_Static_assert(offsetof(solver_settings_type, conv_check) == 2 * sizeof(void*),
               "solver_settings_type: conv_check must follow project");
_Static_assert(offsetof(solver_settings_type, logger) == 3 * sizeof(void*),
               "solver_settings_type: logger must follow conv_check");

_Static_assert(offsetof(stability_settings_type, precond) == 0,
               "stability_settings_type: precond must be the first field");
_Static_assert(offsetof(stability_settings_type, project) == 1 * sizeof(void*),
               "stability_settings_type: project must follow precond");
_Static_assert(offsetof(stability_settings_type, logger)  == 2 * sizeof(void*),
               "stability_settings_type: logger must follow project");

// ---------------------------------------------------------------------------
// Hartmann 6D function
// ---------------------------------------------------------------------------

#define N_PARAM 6
#define N_TERM 4

extern const c_real hartmann6d_alpha[N_TERM];
const c_real *alpha = hartmann6d_alpha;
extern const c_real hartmann6d_A[ N_TERM * N_PARAM ];
#define A(i,j) (hartmann6d_A[(i) + (j)*N_TERM])
extern const c_real hartmann6d_P[ N_TERM * N_PARAM ];
#define P(i,j) (hartmann6d_P[(i) + (j)*N_TERM])
extern const c_real hartmann6d_minimum1[N_PARAM];
const c_real *minimum1 = hartmann6d_minimum1;
extern const c_real hartmann6d_minimum2[N_PARAM];
const c_real *minimum2 = hartmann6d_minimum2;
extern const c_real hartmann6d_saddle_point[N_PARAM];
const c_real *saddle_point = hartmann6d_saddle_point;

// Module-level state that the callbacks read/write, mirroring the `curr_vars` / `hess` 
// globals in the Fortran test.
static c_real curr_vars[N_PARAM];
static c_real hess[N_PARAM][N_PARAM];

static void exp_terms(const c_real x[N_PARAM], c_real out[N_TERM])
{
    for (int i = 0; i < N_TERM; i++) {
        c_real s = 0.0;
        for (int j = 0; j < N_PARAM; j++) {
            c_real d = x[j] - P(i,j);
            s += A(i,j) * d * d;
        }
        out[i] = exp(-s);
    }
}

static c_real hartmann_func(const c_real x[N_PARAM])
{
    c_real e[N_TERM];
    exp_terms(x, e);
    c_real f = 0.0;
    for (int i = 0; i < N_TERM; i++) f -= alpha[i] * e[i];
    return f;
}

static void hartmann_grad(const c_real x[N_PARAM], c_real grad[N_PARAM])
{
    c_real e[N_TERM];
    exp_terms(x, e);
    for (int j = 0; j < N_PARAM; j++) {
        c_real g = 0.0;
        for (int i = 0; i < N_TERM; i++)
            g += 2.0 * alpha[i] * A(i,j) * (x[j] - P(i,j)) * e[i];
        grad[j] = g;
    }
}

static void hartmann_hess(const c_real x[N_PARAM])
{
    c_real e[N_TERM];
    exp_terms(x, e);
    for (int i = 0; i < N_PARAM; i++) {
        c_real h_ii = 0.0;
        for (int k = 0; k < N_TERM; k++) {
            c_real d = x[i] - P(k,i);
            h_ii += alpha[k] * A(k,i) * e[k] * (1.0 - 2.0 * A(k,i) * d * d);
        }
        hess[i][i] = 2.0 * h_ii;
        for (int j = 0; j < i; j++) {
            c_real h_ij = 0.0;
            for (int k = 0; k < N_TERM; k++) {
                h_ij += alpha[k] * A(k,i) * A(k,j)
                        * (x[i] - P(k,i)) * (x[j] - P(k,j)) * e[k];
            }
            hess[i][j] = -4.0 * h_ij;
            hess[j][i] = hess[i][j];
        }
    }
}

// ---------------------------------------------------------------------------
// Callbacks exposed to the Fortran solver via the C ABI.
// ---------------------------------------------------------------------------

static c_int hess_x_fun(const c_real* x, c_real* hx)
{
    for (int i = 0; i < N_PARAM; i++) {
        c_real s = 0.0;
        for (int j = 0; j < N_PARAM; j++) s += hess[i][j] * x[j];
        hx[i] = s;
    }
    return 0;
}

static c_int update_orbs(const c_real* delta_vars, c_real* func, c_real* grad, 
                         c_real* h_diag, hess_x_fp* hess_x_ptr)
{
    for (int i = 0; i < N_PARAM; i++) curr_vars[i] += delta_vars[i];
    *func = hartmann_func(curr_vars);
    hartmann_grad(curr_vars, grad);
    hartmann_hess(curr_vars);
    for (int i = 0; i < N_PARAM; i++) h_diag[i] = hess[i][i];
    *hess_x_ptr = hess_x_fun;
    return 0;
}

static c_int obj_func(const c_real* delta_vars, c_real* func)
{
    c_real x[N_PARAM];
    for (int i = 0; i < N_PARAM; i++) x[i] = curr_vars[i] + delta_vars[i];
    *func = hartmann_func(x);
    return 0;
}

// Identity preconditioner — exercises the precond callback slot without risking a zero 
// vector when mu=0 (which would trip the Gram-Schmidt zero-vector guard).
static c_int precond(const c_real* residual, const c_real* mu, c_real* precond_residual)
{
    (void)mu;
    for (int i = 0; i < N_PARAM; i++) precond_residual[i] = residual[i];
    return 0;
}

static int logger_called = 0;

static void logger(const char* message)
{
    (void)message;
    logger_called = 1;
}

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

// note: not named `near` because that is a legacy keyword macro on MSVC
static int vec_close(const c_real* a, const c_real* b, c_real tol)
{
    for (int i = 0; i < N_PARAM; i++) if (fabs(a[i] - b[i]) > tol) return 0;
    return 1;
}

static int vec_close_either(const c_real* x, const c_real* a, const c_real* b,
                            c_real tol)
{
    return vec_close(x, a, tol) || vec_close(x, b, tol);
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

bool test_solver_settings_init(void)
{
    // get defaults
    void get_default_solver_values(solver_settings_type* settings);
    solver_settings_type defaults = {0};
    get_default_solver_values(&defaults);

    // call function
    solver_settings_type s = solver_settings_init();
    
    // compare values
    bool ok = true;
    if (!s.initialized) { 
        fprintf(stderr, "test_solver_settings_init failed: Settings not initialized.\n"); ok = false; 
    }
    if (s.stability != defaults.stability) { 
        fprintf(stderr, "test_solver_settings_init failed: Stability parameter wrong.\n"); ok = false; 
    }
    if (s.line_search != defaults.line_search) { 
        fprintf(stderr, "test_solver_settings_init failed: Line search parameter wrong.\n"); ok = false; 
    }
    if (fabs(s.conv_tol - defaults.conv_tol) > 1e-15) { 
        fprintf(stderr, "test_solver_settings_init failed: Convergence tolerance parameter wrong.\n"); ok = false; 
    }
    if (fabs(s.start_trust_radius - defaults.start_trust_radius) > 1e-15) { 
        fprintf(stderr, "test_solver_settings_init failed: Starting trust radius parameter wrong.\n"); ok = false; 
    }
    if (s.n_macro != defaults.n_macro) { 
        fprintf(stderr, "test_solver_settings_init failed: Number of macro iterations parameter wrong.\n"); ok = false; 
    }
    if (s.n_micro != defaults.n_micro) { 
        fprintf(stderr, "test_solver_settings_init failed: Number of micro iterations parameter wrong.\n"); ok = false; 
    }
    if (s.jacobi_davidson_start != defaults.jacobi_davidson_start) { 
        fprintf(stderr, "test_solver_settings_init failed: Jacobi-Davidson starting parameter wrong.\n"); ok = false; 
    }
    if (s.n_random_trial_vectors != defaults.n_random_trial_vectors) { 
        fprintf(stderr, "test_solver_settings_init failed: Number of random trial vectors parameter wrong.\n"); ok = false; 
    }
    if (s.seed != defaults.seed) { 
        fprintf(stderr, "test_solver_settings_init failed: Seed parameter wrong.\n"); ok = false; 
    }
    if (s.verbose != defaults.verbose) { 
        fprintf(stderr, "test_solver_settings_init failed: Verbosity parameter wrong.\n"); ok = false; 
    }
    if (strcmp(s.subsystem_solver, defaults.subsystem_solver) != 0) { 
        fprintf(stderr, "test_solver_settings_init failed: Subsystem solver parameter wrong.\n"); ok = false; 
    }
    if (s.precond || s.project || s.conv_check || s.logger) { 
        fprintf(stderr, "test_solver_settings_init failed: Callback pointers should be NULL.\n"); ok = false;
    }

    return ok;
}

bool test_stability_settings_init(void)
{
    // get defaults
    void get_default_stability_values(stability_settings_type* settings);
    stability_settings_type defaults = {0};
    get_default_stability_values(&defaults);

    // call function
    stability_settings_type s = stability_settings_init();
    
    // compare values
    bool ok = true;
    if (!s.initialized) { 
        fprintf(stderr, "test_stability_settings_init failed: Settings not initialized.\n"); ok = false; 
    }
    if (fabs(s.conv_tol - defaults.conv_tol) > 1e-20) { 
        fprintf(stderr, "test_stability_settings_init failed: Convergence tolerance parameter wrong.\n"); ok = false; 
    }
    if (s.n_random_trial_vectors != defaults.n_random_trial_vectors) { 
        fprintf(stderr, "test_stability_settings_init failed: Number of random trial vectors parameter wrong.\n"); ok = false; 
    }
    if (s.n_iter != defaults.n_iter) { 
        fprintf(stderr, "test_stability_settings_init failed: Number of iterations parameter wrong.\n"); ok = false; 
    }
    if (s.jacobi_davidson_start != defaults.jacobi_davidson_start) { 
        fprintf(stderr, "test_stability_settings_init failed: Jacobi-Davidson starting parameter wrong.\n"); ok = false; 
    }
    if (s.seed != defaults.seed) { 
        fprintf(stderr, "test_stability_settings_init failed: Seed parameter wrong.\n"); ok = false; 
    }
    if (s.verbose != defaults.verbose) { 
        fprintf(stderr, "test_stability_settings_init failed: Verbosity parameter wrong.\n"); ok = false; 
    }
    if (strcmp(s.diag_solver, defaults.diag_solver) != 0) { 
        fprintf(stderr, "test_stability_settings_init failed: Diagonal solver parameter wrong.\n"); ok = false; 
    }
    if (s.precond || s.project || s.logger) { 
        fprintf(stderr, "test_stability_settings_init failed: Callback pointers should be NULL.\n"); ok = false; 
    }

    return ok;
}

bool test_solver_c(void)
{
    bool ok = true;

    solver_settings_type settings = solver_settings_init();
    settings.precond = precond;
    settings.logger  = logger;
    settings.verbose = 3;  // ensure the logger callback is exercised
    logger_called = 0;

    // Start in the quadratic region near first minimum
    const c_real start_near_min1[N_PARAM] = {0.20, 0.15, 0.48, 0.28, 0.31, 0.66};
    memcpy(curr_vars, start_near_min1, sizeof(curr_vars));
    c_int error = solver(update_orbs, obj_func, N_PARAM, settings);
    if (error != 0) { fprintf(stderr, "test_solver_c failed: Produced error.\n"); ok = false; }
    if (!vec_close(curr_vars, minimum1, 1e-4)) {
        fprintf(stderr, "test_solver_c failed: Solver did not find minimum.\n");
        ok = false;
    }
    if (!logger_called) { fprintf(stderr, "test_solver_c failed: Logger was not called.\n"); ok = false; }

    // start near a saddle so the solver has to switch to a non-Newton step
    const c_real start_near_saddle[N_PARAM] = {0.35, 0.59, 0.48, 0.40, 0.31, 0.32};
    memcpy(curr_vars, start_near_saddle, sizeof(curr_vars));
    error = solver(update_orbs, obj_func, N_PARAM, settings);
    if (error != 0) { fprintf(stderr, "test_solver_c failed: Produced error when starting near saddle.\n"); ok = false; }
    if (!vec_close_either(curr_vars, minimum1, minimum2, 1e-4)) {
        fprintf(stderr, "test_solver_c failed: Solver did not find minimum when starting near saddle.\n");
        ok = false;
    }

    return ok;
}

bool test_stability_check_c(void)
{
    bool ok = true;

    stability_settings_type settings = stability_settings_init();
    settings.logger = logger;
    settings.verbose = 3;  // ensure the logger callback is exercised
    logger_called = 0;

    // at a minimum, expect stable
    memcpy(curr_vars, minimum1, sizeof(curr_vars));
    hartmann_hess(curr_vars);
    c_real h_diag[N_PARAM];
    for (int i = 0; i < N_PARAM; i++) h_diag[i] = hess[i][i];
    c_real direction[N_PARAM] = {0};
    c_bool stable = false;
    c_int error = stability_check(h_diag, hess_x_fun, N_PARAM, &stable, settings, 
                                  direction);
    if (error != 0) { 
        fprintf(stderr, "test_stability_check_c failed: Produced error.\n"); 
        ok = false; 
    }
    if (!stable) { 
        fprintf(stderr, "test_stability_check_c failed: Stability incorrectly classifies stability of minimum.\n"); 
        ok = false; 
    }
    if (!logger_called) { 
        fprintf(stderr, "test_stability_check_c failed: Logger was not called.\n");
        ok = false; 
    }

    // at a saddle, expect unstable
    memcpy(curr_vars, saddle_point, sizeof(curr_vars));
    hartmann_hess(curr_vars);
    for (int i = 0; i < N_PARAM; i++) h_diag[i] = hess[i][i];

    stable = true;
    error = stability_check(h_diag, hess_x_fun, N_PARAM, &stable, settings, direction);
    if (error != 0) { 
        fprintf(stderr, "test_stability_check_c failed: Produced error near saddle.\n"); 
        ok = false; 
    }
    if (stable) { 
        fprintf(stderr, "test_stability_check_c failed: Stability incorrectly classifies stability of saddle point.\n"); 
        ok = false; 
    }

    // the descent direction at the saddle should align with the known 
    // negative-curvature eigenvector
    static const c_real ref_direction[N_PARAM] = {
        -0.173375920238, -0.518489821791, -6.432848975252e-3,
        -0.340127852882,  3.066460316955e-3, 0.765095650196
    };
    c_real dot = 0.0;
    for (int i = 0; i < N_PARAM; i++) dot += direction[i] * ref_direction[i];
    if (fabs(fabs(dot) - 1.0) > 1e-6) {
        fprintf(stderr, "test_stability_check_c failed: Stability check does not return correct direction for saddle point.\n");
        ok = false;
    }

    // also exercise the no-direction path
    stable = true;
    error = stability_check(h_diag, hess_x_fun, N_PARAM, &stable, settings, NULL);
    if (error != 0) { 
        fprintf(stderr, "test_stability_check_c failed: Produced error when not passing direction.\n"); 
        ok = false; 
    }
    if (stable) { 
        fprintf(stderr, "test_stability_check_c failed: Stability incorrectly classifies stability of saddle point when not passing direction.\n"); 
        ok = false; 
    }

    return ok;
}
