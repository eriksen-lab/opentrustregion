// Copyright (C) 2025- Jonas Greiner
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
//
// Pure-C system test for the public C interface declared in opentrustregion.h.
//
// The Fortran-side c_interface_unit_tests cover the bind(C) wrappers but never
// compile against the C header itself. This test does, so any drift between
// solver_settings_type_c (Fortran) and solver_settings_type (C) — like the
// missing `project` field that slipped through #34 — is caught here.
//
// The numerical test problem is the Hartmann 6D function, mirroring
// tests/opentrustregion_unit_tests.f90 so behaviour can be cross-checked.

#include <math.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "opentrustregion.h"

// ---------------------------------------------------------------------------
// Compile-time layout checks for the C structs.
//
// These only verify that the C header is self-consistent: each field sits
// where the field order claims it does, with no surprise padding before the
// pointer block. Cross-language drift (Fortran vs. C) is caught at runtime by
// the default-value test below — `init_solver_settings` is implemented in
// Fortran, so if the layouts disagree we read garbage instead of the known
// defaults.
// ---------------------------------------------------------------------------

_Static_assert(offsetof(solver_settings_type, precond)    == 0,
               "solver_settings_type: precond must be the first field");
_Static_assert(offsetof(solver_settings_type, project)    == 1 * sizeof(void*),
               "solver_settings_type: project must follow precond");
_Static_assert(offsetof(solver_settings_type, conv_check) == 2 * sizeof(void*),
               "solver_settings_type: conv_check must follow project");
_Static_assert(offsetof(solver_settings_type, logger)     == 3 * sizeof(void*),
               "solver_settings_type: logger must follow conv_check");

_Static_assert(offsetof(stability_settings_type, precond) == 0,
               "stability_settings_type: precond must be the first field");
_Static_assert(offsetof(stability_settings_type, project) == 1 * sizeof(void*),
               "stability_settings_type: project must follow precond");
_Static_assert(offsetof(stability_settings_type, logger)  == 2 * sizeof(void*),
               "stability_settings_type: logger must follow project");

// ---------------------------------------------------------------------------
// Hartmann 6D function (constants copied from
// tests/opentrustregion_unit_tests.f90).
// ---------------------------------------------------------------------------

#define N_PARAM 6

static const c_real alpha[4] = {1.0, 1.2, 3.0, 3.2};
static const c_real A[4][N_PARAM] = {
    {10.0,  3.0, 17.0,  3.5,  1.7,  8.0},
    { 0.05, 10.0, 17.0,  0.1,  8.0, 14.0},
    { 3.0,  3.5,  1.7, 10.0, 17.0,  8.0},
    {17.0,  8.0,  0.05, 10.0,  0.1, 14.0}
};
static const c_real P[4][N_PARAM] = {
    {0.1312, 0.1696, 0.5569, 0.0124, 0.8283, 0.5886},
    {0.2329, 0.4135, 0.8307, 0.3736, 0.1004, 0.9991},
    {0.2348, 0.1451, 0.3522, 0.2883, 0.3047, 0.6650},
    {0.4047, 0.8828, 0.8732, 0.5743, 0.1091, 0.0381}
};

static const c_real minimum1[N_PARAM] = {
    0.20168951, 0.15001069, 0.47687398,
    0.27533243, 0.31165162, 0.65730053
};
static const c_real minimum2[N_PARAM] = {
    0.40465313, 0.88244493, 0.84610160,
    0.57398969, 0.13892673, 0.03849589
};
static const c_real saddle_point[N_PARAM] = {
    0.35278250, 0.59374767, 0.47631257,
    0.40058250, 0.31111531, 0.32397158
};

// Module-level state that the callbacks read/write, mirroring the
// `curr_vars` / `hess` globals in the Fortran test.
static c_real curr_vars[N_PARAM];
static c_real hess[N_PARAM][N_PARAM];

static void exp_terms(const c_real x[N_PARAM], c_real out[4])
{
    for (int i = 0; i < 4; i++) {
        c_real s = 0.0;
        for (int j = 0; j < N_PARAM; j++) {
            c_real d = x[j] - P[i][j];
            s += A[i][j] * d * d;
        }
        out[i] = exp(-s);
    }
}

static c_real hartmann_func(const c_real x[N_PARAM])
{
    c_real e[4];
    exp_terms(x, e);
    c_real f = 0.0;
    for (int i = 0; i < 4; i++) f -= alpha[i] * e[i];
    return f;
}

static void hartmann_grad(const c_real x[N_PARAM], c_real grad[N_PARAM])
{
    c_real e[4];
    exp_terms(x, e);
    for (int j = 0; j < N_PARAM; j++) {
        c_real g = 0.0;
        for (int i = 0; i < 4; i++)
            g += 2.0 * alpha[i] * A[i][j] * (x[j] - P[i][j]) * e[i];
        grad[j] = g;
    }
}

static void hartmann_hess(const c_real x[N_PARAM])
{
    c_real e[4];
    exp_terms(x, e);
    for (int i = 0; i < N_PARAM; i++) {
        c_real h_ii = 0.0;
        for (int k = 0; k < 4; k++) {
            c_real d = x[i] - P[k][i];
            h_ii += alpha[k] * A[k][i] * e[k] * (1.0 - 2.0 * A[k][i] * d * d);
        }
        hess[i][i] = 2.0 * h_ii;
        for (int j = 0; j < i; j++) {
            c_real h_ij = 0.0;
            for (int k = 0; k < 4; k++) {
                h_ij += alpha[k] * A[k][i] * A[k][j]
                        * (x[i] - P[k][i]) * (x[j] - P[k][j]) * e[k];
            }
            hess[i][j] = -4.0 * h_ij;
            hess[j][i] = hess[i][j];
        }
    }
}

// ---------------------------------------------------------------------------
// Callbacks exposed to the Fortran solver via the C ABI.
// ---------------------------------------------------------------------------

static c_int hess_x_cb(const c_real* x, c_real* hx)
{
    for (int i = 0; i < N_PARAM; i++) {
        c_real s = 0.0;
        for (int j = 0; j < N_PARAM; j++) s += hess[i][j] * x[j];
        hx[i] = s;
    }
    return 0;
}

static c_int update_orbs_cb(const c_real* delta_vars,
                            c_real* func, c_real* grad, c_real* h_diag,
                            hess_x_fp* hess_x_ptr)
{
    for (int i = 0; i < N_PARAM; i++) curr_vars[i] += delta_vars[i];
    *func = hartmann_func(curr_vars);
    hartmann_grad(curr_vars, grad);
    hartmann_hess(curr_vars);
    for (int i = 0; i < N_PARAM; i++) h_diag[i] = hess[i][i];
    *hess_x_ptr = hess_x_cb;
    return 0;
}

static c_int obj_func_cb(const c_real* delta_vars, c_real* func)
{
    c_real x[N_PARAM];
    for (int i = 0; i < N_PARAM; i++) x[i] = curr_vars[i] + delta_vars[i];
    *func = hartmann_func(x);
    return 0;
}

// Identity preconditioner — exercises the precond callback slot without
// risking a zero vector when mu=0 (which would trip the Gram-Schmidt
// zero-vector guard).
static c_int precond_cb(const c_real* residual, const c_real* mu,
                        c_real* precond_residual)
{
    (void)mu;
    for (int i = 0; i < N_PARAM; i++) precond_residual[i] = residual[i];
    return 0;
}

static int logger_called = 0;

static void logger_cb(const char* message)
{
    (void)message;
    logger_called = 1;
}

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

static int near(const c_real* a, const c_real* b, c_real tol)
{
    for (int i = 0; i < N_PARAM; i++) if (fabs(a[i] - b[i]) > tol) return 0;
    return 1;
}

static int near_either(const c_real* x, const c_real* a, const c_real* b,
                       c_real tol)
{
    return near(x, a, tol) || near(x, b, tol);
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

static int test_default_settings(void)
{
    // init_solver_settings is implemented in Fortran. If C reads the wrong
    // bytes for any field, the defaults won't match.
    solver_settings_type s = solver_settings_init();
    int ok = 1;
    if (!s.initialized)                            { fprintf(stderr, "solver default initialized != true\n"); ok = 0; }
    if (s.stability)                               { fprintf(stderr, "solver default stability != false\n"); ok = 0; }
    if (s.line_search)                             { fprintf(stderr, "solver default line_search != false\n"); ok = 0; }
    if (fabs(s.conv_tol - 1e-5) > 1e-15)           { fprintf(stderr, "solver default conv_tol wrong: %g\n", s.conv_tol); ok = 0; }
    if (fabs(s.start_trust_radius - 0.4) > 1e-15)  { fprintf(stderr, "solver default start_trust_radius wrong\n"); ok = 0; }
    if (s.n_macro != 150)                          { fprintf(stderr, "solver default n_macro wrong\n"); ok = 0; }
    if (s.n_micro != 50)                           { fprintf(stderr, "solver default n_micro wrong\n"); ok = 0; }
    if (s.jacobi_davidson_start != 30)             { fprintf(stderr, "solver default jacobi_davidson_start wrong\n"); ok = 0; }
    if (s.n_random_trial_vectors != 1)             { fprintf(stderr, "solver default n_random_trial_vectors wrong\n"); ok = 0; }
    if (s.seed != 42)                              { fprintf(stderr, "solver default seed wrong\n"); ok = 0; }
    if (s.verbose != 0)                            { fprintf(stderr, "solver default verbose wrong\n"); ok = 0; }
    if (strcmp(s.subsystem_solver, "davidson") != 0) {
        fprintf(stderr, "solver default subsystem_solver wrong: '%s'\n", s.subsystem_solver);
        ok = 0;
    }
    if (s.precond || s.project || s.conv_check || s.logger) {
        fprintf(stderr, "solver default callback pointers should be NULL\n");
        ok = 0;
    }

    stability_settings_type t = stability_settings_init();
    if (!t.initialized)                            { fprintf(stderr, "stability default initialized != true\n"); ok = 0; }
    if (fabs(t.conv_tol - 1e-8) > 1e-20)           { fprintf(stderr, "stability default conv_tol wrong\n"); ok = 0; }
    if (t.n_random_trial_vectors != 20)            { fprintf(stderr, "stability default n_random_trial_vectors wrong\n"); ok = 0; }
    if (t.n_iter != 100)                           { fprintf(stderr, "stability default n_iter wrong\n"); ok = 0; }
    if (t.jacobi_davidson_start != 50)             { fprintf(stderr, "stability default jacobi_davidson_start wrong\n"); ok = 0; }
    if (t.seed != 42)                              { fprintf(stderr, "stability default seed wrong\n"); ok = 0; }
    if (t.verbose != 0)                            { fprintf(stderr, "stability default verbose wrong\n"); ok = 0; }
    if (strcmp(t.diag_solver, "davidson") != 0) {
        fprintf(stderr, "stability default diag_solver wrong: '%s'\n", t.diag_solver);
        ok = 0;
    }
    if (t.precond || t.project || t.logger) {
        fprintf(stderr, "stability default callback pointers should be NULL\n");
        ok = 0;
    }
    return ok;
}

static int test_solver_runs(void)
{
    int ok = 1;

    // Start in the quadratic region near minimum1. With every optional
    // callback wired up, this also exercises the precond/logger function
    // pointer slots.
    const c_real start_near_min1[N_PARAM] = {
        0.20, 0.15, 0.48, 0.28, 0.31, 0.66
    };
    memcpy(curr_vars, start_near_min1, sizeof(curr_vars));

    solver_settings_type settings = solver_settings_init();
    settings.precond = precond_cb;
    settings.logger  = logger_cb;
    settings.verbose = 3;  // ensure the logger callback is exercised
    logger_called = 0;

    c_int error = solver(update_orbs_cb, obj_func_cb, N_PARAM, settings);
    if (error != 0) { fprintf(stderr, "solver returned error %ld\n", (long)error); ok = 0; }
    if (!logger_called) { fprintf(stderr, "logger callback was never invoked\n"); ok = 0; }
    if (!near(curr_vars, minimum1, 1e-4)) {
        fprintf(stderr, "solver did not converge to the known minimum\n");
        ok = 0;
    }

    // Start near a saddle so the solver has to switch to a non-Newton step.
    const c_real start_near_saddle[N_PARAM] = {
        0.35, 0.59, 0.48, 0.40, 0.31, 0.32
    };
    memcpy(curr_vars, start_near_saddle, sizeof(curr_vars));
    error = solver(update_orbs_cb, obj_func_cb, N_PARAM, settings);
    if (error != 0) { fprintf(stderr, "solver (saddle start) returned error %ld\n", (long)error); ok = 0; }
    if (!near_either(curr_vars, minimum1, minimum2, 1e-4)) {
        fprintf(stderr, "solver from saddle did not reach a known minimum\n");
        ok = 0;
    }

    return ok;
}

static int test_stability_check_runs(void)
{
    int ok = 1;

    // At a true minimum, expect stable.
    memcpy(curr_vars, minimum1, sizeof(curr_vars));
    hartmann_hess(curr_vars);
    c_real h_diag[N_PARAM];
    for (int i = 0; i < N_PARAM; i++) h_diag[i] = hess[i][i];

    stability_settings_type settings = stability_settings_init();
    settings.logger = logger_cb;
    settings.verbose = 3;  // ensure the logger callback is exercised
    logger_called = 0;

    c_real direction[N_PARAM] = {0};
    c_bool stable = false;
    c_int error = stability_check(h_diag, hess_x_cb, N_PARAM, &stable,
                                  settings, direction);
    if (error != 0) { fprintf(stderr, "stability_check returned error %ld\n", (long)error); ok = 0; }
    if (!stable)    { fprintf(stderr, "minimum1 reported as unstable\n"); ok = 0; }
    if (!logger_called) { fprintf(stderr, "stability logger callback was never invoked\n"); ok = 0; }

    // At a saddle, expect unstable.
    memcpy(curr_vars, saddle_point, sizeof(curr_vars));
    hartmann_hess(curr_vars);
    for (int i = 0; i < N_PARAM; i++) h_diag[i] = hess[i][i];

    stable = true;
    error = stability_check(h_diag, hess_x_cb, N_PARAM, &stable, settings,
                            direction);
    if (error != 0) { fprintf(stderr, "stability_check (saddle) returned error %ld\n", (long)error); ok = 0; }
    if (stable)     { fprintf(stderr, "saddle_point reported as stable\n"); ok = 0; }

    // The descent direction at the saddle should align with the known
    // negative-curvature eigenvector (sign-free comparison via |dot|).
    static const c_real ref_direction[N_PARAM] = {
        -0.173375920238, -0.518489821791, -6.432848975252e-3,
        -0.340127852882,  3.066460316955e-3, 0.765095650196
    };
    c_real dot = 0.0;
    for (int i = 0; i < N_PARAM; i++) dot += direction[i] * ref_direction[i];
    if (fabs(fabs(dot) - 1.0) > 1e-6) {
        fprintf(stderr, "stability_check returned wrong descent direction (|dot|=%g)\n",
                fabs(dot));
        ok = 0;
    }

    // Also exercise the no-direction path (kappa NULL).
    stable = true;
    error = stability_check(h_diag, hess_x_cb, N_PARAM, &stable, settings, NULL);
    if (error != 0)  { fprintf(stderr, "stability_check (no kappa) returned error %ld\n", (long)error); ok = 0; }
    if (stable)      { fprintf(stderr, "saddle reported stable on no-kappa path\n"); ok = 0; }

    return ok;
}

int main(void)
{
    int ok = 1;
    ok &= test_default_settings();
    ok &= test_solver_runs();
    ok &= test_stability_check_runs();
    if (ok) {
        printf("c_system_test PASSED\n");
        return 0;
    }
    printf("c_system_test FAILED\n");
    return 1;
}
