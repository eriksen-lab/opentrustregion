// Copyright (C) 2025- Jonas Greiner
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifdef __cplusplus
extern "C" {
#endif

bool solver(
    const bool (*update_orbs)(const double*, double*, double**, double**, bool (**)(const double*, double**)),
    const bool (*obj_func)(const double*, double*),
    const int n_param,
    const bool (*precond)(const double*, const double, double**),
    const bool (*conv_check)(bool*),
    const bool* stability_ptr,
    const bool* line_search_ptr,
    const bool* davidson_ptr,
    const bool* jacobi_davidson_ptr,
    const bool* prefer_jacobi_davidson,
    const double* conv_tol_ptr,
    const int* n_random_trial_vectors_ptr,
    const double* start_trust_radius_ptr,
    const int* n_macro_ptr,
    const int* n_micro_ptr,
    const bool* global_red_factor_ptr,
    const bool* local_red_factor_ptr,
    const int* seed_ptr,
    const int* verbose_ptr,
    const void (*logger)(const char*)
);

bool stability_check(
    const double* h_diag,
    const bool (*hess_x)(const double*, double**),
    const int n_param,
    bool stable,
    double* kappa,
    const bool (*precond)(const double*, const double, double**),
    const bool* jacobi_davidson_ptr,
    const double* conv_tol_ptr,
    const int* n_random_trial_vectors_ptr,
    const int* n_iter_ptr,
    const int* verbose_ptr,
    const void (*logger)(const char*)
);

#ifdef __cplusplus
}
#endif
