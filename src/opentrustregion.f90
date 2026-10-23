! Copyright (C) 2025- Jonas Greiner
!
! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.

module opentrustregion

    use, intrinsic :: iso_fortran_env, only: int32, int64, real64, &
                                             stdout => output_unit, stderr => error_unit

    implicit none

    integer, parameter :: rp = real64
#ifdef USE_ILP64
    integer, parameter :: ip = int64  ! 64-bit integers
#else
    integer, parameter :: ip = int32  ! 32-bit integers
#endif
    integer, parameter :: kw_len = 64

    ! mathematical constants
    real(rp), parameter :: pi = 4.0_rp * atan(1.0_rp)

    ! define trust region parameters
    real(rp), parameter :: &
        trust_radius_shrink_ratio = 0.25_rp, trust_radius_expand_ratio = 0.75_rp, &
        trust_radius_shrink_factor = 0.5_rp, trust_radius_expand_factor = 2.0_rp

    ! define error codes
    integer(ip), parameter :: &
        error_solver = 100, error_solver_max_iter = error_solver + 2, &
        error_stability_check = 200, &
        error_stability_check_max_iter = error_stability_check + 2, &
        error_gram_schmidt_lin_dep = 51, error_obj_func = 1100, &
        error_update_orbs = 1200, error_hess_x = 1300, error_precond = 1400, &
        error_conv_check = 1500, error_project = 1600, &
        error_get_extra_trial_vectors = 1700

    ! define useful parameters
    real(rp), parameter :: numerical_zero = 1e-14_rp, precond_floor = 1e-10_rp, &
                           hess_symm_thres = 1e-12_rp, residual_norm_floor = 1e-12_rp, &
                           level_shift_local_thres = 1e-12_rp, &
                           newton_eigval_thresh = -1e-5_rp

    ! define verbosity levels
    integer(ip), parameter :: verbosity_silent = 0, verbosity_error = 1, &
                              verbosity_warning = 2, verbosity_info = 3, &
                              verbosity_debug = 4

    ! define the options of the subsystem solver and of the diagonalization solver
    character(len=*), parameter :: &
        subsystem_solver_options(5) = [character(len=18) :: "davidson_ls", &
                                       "davidson_ah", "jacobi-davidson_ls", &
                                       "jacobi-davidson_ah", "tcg"], &
        diag_solver_options(2) = [character(len=15) :: "davidson", "jacobi-davidson"]

    ! define warning messages
    character(len=*), parameter :: &
        settings_uninitialized_warning_msg = &
            "Settings were not initialized. All settings are set to default values", &
        stability_settings_uninitialized_warning_msg = &
            "Stability check settings were not initialized. All stability check "// &
            "settings are set to default values", &
        started_at_saddle_point_warning_msg = &
            "Started at saddle point. The algorithm will continue by moving along "// &
            "eigenvector direction corresponding to negative eigenvalue.", &
        reached_saddle_point_warning_msg = &
            "Reached saddle point. This is likely due to symmetry and can be "// &
            "avoided by decreasing the convergence criterion. The algorithm will "// &
            "continue by moving along eigenvector direction corresponding to "// &
            "negative eigenvalue.", &
        reached_saddle_point_no_grad_noise_warning_msg = &
            "Reached saddle point. This is likely due to symmetry and can be "// &
            "avoided by enabling the gradient perturbation grad_noise. The "// &
            "algorithm will continue by moving along eigenvector direction "// &
            "corresponding to negative eigenvalue.", &
        unstable_warning_msg = "Solution not stable. Lowest eigenvalue: ", &
        bisection_imaginary_eigval_warning_msg = &
            "Bisection cannot reach trust region due to imaginary part of lowest "// &
            "eigenvalue. This is purely an artifact of using an asymmetric Hessian.", &
        complex_eigval_warning_msg = &
            "Lowest eigenvalues include a complex conjugate pair, for which "// &
            "normalized real vectors of its invariant subspace are returned "// &
            "instead of eigenvectors. This is purely an artifact of using an "// &
            "asymmetric Hessian.", &
        lin_dep_trial_vector_warning_msg = &
            "Extra trial vector is linearly dependent and is being removed.", &
        minres_max_iter_warning_msg = "MINRES: The iteration limit was reached. "// &
                                      "The approximate solution of the last "// &
                                      "iteration is used.", &
        function_unchanged_warning_msg = "No further improvement of cost function "// &
                                         "possible. Convergence criterion is not "// &
                                         "fulfilled but calculation should be "// &
                                         "converged up to floating point precision.", &
        random_trial_vector_warning_msg = "Number of random trial vectors should "// &
                                          "be smaller than half the number of "// &
                                          "parameters.", &
        project_warning_msg = &
            "Custom projection is provided. To optimize performance, OTR assumes "// &
            "that all other provided routines (update_orbs, hess_x, precond) are "// &
            "already projected onto the relevant orbital rotation subspace. If "// &
            "these routines are not self-projecting, redundant rotations may "// &
            "contaminate the trial space and cause convergence issues."

    ! state of the stability check, which the stability check keeps current in the
    ! record it is given
    type :: stability_state_type
        ! stability check iteration
        integer(ip) :: iteration = 0
        ! current estimate of the lowest eigenvalue and the RMS of its residual
        real(rp) :: min_eigval = 0.0_rp, residual_rms = 0.0_rp
        ! number of Hessian linear transformations
        integer(ip) :: n_hess_x = 0
    end type

    ! state of the optimization, which the solver keeps current in the record it is
    ! given during the optimization
    type :: solver_state_type
        ! macro iteration
        integer(ip) :: macro_iteration = 0
        ! objective function value and RMS gradient at the current point
        real(rp) :: func = 0.0_rp, grad_rms = 0.0_rp
        ! current trust radius
        real(rp) :: trust_radius = 0.0_rp
        ! micro iteration of the subproblem being solved and how many of these were
        ! Jacobi-Davidson iterations
        integer(ip) :: micro_iteration = 0, jacobi_davidson_iteration = 0
        ! Euclidean norm of the residual of the subproblem relative to that of the
        ! gradient it is solved for, as in the stopping test of the micro iterations
        real(rp) :: relative_residual = 0.0_rp
        ! level shift and size of the last step in the trust-region metric
        real(rp) :: level_shift = 0.0_rp, step_size = 0.0_rp
        ! number of orbital updates, objective function evaluations and Hessian linear
        ! transformations, the latter including those of the internal stability checks
        integer(ip) :: n_update_orbs = 0, n_obj_func = 0, n_hess_x = 0
        ! whether the predicted function change of the last trust-region step could
        ! not be resolved while its actual change could not be resolved either or the
        ! step was rejected, or the gradient vanishes exactly, unless the convergence
        ! check passes
        logical :: max_precision_reached = .false.
        ! state of the most recent internal stability check
        type(stability_state_type) :: stability_state = stability_state_type()
    end type

    ! interfaces for callback functions
    abstract interface
        subroutine hess_x_type(x, hess_x, error, context)
            import :: rp, ip

            real(rp), intent(in), target :: x(:)
            real(rp), intent(out), target :: hess_x(:)
            integer(ip), intent(out) :: error
            class(*), intent(in), pointer :: context
        end subroutine hess_x_type
    end interface

    abstract interface
        subroutine update_orbs_type(kappa, func, grad, h_diag, hess_x_funptr, error, &
                                    context)
            import :: rp, hess_x_type, ip

            real(rp), intent(in), target :: kappa(:)
            real(rp), intent(out) :: func
            real(rp), intent(out), target :: grad(:), h_diag(:)
            procedure(hess_x_type), intent(inout), pointer :: hess_x_funptr
            integer(ip), intent(out) :: error
            class(*), intent(in), pointer :: context
        end subroutine update_orbs_type
    end interface

    abstract interface
        function obj_func_type(kappa, error, context) result(func)
            import :: rp, ip

            real(rp), intent(in), target :: kappa(:)
            integer(ip), intent(out) :: error
            class(*), intent(in), pointer :: context
            real(rp) :: func
        end function obj_func_type
    end interface

    abstract interface
        subroutine precond_type(residual, mu, precond_residual, error, context)
            import :: rp, ip

            real(rp), intent(in), target :: residual(:)
            real(rp), intent(in) :: mu
            real(rp), intent(out), target :: precond_residual(:)
            integer(ip), intent(out) :: error
            class(*), intent(in), pointer :: context
        end subroutine precond_type
    end interface

    abstract interface
        subroutine project_type(vector, error, context)
            import :: rp, ip

            real(rp), intent(inout), target :: vector(:)
            integer(ip), intent(out) :: error
            class(*), intent(in), pointer :: context
        end subroutine project_type
    end interface

    abstract interface
        function conv_check_type(error, context) result(converged)
            import :: ip

            integer(ip), intent(out) :: error
            class(*), intent(in), pointer :: context
            logical :: converged
        end function conv_check_type
    end interface

    abstract interface
        subroutine get_extra_trial_vectors_type(trial_vectors, error, context)
            import :: rp, ip

            real(rp), intent(out), target :: trial_vectors(:, :)
            integer(ip), intent(out) :: error
            class(*), intent(in), pointer :: context
        end subroutine get_extra_trial_vectors_type
    end interface

    abstract interface
        subroutine logger_type(message, context)
            character(len=*), intent(in) :: message
            class(*), intent(in), pointer :: context
        end subroutine logger_type
    end interface

    ! derived types for the settings
    type, abstract :: settings_type
        logical :: hess_symm = .true.
        integer(ip) :: n_extra_trial_vectors = 1, seed = 42, verbose = 0
        class(*), pointer :: context => null()
        procedure(precond_type), pointer, nopass :: precond => null()
        procedure(project_type), pointer, nopass :: project => null()
        procedure(get_extra_trial_vectors_type), pointer, nopass :: &
            get_extra_trial_vectors => null()
        procedure(logger_type), pointer, nopass :: logger => null()
    contains
        procedure :: log => print_message
    end type

    type, extends(settings_type) :: stability_settings_type
        logical :: stop_on_instability = .false.
        real(rp) :: conv_tol = 1e-9_rp, stability_thresh = -1e-3_rp
        integer(ip) :: n_random_trial_vectors = -1, n_iter = 100, &
                       jacobi_davidson_start = 50
        character(len=kw_len) :: diag_solver = "davidson"
        procedure(conv_check_type), pointer, nopass :: conv_check => null()
    end type

    type, extends(settings_type) :: solver_settings_type
        logical :: stability = .false., line_search = .false., grad_noise = .true.
        real(rp) :: conv_tol = 1e-5_rp, start_trust_radius = 0.4_rp, &
                    global_red_factor = 1e-3_rp, local_red_factor = 1e-4_rp
        integer(ip) :: n_random_trial_vectors = 0, n_macro = 150, n_micro = 50, &
                       jacobi_davidson_start = 30
        character(len=kw_len) :: subsystem_solver = "davidson_ls"
        type(stability_settings_type) :: stability_settings = stability_settings_type()
        procedure(conv_check_type), pointer, nopass :: conv_check => null()
    contains
        procedure :: print_results
    end type

    ! default settings
    type(stability_settings_type), parameter :: default_stability_settings = &
        stability_settings_type()
    type(solver_settings_type), parameter :: default_solver_settings = &
        solver_settings_type()
    integer(ip), parameter :: standalone_n_random_trial_vectors = 20

contains

    subroutine solver(update_orbs, obj_func, n_param, error, settings, state)
        !
        ! this subroutine is the main solver for orbital optimization
        !
        procedure(update_orbs_type), intent(in), pointer :: update_orbs
        procedure(obj_func_type), intent(in), pointer :: obj_func
        integer(ip), intent(in) :: n_param
        integer(ip), intent(out) :: error
        type(solver_settings_type), intent(in), optional :: settings
        type(solver_state_type), intent(out), target, optional :: state

        real(rp) :: func, grad_norm, grad_rms, mu, n_kappa, kappa_norm
        real(rp), allocatable :: kappa(:), grad(:), h_diag(:)
        logical :: macro_converged, stable, conv_check_passed, local_region
        integer(ip) :: imacro
        character(len=300) :: msg
        procedure(hess_x_type), pointer :: hess_x_funptr
        type(solver_settings_type) :: local_settings
        type(solver_state_type), target :: local_state
        type(solver_state_type), pointer :: current_state
        real(rp), external :: dnrm2

        ! initialize error flag
        error = 0

        ! the first subproblem is perturbed as one in the global region, every later
        ! one as one in the region the previous subproblem was solved in
        local_region = .false.

        ! keep the state in the given record, so that it is current during the
        ! optimization, or in a local one, both of which start from the default values
        if (present(state)) then
            current_state => state
        else
            current_state => local_state
        end if

        ! work on a copy of the settings, or of the default settings if none are
        ! passed, which leaves the caller's settings as they are while the copy is
        ! completed by the defaults that depend on other settings
        if (present(settings)) then
            local_settings = settings
        else
            local_settings = default_solver_settings
        end if

        ! initialize macroiteration convergence
        macro_converged = .false.

        ! initialize stability boolean
        stable = .true.

        ! initialize random number generator
        call init_rng(local_settings%seed)

        ! resolve the defaults that depend on other settings
        call init_defaults(local_settings)

        ! perform sanity check
        call solver_sanity_check(local_settings, n_param, error)
        call add_error_origin(error, error_solver, local_settings)
        if (error /= 0) return

        ! initialize starting trust radius
        current_state%trust_radius = local_settings%start_trust_radius

        ! print header
        call local_settings%log(repeat("-", 109), verbosity_info)
        call local_settings%log( &
            " Iteration |     Objective function     | Gradient RMS |"// &
            " Level shift |   Micro    | Trust radius | Step size ", verbosity_info)
        call local_settings%log( &
            "           |                            |              |"// &
            "             | iterations |              |           ", verbosity_info)
        call local_settings%log(repeat("-", 109), verbosity_info)

        ! allocate arrays
        allocate(kappa(n_param), grad(n_param), h_diag(n_param))

        ! initialize orbital rotation matrix
        kappa = 0.0_rp

        do imacro = 1, local_settings%n_macro
            if (.not. current_state%max_precision_reached) then
                ! calculate cost function, gradient and Hessian diagonal, the Hessian
                ! linear transformation is disassociated first so that an orbital
                ! update which does not provide one is caught
                hess_x_funptr => null()
                call update_orbs(kappa, func, grad, h_diag, hess_x_funptr, error, &
                                 local_settings%context)

                ! increment number of orbital updates
                current_state%n_update_orbs = current_state%n_update_orbs + 1
                call add_error_origin(error, error_update_orbs, local_settings)
                if (error /= 0) return

                ! an orbital update which succeeds without providing a Hessian linear
                ! transformation violates the interface, which the solver reports itself
                if (.not. associated(hess_x_funptr)) then
                    call local_settings%log("Orbital update did not provide a "// &
                                            "Hessian linear transformation.", &
                                            verbosity_error, .true.)
                    error = error_solver + 1
                    return
                end if

                ! calculate gradient norm
                grad_norm = dnrm2(n_param, grad, 1_ip)

                ! calculate RMS gradient
                grad_rms = grad_norm / sqrt(real(n_param, kind=rp))

                ! an exactly vanishing gradient, which a symmetric problem can reach,
                ! leaves the subsystem solvers no subproblem they can solve, so the
                ! point is taken as one at maximum precision
                if (grad_norm <= 0.0_rp) current_state%max_precision_reached = .true.

                ! update the state with the current point, whose row is printed below
                current_state%macro_iteration = imacro - 1
                current_state%func = func
                current_state%grad_rms = grad_rms

                ! log results
                if (imacro == 1) then
                    call local_settings%print_results(imacro - 1, func, grad_rms)
                else if (.not. stable) then
                    call local_settings%print_results(imacro - 1, func, grad_rms, &
                                                      kappa_norm=kappa_norm)
                else if (current_state%jacobi_davidson_iteration > 0) then
                    call local_settings%print_results( &
                        imacro - 1, func, grad_rms, &
                        level_shift=current_state%level_shift, &
                        n_micro=current_state%micro_iteration, &
                        n_jacobi_davidson=current_state%jacobi_davidson_iteration, &
                        trust_radius=current_state%trust_radius, kappa_norm=kappa_norm)
                else
                    call local_settings%print_results( &
                        imacro - 1, func, grad_rms, &
                        level_shift=current_state%level_shift, &
                        n_micro=current_state%micro_iteration, &
                        trust_radius=current_state%trust_radius, kappa_norm=kappa_norm)
                end if
            end if

            ! check for convergence and stability
            if (associated(local_settings%conv_check)) then
                conv_check_passed = &
                    local_settings%conv_check(error, local_settings%context)
                call add_error_origin(error, error_conv_check, local_settings)
                if (error /= 0) return
            else
                conv_check_passed = grad_rms < local_settings%conv_tol
            end if

            ! a point reached by a step along an unstable mode is not taken as
            ! converged, since this step only sets up the trust-region step that leaves
            ! the stationary point and changes the objective function too little for a
            ! convergence check of this change
            if (.not. stable) conv_check_passed = .false.
            stable = .true.
            if (conv_check_passed .or. current_state%max_precision_reached) then
                ! always perform stability check if starting at stationary point or at
                ! an exactly vanishing gradient and continue along the unstable mode if
                ! the stationary point is unstable
                if (local_settings%stability .or. imacro == 1 .or. &
                    grad_norm <= 0.0_rp) then
                    call check_stationary_point(grad, grad_rms, h_diag, hess_x_funptr, &
                                                imacro, local_settings, current_state, &
                                                stable, kappa, kappa_norm, error)
                    if (error /= 0) return
                    if (.not. stable) then
                        current_state%micro_iteration = 0
                        current_state%jacobi_davidson_iteration = 0
                        current_state%relative_residual = 0.0_rp
                        current_state%level_shift = 0.0_rp
                        current_state%step_size = kappa_norm
                        current_state%max_precision_reached = .false.
                        cycle
                    end if
                end if
                ! convergence takes precedence over maximum precision
                if (conv_check_passed) current_state%max_precision_reached = .false.
                macro_converged = .true.
                exit
            end if

            ! stop before solving the trust region subproblem in the last macro
            ! iteration since its step would never be evaluated
            if (imacro == local_settings%n_macro) exit

            if (string_in("davidson", local_settings%subsystem_solver)) then
                ! solve trust region subproblem with (Jacobi-)Davidson
                call level_shifted_davidson(func, grad, grad_norm, h_diag, n_param, &
                                            obj_func, hess_x_funptr, local_settings, &
                                            local_region, current_state, kappa, &
                                            kappa_norm, mu, error)
                call add_error_origin(error, error_solver, local_settings)
                if (error /= 0) return

                ! remember the region the subproblem was solved in, which sets the size
                ! of the gradient perturbation of the next one, and record the level
                ! shift of the step
                local_region = abs(mu) < level_shift_local_thres
                current_state%level_shift = merge(-mu, 0.0_rp, abs(mu) > 0.0_rp)
            else if (local_settings%subsystem_solver == "tcg") then
                ! solve trust region subproblem with truncated conjugate gradient
                call truncated_conjugate_gradient( &
                    func, grad, h_diag, n_param, obj_func, hess_x_funptr, &
                    local_settings, current_state, kappa, kappa_norm, error)
                call add_error_origin(error, error_solver, local_settings)
                if (error /= 0) return
            end if

            ! perform line search
            if (current_state%max_precision_reached) then
                n_kappa = 0.0_rp
            else if (local_settings%line_search) then
                n_kappa = bracket(obj_func, kappa, 0.0_rp, 1.0_rp, local_settings, &
                                  current_state%n_obj_func, error)
                call add_error_origin(error, error_solver, local_settings)
                if (error /= 0) return
            else
                n_kappa = 1.0_rp
            end if

            ! set orbital rotation
            kappa = n_kappa * kappa
            kappa_norm = n_kappa * kappa_norm
            current_state%step_size = kappa_norm

            ! flush output
            flush(stdout)
            flush(stderr)

        end do

        ! deallocate arrays
        deallocate(kappa, grad, h_diag)

        ! stop if no convergence
        if (.not. macro_converged) then
            call local_settings%log("Orbital optimization has not converged!", &
                                    verbosity_error, .true.)
            error = error_solver_max_iter
            return
        end if

        ! finish logging
        call local_settings%log(repeat("-", 109), verbosity_info)
        write(msg, '(A, I0)') "Total number of Hessian linear transformations: ", &
            current_state%n_hess_x
        call local_settings%log(msg, verbosity_info)
        write(msg, '(A, I0)') "Total number of orbital updates: ", &
            current_state%n_update_orbs
        call local_settings%log(msg, verbosity_info)

        ! flush output
        flush(stdout)
        flush(stderr)

    end subroutine solver

    subroutine stability_check(h_diag, hess_x_funptr, stable, error, settings, state, &
                               kappa)
        !
        ! this subroutine performs a stability check
        !
        real(rp), intent(in) :: h_diag(:)
        procedure(hess_x_type), intent(in), pointer :: hess_x_funptr
        logical, intent(out) :: stable
        integer(ip), intent(out) :: error
        type(stability_settings_type), intent(in), optional :: settings
        type(stability_state_type), intent(out), target, optional :: state
        real(rp), intent(out), optional :: kappa(:)

        integer(ip) :: n_param
        real(rp), allocatable :: red_space_basis(:, :), eigvals(:), eigvecs(:, :)
        real(rp) :: eigval
        character(len=300) :: msg
        logical :: stability_converged
        type(stability_settings_type) :: local_settings
        type(stability_state_type), target :: local_state
        type(stability_state_type), pointer :: current_state

        ! initialize error flag
        error = 0

        ! initialize stable
        stable = .false.

        ! keep the state in the given record, so that it is current during the
        ! stability check, or in a local one, both of which start from the default
        ! values
        if (present(state)) then
            current_state => state
        else
            current_state => local_state
        end if

        ! work on a copy of the settings, or of the default settings if none are
        ! passed, which leaves the caller's settings as they are
        if (present(settings)) then
            local_settings = settings
        else
            local_settings = default_stability_settings
        end if

        ! initialize random number generator
        call init_rng(local_settings%seed)

        ! get number of parameters
        n_param = size(h_diag)

        ! perform sanity check
        call stability_sanity_check(local_settings, n_param, error)
        call add_error_origin(error, error_stability_check, local_settings)
        if (error /= 0) return

        ! generate trial vectors
        call get_stability_trial_space(h_diag, red_space_basis, local_settings, error)
        call add_error_origin(error, error_stability_check, local_settings)
        if (error /= 0) return

        ! determine the lowest eigenpair of the Hessian
        allocate(eigvals(1), eigvecs(n_param, 1))
        call block_davidson(hess_x_funptr, h_diag, red_space_basis, &
                            local_settings%hess_symm, 1_ip, local_settings%n_iter, &
                            eigvals, eigvecs, stability_converged, local_settings, &
                            current_state, error)
        call add_error_origin(error, error_stability_check, local_settings)
        if (error /= 0) return
        eigval = eigvals(1)

        ! check if stability check has converged, which it fails to do either within
        ! the given number of iterations or, before these are used up, once its trial
        ! space can not be expanded any further
        if (.not. stability_converged) then
            if (current_state%iteration < local_settings%n_iter) then
                call local_settings%log("Stability check has not converged since "// &
                                        "its trial space can not be expanded any "// &
                                        "further.", verbosity_error, .true.)
            else
                call local_settings%log("Stability check has not converged in the "// &
                                        "given number of iterations.", &
                                        verbosity_error, .true.)
            end if
            error = error_stability_check_max_iter
            return
        end if

        ! return lowest eigenvector, whose eigenvalue the state holds
        if (present(kappa)) kappa = eigvecs(:, 1)

        ! determine if saddle point
        stable = eigval > local_settings%stability_thresh
        if (.not. stable) then
            write(msg, '(A, F0.4)') unstable_warning_msg, eigval
            call local_settings%log(msg, verbosity_warning)
        end if

        ! deallocate quantities from Davidson iterations
        deallocate(red_space_basis, eigvals, eigvecs)

        ! flush output
        flush(stdout)
        flush(stderr)

    end subroutine stability_check

    subroutine check_stationary_point(grad, grad_rms, h_diag, hess_x_funptr, imacro, &
                                      settings, state, stable, kappa, kappa_norm, error)
        !
        ! this subroutine performs the stability check at a stationary point with the
        ! nested stability check settings and, if the stationary point is unstable,
        ! returns the step along the unstable mode
        !
        real(rp), intent(in) :: grad(:), grad_rms, h_diag(:)
        procedure(hess_x_type), intent(in), pointer :: hess_x_funptr
        integer(ip), intent(in) :: imacro
        type(solver_settings_type), intent(in) :: settings
        type(solver_state_type), intent(inout), target :: state
        logical, intent(out) :: stable
        real(rp), intent(out) :: kappa(:), kappa_norm
        integer(ip), intent(out) :: error

        ! initialize error flag
        error = 0

        ! perform stability check and add its Hessian linear transformations to the
        ! solver's
        call stability_check(h_diag, hess_x_funptr, stable, error, &
                             settings%stability_settings, state=state%stability_state, &
                             kappa=kappa)
        state%n_hess_x = state%n_hess_x + state%stability_state%n_hess_x
        if (error /= 0) return
        if (stable) return

        ! move far enough along the unstable mode, whose eigenvalue the state of the
        ! stability check holds, so that gradient is increased by one order of
        ! magnitude
        call saddle_point_step(grad, grad_rms, state%stability_state%min_eigval, &
                               state%trust_radius, settings, kappa, kappa_norm)

        ! warn of the saddle point, advising to perturb the gradient where it is not
        ! perturbed
        if (imacro == 1) then
            call settings%log(started_at_saddle_point_warning_msg, verbosity_warning)
        else if (settings%grad_noise) then
            call settings%log(reached_saddle_point_warning_msg, verbosity_warning)
        else
            call settings%log(reached_saddle_point_no_grad_noise_warning_msg, &
                              verbosity_warning)
        end if

    end subroutine check_stationary_point

    subroutine saddle_point_step(grad, grad_rms, min_eigval, trust_radius, settings, &
                                 kappa, kappa_norm)
        !
        ! this subroutine determines the step along the direction of negative curvature
        ! such that the (linearized) gradient RMS increases to one order of magnitude
        ! above its current value, or above the convergence threshold if the gradient
        ! is already converged, or, if both vanish, the step to the trust radius, and
        ! limits it to the trust radius
        !
        real(rp), intent(in) :: grad(:), grad_rms, min_eigval, trust_radius
        type(solver_settings_type), intent(in) :: settings
        real(rp), intent(inout) :: kappa(:)
        real(rp), intent(out) :: kappa_norm

        integer(ip) :: n_param
        real(rp) :: grad_kappa, n_kappa, target_grad_rms
        real(rp), external :: ddot, dnrm2

        ! number of parameters
        n_param = size(grad)

        ! project gradient onto direction of negative curvature
        grad_kappa = ddot(n_param, kappa, 1_ip, grad, 1_ip)

        ! target gradient RMS is one order of magnitude above the current value, or
        ! above the convergence threshold if the gradient is already converged, to
        ! ensure the step is always large enough to escape the stationary point
        target_grad_rms = 10 * max(settings%conv_tol, grad_rms)

        if (target_grad_rms > 0.0_rp) then
            ! determine step length that reaches the target gradient RMS in the
            ! linearized model along the direction of negative curvature
            n_kappa = (-grad_kappa + sign(1.0_rp, grad_kappa) * sqrt( &
                grad_kappa**2 + real(n_param, kind=rp) * &
                (target_grad_rms**2 - grad_rms**2))) / min_eigval
        else
            ! an exactly vanishing gradient with a vanishing convergence threshold
            ! leaves no target, so the step reaches the trust radius instead, as a
            ! trust region step along the direction of negative curvature would
            n_kappa = trust_radius
        end if

        ! the step is taken without evaluating the objective function, so it is
        ! limited to the trust radius
        n_kappa = sign(min(abs(n_kappa), trust_radius), n_kappa)

        ! scale step
        kappa = n_kappa * kappa
        kappa_norm = dnrm2(n_param, kappa, 1_ip)

    end subroutine saddle_point_step

    subroutine newton_step(grad_norm, red_space_basis, red_space_hess_eigvals, &
                           red_space_hess_right_eigvecs, red_space_hess_left_eigvecs, &
                           solution, red_space_solution)
        !
        ! this subroutine performs a Newton step by solving the Newton equations in
        ! reduced space without a level shift
        !
        real(rp), intent(in) :: grad_norm, red_space_basis(:, :)
        real(rp), intent(in) :: red_space_hess_eigvals(:), &
                                red_space_hess_right_eigvecs(:, :), &
                                red_space_hess_left_eigvecs(:, :)
        real(rp), intent(out) :: solution(:), red_space_solution(:)

        integer(ip) :: n_param, n_red
        real(rp), allocatable :: eigspace_solution(:)
        external :: dgemv

        ! number of parameters
        n_param = size(solution)

        ! reduced space size
        n_red = size(red_space_basis, 2)

        ! construct Newton step in eigenvector basis, in which an eigenvector whose
        ! gradient component is indistinguishable from zero contributes nothing, also
        ! where its eigenvalue vanishes
        allocate(eigspace_solution(n_red))
        eigspace_solution = 0.0_rp
        where (abs(red_space_hess_left_eigvecs(1, :)) > numerical_zero) &
            eigspace_solution = -grad_norm * red_space_hess_left_eigvecs(1, :) / &
                                red_space_hess_eigvals

        ! transform from eigenvector basis to reduced step basis
        call dgemv("N", n_red, n_red, 1.0_rp, red_space_hess_right_eigvecs, n_red, &
                   eigspace_solution, 1_ip, 0.0_rp, red_space_solution, 1_ip)

        ! get solution in full space
        call dgemv("N", n_param, n_red, 1.0_rp, red_space_basis, n_param, &
                   red_space_solution, 1_ip, 0.0_rp, solution, 1_ip)

    end subroutine newton_step

    subroutine bisection_ah(red_space_hess, grad_norm, red_space_basis, &
                            red_space_hess_eigvals, red_space_hess_right_eigvecs, &
                            red_space_hess_left_eigvecs, trust_radius, solution, &
                            red_space_solution, mu, settings, error)
        !
        ! this subroutine performs bisection to find the parameter alpha that matches
        ! the desired trust radius
        !
        real(rp), intent(in) :: red_space_hess(:, :)
        real(rp), intent(in) :: grad_norm, red_space_basis(:, :), &
                                red_space_hess_eigvals(:), &
                                red_space_hess_right_eigvecs(:, :), &
                                red_space_hess_left_eigvecs(:, :), trust_radius
        type(solver_settings_type), intent(in) :: settings
        real(rp), intent(out) :: solution(:), red_space_solution(:), mu
        integer(ip), intent(out) :: error

        real(rp), allocatable :: aug_hess(:, :), eigspace_solution(:)
        real(rp) :: lower_alpha, middle_alpha, upper_alpha, lower_trust_dist, &
                    middle_trust_dist, upper_trust_dist, current_norm
        logical, allocatable :: non_degenerate_mask(:)
        integer(ip) :: n_param, n_red, n_aug, iter, min_idx
        real(rp), parameter :: &
            lower_alpha_bound = 1e-30_rp, upper_alpha_bound = 1e30_rp, &
            alpha_conv_factor = 1e-12_rp, orthogonality_thres = 1e-12_rp
        real(rp), external :: dnrm2

        ! initialize error flag
        error = 0

        ! number of parameters
        n_param = size(solution)

        ! reduced space size
        n_red = size(red_space_basis, 2)

        ! augmented Hessian size
        n_aug = n_red + 1

        ! construct augmented Hessian
        allocate(aug_hess(n_aug, n_aug))
        aug_hess = 0.0_rp
        aug_hess(2:n_aug, 2:n_aug) = red_space_hess

        ! lower and upper bracket for alpha
        lower_alpha = lower_alpha_bound
        upper_alpha = upper_alpha_bound

        ! check if lowest reduced space Hessian eigenvalue is negative and its
        ! (possibly degenerate) eigenspace has no gradient component, i.e. whether this
        ! is the hard case
        min_idx = minloc(red_space_hess_eigvals, dim=1)
        non_degenerate_mask = abs(red_space_hess_eigvals - &
                                  red_space_hess_eigvals(min_idx)) > numerical_zero
        if (red_space_hess_eigvals(min_idx) < 0.0_rp .and. &
            sqrt(sum(red_space_hess_left_eigvecs(1, :)**2, &
                     mask=.not. non_degenerate_mask)) <= orthogonality_thres) then
            ! get crossover point between lowest reduced space Hessian eigenvalue and
            ! second lowest Hessian eigenvalue in augmented Hessian to get lower
            ! boundary for alpha which ensures that the solution has a gradient
            ! component, excluding eigenvalues degenerate with the lowest one to avoid
            ! dividing by zero
            lower_alpha = sqrt(red_space_hess_eigvals(min_idx) / sum( &
                red_space_hess_left_eigvecs(1, :) * &
                red_space_hess_right_eigvecs(1, :) / &
                (red_space_hess_eigvals(min_idx) - red_space_hess_eigvals), &
                mask=non_degenerate_mask)) / grad_norm

            ! construct solution at crossover point in eigenvector basis while avoiding
            ! contributions of lowest Hessian eigenvalue as the gradient component is
            ! vanishing anyways
            eigspace_solution = merge( &
                grad_norm * red_space_hess_left_eigvecs(1, :) / &
                (red_space_hess_eigvals(min_idx) - red_space_hess_eigvals), 0.0_rp, &
                non_degenerate_mask)

            ! transform from eigenvector basis to reduced step basis
            call dgemv("N", n_red, n_red, 1.0_rp, red_space_hess_right_eigvecs, n_red, &
                       eigspace_solution, 1_ip, 0.0_rp, red_space_solution, 1_ip)
            deallocate(eigspace_solution, non_degenerate_mask)

            ! get solution in full space
            call dgemv("N", n_param, n_red, 1.0_rp, red_space_basis, n_param, &
                       red_space_solution, 1_ip, 0.0_rp, solution, 1_ip)
            current_norm = dnrm2(n_param, solution, 1_ip)
            lower_trust_dist = current_norm - trust_radius

            ! check if step lies inside trust region at lower boundary
            if (current_norm < trust_radius) then
                ! fill the rest of the trust radius with the lowest eigenvector
                red_space_solution = &
                    red_space_solution + &
                    sqrt(max(0.0_rp, trust_radius**2 - current_norm**2)) * &
                    red_space_hess_right_eigvecs(:, min_idx)

                ! construct full space solution
                call dgemv("N", n_param, n_red, 1.0_rp, red_space_basis, n_param, &
                           red_space_solution, 1_ip, 0.0_rp, solution, 1_ip)

                ! set level shift
                mu = red_space_hess_eigvals(min_idx)
                deallocate(aug_hess)
                return
            end if
        else
            ! start at Cauchy point to get bracket
            middle_alpha = grad_norm / trust_radius
            call get_ah_lowest_eigenvec(middle_alpha)
            if (error /= 0) return
            current_norm = dnrm2(n_param, solution, 1_ip)
            middle_trust_dist = current_norm - trust_radius
            lower_alpha = middle_alpha
            lower_trust_dist = middle_trust_dist
        end if

        ! start bracketing
        upper_alpha = lower_alpha
        upper_trust_dist = lower_trust_dist

        ! step is too long: alpha is lower bound
        if (current_norm > trust_radius) then
            ! larger alpha is necessary to find upper bound
            do while (current_norm > trust_radius)
                lower_alpha = upper_alpha
                lower_trust_dist = upper_trust_dist
                upper_alpha = upper_alpha * 10.0_rp
                call get_ah_lowest_eigenvec(upper_alpha)
                if (error /= 0) return
                current_norm = dnrm2(n_param, solution, 1_ip)
                upper_trust_dist = current_norm - trust_radius
                if (upper_alpha > upper_alpha_bound) then
                    call settings%log("Unable to find upper bound for alpha in "// &
                                      "bisection.", verbosity_error, .true.)
                    error = 1
                    return
                end if
            end do
        ! step is too short: alpha is upper bound
        else
            ! smaller alpha is necessary to find lower bound
            do while (current_norm < trust_radius)
                upper_alpha = lower_alpha
                upper_trust_dist = lower_trust_dist
                lower_alpha = lower_alpha / 10.0_rp
                call get_ah_lowest_eigenvec(lower_alpha)
                if (error /= 0) return
                current_norm = dnrm2(n_param, solution, 1_ip)
                lower_trust_dist = current_norm - trust_radius
                if (lower_alpha < lower_alpha_bound) then
                    ! the step norm plateaus below the trust radius as alpha approaches
                    ! zero, which means the unconstrained Newton step already lies
                    ! within the trust region
                    call newton_step( &
                        grad_norm, red_space_basis, red_space_hess_eigvals, &
                        red_space_hess_right_eigvecs, red_space_hess_left_eigvecs, &
                        solution, red_space_solution)
                    mu = 0.0_rp
                    return
                end if
            end do
        end if

        ! get middle alpha
        middle_alpha = sqrt(upper_alpha * lower_alpha)
        call get_ah_lowest_eigenvec(middle_alpha)
        if (error /= 0) return
        middle_trust_dist = dnrm2(n_param, solution, 1_ip) - trust_radius

        ! perform bisection to find root, converge to relative threshold to avoid
        ! precision issues
        iter = 0
        do while (upper_alpha - lower_alpha > alpha_conv_factor * upper_alpha)
            ! targeted trust radius is in upper bracket
            if (lower_trust_dist * middle_trust_dist > 0.0_rp) then
                lower_alpha = middle_alpha
                lower_trust_dist = middle_trust_dist
            ! targeted trust radius is in lower bracket
            else
                upper_alpha = middle_alpha
                upper_trust_dist = middle_trust_dist
            end if
            ! get new middle alpha
            middle_alpha = sqrt(upper_alpha * lower_alpha)
            call get_ah_lowest_eigenvec(middle_alpha)
            if (error /= 0) return
            middle_trust_dist = dnrm2(n_param, solution, 1_ip) - trust_radius
            ! check if maximum number of iterations is reached
            iter = iter + 1
            if (iter > 100) then
                call settings%log("Maximum number of bisection iterations reached.", &
                                  verbosity_error, .true.)
                error = 1
                deallocate(aug_hess)
                return
            end if
        end do

        ! deallocate arrays
        deallocate(aug_hess)

    contains

        subroutine get_ah_lowest_eigenvec(alpha)
            !
            ! this subroutine returns the lowest eigenvector for an augmented Hessian
            !
            real(rp), intent(in) :: alpha

            real(rp), allocatable :: eigvec(:)
            external :: dgemv

            ! finish construction of augmented Hessian
            aug_hess(1, 2) = alpha * grad_norm
            aug_hess(2, 1) = alpha * grad_norm

            ! allocate the final eigenvector array
            allocate(eigvec(n_aug))

            ! perform eigendecomposition and get lowest eigenvalue and corresponding
            ! eigenvector
            call mat_min_eig(aug_hess, settings%hess_symm, mu, eigvec, settings, error)
            if (error /= 0) return

            ! check if eigenvector has level-shift component
            if (abs(eigvec(1)) <= tiny(1.0_rp)) then
                call settings%log("Lowest augmented Hessian eigenvector does not "// &
                                  "have level-shift component.", verbosity_error, &
                                  .true.)
                error = 1
                deallocate(eigvec)
                return
            end if

            ! scale eigenvector such that first element is equal to one and divide by
            ! alpha to get solution in reduced space
            red_space_solution = eigvec(2:) / eigvec(1) / alpha
            deallocate(eigvec)

            ! get solution in full space
            call dgemv("N", n_param, n_red, 1.0_rp, red_space_basis, n_param, &
                       red_space_solution, 1_ip, 0.0_rp, solution, 1_ip)

        end subroutine get_ah_lowest_eigenvec

    end subroutine bisection_ah

    subroutine bisection_mu(red_space_hess, grad_norm, red_space_basis, &
                            red_space_hess_eigvals, red_space_hess_right_eigvecs, &
                            red_space_hess_left_eigvecs, trust_radius, solution, &
                            red_space_solution, final_mu, settings, error)
        !
        ! this subroutine performs bisection to find the level shift such that the step
        ! matches the desired trust radius
        !
        real(rp), intent(in) :: red_space_hess(:, :)
        real(rp), intent(in) :: grad_norm, red_space_basis(:, :), &
                                red_space_hess_right_eigvecs(:, :), &
                                red_space_hess_left_eigvecs(:, :), trust_radius
        complex(rp), intent(in) :: red_space_hess_eigvals(:)
        type(solver_settings_type), intent(in) :: settings
        real(rp), intent(out) :: solution(:), red_space_solution(:), final_mu
        integer(ip), intent(out) :: error

        integer(ip) :: n_param, n_red, iter, min_idx
        real(rp) :: current_norm, upper_mu, lower_mu, middle_mu, sign_factor, &
                    eigvec_overlap
        logical :: singular_start
        logical, allocatable :: non_degenerate_mask(:)
        real(rp), parameter :: orthogonality_thres = 1e-12_rp, &
                               positive_definite_shift = 1e-8_rp, &
                               alpha_conv_factor = 1e-12_rp
        external :: dgemv
        real(rp), external :: dnrm2, ddot

        ! initialize error flag
        error = 0

        ! number of parameters
        n_param = size(solution)

        ! reduced space size
        n_red = size(red_space_basis, 2)

        ! initialize singularity flag
        singular_start = .false.

        ! determine maximum level shift, a vanishing lowest eigenvalue is treated like
        ! a negative one since the level-shifted Hessian is singular at a vanishing
        ! level shift
        min_idx = minloc(red_space_hess_eigvals%re, dim=1)
        non_degenerate_mask = abs(red_space_hess_eigvals - &
                                  red_space_hess_eigvals(min_idx)) > numerical_zero
        if (red_space_hess_eigvals(min_idx)%re < numerical_zero) then
            ! check if the (possibly degenerate) eigenspace of the lowest eigenvalue
            ! has a gradient component or if the minimum eigenvalue has a significant
            ! imaginary part since both cases lead to no pole in the solution norm at
            ! the lowest eigenvalue
            if (sqrt(sum(red_space_hess_left_eigvecs(1, :)**2, &
                         mask=.not. non_degenerate_mask)) <= orthogonality_thres .or. &
                abs(red_space_hess_eigvals(min_idx)%im) >= numerical_zero) then
                upper_mu = red_space_hess_eigvals(min_idx)%re
                singular_start = .true.
            ! ensure level-shifted Hessian is positive definite
            else
                upper_mu = red_space_hess_eigvals(min_idx)%re - positive_definite_shift
            end if
        else
            upper_mu = 0.0_rp
        end if

        ! evaluate step length at the lower bound of the shift
        call get_level_shifted_step(upper_mu, singular_start)
        if (error /= 0) return

        ! check if trust radius can be reached in bisection, a step which reaches it
        ! already at this shift is taken as it is, since the bisection could not
        ! converge towards a shift of zero
        if (current_norm <= trust_radius) then
            ! check if reduced space Hessian is non-positive definite otherwise accept
            ! interior solution as Newton step
            if (red_space_hess_eigvals(min_idx)%re <= newton_eigval_thresh) then
                ! hard case occured if lowest eigenvalue is real, trust radius cannot
                ! be reached because the eigenvector corresponding to the lowest
                ! eigenvalue is orthogonal to the gradient, leading to no pole in the
                ! solution norm
                if (abs(red_space_hess_eigvals(min_idx)%im) < numerical_zero) then
                    ! ensure that we are moving in negative gradient direction by
                    ! checking the sign of the gradient component, for an asymmetric
                    ! Hessian this can be non-zero even if the left eigenvector is
                    ! orthogonal to the gradient
                    if (red_space_hess_right_eigvecs(1, min_idx) > 0.0_rp) then
                        sign_factor = -1.0_rp
                    else
                        sign_factor = 1.0_rp
                    end if

                    ! fill the rest of the trust radius with the lowest eigenvector,
                    ! along which the step already has a component unless the
                    ! gradient is orthogonal to its eigenspace, by moving along it to
                    ! the trust region boundary
                    eigvec_overlap = sign_factor * ddot( &
                        n_red, red_space_solution, 1_ip, &
                        red_space_hess_right_eigvecs(:, min_idx), 1_ip)
                    red_space_solution = red_space_solution + sign_factor * ( &
                        sqrt(eigvec_overlap**2 + trust_radius**2 - current_norm**2) - &
                        eigvec_overlap) * red_space_hess_right_eigvecs(:, min_idx)

                    ! construct full space solution
                    call dgemv("N", n_param, n_red, 1.0_rp, red_space_basis, n_param, &
                               red_space_solution, 1_ip, 0.0_rp, solution, 1_ip)

                ! trust radius cannot be reached because imaginary part leads to no
                ! pole in the solution norm
                else
                    call settings%log(bisection_imaginary_eigval_warning_msg, &
                                      verbosity_warning)
                end if

                ! set final level shift to the lowest eigenvalue
                final_mu = red_space_hess_eigvals(min_idx)%re
            else
                ! set level shift either to zero or to a very small shift for which an
                ! interior solution is still accepted
                final_mu = upper_mu
            end if

        else
            ! initialize lower bound for level shift while safe guarding against
            ! level-shift of zero for positive-definite Hessian
            lower_mu = min(-1.0_rp, upper_mu * 2.0_rp)
            call get_level_shifted_step(lower_mu, .false.)
            if (error /= 0) return

            ! find an lower bound for mu where the step falls within the trust region
            do while (current_norm > trust_radius)
                if (lower_mu < -1e15_rp) then
                    call settings%log("Unable to find lower bound for mu in "// &
                                      "bisection.", verbosity_error, .true.)
                    error = 1
                    return
                end if
                lower_mu = lower_mu * 2.0_rp
                call get_level_shifted_step(lower_mu, .false.)
                if (error /= 0) return
            end do

            ! perform direct bisection on the shift parameter mu, which is resolved
            ! relative to itself or, where it converges to zero and the bracket cannot
            ! shrink relative to it, absolutely
            iter = 0
            do while (abs(lower_mu - upper_mu) > &
                      alpha_conv_factor * max(abs(lower_mu), numerical_zero))
                middle_mu = 0.5_rp * (upper_mu + lower_mu)
                call get_level_shifted_step(middle_mu, .false.)
                if (error /= 0) return

                if (current_norm > trust_radius) then
                    upper_mu = middle_mu
                else
                    lower_mu = middle_mu
                end if

                iter = iter + 1
                if (iter > 100) then
                    call settings%log("Maximum number of bisection iterations "// &
                                      "reached.", verbosity_error, .true.)
                    error = 1
                    return
                end if
            end do

            ! get final solution at the converged shift
            call get_level_shifted_step(middle_mu, .false.)
            if (error /= 0) return

            ! set level shift
            final_mu = middle_mu
        end if

    contains

        subroutine get_level_shifted_step(mu, singular)
            !
            ! this subroutine returns the level-shifted step for a given level shift by
            ! solving the linear system, uses LU factorization for non-singular
            ! matrices and the Moore-Penrose inverse constructed from SVD for singular
            ! matrices
            !
            real(rp), intent(in) :: mu
            logical, intent(in) :: singular

            real(rp), allocatable :: level_shifted_hess(:, :), work(:), sing_values(:)
            integer(ip), allocatable :: ipiv(:)
            integer(ip) :: k, info, rank, lwork
            character(len=300) :: msg
            external :: dgesv, dgelss

            ! construct level-shifted Hessian
            level_shifted_hess = red_space_hess
            do k = 1, n_red
                level_shifted_hess(k, k) = level_shifted_hess(k, k) - mu
            end do

            ! initialize right-hand side
            red_space_solution = 0.0_rp
            red_space_solution(1) = -grad_norm

            ! solve level-shifted linear system using LU factorization if matrix is
            ! non-singular
            if (.not. singular) then
                ! solve linear system
                allocate(ipiv(n_red))
                call dgesv(n_red, 1_ip, level_shifted_hess, n_red, ipiv, &
                           red_space_solution, n_red, info)
                deallocate(level_shifted_hess, ipiv)

                ! check for successful execution
                if (info /= 0) then
                    write(msg, '(A, I0)') "Linear solve failed: Error in DGESV, "// &
                        "info = ", info
                    call settings%log(msg, verbosity_error, .true.)
                    error = 1
                    return
                end if

            ! solve level-shifted linear system using SVD if matrix is singular
            else
                ! query optimal workspace size
                allocate(sing_values(n_red), work(1))
                lwork = -1
                call dgelss(n_red, n_red, 1_ip, level_shifted_hess, n_red, &
                            red_space_solution, n_red, sing_values, numerical_zero, &
                            rank, work, lwork, info)
                lwork = int(work(1))
                deallocate(work)
                allocate(work(lwork))

                ! solve linear system
                call dgelss(n_red, n_red, 1_ip, level_shifted_hess, n_red, &
                            red_space_solution, n_red, sing_values, numerical_zero, &
                            rank, work, lwork, info)
                deallocate(level_shifted_hess, work, sing_values)

                ! check for successful execution
                if (info /= 0) then
                    write(msg, '(A, I0)') "Linear solve failed: Error in DGELSS, "// &
                        "info = ", info
                    call settings%log(msg, verbosity_error, .true.)
                    error = 1
                    return
                end if
            end if

            ! construct full space solution
            call dgemv("N", n_param, n_red, 1.0_rp, red_space_basis, n_param, &
                       red_space_solution, 1_ip, 0.0_rp, solution, 1_ip)

            ! calculate norm
            current_norm = dnrm2(n_param, solution, 1_ip)

        end subroutine get_level_shifted_step

    end subroutine bisection_mu

    function bracket(obj_func, kappa, lower, upper, settings, n_obj_func, error) &
        result(n_kappa)
        !
        ! this function brackets a minimum (algorithm from numerical recipes)
        !
        procedure(obj_func_type), intent(in), pointer :: obj_func
        real(rp), intent(in) :: kappa(:), lower, upper
        class(settings_type), intent(in) :: settings
        integer(ip), intent(inout), target :: n_obj_func
        integer(ip), intent(out) :: error

        real(rp) :: n_kappa, f_upper, f_lower, n_a, n_b, n_c, n_u, f_a, f_b, f_c, f_u, &
                    n_u_lim, tmp1, tmp2, val, denom
        integer(ip) :: iter
        real(rp), parameter :: golden_ratio = (1.0_rp + sqrt(5.0_rp)) / 2.0_rp, &
                               grow_limit = 110.0_rp, denom_floor = 1e-21_rp

        ! initialize error flag
        error = 0

        ! initialize step length
        n_kappa = 0.0_rp

        ! evaluate function at upper and lower bounds
        f_lower = obj_func(lower * kappa, error, settings%context)
        n_obj_func = n_obj_func + 1
        call add_error_origin(error, error_obj_func, settings)
        if (error /= 0) return
        f_upper = obj_func(upper * kappa, error, settings%context)
        n_obj_func = n_obj_func + 1
        call add_error_origin(error, error_obj_func, settings)
        if (error /= 0) return

        ! ensure f_a > f_b
        if (f_upper > f_lower) then
            n_a = upper
            n_b = lower
            f_a = f_upper
            f_b = f_lower
        else
            n_a = lower
            n_b = upper
            f_a = f_lower
            f_b = f_upper
        end if

        ! default step
        n_c = n_b + golden_ratio * (n_b - n_a)
        f_c = obj_func(n_c * kappa, error, settings%context)
        n_obj_func = n_obj_func + 1
        call add_error_origin(error, error_obj_func, settings)
        if (error /= 0) return

        ! continue looping until function at middle point is lower than at brackets
        iter = 0
        do while (f_c < f_b)
            ! compute value u by parabolic extrapolation
            tmp1 = (n_b - n_a) * (f_b - f_c)
            tmp2 = (n_b - n_c) * (f_b - f_a)
            val = tmp2 - tmp1
            denom = 2.0_rp * sign(max(abs(val), denom_floor), val)
            n_u = n_b - ((n_b - n_c) * tmp2 - (n_b - n_a) * tmp1) / denom

            ! maximum growth in parabolic fit
            n_u_lim = n_b + grow_limit * (n_c - n_b)

            ! check if u is between n_b and n_c
            if ((n_u - n_c) * (n_b - n_u) > 0.0_rp) then
                ! evaluate function at n_u
                f_u = obj_func(n_u * kappa, error, settings%context)
                n_obj_func = n_obj_func + 1
                call add_error_origin(error, error_obj_func, settings)
                if (error /= 0) return

                ! check if minimum between n_b and n_c
                if (f_u < f_c) then
                    n_a = n_b
                    n_b = n_u
                    f_a = f_b
                    f_b = f_u
                    exit
                ! check if minium between n_a and n_u
                else if (f_u > f_b) then
                    n_c = n_u
                    f_c = f_u
                    exit
                end if

                ! parabolic fit did not help, default step
                n_u = n_c + golden_ratio * (n_c - n_b)
                f_u = obj_func(n_u * kappa, error, settings%context)
                n_obj_func = n_obj_func + 1
                call add_error_origin(error, error_obj_func, settings)
                if (error /= 0) return
            ! limit parabolic fit to its maximum allowed value
            else if ((n_u - n_u_lim) * (n_u_lim - n_c) >= 0.0_rp) then
                n_u = n_u_lim
                f_u = obj_func(n_u * kappa, error, settings%context)
                n_obj_func = n_obj_func + 1
                call add_error_origin(error, error_obj_func, settings)
                if (error /= 0) return
            ! parabolic fit is between n_c and its allowed limit
            else if ((n_u - n_u_lim) * (n_c - n_u) > 0.0_rp) then
                ! evaluate function at n_u
                f_u = obj_func(n_u * kappa, error, settings%context)
                n_obj_func = n_obj_func + 1
                call add_error_origin(error, error_obj_func, settings)
                if (error /= 0) return

                if (f_u < f_c) then
                    n_b = n_c
                    n_c = n_u
                    n_u = n_c + golden_ratio * (n_c - n_b)
                    f_b = f_c
                    f_c = f_u
                    f_u = obj_func(n_u * kappa, error, settings%context)
                    n_obj_func = n_obj_func + 1
                    call add_error_origin(error, error_obj_func, settings)
                    if (error /= 0) return
                end if
            ! reject parabolic fit and use default step
            else
                n_u = n_c + golden_ratio * (n_c - n_b)
                f_u = obj_func(n_u * kappa, error, settings%context)
                n_obj_func = n_obj_func + 1
                call add_error_origin(error, error_obj_func, settings)
                if (error /= 0) return
            end if
            ! remove oldest point
            n_a = n_b
            n_b = n_c
            n_c = n_u
            f_a = f_b
            f_b = f_c
            f_c = f_u
            ! check if maximum number of iterations is reached
            iter = iter + 1
            if (iter > 100) then
                call settings%log("Maximum number of bracketing iterations reached.", &
                                  verbosity_error, .true.)
                error = 1
                return
            end if
        end do

        ! check if miniumum is bracketed
        if (.not. ( &
            ((f_b < f_c .and. f_b <= f_a) .or. (f_b < f_a .and. f_b <= f_c)) .and. &
            ((n_a < n_b .and. n_b < n_c) .or. (n_c < n_b .and. n_b < n_a)))) then
            call settings%log("Line search did not find minimum", verbosity_error, &
                              .true.)
            error = 1
            return
        end if

        ! set new multiplier
        n_kappa = n_b

    end function bracket

    subroutine extend_matrix(matrix, column, row)
        !
        ! this subroutine extends a matrix with a column and a row, whose last elements
        ! have to agree, since both fill the new corner
        !
        real(rp), allocatable, intent(inout) :: matrix(:, :)
        real(rp), intent(in) :: column(:), row(:)

        real(rp), allocatable :: new_matrix(:, :)
        integer(ip) :: nrows, ncols

        nrows = size(matrix, 1, kind=ip)
        ncols = size(matrix, 2, kind=ip)

        ! allocate the larger matrix and copy existing contents only once
        allocate(new_matrix(nrows + 1, ncols + 1))
        new_matrix(:nrows, :ncols) = matrix

        ! add new column and row
        new_matrix(:, ncols + 1) = column
        new_matrix(nrows + 1, :) = row

        ! transfer allocation back into matrix; old storage is freed by move_alloc
        call move_alloc(new_matrix, matrix)

    end subroutine extend_matrix

    subroutine add_column(matrix, new_col)
        !
        ! this subroutine adds a new column to an array
        !
        real(rp), allocatable, intent(inout) :: matrix(:, :)
        real(rp), intent(in) :: new_col(:)

        real(rp), allocatable :: new_matrix(:, :)
        integer(ip) :: nrows, ncols

        nrows = size(matrix, 1, kind=ip)
        ncols = size(matrix, 2, kind=ip)

        ! allocate the wider matrix and copy existing contents only once
        allocate(new_matrix(nrows, ncols + 1))
        new_matrix(:, :ncols) = matrix

        ! add new column
        new_matrix(:, ncols + 1) = new_col

        ! transfer allocation back into matrix; old storage is freed by move_alloc
        call move_alloc(new_matrix, matrix)

    end subroutine add_column

    subroutine mat_diag(matrix, symm_matrix, eigvals, right_eigvecs, left_eigvecs, &
                        settings, error)
        !
        ! this subroutine returns the eigenvalues and eigenvectors of a matrix
        !
        real(rp), intent(in) :: matrix(:, :)
        logical, intent(in) :: symm_matrix
        class(settings_type), intent(in) :: settings
        complex(rp), intent(out) :: eigvals(:)
        real(rp), intent(out) :: right_eigvecs(:, :), left_eigvecs(:, :)
        integer(ip), intent(out) :: error

        real(rp), allocatable :: eigvals_re(:)

        if (symm_matrix) then
            allocate(eigvals_re(size(eigvals)))
            call symm_mat_diag(matrix, eigvals_re, right_eigvecs, settings, error)
            eigvals%re = eigvals_re
            eigvals%im = 0.0_rp
            deallocate(eigvals_re)
            left_eigvecs = right_eigvecs
        else
            call general_mat_diag(matrix, eigvals, right_eigvecs, left_eigvecs, &
                                  settings, error)
        end if

    end subroutine mat_diag

    subroutine symm_mat_diag(symm_matrix, eigvals, eigvecs, settings, error)
        !
        ! this subroutine returns eigenvalues and eigenvectors of a symmetric matrix
        !
        real(rp), intent(in) :: symm_matrix(:, :)
        class(settings_type), intent(in) :: settings
        real(rp), intent(out) :: eigvals(:), eigvecs(:, :)
        integer(ip), intent(out) :: error

        integer(ip) :: n, lwork, info
        real(rp), allocatable :: work(:)
        character(len=300) :: msg
        external :: dsyev

        ! initialize error flag
        error = 0

        ! size of matrix
        n = size(symm_matrix, 1)

        ! copy matrix
        eigvecs = symm_matrix

        ! query optimal workspace size
        lwork = -1
        allocate(work(1))
        call dsyev("V", "U", n, eigvecs, n, eigvals, work, lwork, info)
        lwork = int(work(1), kind=ip)
        deallocate(work)
        allocate(work(lwork))

        ! perform eigendecomposition
        call dsyev("V", "U", n, eigvecs, n, eigvals, work, lwork, info)

        ! deallocate work array
        deallocate(work)

        ! check for successful execution
        if (info /= 0) then
            write(msg, '(A, I0)') "Eigendecomposition failed: Error in DSYEV, "// &
                "info = ", info
            call settings%log(msg, verbosity_error, .true.)
            error = 1
            return
        end if

    end subroutine symm_mat_diag

    subroutine general_mat_diag(matrix, eigvals, right_eigvecs, left_eigvecs, &
                                settings, error)
        !
        ! this subroutine returns the eigenvalues and eigenvectors of a general matrix
        !
        real(rp), intent(in) :: matrix(:, :)
        class(settings_type), intent(in) :: settings
        complex(rp), intent(out) :: eigvals(:)
        real(rp), intent(out) :: right_eigvecs(:, :), left_eigvecs(:, :)
        integer(ip), intent(out) :: error

        integer(ip) :: n, lwork, info, i
        real(rp), allocatable :: eigvals_re(:), eigvals_im(:), work(:), temp(:, :), &
                                 left_re(:)
        real(rp) :: overlap_re, overlap_im, overlap_sq
        character(len=300) :: msg
        real(rp), external :: ddot
        external :: dgeev

        ! initialize error flag
        error = 0

        ! size of matrix
        n = size(matrix, 1)

        ! copy matrix to avoid modification of original matrix
        temp = matrix

        ! query optimal workspace size
        lwork = -1
        allocate(eigvals_re(n), eigvals_im(n), work(1))
        call dgeev("V", "V", n, temp, n, eigvals_re, eigvals_im, left_eigvecs, n, &
                   right_eigvecs, n, work, lwork, info)
        lwork = int(work(1))
        deallocate(work)
        allocate(work(lwork))

        ! perform eigendecomposition
        call dgeev("V", "V", n, temp, n, eigvals_re, eigvals_im, left_eigvecs, n, &
                   right_eigvecs, n, work, lwork, info)

        ! deallocate arrays
        deallocate(temp, work)

        ! check for successful execution
        if (info /= 0) then
            write(msg, '(A, I0)') "Eigendecomposition failed: Error in DGEEV, "// &
                "info = ", info
            call settings%log(msg, verbosity_error, .true.)
            error = 1
            deallocate(eigvals_re, eigvals_im)
            return
        end if

        ! set eigenvalues
        eigvals%re = eigvals_re
        eigvals%im = eigvals_im
        deallocate(eigvals_re, eigvals_im)

        ! scale left eigenvectors such that biorthonormality is fulfilled; dgeev stores
        ! the eigenvectors of a complex conjugate pair as the real and imaginary parts
        ! of one complex eigenvector in two consecutive columns, which are scaled
        ! together; the scaling stops at an eigenvalue whose left and right
        ! eigenvectors are orthogonal
        i = 1
        do while (i <= n)
            if (abs(eigvals(i)%im) > 0.0_rp) then
                ! scale u = l_i + i l_(i+1) such that its overlap u^H v with
                ! v = r_i + i r_(i+1) is one
                overlap_re = &
                    ddot(n, left_eigvecs(:, i), 1_ip, right_eigvecs(:, i), 1_ip) + &
                    ddot(n, left_eigvecs(:, i + 1), 1_ip, right_eigvecs(:, i + 1), 1_ip)
                overlap_im = &
                    ddot(n, left_eigvecs(:, i), 1_ip, right_eigvecs(:, i + 1), 1_ip) - &
                    ddot(n, left_eigvecs(:, i + 1), 1_ip, right_eigvecs(:, i), 1_ip)
                overlap_sq = overlap_re**2 + overlap_im**2
                if (sqrt(overlap_sq) <= numerical_zero) exit
                left_re = left_eigvecs(:, i)
                left_eigvecs(:, i) = (overlap_re * left_re - &
                                      overlap_im * left_eigvecs(:, i + 1)) / overlap_sq
                left_eigvecs(:, i + 1) = (overlap_im * left_re + overlap_re * &
                                          left_eigvecs(:, i + 1)) / overlap_sq
                i = i + 2
            else
                overlap_re = ddot(n, left_eigvecs(:, i), 1_ip, right_eigvecs(:, i), &
                                  1_ip)
                if (abs(overlap_re) <= numerical_zero) exit
                left_eigvecs(:, i) = left_eigvecs(:, i) / overlap_re
                i = i + 1
            end if
        end do

        ! left and right eigenvectors are orthogonal for an eigenvalue of a defective
        ! matrix, whose eigenvectors do not span the space, so that they cannot be
        ! scaled to biorthonormality
        if (i <= n) then
            call settings%log("Eigendecomposition failed: Left and right "// &
                              "eigenvectors are orthogonal for defective matrix.", &
                              verbosity_error, .true.)
            error = 1
        end if

    end subroutine general_mat_diag

    subroutine mat_min_eig(matrix, symm_matrix, lowest_eigval, lowest_eigvec, &
                           settings, error)
        !
        ! this subroutine returns the lowest eigenvalue and corresponding eigenvector
        ! of a matrix
        !
        real(rp), intent(in) :: matrix(:, :)
        logical, intent(in) :: symm_matrix
        class(settings_type), intent(in) :: settings
        real(rp), intent(out) :: lowest_eigval, lowest_eigvec(:)
        integer(ip), intent(out) :: error

        real(rp) :: eigval(1), eigvec(size(lowest_eigvec), 1)

        call mat_lowest_eigpairs(matrix, symm_matrix, eigval, eigvec, settings, error)
        if (error /= 0) return
        lowest_eigval = eigval(1)
        lowest_eigvec = eigvec(:, 1)

    end subroutine mat_min_eig

    subroutine mat_lowest_eigpairs(matrix, symm_matrix, eigvals, eigvecs, settings, &
                                   error)
        !
        ! this subroutine returns the size(eigvals) lowest eigenpairs of a matrix, in
        ! ascending order
        !
        real(rp), intent(in) :: matrix(:, :)
        logical, intent(in) :: symm_matrix
        class(settings_type), intent(in) :: settings
        real(rp), intent(out) :: eigvals(:), eigvecs(:, :)
        integer(ip), intent(out) :: error

        integer(ip) :: n, n_target, i, idx
        real(rp), allocatable :: all_eigvals(:), right_eigvecs(:, :), left_eigvecs(:, :)
        complex(rp), allocatable :: gen_eigvals(:)
        logical, allocatable :: taken(:)
        real(rp), external :: dnrm2

        ! initialize error flag
        error = 0

        ! size of matrix and number of eigenpairs requested
        n = size(matrix, 1)
        n_target = size(eigvals)

        if (symm_matrix) then
            allocate(all_eigvals(n), right_eigvecs(n, n))
            call symm_mat_diag(matrix, all_eigvals, right_eigvecs, settings, error)
            if (error == 0) then
                eigvals(:n_target) = all_eigvals(:n_target)
                eigvecs(:, :n_target) = right_eigvecs(:, :n_target)
            end if
            deallocate(all_eigvals, right_eigvecs)
        else
            allocate(gen_eigvals(n), right_eigvecs(n, n), left_eigvecs(n, n), taken(n))
            call general_mat_diag(matrix, gen_eigvals, right_eigvecs, left_eigvecs, &
                                  settings, error)
            if (error == 0) then
                ! get eigenvalues with lowest real part, the column of a complex
                ! conjugate pair holds the real or imaginary part of its complex
                ! eigenvector, a real vector of the pair's invariant subspace which is
                ! no eigenvector and is therefore normalized and reported
                taken = .false.
                do i = 1, n_target
                    idx = minloc(gen_eigvals%re, dim=1, mask=.not. taken)
                    taken(idx) = .true.
                    eigvals(i) = gen_eigvals(idx)%re
                    eigvecs(:, i) = right_eigvecs(:, idx)
                    if (abs(gen_eigvals(idx)%im) > 0.0_rp) then
                        eigvecs(:, i) = eigvecs(:, i) / dnrm2(n, eigvecs(:, i), 1_ip)
                        call settings%log(complex_eigval_warning_msg, verbosity_warning)
                    end if
                end do
            end if
            deallocate(gen_eigvals, right_eigvecs, left_eigvecs, taken)
        end if

    end subroutine mat_lowest_eigpairs

    subroutine init_rng(seed)
        !
        ! this subroutine initializes the random number generator
        !
        integer(ip), intent(in) :: seed

        integer, allocatable :: seed_arr(:)
        integer :: n

        call random_seed(size=n)
        allocate(seed_arr(n))
        seed_arr = int(seed, kind=kind(seed_arr))
        call random_seed(put=seed_arr)
        deallocate(seed_arr)

    end subroutine init_rng

    subroutine lowest_h_diag_unit_vectors(h_diag, require_negative, unit_vectors)
        !
        ! this subroutine fills the provided columns with unit vectors along the
        ! successively lowest Hessian diagonal elements, which are the cheapest
        ! available guess at the directions of lowest curvature; when only negative
        ! curvature is wanted, the columns for which no negative element is left are
        ! returned vanishing
        !
        real(rp), intent(in) :: h_diag(:)
        logical, intent(in) :: require_negative
        real(rp), intent(out) :: unit_vectors(:, :)

        integer(ip) :: i, min_idx
        real(rp), allocatable :: h_diag_copy(:)

        ! start from vanishing columns so that any column left unfilled below stays so
        unit_vectors = 0.0_rp

        ! place one unit vector per column at the lowest remaining diagonal element
        h_diag_copy = h_diag
        do i = 1, size(unit_vectors, 2, kind=ip)
            min_idx = minloc(h_diag_copy, dim=1)
            if (require_negative .and. h_diag_copy(min_idx) >= 0.0_rp) exit
            unit_vectors(min_idx, i) = 1.0_rp
            h_diag_copy(min_idx) = huge(1.0_rp)
        end do
        deallocate(h_diag_copy)

    end subroutine lowest_h_diag_unit_vectors

    function generate_trial_vectors(grad, grad_norm, h_diag, settings, error) &
        result(red_space_basis)
        !
        ! this function generates trial vectors
        !
        real(rp), intent(in) :: grad(:), grad_norm, h_diag(:)
        type(solver_settings_type), intent(in) :: settings
        integer(ip), intent(out) :: error

        real(rp), allocatable :: red_space_basis(:, :)

        real(rp), allocatable :: leading_block(:, :)
        integer(ip) :: n_param, n_extra, n_leading, n_total, i

        ! initialize error flag
        error = 0

        ! number of parameters
        n_param = size(grad)

        ! number of extra trial vectors that fit alongside the gradient direction
        n_extra = max(min(settings%n_extra_trial_vectors, n_param - 1), 0_ip)

        ! the gradient direction always leads the trial space
        allocate(leading_block(n_param, 1 + n_extra))
        leading_block(:, 1) = grad / grad_norm

        if (n_extra > 0) then
            ! check if extra trial vectors are supplied
            if (associated(settings%get_extra_trial_vectors)) then
                ! extra trial vectors supplied by the caller
                call settings%get_extra_trial_vectors(leading_block(:, 2:), error, &
                                                      settings%context)
                call add_error_origin(error, error_get_extra_trial_vectors, settings)
            else
                ! unit vectors along the lowest Hessian diagonal elements
                call lowest_h_diag_unit_vectors(h_diag, .true., leading_block(:, 2:))
            end if

            ! project the extra trial vectors, which the caller does not have to
            ! project, while the gradient direction is already projected
            if (error == 0 .and. associated(settings%project)) then
                do i = 2, size(leading_block, 2)
                    call settings%project(leading_block(:, i), error, settings%context)
                    call add_error_origin(error, error_project, settings)
                    if (error /= 0) exit
                end do
            end if
        end if

        ! orthonormalize the leading block, dropping vanishing vectors and vectors that
        ! are linearly dependent on the gradient direction
        if (error == 0) call orthogonalize_trial_vectors(leading_block, settings, error)

        ! the caller assigns the function result before it inspects the error flag, so
        ! the result has to be allocated on every path out of here
        if (error /= 0) then
            allocate(red_space_basis(n_param, 0))
            deallocate(leading_block)
            return
        end if

        ! the trial space holds the leading vectors that survived plus the requested
        ! random ones, as many as fit alongside them
        n_leading = size(leading_block, 2)
        n_total = min(n_leading + settings%n_random_trial_vectors, n_param)

        ! fill the trailing columns with random vectors
        allocate(red_space_basis(n_param, n_total))
        red_space_basis(:, :n_leading) = leading_block
        deallocate(leading_block)
        call generate_random_trial_vectors(red_space_basis, n_total - n_leading, &
                                           settings, error)

    end function generate_trial_vectors

    subroutine generate_random_trial_vectors(red_space_basis, n_fill, settings, error)
        !
        ! this subroutine generates random trial vectors, filling the trailing columns
        ! of the provided space, projecting them and orthogonalizing them against the
        ! leading ones
        !
        real(rp), intent(inout) :: red_space_basis(:, :)
        integer(ip), intent(in) :: n_fill
        class(settings_type), intent(in) :: settings
        integer(ip), intent(out) :: error

        integer(ip) :: n_param, n_trial, i, n_attempts
        integer(ip), parameter :: max_rnd_trial_attempts = 100
        real(rp), parameter :: rnd_vector_min_norm = 1e-3_rp
        real(rp), external :: dnrm2

        ! initialize error flag
        error = 0

        ! number of parameters
        n_param = size(red_space_basis, 1)

        ! number of trial vectors
        n_trial = size(red_space_basis, 2)

        do i = n_trial - n_fill + 1, n_trial
            error = error_gram_schmidt_lin_dep
            n_attempts = 0
            do while (error == error_gram_schmidt_lin_dep)
                n_attempts = n_attempts + 1
                if (n_attempts > max_rnd_trial_attempts) then
                    call settings%log("Maximum number of attempts to generate a "// &
                                      "random trial vector linearly independent of "// &
                                      "the existing trial space reached.", &
                                      verbosity_error, .true.)
                    error = 1
                    return
                end if
                call random_number(red_space_basis(:, i))
                red_space_basis(:, i) = 2 * red_space_basis(:, i) - 1
                do while ( &
                    dnrm2(n_param, red_space_basis(:, i), 1_ip) < rnd_vector_min_norm)
                    call random_number(red_space_basis(:, i))
                    red_space_basis(:, i) = 2 * red_space_basis(:, i) - 1
                end do
                if (associated(settings%project)) then
                    call settings%project(red_space_basis(:, i), error, &
                                          settings%context)
                    call add_error_origin(error, error_project, settings)
                    if (error /= 0) return
                end if
                call gram_schmidt(red_space_basis(:, i), red_space_basis(:, :i - 1), &
                                  settings, error)
            end do
            if (error /= 0) return
        end do

    end subroutine generate_random_trial_vectors

    subroutine block_davidson(hess_x_funptr, h_diag, red_space_basis, symm, n_block, &
                              max_iter, eigvals, eigvecs, converged, settings, state, &
                              error, res_tol)
        !
        ! this subroutine determines the lowest eigenpairs of a linear transformation
        ! with a block Davidson procedure
        !
        procedure(hess_x_type), intent(in), pointer :: hess_x_funptr
        real(rp), intent(in) :: h_diag(:)
        real(rp), allocatable, intent(inout) :: red_space_basis(:, :)
        logical, intent(in) :: symm
        integer(ip), intent(in) :: n_block, max_iter
        real(rp), intent(out) :: eigvals(:), eigvecs(:, :)
        logical, intent(out) :: converged
        type(stability_settings_type), intent(in) :: settings
        type(stability_state_type), intent(inout), target :: state
        integer(ip), intent(out) :: error
        real(rp), intent(in), optional :: res_tol

        integer(ip) :: n_param, n_trial, n_target, n_added, n_unconverged, iter, i, &
                       jacobi_davidson_iter
        real(rp), allocatable :: h_basis(:, :), red_space_hess(:, :), red_eigvals(:), &
                                 red_eigvecs(:, :), solution(:, :), h_solution(:, :), &
                                 residual(:, :)
        logical :: pair_converged
        real(rp), external :: dnrm2
        external :: dgemm

        ! initialize error flag
        error = 0

        ! number of parameters and trial vectors
        n_param = size(h_diag)
        n_trial = size(red_space_basis, 2)

        ! calculate linear transformations of basis vectors
        allocate(h_basis(n_param, n_trial))
        do i = 1, n_trial
            call hess_x_funptr(red_space_basis(:, i), h_basis(:, i), error, &
                               settings%context)
            state%n_hess_x = state%n_hess_x + 1
            call add_error_origin(error, error_hess_x, settings)
            if (error /= 0) return
        end do

        ! loop over iterations
        converged = .false.
        do iter = 1, max_iter
            ! count the iteration
            state%iteration = iter

            ! construct reduced space Hessian
            allocate(red_space_hess(n_trial, n_trial))
            call dgemm("T", "N", n_trial, n_trial, n_param, 1.0_rp, red_space_basis, &
                       n_param, h_basis, n_param, 0.0_rp, red_space_hess, n_trial)

            ! solve reduced space problem for the lowest (real) eigenpairs
            n_target = min(n_block, n_trial)
            allocate(red_eigvals(n_target), red_eigvecs(n_trial, n_target))
            call mat_lowest_eigpairs(red_space_hess, symm, red_eigvals, red_eigvecs, &
                                     settings, error)
            call add_error_origin(error, error_stability_check, settings)
            if (error /= 0) return
            deallocate(red_space_hess)

            ! get full space solutions, their linear transformations and residuals
            allocate(solution(n_param, n_target), h_solution(n_param, n_target), &
                     residual(n_param, n_target))
            call dgemm("N", "N", n_param, n_target, n_trial, 1.0_rp, red_space_basis, &
                       n_param, red_eigvecs, n_trial, 0.0_rp, solution, n_param)
            call dgemm("N", "N", n_param, n_target, n_trial, 1.0_rp, h_basis, n_param, &
                       red_eigvecs, n_trial, 0.0_rp, h_solution, n_param)
            do i = 1, n_target
                residual(:, i) = h_solution(:, i) - red_eigvals(i) * solution(:, i)
            end do

            ! the lowest Ritz value estimates the lowest eigenvalue
            state%min_eigval = red_eigvals(1)
            state%residual_rms = dnrm2(n_param, residual(:, 1), 1_ip) / &
                                 sqrt(real(n_param, kind=rp))

            ! store the current Ritz pairs
            eigvals(:n_target) = red_eigvals(:n_target)
            eigvecs(:, :n_target) = solution

            ! check convergence once through the convergence check function, which
            ! reads the lowest eigenpair from the state, unless a residual threshold is
            ! given, and of every pair otherwise
            if (.not. present(res_tol) .and. associated(settings%conv_check)) then
                converged = settings%conv_check(error, settings%context)
                call add_error_origin(error, error_conv_check, settings)
                if (error /= 0) return
            else
                n_unconverged = 0
                do i = 1, n_target
                    if (present(res_tol)) then
                        pair_converged = dnrm2(n_param, residual(:, i), 1_ip) < res_tol
                    else
                        ! a Ritz value of a symmetric Hessian, the only one for which
                        ! the sanity check allows it, is an upper bound on the true
                        ! lowest eigenvalue, so once it has dropped below the
                        ! instability threshold the sign is already certain and
                        ! convergence can be declared regardless of the residual, when
                        ! requested
                        pair_converged = &
                            (settings%stop_on_instability .and. &
                             red_eigvals(i) < settings%stability_thresh) .or. &
                            dnrm2(n_param, residual(:, i), 1_ip) / &
                            sqrt(real(n_param, kind=rp)) < settings%conv_tol
                    end if
                    if (.not. pair_converged) n_unconverged = n_unconverged + 1
                end do
                converged = n_unconverged == 0
            end if
            if (converged) then
                deallocate(red_eigvals, red_eigvecs, solution, h_solution, residual)
                exit
            end if

            ! expand the space by one vector per unconverged Ritz pair
            n_added = 0
            do i = 1, n_target
                ! the reduced space cannot exceed the full space
                if (n_trial >= n_param) exit

                if (present(res_tol)) then
                    if (dnrm2(n_param, residual(:, i), 1_ip) < res_tol) cycle
                end if

                ! solve the Jacobi-Davidson correction equations for a single tracked
                ! eigenpair once the Jacobi-Davidson method is switched on and
                ! precondition the residual otherwise
                if (.not. present(res_tol) .and. n_target == 1 .and. &
                    settings%diag_solver == "jacobi-davidson" .and. &
                    iter > settings%jacobi_davidson_start) then
                    jacobi_davidson_iter = iter - settings%jacobi_davidson_start
                else
                    jacobi_davidson_iter = 0
                end if

                ! add the new trial vector and its linear transformation to the space,
                ! skipping a vector which is linearly dependent on it, with the
                ! preconditioner level-shifted to the stability threshold: without a
                ! shift, a preconditioner which inverts an approximate Hessian
                ! amplifies the directions in which that Hessian is nearly singular,
                ! such as a rotation that leaves the energy unchanged at a
                ! symmetry-broken point, so that the trial space can converge to the
                ! vanishing eigenvalue of such a direction before reaching a negative
                ! one; the shift bounds this amplification by the curvature the check
                ! has to resolve anyway
                call add_trial_vector(residual(:, i), settings%stability_thresh, &
                                      h_diag, jacobi_davidson_iter, solution(:, i), &
                                      red_eigvals(i), hess_x_funptr, red_space_basis, &
                                      h_basis, settings, state%n_hess_x, error)
                if (error == error_gram_schmidt_lin_dep) then
                    error = 0
                    cycle
                end if
                call add_error_origin(error, error_stability_check, settings)
                if (error /= 0) return
                n_trial = n_trial + 1
                n_added = n_added + 1
            end do
            deallocate(red_eigvals, red_eigvecs, solution, h_solution, residual)

            ! the space can not be expanded any further, which leaves the Ritz pairs
            ! unconverged unless it spans the full space, where they are exact
            if (n_added == 0) then
                converged = n_trial >= n_param
                exit
            end if
        end do

        deallocate(h_basis)

    end subroutine block_davidson

    subroutine get_stability_trial_space(h_diag, red_space_basis, settings, error)
        !
        ! this subroutine generates the initial trial space for the stability check
        !
        real(rp), intent(in) :: h_diag(:)
        real(rp), allocatable, intent(out) :: red_space_basis(:, :)
        type(stability_settings_type), intent(in) :: settings
        integer(ip), intent(out) :: error

        integer(ip) :: n_param, n_total, n_extra, i
        real(rp), allocatable :: leading_block(:, :)

        ! initialize error flag
        error = 0

        ! number of parameters
        n_param = size(h_diag)

        ! number of extra trial vectors requested
        n_extra = max(min(settings%n_extra_trial_vectors, n_param), 0_ip)

        allocate(leading_block(n_param, n_extra))
        if (n_extra > 0) then
            ! check if extra trial vectors are supplied
            if (associated(settings%get_extra_trial_vectors)) then
                ! extra trial vectors supplied by the caller take precedence
                call settings%get_extra_trial_vectors(leading_block, error, &
                                                      settings%context)
                call add_error_origin(error, error_get_extra_trial_vectors, settings)
                if (error /= 0) return
            else
                ! unit vectors along the lowest Hessian diagonal elements
                call lowest_h_diag_unit_vectors(h_diag, .false., leading_block)
            end if

            ! project the extra trial vectors, which the caller does not have to project
            if (associated(settings%project)) then
                do i = 1, n_extra
                    call settings%project(leading_block(:, i), error, settings%context)
                    call add_error_origin(error, error_project, settings)
                    if (error /= 0) return
                end do
            end if

            ! orthonormalize the whole leading block, dropping vanishing and linearly
            ! dependent vectors
            call orthogonalize_trial_vectors(leading_block, settings, error)
            call add_error_origin(error, error_stability_check, settings)
            if (error /= 0) return
            n_extra = size(leading_block, 2)
        end if

        ! the trial space holds the leading vectors that survived plus the requested
        ! random ones
        n_total = min(n_extra + settings%n_random_trial_vectors, n_param)

        ! the Davidson iterations need something to start from, which the checks above
        ! no longer guarantee once every leading vector has been dropped
        if (n_total == 0) then
            call settings%log("All extra trial vectors were dropped and no random "// &
                              "trial vectors are requested, leaving no trial space "// &
                              "for the stability check.", verbosity_error, .true.)
            error = 1
            return
        end if

        ! fill the trailing columns with random vectors and remove any linear
        ! dependencies within the resulting space
        allocate(red_space_basis(n_param, n_total))
        red_space_basis(:, :n_extra) = leading_block
        deallocate(leading_block)
        call generate_random_trial_vectors(red_space_basis, n_total - n_extra, &
                                           settings, error)
        call add_error_origin(error, error_stability_check, settings)

    end subroutine get_stability_trial_space

    subroutine orthogonalize_trial_vectors(red_space_basis, settings, error)
        !
        ! this subroutine orthogonalizes the trial vectors using Gram-Schmidt and
        ! removes linearly dependent and vanishing vectors
        !
        real(rp), intent(inout), allocatable :: red_space_basis(:, :)
        class(settings_type), intent(in) :: settings
        integer(ip), intent(out) :: error

        integer(ip) :: n_param, i, num_valid
        real(rp), allocatable :: tmp_red_space_basis(:, :)
        real(rp), external :: dnrm2

        ! initialize error flag
        error = 0

        ! number of parameters
        n_param = size(red_space_basis, 1)

        ! orthonormalize vectors while removing linearly dependent and vanishing vectors
        num_valid = 0
        do i = 1, size(red_space_basis, 2)
            ! a vanishing vector carries no direction and is how a caller signals that
            ! it has nothing to contribute for this slot, so it is dropped before
            ! Gram-Schmidt, whose own zero-vector guard is an unconditional error
            if (dnrm2(n_param, red_space_basis(:, i), 1_ip) < numerical_zero) cycle
            call gram_schmidt(red_space_basis(:, i), red_space_basis(:, :num_valid), &
                              settings, error)
            if (error == error_gram_schmidt_lin_dep) then
                call settings%log(lin_dep_trial_vector_warning_msg, verbosity_warning)
                error = 0
                cycle
            end if
            if (error /= 0) return
            num_valid = num_valid + 1
            if (num_valid /= i) then
                red_space_basis(:, num_valid) = red_space_basis(:, i)
            end if
        end do
        if (num_valid < size(red_space_basis, 2)) then
            tmp_red_space_basis = red_space_basis(:, :num_valid)
            call move_alloc(tmp_red_space_basis, red_space_basis)
        end if

    end subroutine orthogonalize_trial_vectors

    subroutine gram_schmidt(vector, space, settings, error, lin_trans_vector, &
                            lin_trans_space)
        !
        ! this subroutine orthonormalizes a vector with respect to a vector space,
        ! which keeps a vector within the range of the projection if the vector and the
        ! space lie within it, so that neither has to be projected again; this
        ! subroutine can additionally also return a linear transformation of the
        ! orthogonalized vector if the linear transformations of the vector and the
        ! vector space are provided
        !
        real(rp), intent(inout) :: vector(:)
        real(rp), intent(in) :: space(:, :)
        class(settings_type), intent(in) :: settings
        integer(ip), intent(out) :: error
        real(rp), intent(inout), optional :: lin_trans_vector(:)
        real(rp), intent(in), optional :: lin_trans_space(:, :)

        real(rp), allocatable :: orth(:)
        real(rp) :: norm
        integer(ip) :: n_param, n_vectors, iter, i
        real(rp), parameter :: zero_thres = 1e-16_rp, orth_thres = 1e-14_rp
        real(rp), external :: ddot, dnrm2
        external :: dgemv

        ! initialize error flag
        error = 0

        ! vector size
        n_param = size(vector)

        ! number of vectors
        n_vectors = size(space, 2)

        norm = dnrm2(n_param, vector, 1_ip)
        if (norm < zero_thres) then
            call settings%log("Vector passed to Gram-Schmidt procedure is "// &
                              "numerically zero.", verbosity_error, .true.)
            error = 1
            return
        else if (n_vectors > n_param - 1) then
            call settings%log("Number of vectors in Gram-Schmidt procedure larger "// &
                              "than dimension of vector space.", verbosity_error, &
                              .true.)
            error = 1
            return
        end if

        ! normalize the incoming vector so it always enters a pass below with unit norm
        ! (as does every vector produced by the renormalization at the end of a pass);
        ! this keeps the linear dependence check below meaningful regardless of the
        ! incoming vector's original scale
        vector = vector / norm
        if (present(lin_trans_vector)) lin_trans_vector = lin_trans_vector / norm

        ! allocate array for orthogonalities
        allocate(orth(size(space, 2)))

        iter = 0
        do while (.true.)
            do i = 1, n_vectors
                if (present(lin_trans_vector) .and. present(lin_trans_space)) then
                    lin_trans_vector = lin_trans_vector - &
                                       ddot(n_param, vector, 1_ip, space(:, i), &
                                            1_ip) * lin_trans_space(:, i)
                end if
                vector = orthogonal_projection(vector, space(:, i))
            end do
            norm = dnrm2(n_param, vector, 1_ip)
            ! a linearly dependent vector is in all cases handled as a regular outcome
            if (norm < numerical_zero) then
                error = error_gram_schmidt_lin_dep
                return
            end if
            vector = vector / norm
            if (present(lin_trans_vector) .and. present(lin_trans_space)) &
                lin_trans_vector = lin_trans_vector / norm

            call dgemv("T", n_param, n_vectors, 1.0_rp, space, n_param, vector, 1_ip, &
                       0.0_rp, orth, 1_ip)
            if (maxval(abs(orth)) < orth_thres) exit
            iter = iter + 1
            if (iter > 100) then
                call settings%log("Maximum number of Gram-Schmidt iterations "// &
                                  "reached.", verbosity_error, .true.)
                error = 1
                return
            end if
        end do

        ! deallocate arrays
        deallocate(orth)

    end subroutine gram_schmidt

    subroutine level_shifted_diag_precond(vector, mu, h_diag, precond_vector, &
                                          settings, error)
        !
        ! this subroutine defines the default level-shifted diagonal preconditioner
        !
        real(rp), intent(in) :: vector(:), mu, h_diag(:)
        real(rp), intent(out) :: precond_vector(:)
        class(settings_type), intent(in) :: settings
        integer(ip), intent(out) :: error

        ! initialize error flag
        error = 0

        ! check for user-defined preconditioner
        if (associated(settings%precond)) then
            call settings%precond(vector, mu, precond_vector, error, settings%context)
            call add_error_origin(error, error_precond, settings)
            if (error /= 0) return
        ! construct level-shifted preconditioner
        else
            precond_vector = h_diag - mu
            where (abs(precond_vector) < precond_floor)
                precond_vector = precond_floor
            end where
            precond_vector = vector / precond_vector

            ! ensure basis vector stays in subspace
            if (associated(settings%project)) then
                call settings%project(precond_vector, error, settings%context)
                call add_error_origin(error, error_project, settings)
                if (error /= 0) return
            end if
        end if

    end subroutine level_shifted_diag_precond

    subroutine abs_diag_precond(vector, h_diag, precond_vector, settings, error)
        !
        ! this subroutine defines the default absolute diagonal preconditioner
        !
        real(rp), intent(in) :: vector(:), h_diag(:)
        real(rp), intent(out) :: precond_vector(:)
        class(settings_type), intent(in) :: settings
        integer(ip), intent(out) :: error

        ! initialize error flag
        error = 0

        ! check for user-defined preconditioner
        if (associated(settings%precond)) then
            call settings%precond(vector, 0.0_rp, precond_vector, error, &
                                  settings%context)
            call add_error_origin(error, error_precond, settings)
            if (error /= 0) return
        ! construct positive-definite preconditioner
        else
            precond_vector = max(abs(h_diag), precond_floor)
            precond_vector = vector / precond_vector

            ! ensure basis vector stays in subspace
            if (associated(settings%project)) then
                call settings%project(precond_vector, error, settings%context)
                call add_error_origin(error, error_project, settings)
                if (error /= 0) return
            end if
        end if

    end subroutine abs_diag_precond

    function orthogonal_projection(vector, direction) result(complement)
        !
        ! this function removes a certain direction from a vector
        !
        real(rp), intent(in) :: vector(:), direction(:)

        real(rp) :: complement(size(vector))

        real(rp), external :: ddot

        complement = vector - direction * ddot(int(size(vector), kind=ip), vector, &
                                               1_ip, direction, 1_ip)

    end function orthogonal_projection

    subroutine jacobi_davidson_correction(hess_x_funptr, vector, solution, eigval, &
                                          corr_vector, hess_vector, settings, &
                                          n_hess_x, error)
        !
        ! this subroutine performs the Jacobi-Davidson correction but also returns the
        ! Hessian linear transformation since this can be reused
        !
        procedure(hess_x_type), intent(in), pointer :: hess_x_funptr
        real(rp), intent(in) :: vector(:), solution(:), eigval
        real(rp), intent(out) :: corr_vector(:), hess_vector(:)
        class(settings_type), intent(in) :: settings
        integer(ip), intent(inout), target :: n_hess_x
        integer(ip), intent(out) :: error

        ! initialize error flag
        error = 0

        ! project solution out of vector
        corr_vector = orthogonal_projection(vector, solution)

        ! get Hessian linear transformation of projected vector
        call hess_x_funptr(corr_vector, hess_vector, error, settings%context)
        n_hess_x = n_hess_x + 1
        call add_error_origin(error, error_hess_x, settings)
        if (error /= 0) return

        ! finish construction of correction
        corr_vector = &
            orthogonal_projection(hess_vector - eigval * corr_vector, solution)

    end subroutine jacobi_davidson_correction

    subroutine minres(rhs, hess_x_funptr, solution, eigval, r_tol, vec, hvec, &
                      settings, n_hess_x, error, guess, max_iter)
        !
        ! this subroutine uses the minimum residual method to iteratively solve the
        ! linear system for the Jacobi-Davidson correction equation, modified from
        ! SciPy implementation
        !
        real(rp), intent(in) :: rhs(:), r_tol, solution(:), eigval
        procedure(hess_x_type), intent(in), pointer :: hess_x_funptr
        real(rp), intent(out) :: vec(:), hvec(:)
        class(settings_type), intent(in) :: settings
        integer(ip), intent(inout), target :: n_hess_x
        integer(ip), intent(out) :: error
        real(rp), intent(in), optional :: guess(:)
        integer(ip), intent(in), optional :: max_iter

        integer(ip) :: n, max_iterations, iteration
        real(rp), parameter :: eps = epsilon(1.0_rp)
        real(rp) :: beta_start, beta, phi_bar, rhs1, old_beta, alfa, t_norm2, eps_ln, &
                    old_eps, cs, d_bar, sn, delta, g_bar, root, gamma, phi, g_max, &
                    g_min, tmp, rhs2, a_norm, vec_norm, qr_norm, rel_residual, &
                    rel_ls_residual
        real(rp), allocatable :: matvec(:), r1(:), r2(:), y(:), w(:), hw(:), w1(:), &
                                 hw1(:), w2(:), hw2(:), v(:), hv(:)
        logical :: stop_iteration
        real(rp), external :: dnrm2, ddot

        ! initialize error flag
        error = 0

        ! initialze boolean to stop iterations
        stop_iteration = .false.

        ! size of problem
        n = size(rhs)

        ! allocate solution vector
        allocate(matvec(n), r1(n), r2(n), y(n), w(n), hw(n), w1(n), hw1(n), w2(n), &
                 hw2(n), v(n), hv(n))

        ! initial guess
        if (present(guess)) then
            vec = guess
            call jacobi_davidson_correction(hess_x_funptr, vec, solution, eigval, &
                                            matvec, hvec, settings, n_hess_x, error)
            if (error /= 0) return
        else
            vec = 0.0_rp
            hvec = 0.0_rp
            matvec = 0.0_rp
        end if

        ! maximum number of iterations
        if (present(max_iter)) then
            max_iterations = max_iter
        else
            max_iterations = 5 * n
        end if

        r1 = rhs - matvec
        y = r1

        beta_start = dnrm2(n, r1, 1_ip)

        ! check if starting guess already describes solution
        if (beta_start < numerical_zero) return

        ! solution must be zero vector if rhs vanishes
        if (dnrm2(n, rhs, 1_ip) < numerical_zero) then
            vec = 0.0_rp
            hvec = 0.0_rp
            return
        end if

        ! initialize additional quantities
        beta = beta_start
        r2 = r1
        phi_bar = beta_start
        rhs1 = beta_start
        t_norm2 = 0.0_rp
        eps_ln = 0.0_rp
        cs = -1.0_rp
        d_bar = 0.0_rp
        sn = 0.0_rp
        w = 0.0_rp
        hw = 0.0_rp
        w2 = 0.0_rp
        hw2 = 0.0_rp
        g_max = 0.0_rp
        g_min = huge(1.0_rp)
        rhs2 = 0.0_rp

        iteration = 1
        do while (iteration <= max_iterations)
            ! scale trial vector
            v = y / beta

            ! apply Jacobi-Davidson projector to trial vector
            call jacobi_davidson_correction(hess_x_funptr, v, solution, eigval, y, hv, &
                                            settings, n_hess_x, error)
            if (error /= 0) return

            ! get new trial vector
            if (iteration >= 2) y = y - (beta / old_beta) * r1
            alfa = ddot(n, v, 1_ip, y, 1_ip)
            y = y - (alfa / beta) * r2
            r1 = r2
            r2 = y
            old_beta = beta
            beta = dnrm2(n, r2, 1_ip)
            t_norm2 = t_norm2 + alfa**2 + old_beta**2 + beta**2

            if (iteration == 1 .and. beta / beta_start <= 10 * eps) &
                stop_iteration = .true.

            ! apply previous plane rotation
            ! [delta_k epsln_k+1] = [cs  sn][d_bar_k     0   ]
            ! [g_bar_k d_bar_k+1]   [sn -cs][alfa_k  beta_k+1]
            old_eps = eps_ln
            delta = cs * d_bar + sn * alfa
            g_bar = sn * d_bar - cs * alfa
            eps_ln = sn * beta
            d_bar = -cs * beta
            root = sqrt(g_bar**2 + d_bar**2)

            ! compute next plane rotation by calculating cs and sn
            gamma = max(sqrt(g_bar**2 + beta**2), eps)
            cs = g_bar / gamma
            sn = beta / gamma
            phi = cs * phi_bar
            phi_bar = sn * phi_bar

            ! update vec and hvec
            w1 = w2
            w2 = w
            w = (v - old_eps * w1 - delta * w2) / gamma
            hw1 = hw2
            hw2 = hw
            hw = (hv - old_eps * hw1 - delta * hw2) / gamma
            vec = vec + phi * w
            hvec = hvec + phi * hw

            ! round variables
            g_max = max(g_max, gamma)
            g_min = min(g_min, gamma)
            tmp = rhs1 / gamma
            rhs1 = rhs2 - delta * tmp
            rhs2 = -eps_ln * tmp

            ! estimate various norms and test for convergence
            a_norm = sqrt(t_norm2)
            vec_norm = dnrm2(n, vec, 1_ip)
            qr_norm = phi_bar

            ! relative residuals of the linear system and of the least-squares problem,
            ! which are computed separately since Fortran does not guarantee that a
            ! condition guarding a division by a vanishing norm is evaluated first
            if (vec_norm > 0.0_rp .and. a_norm > 0.0_rp) then
                rel_residual = qr_norm / (a_norm * vec_norm)
            else
                rel_residual = huge(1.0_rp)
            end if
            if (a_norm > numerical_zero) then
                rel_ls_residual = root / a_norm
            else
                rel_ls_residual = huge(1.0_rp)
            end if

            ! check if rhs and initial vector are eigenvectors
            if (stop_iteration) then
                call settings%log("MINRES: beta2 = 0. If M = I, b and x are "// &
                                  "eigenvectors.", verbosity_debug)
                exit
            ! ||r||  / (||A|| ||x||)
            else if (rel_residual <= r_tol) then
                call settings%log("MINRES: A solution to Ax = b was found, given "// &
                                  "provided tolerance.", verbosity_debug)
                exit
            ! ||Ar|| / (||A|| ||r||)
            else if (rel_ls_residual <= r_tol) then
                call settings%log("MINRES: A least-squares solution was found, "// &
                                  "given provided tolerance.", verbosity_debug)
                exit
            ! check if reasonable accuracy is achieved with respect to machine precision
            else if (a_norm * vec_norm * eps >= beta_start) then
                call settings%log("MINRES: Reasonable accuracy achieved, given "// &
                                  "machine precision.", verbosity_debug)
                exit
            ! compare estimate of condition number of matrix to machine precision
            else if (g_max / g_min >= 0.1_rp / eps) then
                call settings%log("MINRES: x has converged to an eigenvector.", &
                                  verbosity_debug)
                exit
            ! these tests ensure convergence is still achieved when r_tol approaches
            ! machine precision
            else if (1.0_rp + rel_residual <= 1.0_rp) then
                call settings%log("MINRES: A solution to Ax = b was found, given "// &
                                  "provided tolerance.", verbosity_debug)
                exit
            else if (1.0_rp + rel_ls_residual <= 1.0_rp) then
                call settings%log("MINRES: A least-squares solution was found, "// &
                                  "given provided tolerance.", verbosity_debug)
                exit
            ! check if maximum number of iterations has been reached, after the tests
            ! for convergence, so that an iteration which also converges does not warn
            else if (iteration == max_iterations) then
                call settings%log(minres_max_iter_warning_msg, verbosity_warning)
                exit
            end if

            ! increment iteration counter
            iteration = iteration + 1

        end do

        deallocate(matvec, r1, r2, y, w, hw, w1, hw1, w2, hw2, v, hv)

    end subroutine minres

    subroutine add_trial_vector(residual, level_shift, h_diag, jacobi_davidson_iter, &
                                solution, eigval, hess_x_funptr, red_space_basis, &
                                h_basis, settings, n_hess_x, error)
        !
        ! this subroutine generates a new trial vector from a residual and adds it and
        ! its Hessian linear transformation to the reduced space basis, either by
        ! preconditioning the residual with the level-shifted Hessian diagonal
        ! (Davidson, jacobi_davidson_iter = 0) or by solving the Jacobi-Davidson
        ! correction equations for the solution and its eigenvalue in Jacobi-Davidson
        ! iteration jacobi_davidson_iter, starting from 1, to a tolerance that is
        ! tightened with every iteration, it returns an error without changing the
        ! reduced space basis if the new vector is linearly dependent on it
        !
        real(rp), intent(in) :: residual(:), level_shift, h_diag(:), solution(:), eigval
        integer(ip), intent(in) :: jacobi_davidson_iter
        procedure(hess_x_type), intent(in), pointer :: hess_x_funptr
        real(rp), allocatable, intent(inout) :: red_space_basis(:, :), h_basis(:, :)
        class(settings_type), intent(in) :: settings
        integer(ip), intent(inout), target :: n_hess_x
        integer(ip), intent(out) :: error

        real(rp), allocatable :: basis_vec(:), h_basis_vec(:)
        integer(ip) :: n_param, n_trial
        real(rp), external :: ddot

        ! initialize error flag
        error = 0

        ! number of parameters and trial vectors
        n_param = size(red_space_basis, 1)
        n_trial = size(red_space_basis, 2)

        ! allocate new basis vector and its Hessian linear transformation
        allocate(basis_vec(n_param), h_basis_vec(n_param))

        if (jacobi_davidson_iter == 0) then
            ! precondition residual
            call level_shifted_diag_precond(residual, level_shift, h_diag, basis_vec, &
                                            settings, error)
            if (error /= 0) return

            ! orthonormalize to current orbital space to get new basis vector, a
            ! linearly dependent vector is returned to the caller
            call gram_schmidt(basis_vec, red_space_basis, settings, error)
            if (error /= 0) return

            ! get linear transformation of new basis vector
            call hess_x_funptr(basis_vec, h_basis_vec, error, settings%context)
            n_hess_x = n_hess_x + 1
            call add_error_origin(error, error_hess_x, settings)
            if (error /= 0) return
        else
            ! solve Jacobi-Davidson correction equations
            call minres(-residual, hess_x_funptr, solution, eigval, &
                        3.0_rp**(-(jacobi_davidson_iter - 1)), basis_vec, h_basis_vec, &
                        settings, n_hess_x, error)
            if (error /= 0) return

            ! orthonormalize to current orbital space to get new basis vector, a
            ! linearly dependent vector is returned to the caller
            call gram_schmidt(basis_vec, red_space_basis, settings, error, &
                              lin_trans_vector=h_basis_vec, lin_trans_space=h_basis)
            if (error /= 0) return

            ! check if resulting linear transformation still respects Hessian symmetry
            ! which can happen due to numerical noise accumulation
            if (abs(ddot(n_param, red_space_basis(:, n_trial), 1_ip, h_basis_vec, &
                         1_ip) - ddot(n_param, basis_vec, 1_ip, h_basis(:, n_trial), &
                                      1_ip)) > hess_symm_thres) then
                call hess_x_funptr(basis_vec, h_basis_vec, error, settings%context)
                n_hess_x = n_hess_x + 1
                call add_error_origin(error, error_hess_x, settings)
                if (error /= 0) return
            end if
        end if

        ! add new trial vector and its linear transformation to reduced space basis
        call add_column(red_space_basis, basis_vec)
        call add_column(h_basis, h_basis_vec)

        ! deallocate new basis vector and its Hessian linear transformation
        deallocate(basis_vec, h_basis_vec)

    end subroutine add_trial_vector

    subroutine print_results(self, iteration, func, grad_rms, level_shift, n_micro, &
                             n_jacobi_davidson, trust_radius, kappa_norm)
        !
        ! this subroutine prints rows of the result table
        !
        class(solver_settings_type), intent(in) :: self
        integer(ip), intent(in) :: iteration
        real(rp), intent(in) :: func, grad_rms
        real(rp), intent(in), optional :: level_shift, trust_radius, kappa_norm
        integer(ip), intent(in), optional :: n_micro, n_jacobi_davidson

        character(len=11) :: iteration_str
        character(len=12) :: n_micro_str
        character(len=28) :: func_str
        character(len=14) :: grad_rms_str, trust_radius_str
        character(len=13) :: level_shift_str
        character(len=10) :: kappa_norm_str

        write(iteration_str, '(4X, I4, 3X)') iteration
        write(func_str, '(4X, 1PE21.14, 3X)') func
        write(grad_rms_str, '(2X, 1PE9.2, 3X)') grad_rms

        if (present(level_shift)) then
            write(level_shift_str, '(2X, 1PE9.2, 2X)') level_shift
        else
            write(level_shift_str, '(6X, "-", 6X)')
        end if

        if (present(n_micro)) then
            if (present(n_jacobi_davidson)) then
                write(n_micro_str, '(1X, I3, 1X, "|", 1X, I3, 2X)') n_jacobi_davidson, &
                    n_micro - n_jacobi_davidson
            else
                write(n_micro_str, '(7X, I3, 2X)') n_micro
            end if
        else
            write(n_micro_str, '(8X, "-", 3X)')
        end if

        if (present(trust_radius)) then
            write(trust_radius_str, '(3X, 1PE9.2, 2X)') trust_radius
        else
            write(trust_radius_str, '(8X, "-", 5X)')
        end if

        if (present(kappa_norm)) then
            write(kappa_norm_str, '(1X, 1PE9.2)') kappa_norm
        else
            write(kappa_norm_str, '(6X, "-", 3X)')
        end if

        call self%log(iteration_str//"|"//func_str//"|"//grad_rms_str//"|"// &
                      level_shift_str//"|"//n_micro_str//"|"//trust_radius_str//"|"// &
                      kappa_norm_str, verbosity_info)

    end subroutine print_results

    subroutine print_message(self, message, level, error)
        !
        ! this subroutine performs logging
        !
        class(settings_type), intent(in) :: self
        character(len=*), intent(in) :: message
        integer(ip), intent(in) :: level
        logical, intent(in), optional :: error

        integer(ip), parameter :: max_length = 109
        integer(ip) :: i, out_unit
        character(len=:), dimension(:), allocatable :: substrings

        if (self%verbose >= level) then
            call split_string_by_space(message, max_length, substrings)
            if (associated(self%logger)) then
                do i = 1, size(substrings)
                    call self%logger(" "//substrings(i), self%context)
                end do
            else
                if (.not. present(error)) then
                    out_unit = stdout
                else if (.not. error) then
                    out_unit = stdout
                else
                    out_unit = stderr
                end if
                do i = 1, size(substrings)
                    write(out_unit, '(A)') " "//substrings(i)
                end do
            end if
            deallocate(substrings)
        end if

    end subroutine print_message

    subroutine split_string_by_space(input, max_length, substrings)
        !
        ! this subroutine splits a string by spaces to produce substrings of a maximum
        ! length
        !
        character(len=*), intent(in) :: input
        integer(ip), intent(in) :: max_length
        character(len=:), intent(out), allocatable :: substrings(:)

        integer(ip) :: i, len_input, start_pos, end_pos, space_pos

        character(len=max_length) :: temp_string
        character(len=max_length), allocatable :: temp_substrings(:)

        len_input = len_trim(input)
        start_pos = 1

        ! start from an empty set of substrings
        allocate(character(len=max_length) :: substrings(0))

        do while (start_pos <= len_input)
            end_pos = min(start_pos + max_length - 1, len_input)
            space_pos = 0

            ! find last space before max_length characters
            do i = end_pos, start_pos, -1
                if (input(i:i) == " ") then
                    space_pos = i
                    exit
                end if
            end do

            ! force splitting if no space is found
            if (end_pos == len_input) then
                space_pos = len_input
            else if (space_pos == 0) then
                space_pos = end_pos
            end if

            ! add substring
            temp_string = " "
            temp_string(1:space_pos - start_pos + 1) = input(start_pos:space_pos)
            allocate(temp_substrings(size(substrings) + 1))
            temp_substrings(1:size(substrings)) = substrings
            temp_substrings(size(substrings) + 1) = temp_string
            substrings = temp_substrings
            deallocate(temp_substrings)
            start_pos = space_pos + 1
        end do

    end subroutine split_string_by_space

    subroutine level_shifted_davidson(func, grad, grad_norm, h_diag, n_param, &
                                      obj_func, hess_x_funptr, settings, local_region, &
                                      state, solution, solution_norm, mu, error)
        !
        ! this subroutine performs level-shifted (Jacobi-)Davidson to solve the trust
        ! region subproblem
        !
        real(rp), intent(in) :: func, grad(:), grad_norm, h_diag(:)
        integer(ip), intent(in) :: n_param
        procedure(obj_func_type), pointer, intent(in) :: obj_func
        procedure(hess_x_type), pointer, intent(in) :: hess_x_funptr
        type(solver_settings_type), intent(in) :: settings
        logical, intent(in) :: local_region
        type(solver_state_type), intent(inout), target :: state
        real(rp), intent(out) :: solution(:), solution_norm, mu
        integer(ip), intent(out) :: error

        real(rp), allocatable :: &
            perturbed_grad(:), red_space_basis(:, :), h_basis(:, :), &
            red_space_hess(:, :), red_space_solution(:), h_solution(:), residual(:), &
            solution_normalized(:), last_solution_normalized(:), column_vec(:), &
            row_vec(:), red_space_hess_right_eigvecs(:, :), &
            red_space_hess_left_eigvecs(:, :), red_space_hess_eigvals_re(:)
        complex(rp), allocatable :: red_space_hess_eigvals(:)
        integer(ip) :: n_trial, i, imicro, initial_imicro, min_idx, jacobi_davidson_iter
        logical :: accept_step, micro_converged, newton
        real(rp) :: perturbed_grad_norm, residual_norm, initial_residual_norm, new_func
        real(rp), parameter :: solution_overlap_thresh = 0.5_rp, &
                               residual_norm_max_red_factor = 0.8_rp
        real(rp), external :: dnrm2, ddot
        external :: dgemm, dgemv

        ! initialize error flag
        error = 0

        ! no micro iterations yet, and the residual of the vanishing solution is the
        ! gradient
        state%micro_iteration = 0
        state%jacobi_davidson_iteration = 0
        state%relative_residual = 1.0_rp

        ! solve the subproblem for a gradient perturbed by a small random vector, which
        ! keeps the solution from remaining exactly within the symmetry subspace of the
        ! current point when the Hessian (in particular an approximate one) carries no
        ! negative curvature out of that subspace
        perturbed_grad = grad
        perturbed_grad_norm = grad_norm
        if (settings%grad_noise) then
            call perturb_vector(perturbed_grad, noise_level(settings, local_region), &
                                settings, error)
            if (error /= 0) return
            perturbed_grad_norm = dnrm2(n_param, perturbed_grad, 1_ip)
        end if

        ! generate trial vectors
        red_space_basis = generate_trial_vectors(perturbed_grad, perturbed_grad_norm, &
                                                 h_diag, settings, error)
        if (error /= 0) return

        ! number of trial vectors
        n_trial = size(red_space_basis, 2)

        ! calculate linear transformations of basis vectors
        allocate(h_basis(n_param, n_trial))
        do i = 1, n_trial
            call hess_x_funptr(red_space_basis(:, i), h_basis(:, i), error, &
                               settings%context)
            state%n_hess_x = state%n_hess_x + 1
            call add_error_origin(error, error_hess_x, settings)
            if (error /= 0) return
        end do

        ! construct reduced space Hessian
        allocate(red_space_hess(n_trial, n_trial))
        call dgemm("T", "N", n_trial, n_trial, n_param, 1.0_rp, red_space_basis, &
                   n_param, h_basis, n_param, 0.0_rp, red_space_hess, n_trial)

        ! allocate space for reduced space solution, Hessian linear transformation of
        ! solution, residual and the (last) normalized solution
        allocate(red_space_solution(n_trial), h_solution(n_param), residual(n_param), &
                 solution_normalized(n_param), last_solution_normalized(n_param))

        ! decrease trust radius until micro iterations converge and step is accepted
        last_solution_normalized = 0.0_rp
        accept_step = .false.
        do while (.not. accept_step)
            micro_converged = .false.

            jacobi_davidson_iter = 0
            do imicro = 1, settings%n_micro
                ! count the micro iteration and, once Jacobi-Davidson has started for
                ! this step, the Jacobi-Davidson iteration, the state over all steps of
                ! the subproblem
                state%micro_iteration = state%micro_iteration + 1
                if (jacobi_davidson_iter > 0) then
                    jacobi_davidson_iter = jacobi_davidson_iter + 1
                    state%jacobi_davidson_iteration = &
                        state%jacobi_davidson_iteration + 1
                end if

                ! perform eigendecomposition
                allocate(red_space_hess_right_eigvecs(n_trial, n_trial), &
                         red_space_hess_left_eigvecs(n_trial, n_trial), &
                         red_space_hess_eigvals(n_trial))
                call mat_diag(red_space_hess, settings%hess_symm, &
                              red_space_hess_eigvals, red_space_hess_right_eigvecs, &
                              red_space_hess_left_eigvecs, settings, error)
                if (error /= 0) return

                if (string_in("ls", settings%subsystem_solver)) then
                    ! perform bisection to find the level shift
                    call bisection_mu(red_space_hess, perturbed_grad_norm, &
                                      red_space_basis, red_space_hess_eigvals, &
                                      red_space_hess_right_eigvecs, &
                                      red_space_hess_left_eigvecs, state%trust_radius, &
                                      solution, red_space_solution, mu, settings, error)
                    if (error /= 0) return
                else
                    ! do a Newton step if the model is positive definite and the step
                    ! is within the trust region
                    newton = .false.
                    red_space_hess_eigvals_re = red_space_hess_eigvals%re
                    min_idx = minloc(red_space_hess_eigvals%re, dim=1)
                    if (red_space_hess_eigvals(min_idx)%re > newton_eigval_thresh) then
                        call newton_step( &
                            perturbed_grad_norm, red_space_basis, &
                            red_space_hess_eigvals_re, red_space_hess_right_eigvecs, &
                            red_space_hess_left_eigvecs, solution, red_space_solution)
                        mu = 0.0_rp
                        if (dnrm2(n_param, solution, 1_ip) < state%trust_radius) &
                            newton = .true.
                    end if

                    ! perform bisection on augmented Hessian to find the level shift
                    if (.not. newton) then
                        call bisection_ah( &
                            red_space_hess, perturbed_grad_norm, red_space_basis, &
                            red_space_hess_eigvals_re, red_space_hess_right_eigvecs, &
                            red_space_hess_left_eigvecs, state%trust_radius, solution, &
                            red_space_solution, mu, settings, error)
                        if (error /= 0) return
                    end if
                    deallocate(red_space_hess_eigvals_re)
                end if

                deallocate(red_space_hess_right_eigvecs, red_space_hess_left_eigvecs, &
                           red_space_hess_eigvals)

                ! calculate Hessian linear transformation of solution
                call dgemv("N", n_param, n_trial, 1.0_rp, h_basis, n_param, &
                           red_space_solution, 1_ip, 0.0_rp, h_solution, 1_ip)

                ! calculate residual
                residual = perturbed_grad + h_solution - mu * solution

                ! calculate residual norm
                residual_norm = dnrm2(n_param, residual, 1_ip)
                state%relative_residual = residual_norm / perturbed_grad_norm

                ! get normalized solution vector
                solution_normalized = solution / dnrm2(n_param, solution, 1_ip)

                ! reset initial residual norm if solution changes
                if (ddot(n_param, last_solution_normalized, 1_ip, solution_normalized, &
                         1_ip)**2 < solution_overlap_thresh) then
                    initial_imicro = imicro
                    initial_residual_norm = residual_norm
                end if

                ! check if micro iterations have converged
                if (residual_norm < micro_conv_tol(mu, noise_level( &
                    settings, local_region), perturbed_grad_norm, settings)) then
                    micro_converged = .true.
                    exit
                ! check if Jacobi-Davidson is used and has not been started
                else if (string_in("jacobi", settings%subsystem_solver) .and. &
                         jacobi_davidson_iter == 0) then
                    ! check residual has not decreased sufficiently or if maximum of
                    ! Davidson iterations has been reached
                    if ((imicro - initial_imicro >= 10 .and. residual_norm > &
                         residual_norm_max_red_factor * initial_residual_norm) .or. &
                        imicro > settings%jacobi_davidson_start) then

                        ! switch to Jacobi-Davidson
                        jacobi_davidson_iter = 1
                        state%jacobi_davidson_iteration = &
                            state%jacobi_davidson_iteration + 1
                    end if
                ! check if residual has not decreased sufficiently
                else if (imicro - initial_imicro >= 10 .and. residual_norm > &
                         residual_norm_max_red_factor * initial_residual_norm) then
                    exit
                end if

                ! save current solution
                last_solution_normalized = solution_normalized

                ! the reduced space is not reset when a trust region step is rejected,
                ! so after enough rejected steps it grows past the dimension of the
                ! full parameter space, so in that case stop here since the reduced
                ! space has reached full rank and cannot be expanded further
                if (n_trial >= n_param) then
                    micro_converged = .true.
                    exit
                end if

                ! add new trial vector from preconditioned residual (Davidson) or from
                ! Jacobi-Davidson correction equations in the current Jacobi-Davidson
                ! iteration
                call add_trial_vector(residual, mu, h_diag, jacobi_davidson_iter, &
                                      solution_normalized, mu, hess_x_funptr, &
                                      red_space_basis, h_basis, settings, &
                                      state%n_hess_x, error)
                if (error == error_gram_schmidt_lin_dep) then
                    ! new vector is linearly dependent, so the reduced space cannot
                    ! be usefully expanded further due to degeneracy
                    micro_converged = .true.
                    exit
                else if (error /= 0) then
                    return
                end if

                ! increment trial vector count
                n_trial = n_trial + 1

                ! extend reduced space Hessian
                allocate(column_vec(n_trial))
                call dgemv("T", n_param, n_trial, 1.0_rp, red_space_basis, n_param, &
                           h_basis(:, n_trial), 1_ip, 0.0_rp, column_vec, 1_ip)
                row_vec = column_vec
                if (.not. settings%hess_symm) then
                    call dgemv("T", n_param, n_trial, 1.0_rp, h_basis, n_param, &
                               red_space_basis(:, n_trial), 1_ip, 0.0_rp, row_vec, 1_ip)
                end if
                call extend_matrix(red_space_hess, column_vec, row_vec)
                deallocate(column_vec, row_vec)

                ! reallocate reduced space solution
                deallocate(red_space_solution)
                allocate(red_space_solution(n_trial))
            end do

            ! evaluate function at predicted point
            new_func = obj_func(solution, error, settings%context)
            state%n_obj_func = state%n_obj_func + 1
            call add_error_origin(error, error_obj_func, settings)
            if (error /= 0) return

            ! decide whether to accept step and modify trust radius; the Davidson
            ! trust region is always Euclidean, so this norm must be computed here
            ! explicitly rather than reused from solution_norm, which is only set
            ! once after this loop exits
            accept_step = accept_trust_region_step( &
                solution, dnrm2(n_param, solution, 1_ip), func, new_func - func, &
                ddot(n_param, solution, 1_ip, perturbed_grad + 0.5_rp * h_solution, &
                     1_ip), micro_converged, settings, state)
            if (state%max_precision_reached) exit
        end do

        ! deallocate quantities from microiterations
        deallocate(red_space_solution, red_space_hess, red_space_basis, h_basis, &
                   h_solution, residual, solution_normalized, last_solution_normalized)

        ! get norm of orbital rotation
        solution_norm = dnrm2(n_param, solution, 1_ip)

    end subroutine level_shifted_davidson

    subroutine truncated_conjugate_gradient(func, grad, h_diag, n_param, obj_func, &
                                            hess_x_funptr, settings, state, solution, &
                                            solution_norm, error)
        !
        ! this subroutine performs truncated conjugate gradient to solve the trust
        ! region subproblem
        !
        real(rp), intent(in) :: func, grad(:), h_diag(:)
        integer(ip), intent(in) :: n_param
        procedure(obj_func_type), pointer, intent(in) :: obj_func
        procedure(hess_x_type), pointer, intent(in) :: hess_x_funptr
        type(solver_settings_type), intent(in) :: settings
        type(solver_state_type), intent(inout), target :: state
        real(rp), intent(out) :: solution(:), solution_norm
        integer(ip), intent(out) :: error

        logical :: accept_step, micro_converged
        real(rp) :: model_func, initial_residual_norm, curvature, step_size, &
                    solution_dot, solution_direction_dot, direction_dot, step_length, &
                    model_func_new, new_func
        real(rp), allocatable :: &
            h_solution(:), residual(:), precond_residual(:), direction(:), &
            hess_direction(:), precond_solution(:), precond_direction(:), &
            solution_new(:), h_solution_new(:), residual_new(:), &
            precond_residual_new(:), solutions(:, :), h_solutions(:, :)
        integer(ip) :: i, imicro
        real(rp), parameter :: residual_norm_conv_red_factor = 1e-3_rp
        real(rp), external :: ddot, dnrm2

        ! initialize error flag
        error = 0

        ! no micro iterations yet, and the residual of the vanishing solution is the
        ! gradient
        state%micro_iteration = 0
        state%jacobi_davidson_iteration = 0
        state%relative_residual = 1.0_rp

        ! allocate arrays
        allocate(h_solution(n_param), solutions(n_param, 1), h_solutions(n_param, 1), &
                 precond_residual(n_param))

        ! initialize solution
        solution = 0.0_rp
        h_solution = 0.0_rp
        solutions(:, 1) = solution
        h_solutions(:, 1) = h_solution
        model_func = 0.0_rp

        ! initialize residual, which should include h_solution if not starting at zero,
        ! and perturb it, if the gradient is perturbed, by a random vector whose
        ! relative size is the reduction factor of the local region, which keeps the
        ! solution from remaining exactly within a symmetry subspace of the current
        ! point
        residual = grad
        if (settings%grad_noise) then
            call perturb_vector(residual, noise_level(settings, .true.), settings, &
                                error)
            if (error /= 0) return
        end if

        ! get initial residual norm
        initial_residual_norm = dnrm2(n_param, residual, 1_ip)

        ! initialize preconditioned residual and direction
        call abs_diag_precond(residual, h_diag, precond_residual, settings, error)
        if (error /= 0) return

        direction = -precond_residual

        ! allocate arrays needed to propose a new step
        allocate(hess_direction(n_param), precond_solution(n_param), &
                 precond_direction(n_param), solution_new(n_param), &
                 h_solution_new(n_param), residual_new(n_param), &
                 precond_residual_new(n_param))

        ! start microiteration loop
        micro_converged = .false.
        do imicro = 1, settings%n_micro - 1
            ! count the micro iteration
            state%micro_iteration = state%micro_iteration + 1

            ! get Hessian linear transformation of direction
            call hess_x_funptr(direction, hess_direction, error, settings%context)
            state%n_hess_x = state%n_hess_x + 1
            call add_error_origin(error, error_hess_x, settings)
            if (error /= 0) return

            ! calculate curvature
            curvature = ddot(n_param, direction, 1_ip, hess_direction, 1_ip)

            ! get step size along new direction
            step_size = ddot(n_param, residual, 1_ip, precond_residual, 1_ip) / &
                        curvature

            ! precondition current solution and direction
            call abs_diag_precond(solution, h_diag, precond_solution, settings, error)
            if (error /= 0) return
            call abs_diag_precond(direction, h_diag, precond_direction, settings, error)
            if (error /= 0) return

            ! calculate dot products
            solution_dot = ddot(n_param, solution, 1_ip, precond_solution, 1_ip)
            solution_direction_dot = &
                ddot(n_param, solution, 1_ip, precond_direction, 1_ip)
            direction_dot = ddot(n_param, direction, 1_ip, precond_direction, 1_ip)

            ! calculate total step length
            step_length = solution_dot + 2 * step_size * solution_direction_dot + &
                          step_size**2 * direction_dot

            if (curvature < 0.0_rp .or. step_length >= state%trust_radius**2) then
                ! solve quadratic equation
                step_size = (-solution_direction_dot + sqrt( &
                    solution_direction_dot**2 + direction_dot * &
                    (state%trust_radius**2 - solution_dot))) / direction_dot

                ! get step to boundary and exit
                solution = solution + step_size * direction
                h_solution = h_solution + step_size * hess_direction
                call add_column(solutions, solution)
                call add_column(h_solutions, h_solution)
                micro_converged = .true.
                exit
            end if

            ! get new step
            solution_new = solution + step_size * direction
            h_solution_new = h_solution + step_size * hess_direction

            ! get new model function value
            model_func_new = &
                ddot(n_param, solution_new, 1_ip, grad + 0.5_rp * h_solution_new, 1_ip)

            ! check if model was improved
            if (model_func_new >= model_func) then
                micro_converged = .true.
                exit
            end if

            ! accept step
            solution = solution_new
            h_solution = h_solution_new
            call add_column(solutions, solution)
            call add_column(h_solutions, h_solution)

            ! get residual for model
            residual_new = residual + step_size * hess_direction
            call abs_diag_precond(residual_new, h_diag, precond_residual_new, &
                                  settings, error)
            if (error /= 0) return

            ! update the residual of the subproblem relative to the initial one
            state%relative_residual = dnrm2(n_param, residual_new, 1_ip) / &
                                      initial_residual_norm

            ! check for linear or superlinear (in this case quadratic) convergence,
            ! never below the relative size of the gradient perturbation, which the
            ! micro iterations are not asked to resolve
            if (dnrm2(n_param, residual_new, 1_ip) <= initial_residual_norm * max( &
                min(residual_norm_conv_red_factor, initial_residual_norm), &
                noise_level(settings, .true.))) then
                micro_converged = .true.
                exit
            end if

            ! get new search direction
            direction = -precond_residual_new + ddot( &
                n_param, residual_new, 1_ip, precond_residual_new, 1_ip) / ddot( &
                    n_param, residual, 1_ip, precond_residual, 1_ip) * direction

            ! save new model
            model_func = model_func_new
            residual = residual_new
            precond_residual = precond_residual_new
        end do

        ! deallocate arrays no longer needed
        deallocate(residual, precond_residual, solution_new, h_solution_new, &
                   residual_new, precond_residual_new)

        ! evaluate function at predicted point
        new_func = obj_func(solution, error, settings%context)
        state%n_obj_func = state%n_obj_func + 1
        call add_error_origin(error, error_obj_func, settings)
        if (error /= 0) return

        ! reduce trust region until step is accepted
        do while (.true.)
            ! get norm of the step in the metric of the absolute diagonal
            ! preconditioner, which defines the trust region
            call abs_diag_precond(solution, h_diag, precond_solution, settings, error)
            if (error /= 0) return
            solution_norm = sqrt(ddot(n_param, solution, 1_ip, precond_solution, 1_ip))

            ! decide whether to accept step and modify trust radius
            accept_step = accept_trust_region_step( &
                solution, solution_norm, func, new_func - func, &
                ddot(n_param, solution, 1_ip, grad + 0.5_rp * h_solution, 1_ip), &
                micro_converged, settings, state)
            if (accept_step .or. state%max_precision_reached) exit

            ! check if step exceeds new trust region boundary
            if (solution_norm > state%trust_radius) then
                ! find step that exceeds trust region boundary
                do i = 1, size(solutions, 2)
                    call abs_diag_precond(solutions(:, i), h_diag, precond_solution, &
                                          settings, error)
                    if (error /= 0) return

                    if (ddot(n_param, solutions(:, i), 1_ip, precond_solution, 1_ip) > &
                        state%trust_radius**2) then
                        ! get previous step
                        solution = solutions(:, i - 1)
                        h_solution = h_solutions(:, i - 1)

                        ! get direction
                        direction = solutions(:, i) - solutions(:, i - 1)
                        hess_direction = h_solutions(:, i) - h_solutions(:, i - 1)

                        ! precondition current solution and direction
                        call abs_diag_precond(solution, h_diag, precond_solution, &
                                              settings, error)
                        if (error /= 0) return
                        call abs_diag_precond(direction, h_diag, precond_direction, &
                                              settings, error)
                        if (error /= 0) return

                        ! calculate dot products
                        solution_dot = &
                            ddot(n_param, solution, 1_ip, precond_solution, 1_ip)
                        solution_direction_dot = &
                            ddot(n_param, solution, 1_ip, precond_direction, 1_ip)
                        direction_dot = &
                            ddot(n_param, direction, 1_ip, precond_direction, 1_ip)

                        ! solve quadratic equation
                        step_size = (-solution_direction_dot + sqrt( &
                            solution_direction_dot**2 + direction_dot * &
                            (state%trust_radius**2 - solution_dot))) / direction_dot

                        ! get step to boundary
                        solution = solution + step_size * direction
                        h_solution = h_solution + step_size * hess_direction

                        ! evaluate function at predicted point
                        new_func = obj_func(solution, error, settings%context)
                        state%n_obj_func = state%n_obj_func + 1
                        call add_error_origin(error, error_obj_func, settings)
                        if (error /= 0) return

                        ! exit to retry whether step is accepted
                        micro_converged = .true.
                        exit
                    end if
                end do
            end if
        end do

        ! a step that is not taken has no norm
        if (state%max_precision_reached) solution_norm = 0.0_rp

        ! deallocate arrays
        deallocate(h_solution, solutions, h_solutions, direction, hess_direction, &
                   precond_solution, precond_direction)

    end subroutine truncated_conjugate_gradient

    logical function accept_trust_region_step(solution, solution_norm, func, &
                                              actual_func_diff, pred_func_diff, &
                                              micro_converged, settings, state)
        !
        ! this function checks whether the trust region step is accepted and modifies
        ! the trust radius in the state accordingly, where it also flags whether the
        ! function changes could not be resolved anymore
        !
        real(rp), intent(in) :: solution(:), solution_norm, func, actual_func_diff, &
                                pred_func_diff
        logical, intent(in) :: micro_converged
        type(solver_settings_type), intent(in) :: settings
        type(solver_state_type), intent(inout), target :: state

        real(rp), parameter :: func_resolution_safety_factor = 10.0_rp
        real(rp) :: ratio, func_resolution
        logical :: finite_changes

        ! default to maximum precision not yet reached
        state%max_precision_reached = .false.

        ! check if function value has not decreased up to floating point precision, in
        ! which case the step is not taken
        func_resolution = func_resolution_safety_factor * abs(func) * epsilon(1.0_rp)
        if (abs(actual_func_diff) <= func_resolution .and. &
            abs(pred_func_diff) <= func_resolution) then
            call settings%log(function_unchanged_warning_msg, verbosity_warning)
            state%max_precision_reached = .true.
            accept_trust_region_step = .false.
            return
        end if

        ! calculate ratio of evaluated function and predicted function
        ratio = actual_func_diff / pred_func_diff

        ! decrease trust radius if micro iterations are unable to converge, if function
        ! value has not decreased, if the function changes are not finite, as for an
        ! objective function which overflows at the step, or if individual orbitals
        ! change too much
        finite_changes = abs(actual_func_diff) <= huge(1.0_rp) .and. &
                         abs(pred_func_diff) <= huge(1.0_rp)
        if (.not. micro_converged .or. ratio < 0.0_rp .or. .not. finite_changes .or. &
            any(abs(solution) > pi / 4)) then
            ! a rejected step whose predicted function change can not be resolved can
            ! not be improved by a smaller trust radius, as when the objective function
            ! and the orbital update disagree by more than the resolution, so the step
            ! is not taken and maximum precision is reached
            if (abs(pred_func_diff) <= func_resolution) then
                call settings%log(function_unchanged_warning_msg, verbosity_warning)
                state%max_precision_reached = .true.
                accept_trust_region_step = .false.
                return
            end if
            state%trust_radius = trust_radius_shrink_factor * state%trust_radius
            accept_trust_region_step = .false.
        ! check if step is too long
        else if (ratio < trust_radius_shrink_ratio) then
            state%trust_radius = trust_radius_shrink_factor * state%trust_radius
            accept_trust_region_step = .true.
        ! check if quadratic approximation is valid
        else if (ratio < trust_radius_expand_ratio) then
            accept_trust_region_step = .true.
        ! check if step is potentially too short
        else
            if (solution_norm >= 0.99_rp * state%trust_radius) &
                state%trust_radius = trust_radius_expand_factor * state%trust_radius
            accept_trust_region_step = .true.
        end if

    end function accept_trust_region_step

    subroutine init_defaults(settings)
        !
        ! this subroutine resolves the defaults that depend on other settings: the
        ! nested stability check settings inherit the solver's preconditioning,
        ! projection and logging functions and its host context unless they provide
        ! their own, its verbosity if higher and a non-symmetric Hessian, stop as soon
        ! as an instability is certain for a symmetric Hessian unless they provide a
        ! convergence check and use a single random trial vector unless their number
        ! is set
        !
        type(solver_settings_type), intent(inout) :: settings

        ! let the nested stability check settings inherit preconditioning, projection
        ! and logging functions and the host context if they provide none of their own
        associate (nested => settings%stability_settings)
            if (.not. associated(nested%precond)) nested%precond => settings%precond
            if (.not. associated(nested%project)) nested%project => settings%project
            if (.not. associated(nested%logger)) nested%logger => settings%logger
            if (.not. associated(nested%context)) nested%context => settings%context

            ! inherit the solver's verbosity if higher
            nested%verbose = max(nested%verbose, settings%verbose)

            ! treat the Hessian as non-symmetric if the solver does, since the internal
            ! stability check diagonalizes the solver's Hessian
            nested%hess_symm = nested%hess_symm .and. settings%hess_symm

            ! stop the check as soon as an instability is certain for a symmetric
            ! Hessian unless a convergence check is provided, and use a single random
            ! trial vector unless their number is set
            if (.not. associated(nested%conv_check) .and. nested%hess_symm) &
                nested%stop_on_instability = .true.
            if (nested%n_random_trial_vectors < 0) nested%n_random_trial_vectors = 1
        end associate

    end subroutine init_defaults

    subroutine solver_sanity_check(settings, n_param, error)
        !
        ! this subroutine performs a sanity check for solver input parameters
        !
        type(solver_settings_type), intent(inout) :: settings
        integer(ip), intent(in) :: n_param
        integer(ip), intent(out) :: error

        character(len=300) :: msg

        ! initialize error flag
        error = 0

        ! check that number of parameters is positive
        if (n_param < 1) then
            call settings%log("Number of parameters should be larger than 0.", &
                              verbosity_error, .true.)
            error = 1
            return
        end if

        ! check that number of microiterations is positive since the subsystem solvers
        ! assume that at least one microiteration is performed
        if (settings%n_micro < 1) then
            call settings%log("Number of microiterations should be larger than 0.", &
                              verbosity_error, .true.)
            error = 1
            return
        end if

        ! check that number of random trial vectors is not negative
        if (settings%n_random_trial_vectors < 0) then
            call settings%log("Number of random trial vectors should not be "// &
                              "negative.", verbosity_error, .true.)
            error = 1
            return
        end if

        ! convert strings to lowercase
        settings%subsystem_solver = string_to_lowercase(settings%subsystem_solver)

        ! check that number of random trial vectors is below number of parameters
        if (string_in("davidson", settings%subsystem_solver) .and. &
            settings%n_random_trial_vectors > n_param / 2) then
            settings%n_random_trial_vectors = n_param / 2
            write(msg, '(A, I0, A)') random_trial_vector_warning_msg//" Setting to ", &
                settings%n_random_trial_vectors, "."
            call settings%log(msg, verbosity_warning)
        end if

        ! check whether the extra trial vector function has anything to fill
        if (associated(settings%get_extra_trial_vectors) .and. &
            settings%n_extra_trial_vectors <= 0) then
            call settings%log("Number of extra trial vectors should be larger than "// &
                              "zero when extra trial vector function is passed.", &
                              verbosity_error, .true.)
            error = 1
            return
        end if

        ! check for character options
        if (.not. any(settings%subsystem_solver == subsystem_solver_options)) then
            call settings%log( &
                "Subsystem solver option unknown. Possible values are "// &
                option_list(subsystem_solver_options)//".", verbosity_error, .true.)
            error = 1
            return
        end if

        ! check whether non-symmetric Hessian is requested with a solver that cannot
        ! handle it
        if (.not. settings%hess_symm .and. &
            (string_in("ah", settings%subsystem_solver) .or. &
             string_in("jacobi-davidson", settings%subsystem_solver) .or. &
             settings%subsystem_solver == "tcg")) then
            call settings%log("Non-symmetric Hessian not supported with augmented "// &
                              "Hessian, Jacobi-Davidson or truncated conjugate "// &
                              "gradient solvers.", verbosity_error, .true.)
            error = 1
            return
        end if

        ! check whether projection function is passed, unless the nested stability
        ! check settings inherited it, whose check below prints the same warning
        if (associated(settings%project) .and. &
            .not. associated(settings%stability_settings%project, settings%project)) &
            call settings%log(project_warning_msg, verbosity_warning)

        ! check the nested stability check settings, with the defaults they inherit
        ! from the solver already resolved, since the internal stability check also
        ! runs without a requested stability check, so that settings it would reject
        ! are reported before the optimization starts
        call stability_sanity_check(settings%stability_settings, n_param, error)

    end subroutine solver_sanity_check

    subroutine stability_sanity_check(settings, n_param, error)
        !
        ! this subroutine performs a sanity check for stability input parameters
        !
        type(stability_settings_type), intent(inout) :: settings
        integer(ip), intent(in) :: n_param
        integer(ip), intent(out) :: error
        character(len=300) :: msg

        ! initialize error flag
        error = 0

        ! check that number of parameters is positive
        if (n_param < 1) then
            call settings%log("Number of parameters should be larger than 0.", &
                              verbosity_error, .true.)
            error = 1
            return
        end if

        ! convert strings to lowercase
        settings%diag_solver = string_to_lowercase(settings%diag_solver)

        ! check whether the extra trial vector function has anything to fill
        if (associated(settings%get_extra_trial_vectors) .and. &
            settings%n_extra_trial_vectors <= 0) then
            call settings%log("Number of extra trial vectors should be larger than "// &
                              "zero when extra trial vector function is passed.", &
                              verbosity_error, .true.)
            error = 1
            return
        end if

        ! resolve the automatic number of random trial vectors
        if (settings%n_random_trial_vectors < 0) &
            settings%n_random_trial_vectors = standalone_n_random_trial_vectors

        ! check that number of random trial vectors is below number of parameters
        if (settings%n_random_trial_vectors > n_param / 2) then
            settings%n_random_trial_vectors = n_param / 2
            write(msg, '(A, I0, A)') random_trial_vector_warning_msg//" Setting to ", &
                settings%n_random_trial_vectors, "."
            call settings%log(msg, verbosity_warning)
        end if

        ! check that the trial space is not empty, also for a negative number of extra
        ! trial vectors, which requests none, and once the number of random trial
        ! vectors is reduced
        if (settings%n_extra_trial_vectors <= 0 .and. &
            settings%n_random_trial_vectors == 0) then
            call settings%log("Vanishing number of trial vectors for stability "// &
                              "check.", verbosity_error, .true.)
            error = 1
            return
        end if

        ! check for character options
        if (.not. any(settings%diag_solver == diag_solver_options)) then
            call settings%log( &
                "Diagonalization solver option unknown. Possible values are "// &
                option_list(diag_solver_options)//".", verbosity_error, .true.)
            error = 1
            return
        end if

        ! check whether non-symmetric Hessian is requested with a solver that cannot
        ! handle it
        if (.not. settings%hess_symm .and. settings%diag_solver == "jacobi-davidson") &
            then
            call settings%log("Non-symmetric Hessian not supported with "// &
                              "Jacobi-Davidson solver.", verbosity_error, .true.)
            error = 1
            return
        end if

        ! check whether stopping at an instability is requested for a non-symmetric
        ! Hessian, whose Ritz values are no upper bounds on its eigenvalues
        if (.not. settings%hess_symm .and. settings%stop_on_instability) then
            call settings%log("Stopping at an instability not supported with "// &
                              "non-symmetric Hessian.", verbosity_error, .true.)
            error = 1
            return
        end if

        ! check that the stability threshold does not count a minimum as unstable
        if (settings%stability_thresh > 0.0_rp) then
            call settings%log("Stability threshold should not be positive.", &
                              verbosity_error, .true.)
            error = 1
            return
        end if

        ! check whether projection function is passed
        if (associated(settings%project)) &
            call settings%log(project_warning_msg, verbosity_warning)

    end subroutine stability_sanity_check

    subroutine add_error_origin(error_code, error_origin, settings)
        !
        ! this subroutine modifies the error code by adding the error's origin if it is
        ! not already added
        !
        class(settings_type), intent(in) :: settings
        integer(ip), intent(inout) :: error_code
        integer(ip), intent(in) :: error_origin

        if (error_code > 0) then
            if (error_code < 100) error_code = error_origin + error_code
        else if (error_code < 0) then
            call settings%log("Negative error code encountered.", verbosity_error, &
                              .true.)
            error_code = error_origin + 1
        end if

    end subroutine add_error_origin

    function string_to_lowercase(str) result(lower_str)
        !
        ! this function converts a string to lower case
        !
        character(len=*), intent(in) :: str
        character(len=len(str)) :: lower_str
        integer(ip) :: i, code
        integer(ip), parameter :: ascii_upper_lower_diff = 32

        do i = 1, len(str)
            code = iachar(str(i:i))
            if (code >= iachar("A") .and. code <= iachar("Z")) then
                lower_str(i:i) = achar(code + ascii_upper_lower_diff)
            else
                lower_str(i:i) = str(i:i)
            end if
        end do

    end function string_to_lowercase

    logical function string_in(substr, str)
        !
        ! this function checks whether a string is contained in another string
        !
        character(len=*), intent(in) :: substr, str

        string_in = index(str, substr) > 0

    end function string_in

    subroutine perturb_vector(vector, noise, settings, error)
        !
        ! this subroutine perturbs the input vector by a random vector whose norm is
        ! the given fraction of the norm of the input vector; the random vector is
        ! projected before it is scaled, so that the perturbation keeps its full size
        ! when a projection removes redundant parameters, and the input vector is
        ! left unperturbed in the unlikely case that the projection removes every
        ! random vector drawn
        !
        real(rp), intent(inout) :: vector(:)
        real(rp), intent(in) :: noise
        type(solver_settings_type), intent(in) :: settings
        integer(ip), intent(out) :: error

        integer(ip) :: n_param, i_draw
        real(rp), allocatable :: random_vector(:)
        integer(ip), parameter :: max_draws = 10
        real(rp), external :: dnrm2

        ! initialize error flag
        error = 0

        ! number of parameters
        n_param = size(vector)

        ! draw a random vector until one is left after the projection
        allocate(random_vector(n_param))
        do i_draw = 1, max_draws
            call random_number(random_vector)
            random_vector = 2.0_rp * random_vector - 1.0_rp
            if (associated(settings%project)) then
                call settings%project(random_vector, error, settings%context)
                call add_error_origin(error, error_project, settings)
                if (error /= 0) return
            end if
            if (dnrm2(n_param, random_vector, 1_ip) >= numerical_zero) then
                ! scale the random vector to the given fraction of the input norm
                vector = vector + noise * dnrm2(n_param, vector, 1_ip) * &
                         random_vector / dnrm2(n_param, random_vector, 1_ip)
                exit
            end if
        end do
        deallocate(random_vector)

    end subroutine perturb_vector

    function noise_level(settings, local_region) result(noise)
        !
        ! this function returns the relative size of the random perturbation of the
        ! gradient, the micro-iteration reduction factor of the local or the global
        ! region, or zero if the gradient is not perturbed; the subsystem solvers
        ! never solve to a smaller relative residual than this, so that they are not
        ! asked to resolve the perturbation they were given, whichever region the
        ! subproblem turns out to be solved in
        !
        type(solver_settings_type), intent(in) :: settings
        logical, intent(in) :: local_region
        real(rp) :: noise

        if (.not. settings%grad_noise) then
            noise = 0.0_rp
        else if (local_region) then
            noise = settings%local_red_factor
        else
            noise = settings%global_red_factor
        end if

    end function noise_level

    function micro_conv_tol(level_shift, noise, initial_residual_norm, settings) &
        result(conv_tol)
        !
        ! this function returns the residual norm at which the micro iterations stop,
        ! the reduction factor of the local region, where the level shift vanishes, or
        ! of the global region, never below the gradient perturbation, relative to the
        ! initial residual norm and never below the residual norm floor
        !
        real(rp), intent(in) :: level_shift, noise, initial_residual_norm
        type(solver_settings_type), intent(in) :: settings
        real(rp) :: conv_tol

        real(rp) :: red_factor

        if (abs(level_shift) < level_shift_local_thres) then
            red_factor = settings%local_red_factor
        else
            red_factor = settings%global_red_factor
        end if
        conv_tol = max(max(red_factor, noise) * initial_residual_norm, &
                       residual_norm_floor)

    end function micro_conv_tol

    function option_list(options) result(list)
        !
        ! this function lists options as quoted words separated by commas
        !
        character(len=*), intent(in) :: options(:)
        character(len=:), allocatable :: list

        integer(ip) :: i

        ! append every option in quotes
        list = ""
        do i = 1, size(options)
            if (i > 1) list = list//", "
            list = list//""""//trim(options(i))//""""
        end do

    end function option_list

end module opentrustregion
