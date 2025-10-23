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

    ! mathematical constants
    real(rp), parameter :: pi = 4.0_rp * atan(1.0_rp)

    ! define trust region parameters
    real(rp), parameter :: trust_radius_shrink_ratio = 0.25_rp, &
                           trust_radius_expand_ratio = 0.75_rp, &
                           trust_radius_shrink_factor = 0.7_rp, &
                           trust_radius_expand_factor = 1.2_rp

    ! define error codes
    integer(ip), parameter :: error_solver = 100, error_stability_check = 200, &
                              error_obj_func = 1100, error_update_orbs = 1200, &
                              error_hess_x = 1300, error_precond = 1400, &
                              error_conv_check = 1500

    ! define useful parameters
    real(rp), parameter :: numerical_zero = 1e-14_rp, precond_floor = 1e-10_rp, &
                           hess_symm_thres = 1e-12_rp
    integer(ip), parameter :: jacobi_davidson_start = 30

    ! interfaces for callback functions
    abstract interface
        subroutine hess_x_type(x, hess_x, error)
            import :: rp, ip

            real(rp), intent(in), target :: x(:)
            real(rp), intent(out), target :: hess_x(:)
            integer(ip), intent(out) :: error
        end subroutine hess_x_type
    end interface

    abstract interface
        subroutine update_orbs_type(kappa, func, grad, h_diag, hess_x_funptr, error)
            import :: rp, hess_x_type, ip

            real(rp), intent(in), target :: kappa(:)
            real(rp), intent(out) :: func
            real(rp), intent(out), target :: grad(:), h_diag(:)
            procedure(hess_x_type), intent(out), pointer :: hess_x_funptr
            integer(ip), intent(out) :: error
        end subroutine update_orbs_type
    end interface

    abstract interface
        function obj_func_type(kappa, error) result(func)
            import :: rp, ip

            real(rp), intent(in), target :: kappa(:)
            integer(ip), intent(out) :: error
            real(rp) :: func
        end function obj_func_type
    end interface

    abstract interface
        subroutine precond_type(residual, mu, precond_residual, error)
            import :: rp, ip

            real(rp), intent(in), target :: residual(:)
            real(rp), intent(in) :: mu
            real(rp), intent(out), target :: precond_residual(:)
            integer(ip), intent(out) :: error
        end subroutine precond_type
    end interface

    abstract interface
        function conv_check_type(error) result(converged)
            import :: ip

            integer(ip), intent(out) :: error
            logical :: converged
        end function conv_check_type
    end interface

    abstract interface
        subroutine logger_type(message)
            character(*), intent(in) :: message
        end subroutine logger_type
    end interface

    ! derived type for solver settings
    type, abstract :: settings_type
        logical :: jacobi_davidson, initialized = .false.
        real(rp) :: conv_tol
        integer(ip) :: n_random_trial_vectors, verbose
        procedure(precond_type), pointer, nopass :: precond
        procedure(logger_type), pointer, nopass :: logger
    contains
        procedure :: log
        procedure(init_type), deferred :: init
    end type

    abstract interface
        subroutine init_type(self, error)
            import :: settings_type, ip

            class(settings_type), intent(out) :: self
            integer(ip), intent(out) :: error
        end subroutine init_type
    end interface

    type, extends(settings_type) :: solver_settings_type
        logical :: stability, line_search, davidson, prefer_jacobi_davidson
        real(rp) :: start_trust_radius, global_red_factor, local_red_factor
        integer(ip) :: n_macro, n_micro, seed
        procedure(conv_check_type), pointer, nopass :: conv_check
    contains
        procedure :: init => init_solver_settings, print_results
    end type

    type, extends(settings_type) :: stability_settings_type
        integer(ip) :: n_iter
    contains
        procedure :: init => init_stability_settings
    end type

    ! default settings
    logical, parameter :: solver_stability_default = .false., &
                          solver_line_search_default = .false., &
                          solver_davidson_default = .true., &
                          solver_jacobi_davidson_default = .false., &
                          solver_prefer_jacobi_davidson_default = .false., &
                          stability_jacobi_davidson_default = .true.
    real(rp), parameter :: solver_conv_tol_default = 1e-5_rp, &
                           solver_start_trust_radius_default = 0.4_rp, &
                           solver_global_red_factor_default = 1e-3_rp, &
                           solver_local_red_factor_default = 1e-4_rp, &
                           stability_conv_tol_default = 1e-8_rp
    integer(ip), parameter :: solver_n_random_trial_vectors_default = 1, &
                              solver_n_macro_default = 150, &
                              solver_n_micro_default = 50, &
                              solver_seed_default = 42, &
                              solver_verbose_default = 0, &
                              stability_n_random_trial_vectors_default = 20, &
                              stability_n_iter_default = 100, &
                              stability_verbose_default = 0

    type(solver_settings_type) :: default_solver_settings = &
        solver_settings_type(precond = null(), conv_check = null(), logger = null(), &
                             stability = solver_stability_default, &
                             line_search = solver_line_search_default, &
                             davidson = solver_davidson_default, &
                             jacobi_davidson = solver_jacobi_davidson_default, &
                             prefer_jacobi_davidson = &
                                solver_prefer_jacobi_davidson_default, &
                             initialized = .true., conv_tol = solver_conv_tol_default, &
                             start_trust_radius = solver_start_trust_radius_default, &
                             global_red_factor = solver_global_red_factor_default, &
                             local_red_factor = solver_local_red_factor_default, &
                             n_random_trial_vectors = &
                                solver_n_random_trial_vectors_default, &
                             n_macro = solver_n_macro_default, &
                             n_micro = solver_n_micro_default, &
                             seed = solver_seed_default, &
                             verbose = solver_verbose_default)
    type(stability_settings_type) :: default_stability_settings = &
        stability_settings_type(precond = null(), logger = null(), &
                                jacobi_davidson = stability_jacobi_davidson_default, &
                                initialized = .true., &
                                conv_tol = stability_conv_tol_default, &
                                n_random_trial_vectors = &
                                    stability_n_random_trial_vectors_default, &
                                n_iter = stability_n_iter_default, &
                                verbose = stability_verbose_default)

    ! interfaces for solver and stability_check subroutines
    interface
        subroutine solver_type(update_orbs, obj_func, n_param, error, settings)

            import :: ip, update_orbs_type, obj_func_type, solver_settings_type

            procedure(update_orbs_type), intent(in), pointer :: update_orbs
            procedure(obj_func_type), intent(in), pointer :: obj_func
            integer(ip), intent(in) :: n_param
            integer(ip), intent(out) :: error
            type(solver_settings_type), intent(inout) :: settings

        end subroutine solver_type
    end interface

    interface
        subroutine stability_check_type(h_diag, hess_x, stable, error, settings, kappa)

            import:: rp, ip, hess_x_type, stability_settings_type

            real(rp), intent(in) :: h_diag(:)
            procedure(hess_x_type), intent(in), pointer :: hess_x
            logical, intent(out) :: stable
            integer(ip), intent(out) :: error
            type(stability_settings_type), intent(inout) :: settings
            real(rp), intent(out), optional :: kappa(:)

        end subroutine stability_check_type
    end interface

    ! create function pointers to ensure that routines comply with interface
    procedure(solver_type), pointer :: solver_ptr => solver
    procedure(stability_check_type), pointer :: stability_check_ptr => stability_check

    ! define global variables
    integer(ip) :: tot_orb_update = 0, tot_hess_x = 0

contains

    subroutine solver(update_orbs, obj_func, n_param, error, settings)
        !
        ! this subroutine is the main solver for orbital optimization
        !
        procedure(update_orbs_type), intent(in), pointer :: update_orbs
        procedure(obj_func_type), intent(in), pointer :: obj_func
        integer(ip), intent(in) :: n_param
        integer(ip), intent(out) :: error
        type(solver_settings_type), intent(inout) :: settings

        type(stability_settings_type) :: stability_settings
        real(rp) :: trust_radius, func, grad_norm, grad_rms, mu, new_func, n_kappa, &
                    kappa_norm
        real(rp), allocatable :: kappa(:), grad(:), h_diag(:), solution(:), &
                                 precond_kappa(:)
        logical :: max_precision_reached, macro_converged, stable, &
                   jacobi_davidson_started, conv_check_passed
        integer(ip) :: imacro, imicro, imicro_jacobi_davidson, i
        character(300) :: msg
        integer(ip), parameter :: stability_n_points = 21
        procedure(hess_x_type), pointer :: hess_x_funptr
        real(rp), external :: dnrm2, ddot

        ! initialize error flag
        error = 0

        ! initialize settings
        if (.not. settings%initialized) then
            call settings%init(error)
            call add_error_origin(error, error_solver, settings)
            if (error /= 0) return
            call settings%log("Settings were not initialized. All settings are set "// &
                              "to default values", 2)
        end if

        ! initialize maximum precision convergence
        max_precision_reached = .false.

        ! initialize macroiteration convergence
        macro_converged = .false.

        ! initialize stabilty boolean
        stable = .true.

        ! initialize random number generator
        call init_rng(settings%seed)

        ! initialize starting trust radius
        trust_radius = settings%start_trust_radius

        ! print header
        call settings%log(repeat("-", 109), 3)
        call settings%log(" Iteration |     Objective function     | Gradient RMS |"// &
                          " Level shift |   Micro    | Trust radius | Step size ", 3)
        call settings%log("           |                            |              |"// &
                          "             | iterations |              |           ", 3)
        call settings%log(repeat("-", 109), 3)

        ! allocate arrays
        allocate(kappa(n_param), grad(n_param), h_diag(n_param), solution(n_param), &
                 precond_kappa(n_param))

        ! initialize orbital rotation matrix
        kappa = 0.0_rp

        do imacro = 1, settings%n_macro
            if (.not. max_precision_reached) then
                ! calculate cost function, gradient and Hessian diagonal
                call update_orbs(kappa, func, grad, h_diag, hess_x_funptr, error)
                call add_error_origin(error, error_update_orbs, settings)
                if (error /= 0) return

                ! perform sanity check
                if (imacro == 1) then
                    call sanity_check(settings, n_param, grad, error)
                    call add_error_origin(error, error_solver, settings)
                    if (error /= 0) return
                end if

                ! calculate gradient norm
                grad_norm = dnrm2(n_param, grad, 1)

                ! calculate RMS gradient
                grad_rms = grad_norm / sqrt(real(n_param, kind=rp))

                ! log results
                if (imacro == 1) then
                    call settings%print_results(imacro - 1, func, grad_rms)
                else
                    if (settings%davidson) then
                        kappa_norm = dnrm2(n_param, kappa, 1)
                    else
                        if (associated(settings%precond)) then
                            call settings%precond(kappa, 0.0_rp, precond_kappa, error)
                            call add_error_origin(error, error_precond, settings)
                            if (error /= 0) return
                        else
                            call abs_diag_precond(kappa, h_diag, precond_kappa)
                        end if
                        kappa_norm = sqrt(ddot(n_param, kappa, 1, precond_kappa, 1))
                        mu = 0.0_rp
                        jacobi_davidson_started = .false.
                    end if
                    if (.not. stable) then
                        call settings%print_results(imacro - 1, func, grad_rms, &
                                                    kappa_norm=kappa_norm)
                        stable = .true.
                    else if (jacobi_davidson_started) then
                        call settings%print_results(imacro - 1, func, grad_rms, &
                                                    level_shift=-mu, n_micro=imicro, &
                                                    imicro_jacobi_davidson=&
                                                    imicro_jacobi_davidson, &
                                                    trust_radius=trust_radius, &
                                                    kappa_norm=kappa_norm)
                    else
                        call settings%print_results(imacro - 1, func, grad_rms, &
                                                    level_shift=-mu, n_micro=imicro, &
                                                    trust_radius=trust_radius, &
                                                    kappa_norm=kappa_norm)
                    end if
                end if
            end if

            ! check for convergence and stability
            if (associated(settings%conv_check)) then
                conv_check_passed = settings%conv_check(error)
                call add_error_origin(error, error_conv_check, settings)
                if (error /= 0) return
            else
                conv_check_passed = .false.
            end if
            if (grad_rms < settings%conv_tol .or. max_precision_reached .or. &
                conv_check_passed) then
                ! always perform stability check if starting at stationary point
                if (settings%stability .or. imacro == 1) then
                    call stability_settings%init(error)
                    call add_error_origin(error, error_solver, settings)
                    if (error /= 0) return
                    stability_settings%precond => settings%precond
                    stability_settings%verbose = settings%verbose
                    stability_settings%logger => settings%logger
                    call stability_check(h_diag, hess_x_funptr, stable, error, &
                                         stability_settings, kappa=kappa)
                    call add_error_origin(error, error_stability_check, settings)
                    if (error /= 0) return
                    if (.not. stable) then
                        ! logarithmic line search
                        do i = 1, stability_n_points
                            n_kappa = 10.0_rp**(-(i - 1) / &
                                                real(stability_n_points - 1, rp) * &
                                                10.0_rp)
                            new_func = obj_func(n_kappa*kappa, error)
                            call add_error_origin(error, error_obj_func, settings)
                            if (error /= 0) return
                            if (new_func < func) then
                                kappa = n_kappa*kappa
                                exit
                            end if
                        end do
                        if (new_func >= func) then
                            call settings%log("Line search was unable to find "// &
                                              "lower objective function along "// &
                                              "unstable mode.", 1, .true.)
                            error = error_solver + 1
                            return
                        else if (imacro == 1) then
                            call settings%log("Started at saddle point. The "// &
                                              "algorithm will continue by moving "// &
                                              "along eigenvector direction "// &
                                              "corresponding to negative eigenvalue.", &
                                              1, .true.)
                        else
                            call settings%log("Reached saddle point. This is "// &
                                              "likely due to symmetry and can be "// &
                                              "avoided by increasing the number of "// &
                                              "random trial vectors. The algorithm "// &
                                              "will continue by moving along "// &
                                              "eigenvector direction corresponding "// &
                                              "to negative eigenvalue.", 1, .true.)
                        end if
                        max_precision_reached = .false.
                        cycle
                    else
                        macro_converged = .true.
                        exit
                    end if
                else
                    macro_converged = .true.
                    exit
                end if
            end if

            if (settings%davidson) then
                ! solve trust region subproblem with (Jacobi-)Davidson
                call level_shifted_davidson(func, grad, grad_norm, h_diag, n_param, &
                                            obj_func, hess_x_funptr, settings, &
                                            trust_radius, solution, mu, imicro, &
                                            imicro_jacobi_davidson, &
                                            jacobi_davidson_started, &
                                            max_precision_reached, error)
                call add_error_origin(error, error_solver, settings)
                if (error /= 0) return
            else
                ! solve trust region subproblem with truncated conjugate gradient
                call truncated_conjugate_gradient(func, grad, h_diag, n_param, &
                                                  obj_func, hess_x_funptr, &
                                                  settings, trust_radius, solution, &
                                                  imicro, max_precision_reached, error)
                call add_error_origin(error, error_solver, settings)
                if (error /= 0) return
            end if

            ! perform line search
            if (max_precision_reached) then
                n_kappa = 0.0_rp
            else if (settings%line_search) then
                n_kappa = bracket(obj_func, solution, 0.0_rp, 1.0_rp, settings, error)
                call add_error_origin(error, error_solver, settings)
                if (error /= 0) return
            else
                n_kappa = 1.0_rp
            end if

            ! set orbital rotation
            kappa = n_kappa*solution

            ! flush output
            flush(stdout)
            flush(stderr)

        end do

        ! deallocate arrays
        deallocate(kappa, grad, h_diag, solution, precond_kappa)

        ! increment total number of orbital updates
        tot_orb_update = tot_orb_update + imacro

        ! stop if no convergence
        if (.not. macro_converged) then
            call settings%log("Orbital optimization has not converged!", 1)
            error = error_solver + 1
            return
        end if

        ! finish logging
        call settings%log(repeat("-", 109), 3)
        write (msg, '(A, I0)') "Total number of Hessian linear transformations: ", &
            tot_hess_x
        call settings%log(msg, 3)
        write (msg, '(A, I0)') "Total number of orbital updates: ", tot_orb_update
        call settings%log(msg, 3)

        ! reset global counter variables
        tot_orb_update = 0
        tot_hess_x = 0

        ! flush output
        flush (stdout)
        flush (stderr)

    end subroutine solver

    subroutine stability_check(h_diag, hess_x_funptr, stable, error, settings, kappa)
        !
        ! this subroutine performs a stability check
        !
        real(rp), intent(in) :: h_diag(:)
        procedure(hess_x_type), intent(in), pointer :: hess_x_funptr
        logical, intent(out) :: stable
        integer(ip), intent(out) :: error
        type(stability_settings_type), intent(inout) :: settings
        real(rp), intent(out), optional :: kappa(:)

        integer(ip) :: n_param, n_trial, i, iter
        real(rp), allocatable :: solution(:), h_solution(:), residual(:), &
                                 basis_vec(:), h_basis_vec(:), red_space_basis(:, :), &
                                 h_basis(:, :), red_space_hess(:, :), &
                                 red_space_solution(:), red_space_hess_vec(:)
        real(rp) :: eigval, minres_tol, stability_rms
        character(300) :: msg
        real(rp), parameter :: stability_thresh = -1e-2_rp
        real(rp), external :: dnrm2, ddot
        external :: dgemm, dgemv

        ! initialize error flag
        error = 0

        ! initialize settings
        if (.not. settings%initialized) then
            call settings%init(error)
            call add_error_origin(error, error_stability_check, settings)
            if (error /= 0) return
            call settings%log("Settings were not initialized. All settings are set "// &
                              "to default values", 2)
        end if

        ! get number of parameters
        n_param = size(h_diag)

        ! check that number of random trial vectors is below number of parameters
        if (settings%n_random_trial_vectors > n_param/2) then
            settings%n_random_trial_vectors = n_param/2
            write (msg, '(A, I0, A)') "Number of random trial vectors should be "// &
                "smaller than half the number of parameters. Setting to ", &
                settings%n_random_trial_vectors, "."
            call settings%log(msg, 2)
        end if

        ! generate trial vectors
        allocate(red_space_basis(n_param, 1 + settings%n_random_trial_vectors))
        red_space_basis(:, 1) = 0.0_rp
        red_space_basis(minloc(h_diag), 1) = 1.0_rp
        call generate_random_trial_vectors(red_space_basis, settings, error)
        call add_error_origin(error, error_stability_check, settings)
        if (error /= 0) return

        ! number of trial vectors
        n_trial = size(red_space_basis, 2)

        ! calculate linear transformations of basis vectors
        allocate(h_basis(n_param, n_trial))
        do i = 1, n_trial
            call hess_x_funptr(red_space_basis(:, i), h_basis(:, i), error)
            call add_error_origin(error, error_hess_x, settings)
            if (error /= 0) return
        end do

        ! construct augmented Hessian in reduced space
        allocate(red_space_hess(n_trial, n_trial))
        call dgemm("T", "N", n_trial, n_trial, n_param, 1.0_rp, red_space_basis, &
                   n_param, h_basis, n_param, 0.0_rp, red_space_hess, n_trial)

        ! allocate arrays used throughout Davidson procedure
        allocate(red_space_solution(n_trial), solution(n_param), h_solution(n_param), &
                 residual(n_param), basis_vec(n_param), h_basis_vec(n_param))

        ! loop over iterations
        do iter = 1, settings%n_iter
            ! solve reduced space problem
            call symm_mat_min_eig(red_space_hess, eigval, red_space_solution, &
                                  settings, error)
            call add_error_origin(error, error_stability_check, settings)
            if (error /= 0) return

            ! get full space solution
            call dgemv("N", size(red_space_basis, 1), size(red_space_basis, 2), &
                       1.0_rp, red_space_basis, size(red_space_basis, 1), &
                       red_space_solution, 1, 0.0_rp, solution, 1)

            ! calculate Hessian linear transformation of solution
            call dgemv("N", n_param, size(h_basis, 2), 1.0_rp, h_basis, n_param, &
                       red_space_solution, 1, 0.0_rp, h_solution, 1)

            ! calculate residual
            residual = h_solution - eigval*solution

            ! check convergence
            stability_rms = dnrm2(n_param, residual, 1) / sqrt(real(n_param, kind=rp))
            if (stability_rms < settings%conv_tol) exit

            if (.not. settings%jacobi_davidson .or. iter <= jacobi_davidson_start) then
                ! precondition residual
                if (associated(settings%precond)) then
                    call settings%precond(residual, 0.0_rp, basis_vec, error)
                    call add_error_origin(error, error_precond, settings)
                    if (error /= 0) return
                else
                    call level_shifted_diag_precond(residual, 0.0_rp, h_diag, basis_vec)
                end if

                ! orthonormalize to current orbital space to get new basis vector
                call gram_schmidt(basis_vec, red_space_basis, settings, error)
                call add_error_origin(error, error_stability_check, settings)
                if (error /= 0) return

                ! add linear transformation of new basis vector
                call hess_x_funptr(basis_vec, h_basis_vec, error)
                call add_error_origin(error, error_hess_x, settings)
                if (error /= 0) return

                ! increment Hessian linear transformations
                tot_hess_x = tot_hess_x + 1

            else
                ! solve Jacobi-Davidson correction equations
                minres_tol = 3.0_rp ** (-(iter - jacobi_davidson_start - 1))
                call minres(-residual, hess_x_funptr, solution, eigval, minres_tol, &
                            basis_vec, h_basis_vec, settings, error)
                call add_error_origin(error, error_stability_check, settings)
                if (error /= 0) return

                ! orthonormalize to current orbital space to get new basis 
                ! vector
                call gram_schmidt(basis_vec, red_space_basis, settings, error, &
                                  lin_trans_vector=h_basis_vec, lin_trans_space=h_basis)
                call add_error_origin(error, error_stability_check, settings)
                if (error /= 0) return

                ! check if resulting linear transformation still respects Hessian 
                ! symmetry which can happen due to numerical noise accumulation
                if (abs(ddot(n_param, red_space_basis(:, size(red_space_basis, 2)), 1, &
                             h_basis_vec, 1) - &
                        ddot(n_param, basis_vec, 1, &
                             h_basis(:, size(red_space_basis, 2)), 1)) > &
                    hess_symm_thres) then
                    call hess_x_funptr(basis_vec, h_basis_vec, error)
                    call add_error_origin(error, error_hess_x, settings)
                    if (error /= 0) return
                end if
                
            end if

            ! add new trial vector to orbital space
            call add_column(red_space_basis, basis_vec)

            ! add linear transformation of new basis vector
            call add_column(h_basis, h_basis_vec)

            ! construct new reduced space Hessian
            allocate(red_space_hess_vec(size(red_space_basis, 2)))
            call dgemv("T", n_param, size(red_space_basis, 2), 1.0_rp, &
                       red_space_basis, n_param, &
                       h_basis(:, size(red_space_basis, 2)), 1, 0.0_rp, &
                       red_space_hess_vec, 1)
            call extend_symm_matrix(red_space_hess, red_space_hess_vec)
            deallocate(red_space_hess_vec)

            ! reallocate reduced space solution
            deallocate(red_space_solution)
            allocate(red_space_solution(size(red_space_basis, 2)))

        end do

        ! check if stability check has converged
        if (stability_rms >= settings%conv_tol) &
            call settings%log("Stability check has not converged in the given "// &
                              "number of iterations.", 1, .true.)

        ! increment total number of Hessian linear transformations
        tot_hess_x = tot_hess_x + size(red_space_basis, 2)

        ! determine if saddle point
        stable = eigval > stability_thresh
        
        if (stable) then
            if (present(kappa)) kappa = 0.0_rp
        else
            if (present(kappa)) kappa = solution
            write (msg, '(A, F0.4)') "Solution not stable. Lowest eigenvalue: ", eigval
            call settings%log(msg, 1, .true.)
        end if

        ! deallocate quantities from Davidson iterations
        deallocate(solution, h_solution, residual, basis_vec, h_basis_vec, &
                   red_space_solution, red_space_hess, h_basis, red_space_basis)

        ! flush output
        flush (stdout)
        flush (stderr)

    end subroutine stability_check

    subroutine newton_step(aug_hess, grad_norm, red_space_basis, solution, &
                           red_space_solution, settings, error)
        !
        ! this subroutine performs a Newton step by solving the Newton equations in
        ! reduced space without a level shift
        !
        real(rp), intent(in) :: aug_hess(:, :), grad_norm, red_space_basis(:, :)
        class(settings_type), intent(in) :: settings
        real(rp), intent(out) :: solution(:), red_space_solution(:)
        integer(ip), intent(out) :: error

        integer(ip) :: n_param, n_red, lwork, info
        integer(ip), allocatable :: ipiv(:)
        real(rp), allocatable :: red_hess(:, :), work(:)
        character(300) :: msg
        external :: dsysv, dgemv

        ! initialize error flag
        error = 0

        ! number of parameters
        n_param = size(solution)

        ! reduced space size
        n_red = size(red_space_basis, 2)

        ! reduced space Hessian
        allocate(red_hess(n_red, n_red))
        red_hess = aug_hess(2:, 2:)

        ! set gradient
        red_space_solution = 0.0_rp
        red_space_solution(1) = -grad_norm

        ! query optimal workspace size
        lwork = -1
        allocate(ipiv(n_red), work(1))
        call dsysv("U", n_red, 1, red_hess, n_red, ipiv, red_space_solution, n_red, &
                   work, lwork, info)
        lwork = int(work(1), kind=ip)
        deallocate(work)
        allocate(work(lwork))

        ! solve linear system
        call dsysv("U", n_red, 1, red_hess, n_red, ipiv, red_space_solution, n_red, &
                   work, lwork, info)

        ! deallocate arrays
        deallocate(red_hess, ipiv, work)

        ! check for errors
        if (info /= 0) then
            write (msg, '(A, I0)') "Linear solver failed: Error in DSYSV, info = ", info
            call settings%log(msg, 1, .true.)
            error = 1
            return
        end if

        ! get solution in full space
        call dgemv("N", n_param, n_red, 1.0_rp, red_space_basis, n_param, &
                   red_space_solution, 1, 0.0_rp, solution, 1)

    end subroutine newton_step

    subroutine bisection(aug_hess, grad_norm, red_space_basis, trust_radius, solution, &
                         red_space_solution, mu, bracketed, settings, error)
        !
        ! this subroutine performs bisection to find the parameter alpha that matches
        ! the desired trust radius
        !
        real(rp), intent(inout) :: aug_hess(:, :)
        real(rp), intent(in) :: grad_norm, red_space_basis(:, :), trust_radius
        class(settings_type), intent(in) :: settings
        real(rp), intent(out) :: solution(:), red_space_solution(:), mu
        logical, intent(out) :: bracketed
        integer(ip), intent(out) :: error

        real(rp) :: lower_alpha, middle_alpha, upper_alpha, lower_trust_dist, &
                    middle_trust_dist, upper_trust_dist
        integer(ip) :: n_param, n_red, iter
        real(rp), parameter :: lower_alpha_bound = 1e-4_rp, &
                               upper_alpha_bound = 1e6_rp, &
                               alpha_conv_factor = 1e-12_rp
        real(rp), external :: dnrm2

        ! initialize error flag
        error = 0

        ! initialize bracketing flag
        bracketed = .false.

        ! number of parameters
        n_param = size(solution)

        ! reduced space size
        n_red = size(red_space_basis, 2)

        ! lower and upper bracket for alpha
        lower_alpha = lower_alpha_bound
        upper_alpha = upper_alpha_bound

        ! solve reduced space problem with scaled gradient
        call get_ah_lowest_eigenvec(lower_alpha)
        if (error /= 0) return
        lower_trust_dist = dnrm2(n_param, solution, 1) - trust_radius
        call get_ah_lowest_eigenvec(upper_alpha)
        if (error /= 0) return
        upper_trust_dist = dnrm2(n_param, solution, 1) - trust_radius

        ! check if trust region is within bracketing range
        if ((lower_trust_dist*upper_trust_dist) > 0.0_rp) then
            solution = 0.0_rp
            red_space_solution = 0.0_rp
            mu = 0.0_rp
            return
        end if

        ! get middle alpha
        middle_alpha = sqrt(upper_alpha*lower_alpha)
        call get_ah_lowest_eigenvec(middle_alpha)
        if (error /= 0) return
        middle_trust_dist = dnrm2(n_param, solution, 1) - trust_radius

        ! perform bisection to find root, converge to relative threshold to avoid 
        ! precision issues
        iter = 0
        do while (upper_alpha - lower_alpha > alpha_conv_factor * upper_alpha)
            ! targeted trust radius is in upper bracket
            if (lower_trust_dist*middle_trust_dist > 0.0_rp) then
                lower_alpha = middle_alpha
                lower_trust_dist = middle_trust_dist
                ! targeted trust radius is in lower bracket
            else
                upper_alpha = middle_alpha
                upper_trust_dist = middle_trust_dist
            end if
            ! get new middle alpha
            middle_alpha = sqrt(upper_alpha*lower_alpha)
            call get_ah_lowest_eigenvec(middle_alpha)
            if (error /= 0) return
            middle_trust_dist = dnrm2(n_param, solution, 1) - trust_radius
            ! check if maximum number of iterations is reached
            iter = iter + 1
            if (iter > 100) then
                call settings%log("Maximum number of bisection iterations reached.", &
                                  1, .true.)
                error = 1
                return
            end if
        end do

        bracketed = .true.

    contains

        subroutine get_ah_lowest_eigenvec(alpha)
            !
            ! this subroutine returns the lowest eigenvector for an augmented Hessian
            !
            real(rp), intent(in) :: alpha

            real(rp), allocatable :: eigvec(:)
            external :: dgemv

            ! finish construction of augmented Hessian
            aug_hess(1, 2) = alpha*grad_norm
            aug_hess(2, 1) = alpha*grad_norm

            ! allocate eigenvector
            allocate(eigvec(n_red + 1))

            ! perform eigendecomposition and get lowest eigenvalue and corresponding
            ! eigenvector
            call symm_mat_min_eig(aug_hess, mu, eigvec, settings, error)
            if (error /= 0) return

            ! check if eigenvector has level-shift component
            if (abs(eigvec(1)) <= numerical_zero) then
                call settings%log("Trial subspace too small. Increase "// &
                    "n_random_trial_vectors.", 1, .true.)
                error = 1
                return
            end if

            ! scale eigenvector such that first element is equal to one and divide by
            ! alpha to get solution in reduced space
            red_space_solution = eigvec(2:)/eigvec(1)/alpha
            deallocate(eigvec)

            ! get solution in full space
            call dgemv("N", n_param, n_red, 1.0_rp, red_space_basis, n_param, &
                       red_space_solution, 1, 0.0_rp, solution, 1)

        end subroutine get_ah_lowest_eigenvec

    end subroutine bisection

    function bracket(obj_func, kappa, lower, upper, settings, error) result(n_kappa)
        !
        ! this function brackets a minimum (algorithm from numerical recipes)
        !
        procedure(obj_func_type), intent(in), pointer :: obj_func
        real(rp), intent(in) :: kappa(:), lower, upper
        class(settings_type), intent(in) :: settings
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
        f_lower = obj_func(lower*kappa, error)
        call add_error_origin(error, error_obj_func, settings)
        if (error /= 0) return
        f_upper = obj_func(upper*kappa, error)
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
        n_c = n_b + golden_ratio*(n_b - n_a)
        f_c = obj_func(n_c*kappa, error)
        call add_error_origin(error, error_obj_func, settings)
        if (error /= 0) return

        ! continue looping until function at middle point is lower than at brackets
        iter = 0
        do while (f_c < f_b)
            ! compute value u by parabolic extrapolation
            tmp1 = (n_b - n_a)*(f_b - f_c)
            tmp2 = (n_b - n_c)*(f_b - f_a)
            val = tmp2 - tmp1
            denom = 2.0_rp*sign(max(abs(val), denom_floor), val)
            n_u = n_b - ((n_b - n_c)*tmp2 - (n_b - n_a)*tmp1)/denom

            ! maximum growth in parabolic fit
            n_u_lim = n_b + grow_limit*(n_c - n_b)

            ! check if u is between n_b and n_c
            if ((n_u - n_c)*(n_b - n_u) > 0.0) then
                ! evaluate function at n_u
                f_u = obj_func(n_u*kappa, error)
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
                n_u = n_c + golden_ratio*(n_c - n_b)
                f_u = obj_func(n_u*kappa, error)
                call add_error_origin(error, error_obj_func, settings)
                if (error /= 0) return
            ! limit parabolic fit to its maximum allowed value
            else if ((n_u - n_u_lim)*(n_u_lim - n_c) >= 0.0_rp) then
                n_u = n_u_lim
                f_u = obj_func(n_u*kappa, error)
                call add_error_origin(error, error_obj_func, settings)
                if (error /= 0) return
            ! parabolic fit is between n_c and its allowed limit
            else if ((n_u - n_u_lim)*(n_c - n_u) > 0.0_rp) then
                ! evaluate function at n_u
                f_u = obj_func(n_u*kappa, error)
                call add_error_origin(error, error_obj_func, settings)
                if (error /= 0) return

                if (f_u < f_c) then
                    n_b = n_c
                    n_c = n_u
                    n_u = n_c + golden_ratio*(n_c - n_b)
                    f_b = f_c
                    f_c = f_u
                    f_u = obj_func(n_u*kappa, error)
                    call add_error_origin(error, error_obj_func, settings)
                    if (error /= 0) return
                end if
            ! reject parabolic fit and use default step
            else
                n_u = n_c + golden_ratio*(n_c - n_b)
                f_u = obj_func(n_u*kappa, error)
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
                                  1, .true.)
                error = 1
                return
            end if
        end do

        ! check if miniumum is bracketed
        if (.not. (((f_b < f_c .and. f_b <= f_a) .or. (f_b < f_a .and. f_b <= f_c)) &
                   .and. ((n_a < n_b .and. n_b < n_c) .or. &
                          (n_c < n_b .and. n_b < n_a)))) then
            call settings%log("Line search did not find minimum", 1)
            error = 1
            return
        end if

        ! set new multiplier
        n_kappa = n_b

    end function bracket

    subroutine extend_symm_matrix(matrix, vector)
        !
        ! this subroutine extends a symmetric matrix
        !
        real(rp), allocatable, intent(inout) :: matrix(:, :)
        real(rp), intent(in) :: vector(:)

        real(rp), allocatable :: temp(:, :)

        ! copy data
        temp = matrix

        ! reallocate matrix
        deallocate(matrix)
        allocate(matrix(size(temp, 1) + 1, size(temp, 2) + 1))

        ! copy data back
        matrix(:size(temp, 1), :size(temp, 2)) = temp

        ! add new row and column
        matrix(:, size(temp, 2) + 1) = vector
        matrix(size(temp, 1) + 1, :) = vector

        ! deallocate temporary array
        deallocate(temp)

    end subroutine extend_symm_matrix

    subroutine add_column(matrix, new_col)
        !
        ! this subroutine adds a new column to an array
        !
        real(rp), allocatable, intent(inout) :: matrix(:, :)
        real(rp), intent(in) :: new_col(:)

        real(rp), allocatable :: temp(:, :)

        ! copy data
        temp = matrix

        ! reallocate matrix
        deallocate(matrix)
        allocate(matrix(size(temp, 1), size(temp, 2) + 1))

        ! copy data back
        matrix(:, :size(temp, 2)) = temp

        ! add new column
        matrix(:, size(temp, 2) + 1) = new_col

        ! deallocate temporary array
        deallocate(temp)

    end subroutine add_column

    subroutine symm_mat_min_eig(symm_matrix, lowest_eigval, lowest_eigvec, settings, &
                                error)
        !
        ! this function returns the lowest eigenvalue and corresponding eigenvector of
        ! a symmetric matrix
        !
        real(rp), intent(in) :: symm_matrix(:, :)
        class(settings_type), intent(in) :: settings
        real(rp), intent(out) :: lowest_eigval, lowest_eigvec(:)
        integer(ip), intent(out) :: error

        integer(ip) :: n, lwork, info
        real(rp), allocatable :: work(:), eigvals(:), eigvecs(:, :)
        character(300) :: msg
        external :: dsyev

        ! initialize error flag
        error = 0

        ! size of matrix
        n = size(symm_matrix, 1)

        ! copy matrix
        eigvecs = symm_matrix

        ! query optimal workspace size
        lwork = -1
        allocate(eigvals(n), work(1))
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
            write (msg, '(A, I0)') "Eigendecomposition failed: Error in DSYEV, "// &
                "info = ", info
            call settings%log(msg, 1, .true.)
            error = 1
            return
        end if

        ! get lowest eigenvalue and corresponding eigenvector
        lowest_eigval = eigvals(1)
        lowest_eigvec = eigvecs(:, 1)

        ! deallocate eigenvalues and eigenvectors
        deallocate(eigvals, eigvecs)

    end subroutine symm_mat_min_eig

    real(rp) function min_eigval(matrix, settings, error)
        !
        ! this function calculates the lowest eigenvalue of a symmetric matrix
        !
        real(rp), intent(in) :: matrix(:, :)
        class(settings_type), intent(in) :: settings
        integer(ip), intent(out) :: error

        real(rp), allocatable :: eigvals(:), temp(:, :), work(:)
        integer(ip) :: n, lwork, info
        character(300) :: msg
        external :: dsyev

        ! initialize error flag
        error = 0

        ! size of matrix
        n = size(matrix, 1)

        ! copy matrix to avoid modification of original matrix
        temp = matrix

        ! query optimal workspace size
        lwork = -1
        allocate(eigvals(n), work(1))
        call dsyev("N", "U", n, temp, n, eigvals, work, lwork, info)
        lwork = int(work(1), kind=ip)
        deallocate(work)
        allocate(work(lwork))

        ! compute eigenvalues
        call dsyev("N", "U", n, temp, n, eigvals, work, lwork, info)

        ! deallocate temporary and work array
        deallocate(temp, work)

        ! check for successful execution
        if (info /= 0) then
            write (msg, '(A, I0)') "Eigendecomposition failed: Error in DSYEV, "// &
                "info = ", info
            call settings%log(msg, 1, .true.)
            error = 1
            return
        end if

        ! get lowest eigenvalue
        min_eigval = eigvals(1)

        ! deallocate eigenvalues
        deallocate(eigvals)

    end function min_eigval

    subroutine init_rng(seed)
        !
        ! this subroutine initializes the random number generator
        !
        integer(ip), intent(in) :: seed

        integer, allocatable :: seed_arr(:)
        integer :: n

        call random_seed(size=n)
        allocate(seed_arr(n))
        seed_arr = int(seed, kind=4)
        call random_seed(put=seed_arr)
        deallocate(seed_arr)

    end subroutine init_rng

    function generate_trial_vectors(grad, grad_norm, h_diag, settings, error) &
        result(red_space_basis)
        !
        ! this function generates trial vectors
        !
        real(rp), intent(in) :: grad(:), grad_norm, h_diag(:)
        class(settings_type), intent(in) :: settings
        integer(ip), intent(out) :: error

        real(rp), allocatable :: red_space_basis(:, :)

        integer(ip) :: min_idx, n_vectors
        real(rp), external :: dnrm2

        ! initialize error flag
        error = 0

        ! get minimum Hessian diagonal element
        min_idx = minloc(h_diag, dim=1)

        ! add direction if minimum Hessian diagonal element is negative
        if (h_diag(min_idx) < 0.0_rp) then
            n_vectors = 2
            allocate(red_space_basis(size(grad), n_vectors + &
                     settings%n_random_trial_vectors))
            red_space_basis(:, 1) = grad/grad_norm
            red_space_basis(:, 2) = 0.0_rp
            red_space_basis(min_idx, 2) = 1.0_rp
            call gram_schmidt(red_space_basis(:, 2), &
                              reshape(red_space_basis(:, 1), &
                                      [size(red_space_basis, 1), 1]), &
                              settings, error)
            if (error /= 0) return
        else
            n_vectors = 1
            allocate(red_space_basis(size(grad), n_vectors + &
                     settings%n_random_trial_vectors))
            red_space_basis(:, 1) = grad/grad_norm
        end if

        call generate_random_trial_vectors(red_space_basis, settings, error)

    end function generate_trial_vectors

    subroutine generate_random_trial_vectors(red_space_basis, settings, error)
        !
        ! this subroutine generates random trial vectors
        !
        real(rp), intent(inout) :: red_space_basis(:, :)
        class(settings_type), intent(in) :: settings
        integer(ip), intent(out) :: error

        integer(ip) :: n_trial, i
        real(rp), parameter :: rnd_vector_min_norm = 1e-3_rp
        real(rp), external :: dnrm2

        ! initialize error flag
        error = 0

        ! number of trial vectors
        n_trial = size(red_space_basis, 2)

        do i = n_trial - settings%n_random_trial_vectors + 1, n_trial
            call random_number(red_space_basis(:, i))
            red_space_basis(:, i) = 2*red_space_basis(:, i) - 1
            do while (dnrm2(size(red_space_basis(:, i)), red_space_basis(:, i), 1) &
                      < rnd_vector_min_norm)
                call random_number(red_space_basis(:, i))
                red_space_basis(:, i) = 2*red_space_basis(:, i) - 1
            end do
            call gram_schmidt(red_space_basis(:, i), red_space_basis(:, :i - 1), &
                              settings, error)
            if (error /= 0) return
        end do

    end subroutine generate_random_trial_vectors

    subroutine gram_schmidt(vector, space, settings, error, lin_trans_vector, &
                            lin_trans_space)
        !
        ! this function orthonormalizes a vector with respect to a vector space
        ! this function can additionally also return a linear transformation of the 
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
        real(rp) :: dot, norm
        integer(ip) :: iter, i
        real(rp), parameter :: zero_thres = 1e-16_rp, orth_thres = 1e-14_rp
        real(rp), external :: ddot, dnrm2

        ! initialize error flag
        error = 0

        if (dnrm2(size(vector), vector, 1) < zero_thres) then
            call settings%log("Vector passed to Gram-Schmidt procedure is "// &
                "numerically zero.", 1, .true.)
            error = 1
            return
        else if (size(space, 2) > size(space, 1) - 1) then
            call settings%log("Number of vectors in Gram-Schmidt procedure larger "// &
                "than dimension of vector space.", 1, .true.)
            error = 1
            return
        end if
        
        ! allocate array for orthogonalities
        allocate(orth(size(space, 2)))

        iter = 0
        if (.not. (present(lin_trans_vector) .and. present(lin_trans_space))) then
            do while (.true.)
                do i = 1, size(space, 2)
                    vector = orthogonal_projection(vector, space(:, i))
                end do
                vector = vector / dnrm2(size(vector), vector, 1)

                call dgemv("T", size(vector), size(space, 2), 1.0_rp, space, &
                           size(vector), vector, 1, 0.0_rp, orth, 1)
                if (maxval(abs(orth)) < orth_thres) exit
                iter = iter + 1
                if (iter > 100) then
                    call settings%log("Maximum number of Gram-Schmidt iterations "// &
                                      "reached.", 1, .true.)
                    error = 1
                    return
                end if
            end do
        else
            do while (.true.)
                do i = 1, size(space, 2)
                    dot = ddot(size(vector), vector, 1, space(:, i), 1)
                    vector = vector - dot * space(:, i)
                    lin_trans_vector = lin_trans_vector - dot * lin_trans_space(:, i)
                end do
                norm = dnrm2(size(vector), vector, 1)
                vector = vector / norm
                lin_trans_vector = lin_trans_vector / norm

                call dgemv("T", size(vector), size(space, 2), 1.0_rp, space, &
                           size(vector), vector, 1, 0.0_rp, orth, 1)
                if (maxval(abs(orth)) < orth_thres) exit
                iter = iter + 1
                if (iter > 100) then
                    call settings%log("Maximum number of Gram-Schmidt iterations "// &
                                      "reached.", 1, .true.)
                    error = 1
                    return
                end if
            end do
        end if

        ! allocate array for orthogonalities
        deallocate(orth)

    end subroutine gram_schmidt

    subroutine init_solver_settings(self, error)
        !
        ! this subroutine sets the optional settings to their default values
        !
        class(solver_settings_type), intent(out) :: self
        integer(ip), intent(out) :: error

        ! initialize error flag
        error = 0

        ! ensure that class is actually solver_settings_type and not a subclass
        select type(settings => self)
        type is (solver_settings_type)
            settings = default_solver_settings
        class default
            call settings%log("Solver settings could not be initialized because "// &
                              "initialization routine received the wrong type. The "// &
                              "type solver_settings_type was likely subclassed "// &
                              "without providing an initialization routine.", 1, .true.)
            error = 1
        end select

    end subroutine init_solver_settings

    subroutine init_stability_settings(self, error)
        !
        ! this subroutine sets the optional settings to their default values
        !
        class(stability_settings_type), intent(out) :: self
        integer(ip), intent(out) :: error

        ! initialize error flag
        error = 0

        ! ensure that class is actually stability_settings_type and not a subclass
        select type(settings => self)
        type is (stability_settings_type)
            settings = default_stability_settings
        class default
            call settings%log("Stability settings could not be initialized because "// &
                              "initialization routine received the wrong type. The "// &
                              "type stability_settings_type was likely subclassed "// &
                              "without providing an initialization routine.", 1, .true.)
            error = 1
        end select

    end subroutine init_stability_settings

    subroutine level_shifted_diag_precond(vector, mu, h_diag, precond_vector)
        !
        ! this function defines the default level-shifted diagonal preconditioner
        !
        real(rp), intent(in) :: vector(:), mu, h_diag(:)
        real(rp), intent(out) :: precond_vector(:)
        
        ! construct level-shifted preconditioner
        precond_vector = h_diag - mu
        where (abs(precond_vector) < precond_floor)
            precond_vector = precond_floor
        end where

        ! precondition vector
        precond_vector = vector / precond_vector
        
    end subroutine level_shifted_diag_precond

    subroutine abs_diag_precond(vector, h_diag, precond_vector)
        !
        ! this function defines the default absolute diagonal preconditioner
        !
        real(rp), intent(in) :: vector(:), h_diag(:)
        real(rp), intent(out) :: precond_vector(:)

        ! construct positive-definite preconditioner
        precond_vector = max(abs(h_diag), precond_floor)
        
        ! precondition vector
        precond_vector = vector / precond_vector
        
    end subroutine abs_diag_precond

    function orthogonal_projection(vector, direction) result(complement)
        !
        ! this function removes a certain direction from a vector
        !
        real(rp), intent(in) :: vector(:), direction(:)

        real(rp) :: complement(size(vector))

        real(rp), external :: ddot

        complement = vector - direction * ddot(size(vector), vector, 1, direction, 1)

    end function orthogonal_projection

    subroutine jacobi_davidson_correction(hess_x_funptr, vector, solution, eigval, &
                                          corr_vector, hess_vector, settings, error)
        !
        ! this subroutine performs the Jacobi-Davidson correction but also returns the 
        ! Hessian linear transformation since this can be reused
        !
        procedure(hess_x_type), intent(in), pointer :: hess_x_funptr
        real(rp), intent(in) :: vector(:), solution(:), eigval
        real(rp), intent(out) :: corr_vector(:), hess_vector(:)
        class(settings_type), intent(in) :: settings
        integer(ip), intent(out) :: error

        ! initialize error flag
        error = 0
        
        ! project solution out of vector
        corr_vector = orthogonal_projection(vector, solution)

        ! get Hessian linear transformation of projected vector
        call hess_x_funptr(corr_vector, hess_vector, error)
        call add_error_origin(error, error_hess_x, settings)
        if (error /= 0) return

        ! finish construction of correction
        corr_vector = orthogonal_projection(hess_vector - eigval * corr_vector, &
                                            solution)
    
    end subroutine jacobi_davidson_correction

    subroutine minres(rhs, hess_x_funptr, solution, eigval, r_tol, vec, hvec, &
                      settings, error, guess, max_iter)
        !
        ! this function uses the minimum residual method to iteratively solve the 
        ! linear system for the Jacobi-Davidson correction equation, modified from 
        ! SciPy implementation
        !
        real(rp), intent(in) :: rhs(:), r_tol, solution(:), eigval
        procedure(hess_x_type), intent(in), pointer :: hess_x_funptr
        real(rp), intent(out) :: vec(:), hvec(:)
        class(settings_type), intent(in) :: settings
        integer(ip), intent(out) :: error
        real(rp), intent(in), optional :: guess(:)
        integer(ip), intent(in), optional :: max_iter

        integer(ip) :: n, max_iterations, iteration
        real(rp), parameter :: eps = epsilon(1.0_rp)
        real(rp) :: beta_start, beta, phi_bar, rhs1, old_beta, alfa, t_norm2, eps_ln, &
                    old_eps, cs, d_bar, sn, delta, g_bar, root, gamma, phi, g_max, &
                    g_min, tmp, rhs2, a_norm, vec_norm, qr_norm
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
                                            matvec, hvec, settings, error)
            if (error /= 0) return
            tot_hess_x = tot_hess_x + 1
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

        beta_start = dnrm2(n, r1, 1)

        ! check if starting guess already describes solution
        if (beta_start < numerical_zero) return

        ! solution must be zero vector if rhs vanishes
        if (dnrm2(n, rhs, 1) < numerical_zero) then
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
                                            settings, error)
            if (error /= 0) return
            tot_hess_x = tot_hess_x + 1

            ! get new trial vector
            if (iteration >= 2) y = y - (beta / old_beta) * r1
            alfa = ddot(n, v, 1, y, 1)
            y = y - (alfa / beta) * r2
            r1 = r2
            r2 = y
            old_beta = beta
            beta = dnrm2(n, r2, 1)
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
            root = sqrt(g_bar ** 2 + d_bar ** 2)

            ! compute next plane rotation by calculating cs and sn
            gamma = max(sqrt(g_bar ** 2 + beta ** 2), eps)
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
            vec_norm = dnrm2(n, vec, 1)
            qr_norm = phi_bar

            ! check if rhs and initial vector are eigenvectors
            if (stop_iteration) then
                call settings%log("MINRES: beta2 = 0. If M = I, b and x are "// &
                                  "eigenvectors.", 4)
                exit
            ! ||r||  / (||A|| ||x||)
            else if (vec_norm > 0.0_rp .and. a_norm > 0.0_rp .and. &
                qr_norm / (a_norm * vec_norm) <= r_tol) then
                    call settings%log("MINRES: A solution to Ax = b was found, "// &
                                      "given provided tolerance.", 4)
                exit
            ! ||Ar|| / (||A|| ||r||)
            else if (a_norm < numerical_zero .and. root / a_norm <= r_tol) then
                call settings%log("MINRES: A least-squares solution was found, "// &
                                  "given provided tolerance.", 4)
                exit
            ! check if reasonable accuracy is achieved with respect to machine precision
            else if (a_norm * vec_norm * eps >= beta_start) then
                call settings%log("MINRES: Reasonable accuracy achieved, given "// &
                                  "machine precision.", 4)
                exit
            ! compare estimate of condition number of matrix to machine precision
            else if (g_max / g_min >= 0.1 / eps) then
                call settings%log("MINRES: x has converged to an eigenvector.", 4)
                exit
            ! check if maximum number of iterations has been reached
            else if (iteration == max_iterations) then
                call settings%log("MINRES: The iteration limit was reached.", 1, .true.)
                error = 1
                return
            ! these tests ensure convergence is still achieved when r_tol 
            ! approaches machine precision
            else if (vec_norm > 0.0_rp .and. a_norm > 0.0_rp .and. &
                     1.0_rp + qr_norm / (a_norm * vec_norm) <= 1.0_rp) then
                call settings%log("MINRES: A solution to Ax = b was found, given "// &
                                  "provided tolerance.", 4)
                exit
            else if (a_norm < numerical_zero .and. 1.0_rp + root / a_norm <= 1.0_rp) &
                then
                call settings%log("MINRES: A least-squares solution was found, "// &
                                  "given provided tolerance.", 4)
                exit
            end if

            ! increment iteration counter
            iteration = iteration + 1

        end do

        deallocate(matvec, r1, r2, y, w, hw, w1, hw1, w2, hw2, v, hv)

    end subroutine minres

    subroutine print_results(self, iteration, func, grad_rms, level_shift, n_micro, &
                             imicro_jacobi_davidson, trust_radius, kappa_norm)
        !
        ! this function prints rows of the result table
        !
        class(solver_settings_type), intent(in) :: self
        integer(ip), intent(in) :: iteration
        real(rp), intent(in) :: func, grad_rms
        real(rp), intent(in), optional :: level_shift, trust_radius, kappa_norm
        integer(ip), intent(in), optional :: n_micro, imicro_jacobi_davidson

        character(11) :: iteration_str
        character(12) :: n_micro_str
        character(28) :: func_str
        character(14) :: grad_rms_str, trust_radius_str
        character(13) :: level_shift_str
        character(10) :: kappa_norm_str

        write(iteration_str, '(4X, I4, 3X)') iteration
        write(func_str, '(4X, 1PE21.14, 3X)') func
        write(grad_rms_str, '(2X, 1PE9.2, 3X)') grad_rms

        if (present(level_shift)) then
            write(level_shift_str, '(2X, 1PE9.2, 2X)') level_shift
        else
            write(level_shift_str, '(6X, "-", 6X)')
        end if

        if (present(n_micro)) then
            if (present(imicro_jacobi_davidson)) then
                write(n_micro_str, '(X, I3, X, "|", X, I3, 2X)') &
                    n_micro - imicro_jacobi_davidson + 1, imicro_jacobi_davidson - 1
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
            write(kappa_norm_str, '(X, 1PE9.2)') kappa_norm
        else
            write(kappa_norm_str, '(6X, "-", 3X)')
        end if

        call self%log(iteration_str // "|" // func_str // "|" // grad_rms_str &
                      // "|" // level_shift_str // "|" // n_micro_str // "|" // &
                      trust_radius_str // "|" // kappa_norm_str, 3)

    end subroutine print_results

    subroutine log(self, message, level, error)
        !
        ! this function performs logging
        !
        class(settings_type), intent(in) :: self
        character(*), intent(in) :: message
        integer(ip), intent(in) :: level
        logical, intent(in), optional :: error

        integer(ip), parameter :: max_length = 109
        integer(ip) :: i, out_unit
        character(:), dimension(:), allocatable :: substrings

        if (self%verbose >= level) then
            call split_string_by_space(message, max_length, substrings)
            if (associated(self%logger)) then
                do i = 1, size(substrings)
                    call self%logger(" " // substrings(i))
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
                    write(out_unit, '(A)') " " // substrings(i)
                end do
            end if
            deallocate(substrings)
        end if

    end subroutine log

    subroutine split_string_by_space(input, max_length, substrings)
        !
        ! this function splits a string by spaces to produce substrings of a maximum 
        ! length
        !
        character(*), intent(in) :: input
        integer(ip), intent(in) :: max_length
        character(:), intent(out), allocatable :: substrings(:)
    
        integer(ip) :: i, len_input, start_pos, end_pos, space_pos
        
        character(max_length) :: temp_string
        character(max_length), allocatable :: temp_substrings(:)

        len_input = len_trim(input)
        start_pos = 1
    
        do while (start_pos <= len_input)
            end_pos = min(start_pos + max_length - 1, len_input)
            space_pos = 0
    
            ! find last space before max_length characters
            do i = end_pos, start_pos, -1
                if (input(i:i) == ' ') then
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
            temp_string(1:space_pos-start_pos+1) = input(start_pos:space_pos)
            if (allocated(substrings)) then
                allocate(temp_substrings(size(substrings) + 1))
                temp_substrings(1:size(substrings)) = substrings
                substrings = temp_substrings
                deallocate(temp_substrings)
            else
                allocate(character(max_length) :: substrings(1))
            end if
            substrings(size(substrings)) = temp_string
            start_pos = space_pos + 1
        end do

    end subroutine split_string_by_space

    subroutine level_shifted_davidson(func, grad, grad_norm, h_diag, n_param, &
                                      obj_func, hess_x_funptr, settings, trust_radius, &
                                      solution, mu, imicro, imicro_jacobi_davidson, &
                                      jacobi_davidson_started, max_precision_reached, &
                                      error)
        !
        ! this subroutine performs level-shifted (Jacobi-)Davidson to solve the trust 
        ! region subproblem
        !
        real(rp), intent(in) :: func, grad(:), grad_norm, h_diag(:)
        integer(ip), intent(in) :: n_param
        procedure(obj_func_type), pointer, intent(in) :: obj_func
        procedure(hess_x_type), pointer, intent(in) :: hess_x_funptr
        type(solver_settings_type), intent(in) :: settings
        real(rp), intent(inout) :: trust_radius
        real(rp), intent(out) :: solution(:), mu
        integer(ip), intent(out) :: imicro, imicro_jacobi_davidson, error
        logical, intent(out) :: jacobi_davidson_started, max_precision_reached

        real(rp), allocatable :: red_space_basis(:, :), h_basis(:, :), aug_hess(:, :), &
                                 red_space_solution(:), red_hess_vec(:), basis_vec(:), &
                                 h_basis_vec(:), h_solution(:), residual(:), &
                                 solution_normalized(:), last_solution_normalized(:)
        integer(ip) :: n_trial, i, initial_imicro                          
        logical :: accept_step, micro_converged, func_evaluated, newton, bracketed
        real(rp) :: aug_hess_min_eigval, residual_norm, red_factor, &
                    initial_residual_norm, new_func, ratio, minres_tol
        real(rp), parameter :: newton_eigval_thresh = -1e-5_rp, &
                               level_shift_local_thres = 1e-12_rp, &
                               solution_overlap_thresh = 0.5_rp, &
                               residual_norm_floor = 1e-12_rp, &
                               residual_norm_max_red_factor = 0.8_rp, &
                               n_micro_davidson_factor = 0.8_rp, &
                               jacobi_davidson_ratio = 0.75_rp, &
                               jacobi_davidson_boundary_factor = 0.99_rp
        real(rp), external :: dnrm2, ddot

        ! initialize error flag
        error = 0

        ! generate trial vectors
        red_space_basis = generate_trial_vectors(grad, grad_norm, h_diag, settings, &
                                                 error)
        if (error /= 0) return

        ! number of trial vectors
        n_trial = size(red_space_basis, 2)

        ! increment number of Hessian linear transformations
        tot_hess_x = tot_hess_x + n_trial

        ! calculate linear transformations of basis vectors
        allocate(h_basis(n_param, n_trial))
        do i = 1, n_trial
            call hess_x_funptr(red_space_basis(:, i), h_basis(:, i), error)
            call add_error_origin(error, error_hess_x, settings)
            if (error /= 0) return
        end do

        ! construct augmented Hessian in reduced space
        allocate(aug_hess(n_trial + 1, n_trial + 1))
        aug_hess = 0.0_rp
        call dgemm("T", "N", n_trial, n_trial, n_param, 1.0_rp, red_space_basis, &
                   n_param, h_basis, n_param, 0.0_rp, aug_hess(2, 2), n_trial + 1)

        ! allocate space for reduced space solution and Hessian linear transformation
        ! of basis vector
        allocate(red_space_solution(n_trial), h_solution(n_param), basis_vec(n_param), &
                 h_basis_vec(n_param), last_solution_normalized(n_param))

        ! decrease trust radius until micro iterations converge and step is accepted
        last_solution_normalized = 0.0_rp
        accept_step = .false.
        do while (.not. accept_step)
            micro_converged = .false.
            func_evaluated = .false.

            jacobi_davidson_started = .false.
            do imicro = 1, settings%n_micro
                ! do a Newton step if the model is positive definite and the step is 
                ! within the trust region
                newton = .false.
                aug_hess_min_eigval = min_eigval(aug_hess(2:, 2:), settings, error)
                if (error /= 0) return
                if (aug_hess_min_eigval > newton_eigval_thresh) then
                    call newton_step(aug_hess, grad_norm, red_space_basis, &
                                     solution, red_space_solution, settings, error)
                    if (error /= 0) return
                    mu = 0.0_rp
                    if (dnrm2(n_param, solution, 1) < trust_radius) newton = .true.
                end if

                ! otherwise perform bisection to find the level shift
                if (.not. newton) then
                    call bisection(aug_hess, grad_norm, red_space_basis, trust_radius, &
                                   solution, red_space_solution, mu, bracketed, &
                                   settings, error)
                    if (error /= 0) return
                    if (.not. bracketed) exit
                end if

                ! calculate Hessian linear transformation of solution
                call dgemv("N", n_param, size(h_basis, 2), 1.0_rp, h_basis, n_param, &
                           red_space_solution, 1, 0.0_rp, h_solution, 1)

                ! calculate residual
                residual = grad + h_solution - mu*solution

                ! calculate residual norm
                residual_norm = dnrm2(n_param, residual, 1)

                ! determine reduction factor depending on whether local region is
                ! reached
                if (abs(mu) < level_shift_local_thres) then
                    red_factor = settings%local_red_factor
                else
                    red_factor = settings%global_red_factor
                end if

                ! get normalized solution vector
                solution_normalized = solution / dnrm2(n_param, solution, 1)

                ! reset initial residual norm if solution changes
                if (ddot(n_param, last_solution_normalized, 1, solution_normalized, &
                         1)**2 < solution_overlap_thresh) then
                    initial_imicro = imicro
                    initial_residual_norm = residual_norm
                end if

                ! check if micro iterations have converged
                if (residual_norm < max(red_factor * grad_norm, residual_norm_floor)) &
                    then
                    micro_converged = .true.
                    exit
                ! check if Jacobi-Davidson is used and has not been started
                else if (settings%jacobi_davidson .and. .not. jacobi_davidson_started) &
                    then
                    ! check residual has not decreased sufficiently or if maximum of 
                    ! Davidson iterations has been reached
                    if ((imicro - initial_imicro >= 10 .and. residual_norm > &
                         residual_norm_max_red_factor * initial_residual_norm) .or. &
                        imicro > n_micro_davidson_factor * settings%n_micro) then
                        ! evaluate function at approximate point
                        new_func = obj_func(solution, error)
                        call add_error_origin(error, error_obj_func, settings)
                        if (error /= 0) return

                        ! calculate ratio of evaluated function and predicted function
                        ratio = (new_func - func) / ddot(n_param, solution, 1, &
                                                         grad + 0.5_rp * h_solution, 1)

                        ! switch to Jacobi-Davidson only if current solution would lead 
                        ! to trust radius increase when the solution is already at the 
                        ! trust region boundary
                        if (settings%prefer_jacobi_davidson .or. &
                            (ratio > jacobi_davidson_ratio .and. &
                             dnrm2(n_param, solution, 1) > &
                             jacobi_davidson_boundary_factor * trust_radius)) then
                            jacobi_davidson_started = .true.
                            imicro_jacobi_davidson = imicro
                        ! decrease trust radius
                        else
                            func_evaluated = .true.
                            exit
                        end if
                    end if
                else if (imicro - initial_imicro >= 10 .and. residual_norm > &
                         residual_norm_max_red_factor * initial_residual_norm) then
                    exit
                end if

                ! save current solution
                last_solution_normalized = solution_normalized

                if (.not. jacobi_davidson_started) then
                    ! precondition residual
                    if (associated(settings%precond)) then
                        call settings%precond(residual, mu, basis_vec, error)
                        call add_error_origin(error, error_precond, settings)
                        if (error /= 0) return
                    else
                        call level_shifted_diag_precond(residual, mu, h_diag, basis_vec)
                    end if

                    ! orthonormalize to current orbital space to get new basis vector
                    call gram_schmidt(basis_vec, red_space_basis, settings, error)
                    if (error /= 0) return

                    ! add linear transformation of new basis vector
                    call hess_x_funptr(basis_vec, h_basis_vec, error)
                    call add_error_origin(error, error_hess_x, settings)
                    if (error /= 0) return

                    ! increment Hessian linear transformations
                    tot_hess_x = tot_hess_x + 1

                else
                    ! solve Jacobi-Davidson correction equations
                    minres_tol = 3.0_rp ** (-(imicro - imicro_jacobi_davidson))
                    call minres(-residual, hess_x_funptr, solution_normalized, mu, &
                                minres_tol, basis_vec, h_basis_vec, settings, error)
                    if (error /= 0) return

                    ! orthonormalize to current orbital space to get new basis vector
                    call gram_schmidt(basis_vec, red_space_basis, settings, error, &
                                      lin_trans_vector=h_basis_vec, &
                                      lin_trans_space=h_basis)
                    if (error /= 0) return

                    ! check if resulting linear transformation still respects Hessian 
                    ! symmetry which can happen due to numerical noise accumulation
                    if (abs(ddot(n_param, &
                                 red_space_basis(:, size(red_space_basis, 2)), 1, &
                                 h_basis_vec, 1) - &
                            ddot(n_param, basis_vec, 1, &
                                 h_basis(:, size(red_space_basis, 2)), 1)) > &
                                 hess_symm_thres) then
                        call hess_x_funptr(basis_vec, h_basis_vec, error)
                        call add_error_origin(error, error_hess_x, settings)
                        if (error /= 0) return
                    end if

                end if

                ! add new trial vector to orbital space
                call add_column(red_space_basis, basis_vec)

                ! add linear transformation of new basis vector
                call add_column(h_basis, h_basis_vec)

                ! construct new augmented Hessian
                allocate(red_hess_vec(size(red_space_basis, 2) + 1))
                red_hess_vec(1) = 0.0_rp
                call dgemv("T", n_param, size(red_space_basis, 2), 1.0_rp, &
                           red_space_basis, n_param, &
                           h_basis(:, size(red_space_basis, 2)), 1, 0.0_rp, &
                           red_hess_vec(2:), 1)
                call extend_symm_matrix(aug_hess, red_hess_vec)
                deallocate(red_hess_vec)

                ! reallocate reduced space solution
                deallocate(red_space_solution)
                allocate(red_space_solution(size(red_space_basis, 2)))
            end do

            if (.not. func_evaluated) then
                ! evaluate function at predicted point
                new_func = obj_func(solution, error)
                call add_error_origin(error, error_obj_func, settings)
                if (error /= 0) return

                ! calculate ratio of evaluated function and predicted function
                ratio = (new_func - func) / ddot(n_param, solution, 1, &
                                                 grad + 0.5_rp * h_solution, 1)
            end if

            ! decide whether to accept step and modify trust radius
            accept_step = accept_trust_region_step(solution, ratio, micro_converged, &
                                                   settings, trust_radius, &
                                                   max_precision_reached)
            if (max_precision_reached) exit
        end do

        ! deallocate quantities from microiterations
        deallocate(red_space_solution, aug_hess, red_space_basis, h_basis, h_solution, &
                   residual, basis_vec, h_basis_vec, solution_normalized, &
                   last_solution_normalized)

    end subroutine level_shifted_davidson

    subroutine truncated_conjugate_gradient(func, grad, h_diag, n_param, obj_func, &
                                            hess_x_funptr, settings, trust_radius, &
                                            solution, imicro, max_precision_reached, &
                                            error)
        !
        ! this subroutine performs truncated conjugate gradient to solve the trust 
        ! region subproblem
        !
        real(rp), intent(in) :: func, grad(:), h_diag(:)
        integer(ip), intent(in) :: n_param
        procedure(obj_func_type), pointer, intent(in) :: obj_func
        procedure(hess_x_type), pointer, intent(in) :: hess_x_funptr
        type(solver_settings_type), intent(in) :: settings
        real(rp), intent(inout) :: trust_radius
        real(rp), intent(out) :: solution(:)
        integer(ip), intent(out) :: imicro, error
        logical, intent(out) :: max_precision_reached

        logical :: accept_step, micro_converged
        real(rp) :: model_func, initial_residual_norm, curvature, step_size, &
                    solution_dot, solution_direction_dot, direction_dot, step_length, &
                    model_func_new, new_func, ratio
        real(rp), allocatable :: h_solution(:), residual(:), precond_residual(:), &
                                 direction(:), hess_direction(:), precond_solution(:), &
                                 precond_direction(:), solution_new(:), &
                                 h_solution_new(:), residual_new(:), &
                                 precond_residual_new(:), random_vector(:), &
                                 solutions(:, :), h_solutions(:, :)
        integer(ip) :: i
        real(rp), parameter :: random_noise_scale = 1e-4_rp, &
                               residual_norm_conv_red_factor = 1e-3_rp
        real(rp), external :: ddot, dnrm2

        ! initialize error flag
        error = 0

        ! allocate arrays
        allocate(h_solution(n_param), solutions(n_param, 1), h_solutions(n_param, 1), &
                 precond_residual(n_param))

        ! initialize solution
        solution = 0.0_rp
        h_solution = 0.0_rp
        solutions(:, 1) = solution
        h_solutions(:, 1) = h_solution
        model_func = 0.0_rp
        
        ! initialize residual and add random noise, residual should include h_solution 
        ! if not starting at zero
        residual = grad
        allocate(random_vector(n_param))
        random_vector = 0.0_rp
        do while (dnrm2(n_param, random_vector, 1) < numerical_zero)
            call random_number(random_vector)
            random_vector = 2 * random_vector - 1
        end do
        residual = residual + random_noise_scale * dnrm2(n_param, residual, 1) * &
                   random_vector / dnrm2(n_param, random_vector, 1)
        deallocate(random_vector)

        ! get initial residual norm
        initial_residual_norm = dnrm2(n_param, residual, 1)

        ! initialize preconditioned residual and direction,
        if (associated(settings%precond)) then
            call settings%precond(residual, 0.0_rp, precond_residual, error)
            call add_error_origin(error, error_precond, settings)
            if (error /= 0) return
        else
            call abs_diag_precond(residual, h_diag, precond_residual)
        end if
        direction = -precond_residual

        ! allocate arrays needed to propose a new step
        allocate(hess_direction(n_param), precond_solution(n_param), &
                 precond_direction(n_param), solution_new(n_param), &
                 h_solution_new(n_param), residual_new(n_param), &
                 precond_residual_new(n_param))
        
        ! start microiteration loop
        micro_converged = .false.
        do imicro = 1, settings%n_micro - 1
            ! get Hessian linear transformation of direction
            call hess_x_funptr(direction, hess_direction, error)
            call add_error_origin(error, error_hess_x, settings)
            if (error /= 0) return

            ! increment Hessian linear transformations
            tot_hess_x = tot_hess_x + 1

            ! calculate curvature
            curvature = ddot(n_param, direction, 1, hess_direction, 1)

            ! get step size along new direction
            step_size = ddot(n_param, residual, 1, precond_residual, 1) / curvature

            ! precondition current solution and direction
            if (associated(settings%precond)) then
                call settings%precond(solution, 0.0_rp, precond_solution, error)
                call add_error_origin(error, error_precond, settings)
                if (error /= 0) return
                call settings%precond(direction, 0.0_rp, precond_direction, error)
                call add_error_origin(error, error_precond, settings)
                if (error /= 0) return
            else
                call abs_diag_precond(solution, h_diag, precond_solution)
                call abs_diag_precond(direction, h_diag, precond_direction)
            end if

            ! calculate dot products
            solution_dot = ddot(n_param, solution, 1, precond_solution, 1)
            solution_direction_dot = ddot(n_param, solution, 1, precond_direction, 1)
            direction_dot = ddot(n_param, direction, 1, precond_direction, 1)

            ! calculate total step length
            step_length = solution_dot + 2 * step_size * solution_direction_dot + &
                          step_size ** 2 * direction_dot

            if (curvature < 0.0_rp .or. step_length >= trust_radius ** 2) then
                ! solve quadratic equation
                step_size = (-solution_direction_dot + &
                                sqrt(solution_direction_dot ** 2 + direction_dot * &
                                    (trust_radius ** 2 - solution_dot))) / &
                            direction_dot

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
            model_func_new = ddot(n_param, solution_new, 1, &
                                  grad + 0.5_rp * h_solution_new, 1)

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
            if (associated(settings%precond)) then
                call settings%precond(residual_new, 0.0_rp, precond_residual_new, error)
                call add_error_origin(error, error_precond, settings)
                if (error /= 0) return
            else
                call abs_diag_precond(residual_new, h_diag, precond_residual_new)
            end if

            ! check for linear or superlinear (in this case quadratic) convergence
            if (dnrm2(n_param, residual_new, 1) <= initial_residual_norm * &
                min(residual_norm_conv_red_factor, initial_residual_norm)) then
                micro_converged = .true.
                exit
            end if

            ! get new search direction
            direction = -precond_residual_new + &
                        ddot(n_param, residual_new, 1, precond_residual_new, 1) / &
                        ddot(n_param, residual, 1, precond_residual, 1) * direction
            
            ! save new model
            model_func = model_func_new
            residual = residual_new
            precond_residual = precond_residual_new
        end do

        ! deallocate arrays no longer needed
        deallocate(residual, precond_residual, solution_new, h_solution_new, &
                   residual_new, precond_residual_new)

        ! evaluate function at predicted point
        new_func = obj_func(solution, error)
        call add_error_origin(error, error_obj_func, settings)
        if (error /= 0) return

        if (abs(new_func - func) / max(abs(new_func), abs(func)) > numerical_zero) then
            ! calculate ratio of evaluated function and predicted function
            ratio = (new_func - func) / ddot(n_param, solution, 1, &
                                             grad + 0.5_rp * h_solution, 1)

            ! reduce trust region until step is accepted
            do while (.true.)
                ! decide whether to accept step and modify trust radius
                accept_step = accept_trust_region_step(solution, ratio, &
                                                       micro_converged, settings, &
                                                       trust_radius, &
                                                       max_precision_reached)
                if (accept_step .or. max_precision_reached) exit

                ! check if step exceeds new trust region boundary
                if (associated(settings%precond)) then
                    call settings%precond(solution, 0.0_rp, precond_solution, error)
                    if (error /= 0) return
                else
                    call abs_diag_precond(solution, h_diag, precond_solution)
                end if
                if (ddot(n_param, solution, 1, precond_solution, 1) > &
                    trust_radius ** 2) then
                    ! find step that exceeds trust region boundary
                    do i = 1, size(solutions, 2)
                        if (associated(settings%precond)) then
                            call settings%precond(solutions(:, i), 0.0_rp, &
                                                  precond_solution, error)
                            if (error /= 0) return
                        else
                            call abs_diag_precond(solutions(:, i), h_diag, &
                                                  precond_solution)
                        end if
                        if (ddot(n_param, solutions(:, i), 1, precond_solution, 1) > &
                            trust_radius ** 2) then
                            ! get previous step
                            solution = solutions(:, i - 1)
                            h_solution = h_solutions(:, i - 1)

                            ! get direction
                            direction = solutions(:, i) - solutions(:, i - 1)
                            hess_direction = h_solutions(:, i) - h_solutions(:, i - 1)

                            ! precondition current solution and direction
                            if (associated(settings%precond)) then
                                call settings%precond(solution, 0.0_rp, &
                                                      precond_solution, error)
                                call add_error_origin(error, error_precond, settings)
                                if (error /= 0) return
                                call settings%precond(direction, 0.0_rp, &
                                                      precond_direction, error)
                                call add_error_origin(error, error_precond, settings)
                                if (error /= 0) return
                            else
                                call abs_diag_precond(solution, h_diag, &
                                                      precond_solution)
                                call abs_diag_precond(direction, h_diag, &
                                                      precond_direction)
                            end if

                            ! calculate dot products
                            solution_dot = ddot(n_param, solution, 1, &
                                                precond_solution, 1)
                            solution_direction_dot = ddot(n_param, solution, 1, &
                                                          precond_direction, 1)
                            direction_dot = ddot(n_param, direction, 1, &
                                                 precond_direction, 1)

                            ! solve quadratic equation
                            step_size = (-solution_direction_dot + &
                                         sqrt(solution_direction_dot ** 2 + &
                                              direction_dot * &
                                              (trust_radius ** 2 - solution_dot))) / &
                                        direction_dot

                            ! get step to boundary
                            solution = solution + step_size * direction
                            h_solution = h_solution + step_size * hess_direction

                            ! evaluate function at predicted point
                            new_func = obj_func(solution, error)
                            call add_error_origin(error, error_obj_func, settings)
                            if (error /= 0) return

                            ! calculate ratio of evaluated function and predicted 
                            ! function
                            ratio = (new_func - func) / ddot(n_param, solution, 1, &
                                                             grad + 0.5_rp * &
                                                             h_solution, 1)

                            ! exit to retry whether step is accepted
                            micro_converged = .true.
                            exit
                        end if
                    end do
                end if
            end do
        else
            call settings%log("Function value barely changed. Convergence "// &
                              "criterion is not fulfilled but calculation should "// &
                              "be converged up to floating point precision.", 1, .true.)
            max_precision_reached = .true.
        end if

        ! deallocate arrays
        deallocate(h_solution, solutions, h_solutions, direction, hess_direction, &
                   precond_solution, precond_direction)

    end subroutine truncated_conjugate_gradient

    logical function accept_trust_region_step(solution, ratio, micro_converged, &
                                              settings, trust_radius, &
                                              max_precision_reached)
        !
        ! this function checks whether the trust region step is accepted and modified 
        ! the trust region accordingly
        !
        real(rp), intent(in) :: solution(:), ratio
        logical, intent(in) :: micro_converged
        type(solver_settings_type), intent(in) :: settings
        real(rp), intent(inout) :: trust_radius
        logical, intent(out) :: max_precision_reached

        ! default to maximum precision not yet reached
        max_precision_reached = .false.

        ! decrease trust radius if micro iterations are unable to converge, if function 
        ! value has not decreased or if individual orbitals change too much
        if (.not. micro_converged .or. ratio < 0.0_rp .or. any(abs(solution) > pi/4)) &
            then
            trust_radius = trust_radius_shrink_factor * trust_radius
            accept_trust_region_step = .false.
            if (trust_radius < numerical_zero) then
                call settings%log("Trust radius too small. Convergence criterion "// &
                                  "is not fulfilled but calculation should be "// &
                                  "converged up to floating point precision.", 1, &
                                  .true.)
                max_precision_reached = .true.
                return
            end if
        ! check if step is too long
        else if (ratio < trust_radius_shrink_ratio) then
            trust_radius = trust_radius_shrink_factor * trust_radius
            accept_trust_region_step = .true.
        ! check if quadratic approximation is valid
        else if (ratio < trust_radius_expand_ratio) then
            accept_trust_region_step = .true.
        ! check if step is potentially too short
        else
            trust_radius = trust_radius_expand_factor * trust_radius
            accept_trust_region_step = .true.
        end if

    end function accept_trust_region_step

    subroutine sanity_check(settings, n_param, grad, error)
        !
        ! this subroutine performs a sanity check for input parameters
        !
        type(solver_settings_type), intent(inout) :: settings
        integer(ip), intent(in) :: n_param
        real(rp), intent(in) :: grad(:)
        integer(ip), intent(out) :: error

        character(300) :: msg

        ! initialize error flag
        error = 0

        ! check that number of parameters is positive
        if (n_param < 1) then
            call settings%log("Number of parameters should be larger than 0.", 1, &
                              .true.)
            error = 1
            return
        end if

        ! check that number of random trial vectors is below number of parameters
        if (settings%davidson .and. settings%n_random_trial_vectors > n_param/2) then
            settings%n_random_trial_vectors = n_param/2
            write (msg, '(A, I0, A)') "Number of random trial vectors should be "// &
                "smaller than half the number of parameters. Setting to ", &
                settings%n_random_trial_vectors, "."
            call settings%log(msg, 2)
        end if

        ! sanity check for gradient size
        if (size(grad) /= n_param) then
            call settings%log("Size of gradient array returned by subroutine "// &
                              "update_orbs does not equal number of parameters.", 1, &
                              .true.)
            error = 1
            return
        end if

    end subroutine sanity_check

    subroutine add_error_origin(error_code, error_origin, settings)
        !
        ! this function modifies the error code by adding the error's origin if it is 
        ! not already added
        !
        class(settings_type), intent(in) :: settings
        integer(ip), intent(inout) :: error_code
        integer(ip), intent(in) :: error_origin

        if (error_code > 0) then
            if (error_code < 100) error_code = error_origin + error_code
        else if (error_code < 0) then
            call settings%log("Negative error code encountered.", 1, .true.)
            error_code = error_origin + 1
        end if

    end subroutine add_error_origin

end module opentrustregion
