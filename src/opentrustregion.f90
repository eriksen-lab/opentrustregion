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
    real(rp), parameter :: pi = 4.d0*atan(1.d0)

    ! define default optional arguments
    logical, parameter :: solver_stability_default = .false., &
                          solver_line_search_default = .false., &
                          solver_davidson_default = .true., &
                          solver_jacobi_davidson_default = .false., &
                          solver_prefer_jacobi_davidson_default = .false., &
                          stability_jacobi_davidson_default = .true.
    real(rp), parameter :: solver_conv_tol_default = 1d-5, &
                           solver_start_trust_radius_default = 0.4d0, &
                           solver_global_red_factor_default = 1d-3, &
                           solver_local_red_factor_default = 1d-4, &
                           stability_conv_tol_default = 1d-8
    integer(ip), parameter :: solver_n_random_trial_vectors_default = 1, &
                              solver_n_macro_default = 150, &
                              solver_n_micro_default = 50, &
                              solver_seed_default = 42, solver_verbose_default = 0, &
                              stability_n_random_trial_vectors_default = 20, &
                              stability_n_iter_default = 100, &
                              stability_verbose_default = 0

    ! define trust region parameters
    real(rp), parameter :: trust_radius_shrink_ratio = 0.25d0, &
                           trust_radius_expand_ratio = 0.75d0, &
                           trust_radius_shrink_factor = 0.7d0, &
                           trust_radius_expand_factor = 1.2d0

    ! derived type for solver settings
    type :: settings_type
        logical :: jacobi_davidson
        real(rp) :: conv_tol
        integer(ip) :: n_random_trial_vectors, verbose
        procedure(logger_type), pointer, nopass :: logger => null()
    contains
        procedure :: log
    end type

    type, extends(settings_type) :: solver_settings_type
        logical :: stability, line_search, davidson, prefer_jacobi_davidson
        real(rp) :: start_trust_radius, global_red_factor, local_red_factor
        integer(ip) :: n_macro, n_micro, seed
    contains
        procedure :: init_solver_settings, print_results
    end type

    type, extends(settings_type) :: stability_settings_type
        integer(ip) :: n_iter
    contains
        procedure :: init_stability_settings
    end type

    ! interface for setting default values
    interface set_default
        module procedure set_default_real, set_default_integer, set_default_logical
    end interface

    ! define global variables
    integer(ip) :: tot_orb_update = 0, tot_hess_x = 0

    ! interfaces for callback functions
    abstract interface
        function hess_x_type(x, error) result(hess_x)
            import :: rp

            real(rp), intent(in) :: x(:)
            logical, intent(out) :: error
            real(rp) :: hess_x(size(x))
        end function hess_x_type
    end interface

    abstract interface
        subroutine update_orbs_type(kappa, func, grad, h_diag, hess_x_funptr, error)
            import :: rp, hess_x_type

            real(rp), intent(in) :: kappa(:)
            real(rp), intent(out) :: func, grad(:), h_diag(:)
            procedure(hess_x_type), intent(out), pointer :: hess_x_funptr
            logical, intent(out) :: error
        end subroutine update_orbs_type
    end interface

    abstract interface
        function obj_func_type(kappa, error) result(func)
            import :: rp

            real(rp), intent(in) :: kappa(:)
            logical, intent(out) :: error
            real(rp) :: func
        end function obj_func_type
    end interface

    abstract interface
        function precond_type(residual, mu, error) result(precond_residual)
            import :: rp

            real(rp), intent(in) :: residual(:), mu
            logical, intent(out) :: error
            real(rp) :: precond_residual(size(residual))
        end function precond_type
    end interface

    abstract interface
        function conv_check_type(error) result(converged)
            logical, intent(out) :: error
            logical :: converged
        end function conv_check_type
    end interface

    abstract interface
        subroutine logger_type(message)
            character(*), intent(in) :: message
        end subroutine logger_type
    end interface

contains

    subroutine solver(update_orbs, obj_func, n_param, error, precond, conv_check, &
                      stability, line_search, davidson, jacobi_davidson, &
                      prefer_jacobi_davidson, conv_tol, n_random_trial_vectors, &
                      start_trust_radius, n_macro, n_micro, global_red_factor, &
                      local_red_factor, seed, verbose, logger)
        !
        ! this subroutine is the main solver for orbital optimization
        !
        procedure(update_orbs_type), intent(in), pointer :: update_orbs
        procedure(obj_func_type), intent(in), pointer :: obj_func
        integer(ip), intent(in) :: n_param
        procedure(precond_type), intent(in), pointer, optional :: precond
        procedure(conv_check_type), intent(in), pointer, optional :: conv_check
        logical, intent(out) :: error
        logical, intent(in), optional :: stability, line_search, davidson, &
                                         jacobi_davidson, prefer_jacobi_davidson
        real(rp), intent(in), optional :: conv_tol, start_trust_radius, &
                                          global_red_factor, local_red_factor
        integer(ip), intent(in), optional :: n_random_trial_vectors, n_macro, n_micro, &
                                             seed, verbose
        procedure(logger_type), intent(in), pointer, optional :: logger

        type(solver_settings_type) :: settings
        real(rp) :: trust_radius, func, grad_norm, grad_rms, mu, new_func, n_kappa, &
                    kappa_norm
        real(rp), dimension(n_param) :: kappa, grad, h_diag, solution
        logical :: max_precision_reached, macro_converged, stable, &
                   jacobi_davidson_started, use_precond, conv_check_passed
        integer(ip) :: imacro, imicro, imicro_jacobi_davidson, i
        character(300) :: msg
        integer(ip), parameter :: stability_n_points = 21
        procedure(hess_x_type), pointer :: hess_x_funptr
        real(rp), external :: dnrm2, ddot

        ! initialize error flag
        error = .false.

        ! initialize maximum precision convergence
        max_precision_reached = .false.

        ! initialize macroiteration convergence
        macro_converged = .false.

        ! initialize stabilty boolean
        stable = .true.

        ! initialize settings
        call settings%init_solver_settings(stability, line_search, davidson, &
                                           jacobi_davidson, prefer_jacobi_davidson, &
                                           conv_tol, n_random_trial_vectors, &
                                           start_trust_radius, n_macro, n_micro, &
                                           global_red_factor, local_red_factor, seed, &
                                           verbose, logger)

        ! initialize random number generator
        call init_rng(settings%seed)

        ! initialize orbital rotation matrix
        kappa = 0.d0

        ! initialize starting trust radius
        trust_radius = settings%start_trust_radius

        ! check if preconditioner is passed
        if (present(precond)) then
            use_precond = associated(precond)
        else
            use_precond = .false.
        end if

        ! print header
        call settings%log(repeat("-", 109), 3)
        call settings%log(" Iteration |     Objective function     | Gradient RMS |"// &
                          " Level shift |   Micro    | Trust radius | Step size ", 3)
        call settings%log("           |                            |              |"// &
                          "             | iterations |              |           ", 3)
        call settings%log(repeat("-", 109), 3)

        do imacro = 1, settings%n_macro
            if (.not. max_precision_reached) then
                ! calculate cost function, gradient and Hessian diagonal
                call update_orbs(kappa, func, grad, h_diag, hess_x_funptr, error)
                if (error) return

                ! perform sanity check
                if (imacro == 1) then
                    call sanity_check(settings, n_param, grad, error)
                    if (error) return
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
                        if (use_precond) then
                            kappa_norm = sqrt(ddot(n_param, kappa, 1, &
                                                   precond(kappa, 0.d0, error), 1))
                            if (error) return
                            
                        else
                            kappa_norm = sqrt(ddot(n_param, kappa, 1, &
                                                   abs_diag_precond(kappa, h_diag), 1))
                        end if
                        mu = 0.d0
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
            if (present(conv_check)) then
                if (associated(conv_check)) then
                    conv_check_passed = conv_check(error)
                    if (error) return
                else
                    conv_check_passed = .false.
                end if
            else
                conv_check_passed = .false.
            end if
            if (grad_rms < settings%conv_tol .or. max_precision_reached .or. &
                conv_check_passed) then
                ! always perform stability check if starting at stationary point
                if (settings%stability .or. imacro == 1) then
                    call stability_check(h_diag, hess_x_funptr, stable, kappa, error, &
                                         precond=precond, verbose=settings%verbose, &
                                         logger=logger)
                    if (error) return
                    if (.not. stable) then
                        ! logarithmic line search
                        do i = 1, stability_n_points
                            n_kappa = 10.d0**(-(i - 1)/ &
                                              real(stability_n_points - 1, rp)*10.d0)
                            new_func = obj_func(n_kappa*kappa, error)
                            if (error) return
                            if (new_func < func) then
                                kappa = n_kappa*kappa
                                exit
                            end if
                        end do
                        if (new_func >= func) then
                            call settings%log("Line search was unable to find "// &
                                              "lower objective function along "// &
                                              "unstable mode.", 1, .true.)
                            error = .true.
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
                                            use_precond, precond, trust_radius, &
                                            solution, mu, imicro, &
                                            imicro_jacobi_davidson, &
                                            jacobi_davidson_started, &
                                            max_precision_reached, &
                                            error)
                if (error) return
            else
                ! solve trust region subproblem with truncated conjugate gradient
                call truncated_conjugate_gradient(func, grad, h_diag, n_param, &
                                                  obj_func, hess_x_funptr, &
                                                  use_precond, precond, settings, &
                                                  trust_radius, solution, imicro, &
                                                  max_precision_reached, error)
                if (error) return
            end if

            ! perform line search
            if (max_precision_reached) then
                n_kappa = 0.d0
            else if (settings%line_search) then
                n_kappa = bracket(obj_func, solution, 0.d0, 1.d0, settings, error)
                if (error) return
            else
                n_kappa = 1.d0
            end if

            ! set orbital rotation
            kappa = n_kappa*solution

            ! flush output
            flush(stdout)
            flush(stderr)

        end do

        ! increment total number of orbital updates
        tot_orb_update = tot_orb_update + imacro

        ! stop if no convergence
        if (.not. macro_converged) then
            call settings%log("Orbital optimization has not converged!", 1)
            error = .true.
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

    subroutine stability_check(h_diag, hess_x_funptr, stable, kappa, error, precond, &
                               jacobi_davidson, conv_tol, n_random_trial_vectors, &
                               n_iter, verbose, logger)
        !
        ! this subroutine performs a stability check
        !
        real(rp), intent(in) :: h_diag(:)
        procedure(hess_x_type), intent(in), pointer :: hess_x_funptr
        logical, intent(out) :: stable, error
        real(rp), intent(out) :: kappa(:)
        procedure(precond_type), intent(in), pointer, optional :: precond
        logical, intent(in), optional :: jacobi_davidson
        real(rp), intent(in), optional :: conv_tol
        integer(ip), intent(in), optional :: n_random_trial_vectors, n_iter, verbose
        procedure(logger_type), intent(in), pointer, optional :: logger

        type(stability_settings_type) :: settings
        integer(ip) :: n_param, ntrial, i, iter
        real(rp), dimension(size(h_diag)) :: solution, h_solution, residual, &
                                             basis_vec, h_basis_vec
        real(rp) :: eigval, minres_tol, stability_rms
        real(rp), allocatable :: red_space_basis(:, :), h_basis(:, :), &
                                 red_space_hess(:, :), red_space_solution(:), &
                                 red_space_hess_vec(:)
        logical :: use_precond
        character(300) :: msg
        real(rp), external :: dnrm2, ddot
        external :: dgemm, dgemv

        ! initialize error flag
        error = .false.

        ! initialize settings
        call settings%init_stability_settings(jacobi_davidson, conv_tol, &
                                              n_random_trial_vectors, n_iter, verbose, &
                                              logger)

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

        ! check if preconditioner is passed
        if (present(precond)) then
            use_precond = associated(precond)
        else
            use_precond = .false.
        end if

        ! generate trial vectors
        allocate (red_space_basis(size(h_diag), 1 + settings%n_random_trial_vectors))
        red_space_basis(:, 1) = 0.d0
        red_space_basis(minloc(h_diag), 1) = 1.d0
        call generate_random_trial_vectors(red_space_basis, settings, error)
        if (error) return

        ! number of trial vectors
        ntrial = size(red_space_basis, 2)

        ! calculate linear transformations of basis vectors
        allocate (h_basis(n_param, ntrial))
        do i = 1, ntrial
            h_basis(:, i) = hess_x_funptr(red_space_basis(:, i), error)
            if (error) return
        end do

        ! construct augmented Hessian in reduced space
        allocate (red_space_hess(ntrial, ntrial))
        call dgemm("T", "N", ntrial, ntrial, n_param, 1.d0, red_space_basis, &
                   n_param, h_basis, n_param, 0.d0, red_space_hess, ntrial)

        ! allocate space for reduced space solution
        allocate (red_space_solution(size(red_space_basis, 2)))

        ! loop over iterations
        do iter = 1, settings%n_iter
            ! solve reduced space problem
            call symm_mat_min_eig(red_space_hess, eigval, red_space_solution, &
                                  settings, error)
            if (error) return

            ! get full space solution
            call dgemv("N", size(red_space_basis, 1), size(red_space_basis, 2), 1.d0, &
                       red_space_basis, size(red_space_basis, 1), red_space_solution, &
                       1, 0.d0, solution, 1)

            ! calculate Hessian linear transformation of solution
            call dgemv("N", n_param, size(h_basis, 2), 1.d0, h_basis, n_param, &
                       red_space_solution, 1, 0.d0, h_solution, 1)

            ! calculate residual
            residual = h_solution - eigval*solution

            ! check convergence
            stability_rms = dnrm2(n_param, residual, 1) / sqrt(real(n_param, kind=rp))
            if (stability_rms < settings%conv_tol) exit

            if (.not. settings%jacobi_davidson .or. iter <= 30) then
                ! precondition residual
                if (use_precond) then
                    basis_vec = precond(residual, 0.d0, error)
                    if (error) return
                else
                    basis_vec = level_shifted_diag_precond(residual, 0.d0, h_diag)
                end if

                ! orthonormalize to current orbital space to get new basis vector
                call gram_schmidt(basis_vec, red_space_basis, settings, error)
                if (error) return

                ! add linear transformation of new basis vector
                h_basis_vec = hess_x_funptr(basis_vec, error)
                if (error) return

                ! increment Hessian linear transformations
                tot_hess_x = tot_hess_x + 1

            else
                ! solve Jacobi-Davidson correction equations
                minres_tol = 3.d0 ** (-(iter - 31))
                call minres(-residual, hess_x_funptr, solution, eigval, minres_tol, &
                            basis_vec, h_basis_vec, settings, error)
                if (error) return

                ! orthonormalize to current orbital space to get new basis 
                ! vector
                call gram_schmidt(basis_vec, red_space_basis, settings, error, &
                                  lin_trans_vector=h_basis_vec, lin_trans_space=h_basis)
                if (error) return

                ! check if resulting linear transformation still respects Hessian 
                ! symmetry which can happen due to numerical noise accumulation
                if (abs(ddot(n_param, red_space_basis(:, size(red_space_basis, 2)), 1, &
                             h_basis_vec, 1) - &
                        ddot(n_param, basis_vec, 1, &
                             h_basis(:, size(red_space_basis, 2)), 1)) > 1d-12) then
                    h_basis_vec = hess_x_funptr(basis_vec, error)
                    if (error) return
                end if
                
            end if

            ! add new trial vector to orbital space
            call add_column(red_space_basis, basis_vec)

            ! add linear transformation of new basis vector
            call add_column(h_basis, h_basis_vec)

            ! construct new reduced space Hessian
            allocate (red_space_hess_vec(size(red_space_basis, 2)))
            call dgemv("T", n_param, size(red_space_basis, 2), 1.d0, &
                       red_space_basis, n_param, &
                       h_basis(:, size(red_space_basis, 2)), 1, 0.d0, &
                       red_space_hess_vec, 1)
            call extend_symm_matrix(red_space_hess, red_space_hess_vec)
            deallocate (red_space_hess_vec)

            ! reallocate reduced space solution
            deallocate (red_space_solution)
            allocate (red_space_solution(size(red_space_basis, 2)))

        end do

        ! check if stability check has converged
        if (stability_rms >= settings%conv_tol) &
            call settings%log("Stability check has not converged in the given "// &
                              "number of iterations.", 1, .true.)

        ! increment total number of Hessian linear transformations
        tot_hess_x = tot_hess_x + size(red_space_basis, 2)

        ! deallocate quantities from Davidson iterations
        deallocate (red_space_solution)
        deallocate (red_space_hess)
        deallocate (h_basis)
        deallocate (red_space_basis)

        ! determine if saddle point
        stable = eigval > -1.d-2
        if (stable) then
            kappa = 0.d0
        else
            kappa = solution
            write (msg, '(A, F0.4)') "Solution not stable. Lowest eigenvalue: ", eigval
            call settings%log(msg, 1, .true.)
        end if

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
        logical, intent(out) :: error

        integer(ip) :: nred, lwork, info, ipiv(size(red_space_basis, 2))
        real(rp) :: red_hess(size(red_space_basis, 2), size(red_space_basis, 2))
        real(rp), allocatable :: work(:)
        character(300) :: msg
        external :: dsysv, dgemv

        ! initialize error flag
        error = .false.

        ! reduced space size
        nred = size(red_space_basis, 2)

        ! reduced space Hessian
        red_hess = aug_hess(2:, 2:)

        ! set gradient
        red_space_solution = 0.d0
        red_space_solution(1) = -grad_norm

        ! query optimal workspace size
        lwork = -1
        allocate (work(1))
        call dsysv("U", nred, 1, red_hess, nred, ipiv, red_space_solution, nred, work, &
                   lwork, info)
        lwork = int(work(1))
        deallocate (work)
        allocate (work(lwork))

        ! solve linear system
        call dsysv("U", nred, 1, red_hess, nred, ipiv, red_space_solution, nred, work, &
                   lwork, info)

        ! deallocate work array
        deallocate (work)

        ! check for errors
        if (info /= 0) then
            write (msg, '(A, I0)') "Linear solver failed: Error in DSYSV, info = ", info
            call settings%log(msg, 1, .true.)
            error = .true.
            return
        end if

        ! get solution in full space
        call dgemv("N", size(red_space_basis, 1), size(red_space_basis, 2), 1.d0, &
                   red_space_basis, size(red_space_basis, 1), red_space_solution, 1, &
                   0.d0, solution, 1)

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
        logical, intent(out) :: bracketed, error

        real(rp) :: lower_alpha, middle_alpha, upper_alpha, lower_trust_dist, &
                    middle_trust_dist, upper_trust_dist
        integer(ip) :: iter
        real(rp), external :: dnrm2

        ! initialize error flag
        error = .false.

        ! initialize bracketing flag
        bracketed = .false.

        ! lower and upper bracket for alpha
        lower_alpha = 1.d-4
        upper_alpha = 1.d6

        ! solve reduced space problem with scaled gradient
        call get_ah_lowest_eigenvec(lower_alpha)
        if (error) return
        lower_trust_dist = dnrm2(size(solution), solution, 1) - trust_radius
        call get_ah_lowest_eigenvec(upper_alpha)
        if (error) return
        upper_trust_dist = dnrm2(size(solution), solution, 1) - trust_radius

        ! check if trust region is within bracketing range
        if ((lower_trust_dist*upper_trust_dist) > 0.d0) then
            solution = 0.d0
            red_space_solution = 0.d0
            mu = 0.d0
            return
        end if

        ! get middle alpha
        middle_alpha = sqrt(upper_alpha*lower_alpha)
        call get_ah_lowest_eigenvec(middle_alpha)
        if (error) return
        middle_trust_dist = dnrm2(size(solution), solution, 1) - trust_radius

        ! perform bisection to find root, converge to relative threshold to avoid 
        ! precision issues
        iter = 0
        do while (upper_alpha - lower_alpha > 1.d-12 * upper_alpha)
            ! targeted trust radius is in upper bracket
            if (lower_trust_dist*middle_trust_dist > 0.d0) then
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
            if (error) return
            middle_trust_dist = dnrm2(size(solution), solution, 1) - trust_radius
            ! check if maximum number of iterations is reached
            iter = iter + 1
            if (iter > 100) then
                call settings%log("Maximum number of bisection iterations reached.", &
                                  1, .true.)
                error = .true.
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

            real(rp) :: eigvec(size(aug_hess, 1))
            external :: dgemv

            ! finish construction of augmented Hessian
            aug_hess(1, 2) = alpha*grad_norm
            aug_hess(2, 1) = alpha*grad_norm

            ! perform eigendecomposition and get lowest eigenvalue and corresponding
            ! eigenvector
            call symm_mat_min_eig(aug_hess, mu, eigvec, settings, error)
            if (error) return

            ! check if eigenvector has level-shift component
            if (abs(eigvec(1)) == 0.d0) then
                call settings%log("Trial subspace too small. Increase "// &
                    "n_random_trial_vectors.", 1, .true.)
                error = .true.
                return
            end if

            ! scale eigenvector such that first element is equal to one and divide by
            ! alpha to get solution in reduced space
            red_space_solution = eigvec(2:)/eigvec(1)/alpha

            ! get solution in full space
            call dgemv("N", size(red_space_basis, 1), size(red_space_basis, 2), 1.d0, &
                       red_space_basis, size(red_space_basis, 1), red_space_solution, &
                       1, 0.d0, solution, 1)

        end subroutine get_ah_lowest_eigenvec

    end subroutine bisection

    function bracket(obj_func, kappa, lower, upper, settings, error) result(n_kappa)
        !
        ! this function brackets a minimum (algorithm from numerical recipes)
        !
        procedure(obj_func_type), intent(in), pointer :: obj_func
        real(rp), intent(in) :: kappa(:), lower, upper
        class(settings_type), intent(in) :: settings
        logical, intent(out) :: error

        real(rp) :: n_kappa, f_upper, f_lower, n_a, n_b, n_c, n_u, f_a, f_b, f_c, f_u, &
                    n_u_lim, tmp1, tmp2, val, denom
        integer(ip) :: iter
        real(rp), parameter :: golden_ratio = (1.d0 + sqrt(5.d0))/2.d0, &
                               grow_limit = 110.d0

        ! initialize error flag
        error = .false.

        ! evaluate function at upper and lower bounds
        f_lower = obj_func(lower*kappa, error)
        if (error) return
        f_upper = obj_func(upper*kappa, error)
        if (error) return

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
        if (error) return

        ! continue looping until function at middle point is lower than at brackets
        iter = 0
        do while (f_c < f_b)
            ! compute value u by parabolic extrapolation
            tmp1 = (n_b - n_a)*(f_b - f_c)
            tmp2 = (n_b - n_c)*(f_b - f_a)
            val = tmp2 - tmp1
            denom = 2.d0*sign(max(abs(val), 1.d-21), val)
            n_u = n_b - ((n_b - n_c)*tmp2 - (n_b - n_a)*tmp1)/denom

            ! maximum growth in parabolic fit
            n_u_lim = n_b + grow_limit*(n_c - n_b)

            ! check if u is between n_b and n_c
            if ((n_u - n_c)*(n_b - n_u) > 0.0) then
                ! evaluate function at n_u
                f_u = obj_func(n_u*kappa, error)
                if (error) return

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
                if (error) return
            ! limit parabolic fit to its maximum allowed value
            else if ((n_u - n_u_lim)*(n_u_lim - n_c) >= 0.d0) then
                n_u = n_u_lim
                f_u = obj_func(n_u*kappa, error)
                if (error) return
            ! parabolic fit is between n_c and its allowed limit
            else if ((n_u - n_u_lim)*(n_c - n_u) > 0.d0) then
                ! evaluate function at n_u
                f_u = obj_func(n_u*kappa, error)
                if (error) return

                if (f_u < f_c) then
                    n_b = n_c
                    n_c = n_u
                    n_u = n_c + golden_ratio*(n_c - n_b)
                    f_b = f_c
                    f_c = f_u
                    f_u = obj_func(n_u*kappa, error)
                    if (error) return
                end if
            ! reject parabolic fit and use default step
            else
                n_u = n_c + golden_ratio*(n_c - n_b)
                f_u = obj_func(n_u*kappa, error)
                if (error) return
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
                error = .true.
                return
            end if
        end do

        ! check if miniumum is bracketed
        if (.not. (((f_b < f_c .and. f_b <= f_a) .or. (f_b < f_a .and. f_b <= f_c)) &
                   .and. ((n_a < n_b .and. n_b < n_c) .or. &
                          (n_c < n_b .and. n_b < n_a)))) then
            call settings%log("Line search did not find minimum", 1)
            error = .true.
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

        ! allocate temporary array
        allocate (temp(size(matrix, 1), size(matrix, 2)))

        ! copy data
        temp = matrix

        ! reallocate matrix
        deallocate (matrix)
        allocate (matrix(size(temp, 1) + 1, size(temp, 2) + 1))

        ! copy data back
        matrix(:size(temp, 1), :size(temp, 2)) = temp

        ! add new row and column
        matrix(:, size(temp, 2) + 1) = vector
        matrix(size(temp, 1) + 1, :) = vector

        ! deallocate temporary array
        deallocate (temp)

    end subroutine extend_symm_matrix

    subroutine add_column(matrix, new_col)
        !
        ! this subroutine adds a new column to an array
        !
        real(rp), allocatable, intent(inout) :: matrix(:, :)
        real(rp), intent(in) :: new_col(:)

        real(rp), allocatable :: temp(:, :)

        ! allocate temporary array
        allocate (temp(size(matrix, 1), size(matrix, 2)))

        ! copy data
        temp = matrix

        ! reallocate matrix
        deallocate (matrix)
        allocate (matrix(size(temp, 1), size(temp, 2) + 1))

        ! copy data back
        matrix(:, :size(temp, 2)) = temp

        ! add new column
        matrix(:, size(temp, 2) + 1) = new_col

        ! deallocate temporary array
        deallocate (temp)

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
        logical, intent(out) :: error

        integer(ip) :: n, lwork, info
        real(rp), allocatable :: work(:)
        real(rp) :: eigvals(size(symm_matrix, 1)), &
                    eigvecs(size(symm_matrix, 1), size(symm_matrix, 2))
        character(300) :: msg
        external :: dsyev

        ! initialize error flag
        error = .false.

        ! size of matrix
        n = size(symm_matrix, 1)

        ! copy matrix
        eigvecs = symm_matrix

        ! query optimal workspace size
        lwork = -1
        allocate (work(1))
        call dsyev("V", "U", n, eigvecs, n, eigvals, work, lwork, info)
        lwork = int(work(1))
        deallocate (work)
        allocate (work(lwork))

        ! perform eigendecomposition
        call dsyev("V", "U", n, eigvecs, n, eigvals, work, lwork, info)

        ! deallocate work array
        deallocate (work)

        ! check for successful execution
        if (info /= 0) then
            write (msg, '(A, I0)') "Eigendecomposition failed: Error in DSYEV, "// &
                "info = ", info
            call settings%log(msg, 1, .true.)
            error = .true.
            return
        end if

        ! get lowest eigenvalue and corresponding eigenvector
        lowest_eigval = eigvals(1)
        lowest_eigvec = eigvecs(:, 1)

    end subroutine symm_mat_min_eig

    real(rp) function min_eigval(matrix, settings, error)
        !
        ! this function calculates the lowest eigenvalue of a symmetric matrix
        !
        real(rp), intent(in) :: matrix(:, :)
        class(settings_type), intent(in) :: settings
        logical, intent(out) :: error

        real(rp) :: eigvals(size(matrix, 1)), temp(size(matrix, 1), size(matrix, 2))
        integer(ip) :: n, lwork, info
        real(rp), allocatable :: work(:)
        character(300) :: msg
        external :: dsyev

        ! initialize error flag
        error = .false.

        ! size of matrix
        n = size(matrix, 1)

        ! copy matrix to avoid modification of original matrix
        temp = matrix

        ! query optimal workspace size
        lwork = -1
        allocate (work(1))
        call dsyev("N", "U", n, temp, n, eigvals, work, lwork, info)
        lwork = int(work(1))
        deallocate (work)
        allocate (work(lwork))

        ! compute eigenvalues
        call dsyev("N", "U", n, temp, n, eigvals, work, lwork, info)

        ! deallocate work array
        deallocate (work)

        ! check for successful execution
        if (info /= 0) then
            write (msg, '(A, I0)') "Eigendecomposition failed: Error in DSYEV, "// &
                "info = ", info
            call settings%log(msg, 1, .true.)
            error = .true.
            return
        end if

        ! get lowest eigenvalue
        min_eigval = eigvals(1)

    end function min_eigval

    subroutine init_rng(seed)
        !
        ! this subroutine initializes the random number generator
        !
        integer(ip), intent(in) :: seed

        integer, allocatable :: seed_arr(:)
        integer :: n

        call random_seed(size=n)
        allocate (seed_arr(n))
        seed_arr = int(seed, kind=4)
        call random_seed(put=seed_arr)
        deallocate (seed_arr)

    end subroutine init_rng

    function generate_trial_vectors(grad, grad_norm, h_diag, settings, error) &
        result(red_space_basis)
        !
        ! this function generates trial vectors
        !
        real(rp), intent(in) :: grad(:), grad_norm, h_diag(:)
        class(settings_type), intent(in) :: settings
        logical, intent(out) :: error

        real(rp), allocatable :: red_space_basis(:, :)

        integer(ip) :: min_idx, n_vectors
        real(rp) :: trial(size(grad))
        real(rp), external :: dnrm2

        ! initialize error flag
        error = .false.

        min_idx = minloc(h_diag, dim=1)

        if (h_diag(min_idx) < 0.d0) then
            n_vectors = 2
            allocate (red_space_basis(size(grad), n_vectors + &
                      settings%n_random_trial_vectors))
            red_space_basis(:, 1) = grad/grad_norm
            trial = 0.d0
            trial(min_idx) = 1.d0
            call gram_schmidt(trial, reshape(red_space_basis(:, 1), &
                                             [size(red_space_basis, 1), 1]), settings, &
                              error)
            if (error) return
            red_space_basis(:, 2) = trial
        else
            n_vectors = 1
            allocate (red_space_basis(size(grad), n_vectors + &
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
        logical, intent(out) :: error

        integer(ip) :: i
        real(rp) :: trial(size(red_space_basis, 1))
        real(rp), external :: dnrm2

        ! initialize error flag
        error = .false.

        do i = size(red_space_basis, 2) - settings%n_random_trial_vectors + 1, &
            size(red_space_basis, 2)
            call random_number(trial)
            trial = 2*trial - 1
            do while (dnrm2(size(trial), trial, 1) < 1.d-3)
                call random_number(trial)
                trial = 2*trial - 1
            end do
            call gram_schmidt(trial, red_space_basis(:, :i - 1), settings, &
                              error)
            if (error) return
            red_space_basis(:, i) = trial
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
        logical, intent(out) :: error
        real(rp), intent(inout), optional :: lin_trans_vector(:)
        real(rp), intent(in), optional :: lin_trans_space(:, :)

        real(rp) :: orth(size(space, 2)), dot, norm

        integer(ip) :: iter, i
        real(rp), external :: ddot, dnrm2

        ! initialize error flag
        error = .false.

        if (dnrm2(size(vector), vector, 1) < 1.d-16) then
            call settings%log("Vector passed to Gram-Schmidt procedure is "// &
                "numerically zero.", 1, .true.)
            error = .true.
            return
        else if (size(space, 2) > size(space, 1) - 1) then
            call settings%log("Number of vectors in Gram-Schmidt procedure larger "// &
                "than dimension of vector space.", 1, .true.)
            error = .true.
            return
        end if

        iter = 0
        if (.not. (present(lin_trans_vector) .and. present(lin_trans_space))) then
            do while (.true.)
                do i = 1, size(space, 2)
                    vector = orthogonal_projection(vector, space(:, i))
                end do
                vector = vector / dnrm2(size(vector), vector, 1)

                call dgemv("T", size(vector), size(space, 2), 1.d0, space, &
                           size(vector), vector, 1, 0.d0, orth, 1)
                if (maxval(abs(orth)) < 1.d-14) exit
                iter = iter + 1
                if (iter > 100) then
                    call settings%log("Maximum number of Gram-Schmidt iterations "// &
                                      "reached.", 1, .true.)
                    error = .true.
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

                call dgemv("T", size(vector), size(space, 2), 1.d0, space, &
                           size(vector), vector, 1, 0.d0, orth, 1)
                if (maxval(abs(orth)) < 1.d-14) exit
                iter = iter + 1
                if (iter > 100) then
                    call settings%log("Maximum number of Gram-Schmidt iterations "// &
                                      "reached.", 1, .true.)
                    error = .true.
                    return
                end if
            end do
        end if

    end subroutine gram_schmidt

    subroutine init_solver_settings(self, stability, line_search, davidson, &
                                    jacobi_davidson, prefer_jacobi_davidson, conv_tol, &
                                    n_random_trial_vectors, start_trust_radius, &
                                    n_macro, n_micro, global_red_factor, &
                                    local_red_factor, seed, verbose, logger)
        !
        ! this subroutine sets the optional settings to their default values
        !
        class(solver_settings_type), intent(inout) :: self
        logical, intent(in), optional :: stability, line_search, davidson, &
                                         jacobi_davidson, prefer_jacobi_davidson
        real(rp), intent(in), optional :: conv_tol, start_trust_radius, &
                                          global_red_factor, local_red_factor
        integer(ip), intent(in), optional :: n_random_trial_vectors, n_macro, n_micro, &
                                             seed, verbose
        procedure(logger_type), intent(in), pointer, optional :: logger

        self%stability = set_default(stability, solver_stability_default)
        self%line_search = set_default(line_search, solver_line_search_default)
        self%davidson = set_default(davidson, solver_davidson_default)
        self%jacobi_davidson = set_default(jacobi_davidson, &
                                           solver_jacobi_davidson_default)
        self%prefer_jacobi_davidson = set_default(prefer_jacobi_davidson, &
                                                  solver_prefer_jacobi_davidson_default)
        self%conv_tol = set_default(conv_tol, solver_conv_tol_default)
        self%n_random_trial_vectors = set_default(n_random_trial_vectors, &
                                                  solver_n_random_trial_vectors_default)
        self%start_trust_radius = set_default(start_trust_radius, &
                                              solver_start_trust_radius_default)
        self%n_macro = set_default(n_macro, solver_n_macro_default)
        self%n_micro = set_default(n_micro, solver_n_micro_default)
        self%global_red_factor = set_default(global_red_factor, &
                                             solver_global_red_factor_default)
        self%local_red_factor = set_default(local_red_factor, &
                                            solver_local_red_factor_default)
        self%seed = set_default(seed, solver_seed_default)
        self%verbose = set_default(verbose, solver_verbose_default)
        if (present(logger)) self%logger => logger

    end subroutine init_solver_settings

    subroutine init_stability_settings(self, jacobi_davidson, conv_tol, &
                                       n_random_trial_vectors, n_iter, verbose, &
                                       logger)
        !
        ! this subroutine sets the optional settings to their default values
        !
        class(stability_settings_type), intent(inout) :: self
        logical, intent(in), optional :: jacobi_davidson
        real(rp), intent(in), optional :: conv_tol
        integer(ip), intent(in), optional :: n_random_trial_vectors, n_iter, verbose
        procedure(logger_type), intent(in), pointer, optional :: logger

        self%jacobi_davidson = set_default(jacobi_davidson, &
                                           stability_jacobi_davidson_default)
        self%conv_tol = set_default(conv_tol, stability_conv_tol_default)
        self%n_random_trial_vectors = set_default(n_random_trial_vectors, &
                                               stability_n_random_trial_vectors_default)
        self%n_iter = set_default(n_iter, stability_n_iter_default)
        self%verbose = set_default(verbose, stability_verbose_default)
        if (present(logger)) self%logger => logger

    end subroutine init_stability_settings

    function set_default_real(optional_value, default_value) result(variable)
        !
        ! this function sets a default value for reals
        !
        real(rp), intent(in)  :: default_value
        real(rp), intent(in), optional :: optional_value

        real(rp) :: variable

        if (present(optional_value)) then
            variable = optional_value
        else
            variable = default_value
        end if

    end function set_default_real

    function set_default_integer(optional_value, default_value) result(variable)
        !
        ! this function sets a default value for integers
        !
        integer(ip), intent(in)  :: default_value
        integer(ip), intent(in), optional :: optional_value

        integer(ip) :: variable

        if (present(optional_value)) then
            variable = optional_value
        else
            variable = default_value
        end if

    end function set_default_integer

    function set_default_logical(optional_value, default_value) result(variable)
        !
        ! this function sets a default value for logicals
        !
        logical, intent(in)  :: default_value
        logical, intent(in), optional :: optional_value

        logical :: variable

        if (present(optional_value)) then
            variable = optional_value
        else
            variable = default_value
        end if

    end function set_default_logical

    function level_shifted_diag_precond(vector, mu, h_diag) result(precond_vector)
        !
        ! this function defines the default level-shifted diagonal preconditioner
        !
        real(rp), intent(in) :: vector(:), mu, h_diag(:)
        real(rp) :: precond_vector(size(vector))

        real(rp) :: precond(size(vector))
        
        ! construct level-shifted preconditioner
        precond = h_diag - mu
        where (abs(precond) < 1d-10)
            precond = 1d-10
        end where

        ! precondition vector
        precond_vector = vector / precond
        
    end function level_shifted_diag_precond

    function abs_diag_precond(vector, h_diag) result(precond_vector)
        !
        ! this function defines the default absolute diagonal preconditioner
        !
        real(rp), intent(in) :: vector(:), h_diag(:)
        real(rp) :: precond_vector(size(vector))

        real(rp) :: precond(size(vector))

        ! construct positive-definite preconditioner
        precond = max(abs(h_diag), 1.d-10)
        
        ! precondition vector
        precond_vector = vector / precond
        
    end function abs_diag_precond

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
                                          corr_vector, hess_vector, error)
        !
        ! this subroutine performs the Jacobi-Davidson correction but also returns the 
        ! Hessian linear transformation since this can be reused
        !
        procedure(hess_x_type), intent(in), pointer :: hess_x_funptr
        real(rp), intent(in) :: vector(:), solution(:), eigval
        real(rp), intent(out) :: corr_vector(:), hess_vector(:)
        logical, intent(out) :: error

        real(rp), allocatable :: proj_vector(:)

        ! initialize error flag
        error = .false.
        
        ! project solution out of vector
        proj_vector = orthogonal_projection(vector, solution)

        ! get Hessian linear transformation of projected vector
        hess_vector = hess_x_funptr(proj_vector, error)
        if (error) return

        ! finish construction of correction
        corr_vector = orthogonal_projection(hess_vector - eigval * proj_vector, &
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
        logical, intent(out) :: error
        real(rp), intent(in), optional :: guess(:)
        integer(ip), intent(in), optional :: max_iter

        integer(ip) :: n, max_iterations, iteration
        real(rp), parameter :: eps = epsilon(1.d0)
        real(rp) :: beta_start, beta, phi_bar, rhs1, old_beta, alfa, t_norm2, eps_ln, &
                    old_eps, cs, d_bar, sn, delta, g_bar, root, gamma, phi, g_max, &
                    g_min, tmp, rhs2, a_norm, vec_norm, qr_norm
        real(rp), allocatable :: matvec(:), r1(:), r2(:), y(:), w(:), hw(:), w1(:), &
                                 hw1(:), w2(:), hw2(:), v(:), hv(:)
        logical :: stop_iteration
        real(rp), external :: dnrm2, ddot

        ! initialize error flag
        error = .false.

        ! initialze boolean to stop iterations
        stop_iteration = .false.

        ! size of problem
        n = size(rhs)

        ! allocate solution vector
        allocate (matvec(n))
        allocate (r1(n))
        allocate (r2(n))
        allocate (y(n))
        allocate (w(n))
        allocate (hw(n))
        allocate (w1(n))
        allocate (hw1(n))
        allocate (w2(n))
        allocate (hw2(n))
        allocate (v(n))
        allocate (hv(n))

        ! initial guess
        if (present(guess)) then
            vec = guess
            call jacobi_davidson_correction(hess_x_funptr, vec, solution, eigval, &
                                            matvec, hvec, error)
            if (error) return
            tot_hess_x = tot_hess_x + 1
        else
            vec = 0.d0
            hvec = 0.d0
            matvec = 0.d0
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
        if (beta_start < 1d-14) return

        ! solution must be zero vector if rhs vanishes
        if (dnrm2(n, rhs, 1) < 1d-14) then
            vec = 0.d0
            hvec = 0.d0
            return
        end if

        ! initialize additional quantities
        beta = beta_start
        r2 = r1
        phi_bar = beta_start
        rhs1 = beta_start
        t_norm2 = 0.d0
        eps_ln = 0.d0
        cs = -1.d0
        d_bar = 0.d0
        sn = 0.d0
        w = 0.d0
        hw = 0.d0
        w2 = 0.d0
        hw2 = 0.d0
        g_max = 0.d0
        g_min = huge(1.d0)
        rhs2 = 0.d0

        iteration = 1
        do while (iteration <= max_iterations)
            ! scale trial vector
            v = y / beta

            ! apply Jacobi-Davidson projector to trial vector
            call jacobi_davidson_correction(hess_x_funptr, v, solution, eigval, y, hv, &
                                            error)
            if (error) return
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
            else if (vec_norm > 0.d0 .and. a_norm > 0.d0 .and. &
                qr_norm / (a_norm * vec_norm) <= r_tol) then
                    call settings%log("MINRES: A solution to Ax = b was found, "// &
                                      "given provided tolerance.", 4)
                exit
            ! ||Ar|| / (||A|| ||r||)
            else if (a_norm < 1d-14 .and. root / a_norm <= r_tol) then
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
                error = .true.
                return
            ! these tests ensure convergence is still achieved when r_tol 
            ! approaches machine precision
            else if (vec_norm > 0.d0 .and. a_norm > 0.d0 .and. &
                1.d0 + qr_norm / (a_norm * vec_norm) <= 1.d0) then
                call settings%log("MINRES: A solution to Ax = b was found, given "// &
                                  "provided tolerance.", 4)
                exit
            else if (a_norm < 1d-14 .and. 1.d0 + root / a_norm <= 1.d0) then
                call settings%log("MINRES: A least-squares solution was found, "// &
                                  "given provided tolerance.", 4)
                exit
            end if

            ! increment iteration counter
            iteration = iteration + 1

        end do

        deallocate (matvec)
        deallocate (r1)
        deallocate (r2)
        deallocate (y)
        deallocate (w)
        deallocate (hw)
        deallocate (w1)
        deallocate (hw1)
        deallocate (w2)
        deallocate (hw2)
        deallocate (v)
        deallocate (hv)

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
        end if

    end subroutine log

    subroutine split_string_by_space(input, max_length, substrings)
        !
        ! this function splits a string by spaces to produce substrings of a maximum 
        ! length
        !
        character(*), intent(in) :: input
        integer(ip),intent(in) :: max_length
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
                                      obj_func, hess_x_funptr, settings, use_precond, &
                                      precond, trust_radius, solution, mu, imicro, &
                                      imicro_jacobi_davidson, jacobi_davidson_started, &
                                      max_precision_reached, error)
        !
        ! this subroutine performs level-shifted (Jacobi-)Davidson to solve the trust 
        ! region subproblem
        !
        real(rp), intent(in) :: func, grad(:), grad_norm, h_diag(:)
        integer(ip), intent(in) :: n_param
        procedure(obj_func_type), pointer, intent(in) :: obj_func
        procedure(hess_x_type), pointer, intent(in) :: hess_x_funptr
        class(solver_settings_type), intent(in) :: settings
        logical, intent(in) :: use_precond
        procedure(precond_type), intent(in), pointer, optional :: precond
        real(rp), intent(inout) :: trust_radius
        real(rp), intent(out) :: solution(:), mu
        integer(ip), intent(out) :: imicro, imicro_jacobi_davidson
        logical, intent(out) :: jacobi_davidson_started, max_precision_reached, error

        real(rp), allocatable :: red_space_basis(:, :), h_basis(:, :), aug_hess(:, :), &
                                 red_space_solution(:), red_hess_vec(:)
        integer(ip) :: ntrial, i, initial_imicro
        real(rp), dimension(n_param) :: last_solution_normalized, h_solution, &
                                        residual, solution_normalized, basis_vec, &
                                        h_basis_vec
        logical :: accept_step, micro_converged, func_evaluated, &
                   newton, bracketed
        real(rp) :: aug_hess_min_eigval, residual_norm, red_factor, &
                    initial_residual_norm, new_func, ratio, minres_tol
        real(rp), external :: dnrm2, ddot

        ! initialize error flag
        error = .false.

        ! generate trial vectors
        red_space_basis = generate_trial_vectors(grad, grad_norm, h_diag, settings, &
                                                 error)
        if (error) return

        ! number of trial vectors
        ntrial = size(red_space_basis, 2)

        ! increment number of Hessian linear transformations
        tot_hess_x = tot_hess_x + ntrial

        ! calculate linear transformations of basis vectors
        allocate (h_basis(n_param, ntrial))
        do i = 1, ntrial
            h_basis(:, i) = hess_x_funptr(red_space_basis(:, i), error)
            if (error) return
        end do

        ! construct augmented Hessian in reduced space
        allocate (aug_hess(ntrial + 1, ntrial + 1))
        aug_hess = 0.d0
        call dgemm("T", "N", ntrial, ntrial, n_param, 1.d0, red_space_basis, &
                   n_param, h_basis, n_param, 0.d0, aug_hess(2:, 2:), ntrial)

        ! allocate space for reduced space solution
        allocate (red_space_solution(size(red_space_basis, 2)))

        ! decrease trust radius until micro iterations converge and step is accepted
        last_solution_normalized = 0.d0
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
                if (error) return
                if (aug_hess_min_eigval > -1.d-5) then
                    call newton_step(aug_hess, grad_norm, red_space_basis, &
                                     solution, red_space_solution, settings, error)
                    if (error) return
                    mu = 0.d0
                    if (dnrm2(n_param, solution, 1) < trust_radius) newton = .true.
                end if

                ! otherwise perform bisection to find the level shift
                if (.not. newton) then
                    call bisection(aug_hess, grad_norm, red_space_basis, trust_radius, &
                                   solution, red_space_solution, mu, bracketed, &
                                   settings, error)
                    if (error) return
                    if (.not. bracketed) exit
                end if

                ! calculate Hessian linear transformation of solution
                call dgemv("N", n_param, size(h_basis, 2), 1.d0, h_basis, n_param, &
                           red_space_solution, 1, 0.d0, h_solution, 1)

                ! calculate residual
                residual = grad + h_solution - mu*solution

                ! calculate residual norm
                residual_norm = dnrm2(n_param, residual, 1)

                ! determine reduction factor depending on whether local region is
                ! reached
                if (abs(mu) < 1.d-12) then
                    red_factor = settings%local_red_factor
                else
                    red_factor = settings%global_red_factor
                end if

                ! get normalized solution vector
                solution_normalized = solution / dnrm2(n_param, solution, 1)

                ! reset initial residual norm if solution changes
                if (ddot(n_param, last_solution_normalized, 1, solution_normalized, &
                         1)**2 < 0.5d0) then
                    initial_imicro = imicro
                    initial_residual_norm = residual_norm
                end if

                ! check if micro iterations have converged
                if (residual_norm < max(red_factor*grad_norm, 1d-12)) then
                    micro_converged = .true.
                    exit
                ! check if Jacobi-Davidson is used and has not been started
                else if (settings%jacobi_davidson .and. .not. jacobi_davidson_started) then
                    ! check residual has not decreased sufficiently or if maximum of 
                    ! Davidson iterations has been reached
                    if ((imicro - initial_imicro >= 10 .and. &
                         residual_norm > 0.8*initial_residual_norm) .or. &
                        imicro > 0.8 * settings%n_micro) then
                        ! evaluate function at approximate point
                        new_func = obj_func(solution, error)
                        if (error) return

                        ! calculate ratio of evaluated function and predicted function
                        ratio = (new_func - func) / ddot(n_param, solution, 1, &
                                                        grad + 0.5*h_solution, 1)

                        ! switch to Jacobi-Davidson only if current solution would lead 
                        ! to trust radius increase when the solution is already at the 
                        ! trust region boundary
                        if (settings%prefer_jacobi_davidson .or. &
                            (ratio > 0.75d0 .and. dnrm2(n_param, solution, 1) &
                             > 0.99d0 * trust_radius)) then
                            jacobi_davidson_started = .true.
                            imicro_jacobi_davidson = imicro
                        ! decrease trust radius
                        else
                            func_evaluated = .true.
                            exit
                        end if
                    end if
                else if (imicro - initial_imicro >= 10 .and. &
                         residual_norm > 0.8*initial_residual_norm) then
                    exit
                end if

                ! save current solution
                last_solution_normalized = solution_normalized

                if (.not. jacobi_davidson_started) then
                    ! precondition residual
                    if (use_precond) then
                        basis_vec = precond(residual, mu, error)
                        if (error) return
                    else
                        basis_vec = level_shifted_diag_precond(residual, mu, h_diag)
                    end if

                    ! orthonormalize to current orbital space to get new basis vector
                    call gram_schmidt(basis_vec, red_space_basis, settings, error)
                    if (error) return

                    ! add linear transformation of new basis vector
                    h_basis_vec = hess_x_funptr(basis_vec, error)
                    if (error) return

                    ! increment Hessian linear transformations
                    tot_hess_x = tot_hess_x + 1

                else
                    ! solve Jacobi-Davidson correction equations
                    minres_tol = 3.d0 ** (-(imicro - imicro_jacobi_davidson))
                    call minres(-residual, hess_x_funptr, solution_normalized, mu, &
                                minres_tol, basis_vec, h_basis_vec, settings, error)
                    if (error) return

                    ! orthonormalize to current orbital space to get new basis vector
                    call gram_schmidt(basis_vec, red_space_basis, settings, error, &
                                      lin_trans_vector=h_basis_vec, &
                                      lin_trans_space=h_basis)
                    if (error) return

                    ! check if resulting linear transformation still respects Hessian 
                    ! symmetry which can happen due to numerical noise accumulation
                    if (abs(ddot(n_param, &
                                 red_space_basis(:, size(red_space_basis, 2)), 1, &
                                 h_basis_vec, 1) - &
                            ddot(n_param, basis_vec, 1, &
                                 h_basis(:, size(red_space_basis, 2)), 1)) > 1d-12) then
                        h_basis_vec = hess_x_funptr(basis_vec, error)
                        if (error) return
                    end if

                end if

                ! add new trial vector to orbital space
                call add_column(red_space_basis, basis_vec)

                ! add linear transformation of new basis vector
                call add_column(h_basis, h_basis_vec)

                ! construct new augmented Hessian
                allocate (red_hess_vec(size(red_space_basis, 2) + 1))
                red_hess_vec(1) = 0.d0
                call dgemv("T", n_param, size(red_space_basis, 2), 1.d0, &
                           red_space_basis, n_param, &
                           h_basis(:, size(red_space_basis, 2)), 1, 0.d0, &
                           red_hess_vec(2:), 1)
                call extend_symm_matrix(aug_hess, red_hess_vec)
                deallocate (red_hess_vec)

                ! reallocate reduced space solution
                deallocate (red_space_solution)
                allocate (red_space_solution(size(red_space_basis, 2)))
            end do

            if (.not. func_evaluated) then
                ! evaluate function at predicted point
                new_func = obj_func(solution, error)
                if (error) return

                ! calculate ratio of evaluated function and predicted function
                ratio = (new_func - func) / ddot(n_param, solution, 1, &
                                                 grad + 0.5d0*h_solution, 1)
            end if

            ! decide whether to accept step and modify trust radius
            accept_step = accept_trust_region_step(solution, ratio, micro_converged, &
                                                   settings, trust_radius, &
                                                   max_precision_reached)
            if (max_precision_reached) exit
        end do

        ! deallocate quantities from microiterations
        deallocate (red_space_solution)
        deallocate (aug_hess)
        deallocate (h_basis)
        deallocate (red_space_basis)

    end subroutine level_shifted_davidson

    subroutine truncated_conjugate_gradient(func, grad, h_diag, n_param, obj_func, &
                                            hess_x_funptr, use_precond, precond, &
                                            settings, trust_radius, solution, imicro, &
                                            max_precision_reached, error)
        !
        ! this subroutine performs truncated conjugate gradient to solve the trust 
        ! region subproblem
        !
        real(rp), intent(in) :: func, grad(:), h_diag(:)
        integer(ip), intent(in) :: n_param
        procedure(obj_func_type), pointer, intent(in) :: obj_func
        procedure(hess_x_type), pointer, intent(in) :: hess_x_funptr
        class(solver_settings_type), intent(in) :: settings
        logical, intent(in) :: use_precond
        procedure(precond_type), intent(in), pointer, optional :: precond
        real(rp), intent(inout) :: trust_radius
        real(rp), intent(out) :: solution(:)
        integer(ip), intent(out) :: imicro
        logical, intent(out) :: max_precision_reached, error

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
        real(rp), external :: ddot, dnrm2

        ! initialize error flag
        error = .false.

        ! allocate arrays
        allocate (h_solution(n_param))
        allocate (random_vector(n_param))
        allocate (solutions(n_param, 1))
        allocate (h_solutions(n_param, 1))

        ! initialize solution
        solution = 0.d0
        h_solution = 0.d0
        solutions(:, 1) = solution
        h_solutions(:, 1) = h_solution
        model_func = 0.d0
        
        ! initialize residual and add random noise, residual should include h_solution 
        ! if not starting at zero
        residual = grad
        random_vector = 0.d0
        do while (dnrm2(size(random_vector), random_vector, 1) == 0.d0)
            call random_number(random_vector)
            random_vector = 2 * random_vector - 1
        end do
        residual = residual + 1.d-4 * dnrm2(n_param, residual, 1) * random_vector / &
                   dnrm2(n_param, random_vector, 1)
        initial_residual_norm = dnrm2(n_param, residual, 1)

        ! initialize preconditioned residual and direction,
        if (use_precond) then
            precond_residual = precond(residual, 0.d0, error)
            if (error) return
        else
            precond_residual = abs_diag_precond(residual, h_diag)
        end if
        direction = -precond_residual

        ! start microiteration loop
        micro_converged = .false.
        do imicro = 1, settings%n_micro - 1
            ! get Hessian linear transformation of direction
            hess_direction = hess_x_funptr(direction, error)
            if (error) return

            ! increment Hessian linear transformations
            tot_hess_x = tot_hess_x + 1

            ! calculate curvature
            curvature = ddot(n_param, direction, 1, hess_direction, 1)

            ! get step size along new direction
            step_size = ddot(n_param, residual, 1, precond_residual, 1) / curvature

            ! precondition current solution and direction
            if (use_precond) then
                precond_solution = precond(solution, 0.d0, error)
                if (error) return
                precond_direction = precond(direction, 0.d0, error)
                if (error) return
            else
                precond_solution = abs_diag_precond(solution, h_diag)
                precond_direction = abs_diag_precond(direction, h_diag)
            end if

            ! calculate dot products
            solution_dot = ddot(n_param, solution, 1, precond_solution, 1)
            solution_direction_dot = ddot(n_param, solution, 1, precond_direction, 1)
            direction_dot = ddot(n_param, direction, 1, precond_direction, 1)

            ! calculate total step length
            step_length = solution_dot + 2 * step_size * solution_direction_dot + &
                          step_size ** 2 * direction_dot

            if (curvature < 0.d0 .or. step_length >= trust_radius ** 2) then
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
            model_func_new = ddot(size(solution_new), solution_new, 1, &
                                  grad + 0.5d0 * h_solution_new, 1)

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
            if (use_precond) then
                precond_residual_new = precond(residual_new, 0.d0, error)
                if (error) return
            else
                precond_residual_new = abs_diag_precond(residual_new, h_diag)
            end if

            ! check for linear or superlinear (in this case quadratic) convergence
            if (dnrm2(n_param, residual_new, 1) <= initial_residual_norm * &
                min(1.d-3, initial_residual_norm)) then
                micro_converged = .true.
                exit
            end if

            ! get new search direction
            direction = -precond_residual_new + &
                        ddot(size(residual_new), residual_new, 1, &
                                precond_residual_new, 1) / &
                        ddot(size(residual), residual, 1, precond_residual, 1) * &
                        direction
            
            ! save new model
            model_func = model_func_new
            residual = residual_new
            precond_residual = precond_residual_new
        end do

        ! evaluate function at predicted point
        new_func = obj_func(solution, error)
        if (error) return

        if (abs(new_func - func) / max(abs(new_func), abs(func)) > 1.d-14) then
            ! calculate ratio of evaluated function and predicted function
            ratio = (new_func - func) / ddot(n_param, solution, 1, &
                                             grad + 0.5d0*h_solution, 1)

            ! reduce trust region until step is accepted
            do while (.true.)
                ! decide whether to accept step and modify trust radius
                accept_step = accept_trust_region_step(solution, ratio, &
                                                       micro_converged, settings, &
                                                       trust_radius, &
                                                       max_precision_reached)
                if (accept_step .or. max_precision_reached) exit

                ! check if step exceeds new trust region boundary
                if (ddot(n_param, solution, 1, abs_diag_precond(solution, h_diag), 1) &
                    > trust_radius ** 2) then
                    ! find step that exceeds trust region boundary
                    do i = 1, size(solutions, 2)
                        if (ddot(n_param, solutions(:, i), 1, &
                                 abs_diag_precond(solutions(:, i), h_diag), 1) > &
                            trust_radius ** 2) then
                            ! get previous step
                            solution = solutions(:, i - 1)
                            h_solution = h_solutions(:, i - 1)

                            ! get direction
                            direction = solutions(:, i) - solutions(:, i - 1)
                            hess_direction = h_solutions(:, i) - h_solutions(:, i - 1)

                            ! precondition current solution and direction
                            if (use_precond) then
                                precond_solution = precond(solution, 0.d0, error)
                                if (error) return
                                precond_direction = precond(direction, 0.d0, error)
                                if (error) return
                            else
                                precond_solution = abs_diag_precond(solution, h_diag)
                                precond_direction = abs_diag_precond(direction, h_diag)
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
                            if (error) return

                            ! calculate ratio of evaluated function and predicted 
                            ! function
                            ratio = (new_func - func) / ddot(n_param, solution, 1, &
                                                             grad + 0.5d0*h_solution, 1)

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
        deallocate (h_solution)
        deallocate (random_vector)
        deallocate (solutions)
        deallocate (h_solutions)

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
        class(solver_settings_type), intent(in) :: settings
        real(rp), intent(inout) :: trust_radius
        logical, intent(out) :: max_precision_reached

        ! default to maximum precision not yet reached
        max_precision_reached = .false.

        ! decrease trust radius if micro iterations are unable to converge, if function 
        ! value has not decreased or if individual orbitals change too much
        if (.not. micro_converged .or. ratio < 0.d0 .or. &
            any(abs(solution) > pi/4)) then
            trust_radius = 0.7d0*trust_radius
            accept_trust_region_step = .false.
            if (trust_radius < 1.d-14) then
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
        class(solver_settings_type), intent(inout) :: settings
        integer(ip), intent(in) :: n_param
        real(rp), intent(in) :: grad(:)
        logical, intent(out) :: error

        character(300) :: msg

        ! initialize error flag
        error = .false.

        ! check that number of parameters is positive
        if (n_param < 1) then
            call settings%log("Number of parameters should be larger than 0.", 1, &
                              .true.)
            error = .true.
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
                              "update_orbs does not equal number of parameters ", 1, &
                              .true.)
            error = .true.
            return
        end if

    end subroutine sanity_check

end module opentrustregion
