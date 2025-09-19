! Copyright (C) 2025- Jonas Greiner
!
! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.

module c_interface_unit_tests

    use opentrustregion, only: rp, stderr
    use, intrinsic :: iso_c_binding, only: c_long, c_double, c_bool, c_ptr, c_loc, &
                                           c_null_ptr, c_funptr, c_funloc, &
                                           c_null_funptr, c_char

    implicit none

    real(rp), parameter :: tol = 1.d-10

    ! these pointer targets need to be defined out here to ensure that these arrays
    ! do not go out of scope when mock_update_orbs or mock_hess_x exit
    integer(c_long), parameter :: n_param = 3_c_long
    real(c_double), dimension(n_param), target :: grad_c, h_diag_c, hess_x_c, &
                                                  precond_residual_c

    ! logical to test logging function
    logical :: test_logger

contains

    function mock_update_orbs(kappa, func, grad_c_ptr, h_diag_c_ptr, &
                                hess_x_c_funptr) result(error) bind(C)
        !
        ! this subroutine is a test subroutine for the orbital update C function
        !
        real(c_double), intent(in) :: kappa(*)
        real(c_double), intent(out) :: func
        type(c_ptr), intent(out) :: grad_c_ptr, h_diag_c_ptr
        type(c_funptr), intent(out) :: hess_x_c_funptr
        logical(c_bool) :: error

        func = sum(kappa(:n_param))

        grad_c = 2*kappa(:n_param)
        grad_c_ptr = c_loc(grad_c)

        h_diag_c = 3*kappa(:n_param)
        h_diag_c_ptr = c_loc(h_diag_c)

        hess_x_c_funptr = c_funloc(mock_hess_x)

        error = .false.

    end function mock_update_orbs

    function mock_hess_x(x, hess_x_c_ptr) result(error) bind(C)
        !
        ! this subroutine is a test subroutine for the Hessian linear transformation
        ! C function
        !
        real(c_double), intent(in) :: x(*)
        type(c_ptr), intent(out) :: hess_x_c_ptr
        logical(c_bool) :: error

        hess_x_c = 4*x(:n_param)
        hess_x_c_ptr = c_loc(hess_x_c)

        error = .false.

    end function mock_hess_x

    function mock_obj_func(kappa, func) result(error) bind(C)
        !
        ! this function is a test function for the C objective function
        !
        real(c_double), intent(in) :: kappa(*)
        real(c_double), intent(out) :: func
        logical(c_bool) :: error

        func = sum(kappa(:n_param))

        error = .false.

    end function mock_obj_func

    function mock_precond(residual, mu, precond_residual_c_ptr) result(error) bind(C)
        !
        ! this function is a test function for the C preconditioner function
        !
        real(c_double), intent(in) :: residual(*), mu
        type(c_ptr), intent(out) :: precond_residual_c_ptr
        logical(c_bool) :: error

        precond_residual_c = mu*residual(:n_param)
        precond_residual_c_ptr = c_loc(precond_residual_c)

        error = .false.

    end function mock_precond

    function mock_conv_check(converged) result(error) bind(C)
        !
        ! this function is a test function for the convergence check function
        !
        logical(c_bool), intent(out) :: converged
        logical(c_bool) :: error

        converged = .true.

        error = .false.
        
    end function mock_conv_check

    subroutine mock_logger(message_c) bind(C)
        !
        ! this function is a test function for the C logging function
        !
        character(kind=c_char), intent(in) :: message_c(*)
        character(4) :: message

        message = transfer(message_c(1:4), message)
        if (message == "test") test_logger = .true.

    end subroutine mock_logger

    logical(c_bool) function test_solver_c_wrapper() bind(C)
        !
        ! this function tests the C wrapper for the solver
        !
        use c_interface, only: solver, solver_c_wrapper
        use opentrustregion_mock, only: mock_solver, test_passed

        type(c_funptr) :: update_orbs_c_funptr, obj_func_c_funptr, &
                          precond_c_funptr = c_null_funptr, &
                          conv_check_c_funptr = c_null_funptr, &
                          logger_c_funptr = c_null_funptr
        logical(c_bool) :: error
        integer(c_long), target :: n_random_trial_vectors = 5_c_long, &
                                   n_macro = 300_c_long, n_micro = 200_c_long, &
                                   seed = 33_c_long, verbose = 3_c_long
        logical(c_bool), target :: stability = .false., line_search = .true., &
                                   davidson = .false., jacobi_davidson = .false., &
                                   prefer_jacobi_davidson = .true.
        real(c_double), target :: conv_tol = 1e-3_c_double, &
                                  start_trust_radius = 0.2_c_double, &
                                  global_red_factor = 1e-2_c_double, &
                                  local_red_factor = 1e-3_c_double
        type(c_ptr) :: stability_c_ptr = c_null_ptr, line_search_c_ptr = c_null_ptr, &
                       davidson_c_ptr = c_null_ptr, &
                       jacobi_davidson_c_ptr = c_null_ptr, &
                       prefer_jacobi_davidson_c_ptr = c_null_ptr, &
                       conv_tol_c_ptr = c_null_ptr, &
                       n_random_trial_vectors_c_ptr = c_null_ptr, &
                       start_trust_radius_c_ptr = c_null_ptr, &
                       n_macro_c_ptr = c_null_ptr, &
                       n_micro_c_ptr = c_null_ptr, &
                       global_red_factor_c_ptr = c_null_ptr, &
                       local_red_factor_c_ptr = c_null_ptr, &
                       seed_c_ptr = c_null_ptr, verbose_c_ptr = c_null_ptr

        ! assume tests pass
        test_solver_c_wrapper = .true.

        ! inject mock function
        solver => mock_solver

        ! get C function pointers to Fortran functions
        update_orbs_c_funptr = c_funloc(mock_update_orbs)
        obj_func_c_funptr = c_funloc(mock_obj_func)

        ! call solver first without associated optional arguments which should produce
        ! default values
        error = solver_c_wrapper(update_orbs_c_funptr, obj_func_c_funptr, n_param, &
                                 precond_c_funptr, conv_check_c_funptr, &
                                 stability_c_ptr, line_search_c_ptr, davidson_c_ptr, &
                                 jacobi_davidson_c_ptr, prefer_jacobi_davidson_c_ptr, &
                                 conv_tol_c_ptr, n_random_trial_vectors_c_ptr, &
                                 start_trust_radius_c_ptr, n_macro_c_ptr, &
                                 n_micro_c_ptr, global_red_factor_c_ptr, &
                                 local_red_factor_c_ptr, seed_c_ptr, verbose_c_ptr, &
                                 logger_c_funptr)
                              

        ! check if test has passed
        test_solver_c_wrapper = test_passed

        ! check if output variables are as expected
        if (error) then
            test_solver_c_wrapper = .false.
            write (stderr, *) "test_solver_c_wrapper failed: Returned error "// &
                "boolean wrong."
        end if

        ! associate optional arguments with values
        precond_c_funptr = c_funloc(mock_precond)
        conv_check_c_funptr = c_funloc(mock_conv_check)
        stability_c_ptr = c_loc(stability)
        line_search_c_ptr = c_loc(line_search)
        davidson_c_ptr = c_loc(davidson)
        jacobi_davidson_c_ptr = c_loc(jacobi_davidson)
        prefer_jacobi_davidson_c_ptr = c_loc(prefer_jacobi_davidson)
        conv_tol_c_ptr = c_loc(conv_tol)
        n_random_trial_vectors_c_ptr = c_loc(n_random_trial_vectors)
        start_trust_radius_c_ptr = c_loc(start_trust_radius)
        n_macro_c_ptr = c_loc(n_macro)
        n_micro_c_ptr = c_loc(n_micro)
        global_red_factor_c_ptr = c_loc(global_red_factor)
        local_red_factor_c_ptr = c_loc(local_red_factor)
        seed_c_ptr = c_loc(seed)
        verbose_c_ptr = c_loc(verbose)
        logger_c_funptr = c_funloc(mock_logger)

        ! initialize logger logical
        test_logger = .true.     

        ! call solver with associated optional arguments
        error = solver_c_wrapper(update_orbs_c_funptr, obj_func_c_funptr, n_param, &
                                 precond_c_funptr, conv_check_c_funptr, &
                                 stability_c_ptr, line_search_c_ptr, davidson_c_ptr, &
                                 jacobi_davidson_c_ptr, prefer_jacobi_davidson_c_ptr, &
                                 conv_tol_c_ptr, n_random_trial_vectors_c_ptr, &
                                 start_trust_radius_c_ptr, n_macro_c_ptr, &
                                 n_micro_c_ptr, global_red_factor_c_ptr, &
                                 local_red_factor_c_ptr, seed_c_ptr, verbose_c_ptr, &
                                 logger_c_funptr)

        ! check if logging subroutine was correctly called
        if (.not. test_logger) then
            test_solver_c_wrapper = .false.
            write (stderr, *) "test_solver_c_wrapper failed: Called logging "// &
                "subroutine wrong."
        end if

        ! check if output variables are as expected
        if (error) then
            test_solver_c_wrapper = .false.
            write (stderr, *) "test_solver_c_wrapper failed: Returned error "// &
                "boolean wrong."
        end if

        ! check if test has passed
        test_solver_c_wrapper = test_solver_c_wrapper .and. test_passed

    end function test_solver_c_wrapper

    logical(c_bool) function test_stability_check_c_wrapper() bind(C)
        !
        ! this function tests the C wrapper for the stability check
        !
        use c_interface, only: stability_check, stability_check_c_wrapper
        use opentrustregion_mock, only: mock_stability_check, test_passed

        real(c_double), dimension(n_param) :: kappa
        type(c_funptr) :: hess_x_c_funptr, precond_c_funptr = c_null_funptr, &
                          logger_c_funptr = c_null_funptr
        logical(c_bool) :: stable, error
        logical(c_bool), target :: jacobi_davidson = .false.
        real(c_double), target :: conv_tol = 1e-3_c_double
        integer(c_long), target :: n_random_trial_vectors = 3_c_long, &
                                   n_iter = 50_c_long, verbose = 3_c_long
        type(c_ptr) :: jacobi_davidson_c_ptr = c_null_ptr, &
                       conv_tol_c_ptr = c_null_ptr, &
                       n_random_trial_vectors_c_ptr = c_null_ptr, &
                       n_iter_c_ptr = c_null_ptr, verbose_c_ptr = c_null_ptr

        ! assume tests pass
        test_stability_check_c_wrapper = .true.

        ! inject mock function
        stability_check => mock_stability_check

        ! get C function pointers to Fortran functions
        hess_x_c_funptr = c_funloc(mock_hess_x)

        ! initialize Hessian diagonal and get C pointers
        h_diag_c = 3.d0

        ! call stability check first without associated optional arguments which should 
        ! produce default values
        error = stability_check_c_wrapper(h_diag_c, hess_x_c_funptr, n_param, stable, &
                                          kappa, precond_c_funptr, &
                                          jacobi_davidson_c_ptr, conv_tol_c_ptr, &
                                          n_random_trial_vectors_c_ptr, n_iter_c_ptr, &
                                          verbose_c_ptr, logger_c_funptr)

        ! check if test has passed
        test_stability_check_c_wrapper = test_passed

        ! check if output variables are as expected
        if (stable) then
            test_stability_check_c_wrapper = .false.
            write (stderr, *) "test_stability_check_c_wrapper failed: Returned "// &
                "stability boolean wrong."
        end if

        if (any(abs(kappa - 1.d0) > tol)) then
            test_stability_check_c_wrapper = .false.
            write (stderr, *) "test_stability_check_c_wrapper failed: Returned "// &
                "direction wrong."
        end if

        if (error) then
            test_stability_check_c_wrapper = .false.
            write (stderr, *) "test_stability_check_c_wrapper failed: Returned "// &
                "error boolean wrong."
        end if

        ! associate optional arguments with values
        precond_c_funptr = c_funloc(mock_precond)
        jacobi_davidson_c_ptr = c_loc(jacobi_davidson)
        conv_tol_c_ptr = c_loc(conv_tol)
        n_random_trial_vectors_c_ptr = c_loc(n_random_trial_vectors)
        n_iter_c_ptr = c_loc(n_iter)
        verbose_c_ptr = c_loc(verbose)
        logger_c_funptr = c_funloc(mock_logger)

        ! initialize logger logical
        test_logger = .true.

        ! call stability check with associated optional arguments
        error = stability_check_c_wrapper(h_diag_c, hess_x_c_funptr, n_param, stable, &
                                          kappa, precond_c_funptr, &
                                          jacobi_davidson_c_ptr, conv_tol_c_ptr, &
                                          n_random_trial_vectors_c_ptr, n_iter_c_ptr, &
                                          verbose_c_ptr, logger_c_funptr)

        ! check if logging subroutine was correctly called
        if (.not. test_logger) then
            test_stability_check_c_wrapper = .false.
            write (stderr, *) "test_stability_check_c_wrapper failed: Called "// &
                "logging subroutine wrong."
        end if

        ! check if output variables are as expected
        if (stable) then
            test_stability_check_c_wrapper = .false.
            write (stderr, *) "test_stability_check_c_wrapper failed: Returned "// &
                "stability boolean wrong."
        end if

        if (any(abs(kappa - 1.d0) > tol)) then
            test_stability_check_c_wrapper = .false.
            write (stderr, *) "test_stability_check_c_wrapper failed: Returned "// &
                "direction wrong."
        end if

        if (error) then
            test_stability_check_c_wrapper = .false.
            write (stderr, *) "test_stability_check_c_wrapper failed: Returned "// &
                "error boolean wrong."
        end if

        ! check if test has passed
        test_stability_check_c_wrapper = test_passed .and. &
                                         test_stability_check_c_wrapper

    end function test_stability_check_c_wrapper

    logical(c_bool) function test_update_orbs_c_wrapper() bind(C)
        !
        ! this function tests the C wrapper for the orbital update
        !
        use opentrustregion, only: hess_x_type
        use c_interface, only: update_orbs_before_wrapping, update_orbs_c_wrapper

        real(rp), dimension(n_param) :: kappa, grad, h_diag, x, hess_x
        real(rp) :: func
        procedure(hess_x_type), pointer :: hess_x_funptr
        logical :: error

        ! assume tests pass
        test_update_orbs_c_wrapper = .true.

        ! initialize kappa
        kappa = 1.d0

        ! inject mock subroutine
        update_orbs_before_wrapping => mock_update_orbs

        ! call orbital updating subroutine
        call update_orbs_c_wrapper(kappa, func, grad, h_diag, hess_x_funptr, error)

        ! check if error is as expected
        if (error) then
            test_update_orbs_c_wrapper = .false.
            write (stderr, *) "test_update_orbs_c_wrapper failed: Returned error."
        end if

        ! check if function value is as expected
        if (abs(func - 3.d0) > tol) then
            test_update_orbs_c_wrapper = .false.
            write (stderr, *) "test_update_orbs_c_wrapper failed: Returned "// &
                "objective function wrong."
        end if

        ! check if gradient is as expected
        if (any(abs(grad - 2.d0) > tol)) then
            test_update_orbs_c_wrapper = .false.
            write (stderr, *) "test_update_orbs_c_wrapper failed: Returned "// &
                "gradient wrong."
        end if

        ! check if Hessian diagonal is as expected
        if (any(abs(h_diag - 3.d0) > tol)) then
            test_update_orbs_c_wrapper = .false.
            write (stderr, *) "test_update_orbs_c_wrapper failed: Returned Hessian "// &
                "diagonal wrong."
        end if

        ! call Hessian linear transformation function
        x = 1.d0
        hess_x = hess_x_funptr(x, error)

        ! check if error is as expected
        if (error) then
            test_update_orbs_c_wrapper = .false.
            write (stderr, *) "test_update_orbs_c_wrapper failed: Returned Hessian "// &
                "linear transformation function produced error."
        end if

        ! check if Hessian linear transformation is as expected
        if (any(abs(hess_x - 4.d0) > tol)) then
            test_update_orbs_c_wrapper = .false.
            write (stderr, *) "test_update_orbs_c_wrapper failed: Returned Hessian "// &
                "linear transformation of returned Hessian linear transformation "// &
                "function wrong."
        end if

    end function test_update_orbs_c_wrapper

    logical(c_bool) function test_hess_x_c_wrapper() bind(C)
        !
        ! this function tests the C wrapper for the Hessian linear transformation
        !
        use c_interface, only: hess_x_before_wrapping, hess_x_c_wrapper

        real(rp), dimension(n_param) :: x, hess_x
        logical :: error

        ! assume tests pass
        test_hess_x_c_wrapper = .true.

        ! inject mock subroutine
        hess_x_before_wrapping => mock_hess_x

        ! call function
        x = 1.d0
        hess_x = hess_x_c_wrapper(x, error)

        ! check if error is as expected
        if (error) then
            test_hess_x_c_wrapper = .false.
            write (stderr, *) "test_hess_x_c_wrapper failed: Returned error."
        end if

        ! check if Hessian linear transformation is as expected
        if (any(abs(hess_x - 4.d0) > tol)) then
            test_hess_x_c_wrapper = .false.
            write (stderr, *) "test_hess_x_c_wrapper failed: Returned Hessian "// &
                "linear transformation of returned Hessian linear transformation "// &
                "function wrong."
        end if

    end function test_hess_x_c_wrapper

    logical(c_bool) function test_obj_func_c_wrapper() bind(C)
        !
        ! this function tests the C wrapper for the objective function
        !
        use c_interface, only: obj_func_before_wrapping, obj_func_c_wrapper

        real(rp) :: kappa(n_param), obj_func
        logical :: error

        ! assume tests pass
        test_obj_func_c_wrapper = .true.

        ! initialize kappa
        kappa = 1.d0

        ! inject mock function
        obj_func_before_wrapping => mock_obj_func

        ! call function
        obj_func = obj_func_c_wrapper(kappa, error)

        ! check if error is as expected
        if (error) then
            test_obj_func_c_wrapper = .false.
            write (stderr, *) "test_obj_func_c_wrapper failed: Returned error."
        end if

        ! check if function value is as expected
        if (abs(obj_func - 3.d0) > tol) then
            test_obj_func_c_wrapper = .false.
            write (stderr, *) "test_obj_func_c_wrapper failed: Returned objective "// &
                "function wrong."
        end if

    end function test_obj_func_c_wrapper

    logical(c_bool) function test_precond_c_wrapper() bind(C)
        !
        ! this function tests the C wrapper for the preconditioner function
        !
        use c_interface, only: precond_before_wrapping, precond_c_wrapper

        real(rp), dimension(n_param) :: residual, precond_residual
        logical :: error

        ! assume tests pass
        test_precond_c_wrapper = .true.

        ! initialize kappa
        residual = 1.d0

        ! inject mock subroutine
        precond_before_wrapping => mock_precond

        ! call function
        precond_residual = precond_c_wrapper(residual, 5.d0, error)

        ! check if error is as expected
        if (error) then
            test_precond_c_wrapper = .false.
            write (stderr, *) "test_precond_c_wrapper failed: Returned error."
        end if

        ! check if preconditioned residual is as expected
        if (any(abs(precond_residual - 5.d0) > tol)) then
            test_precond_c_wrapper = .false.
            write (stderr, *) "test_precond_c_wrapper failed: Returned "// &
                "preconditioner wrong."
        end if

    end function test_precond_c_wrapper

    logical(c_bool) function test_conv_check_c_wrapper() bind(C)
        !
        ! this function tests the C wrapper for the convergence check function
        !
        use c_interface, only: conv_check_before_wrapping, conv_check_c_wrapper

        logical :: converged, error

        ! assume tests pass
        test_conv_check_c_wrapper = .true.

        ! inject mock function
        conv_check_before_wrapping => mock_conv_check

        ! call function
        converged = conv_check_c_wrapper(error)

        ! check if error is as expected
        if (error) then
            test_conv_check_c_wrapper = .false.
            write (stderr, *) "test_conv_check_c_wrapper failed: Returned error."
        end if

        ! check if convergence boolean is as expected
        if (.not. converged) then
            test_conv_check_c_wrapper = .false.
            write (stderr, *) "test_conv_check_c_wrapper failed: Returned "// &
                "convergence logical wrong."
        end if

    end function test_conv_check_c_wrapper

    logical(c_bool) function test_logger_c_wrapper() bind(C)
        !
        ! this function tests the C wrapper for the logging function
        !
        use c_interface, only: logger_before_wrapping, logger_c_wrapper

        logical :: error

        ! assume tests pass
        test_logger_c_wrapper = .true.

        ! inject mock subroutine
        logger_before_wrapping => mock_logger

        ! call subroutine
        test_logger = .false.
        call logger_c_wrapper("test")

        ! check if logging test boolean is as expected
        if (.not. test_logger) then
            test_logger_c_wrapper = .false.
            write (stderr, *) "test_logger_c_wrapper failed: Returned logging "// &
                "subroutine wrong."
        end if

    end function test_logger_c_wrapper

    logical(c_bool) function test_set_default_c_ptr() bind(C)
        !
        ! this function tests the function that sets default values for variables
        !
        use c_interface, only: set_default_c_ptr

        real(c_double), target :: real_target = 2.0_c_double
        logical(c_bool), target :: bool_target = .true.
        integer(c_long), target :: int_target = 2_c_long

        ! assume tests pass
        test_set_default_c_ptr = .true.

        ! test different combinations and types of optional and default arguments
        if (abs(set_default_c_ptr(c_loc(real_target), 1.d0) - 2.d0) > tol) then
            write (stderr, *) "test_set_default_c_ptr failed: Optional real "// &
                "argument not set correctly."
            test_set_default_c_ptr = .false.
        end if

        if (abs(set_default_c_ptr(c_null_ptr, 1.d0) - 1.d0) > tol) then
            write (stderr, *) "test_set_default_c_ptr failed: Default real "// &
                "argument not set correctly."
            test_set_default_c_ptr = .false.
        end if

        if (set_default_c_ptr(c_loc(bool_target), .false.) .neqv. .true.) then
            write (stderr, *) "test_set_default_c_ptr failed: Optional logical "// &
                "argument not set correctly."
            test_set_default_c_ptr = .false.
        end if

        if (set_default_c_ptr(c_null_ptr, .false.) .neqv. .false.) then
            write (stderr, *) "test_set_default_c_ptr failed: Default logical "// &
                "argument not set correctly."
            test_set_default_c_ptr = .false.
        end if

        if (set_default_c_ptr(c_loc(int_target), 1) /= 2) then
            write (stderr, *) "test_set_default_c_ptr failed: Optional integer "// &
                "argument not set correctly."
            test_set_default_c_ptr = .false.
        end if

        if (set_default_c_ptr(c_null_ptr, 1) /= 1) then
            write (stderr, *) "test_set_default_c_ptr failed: Default integer "// &
                "argument not set correctly."
            test_set_default_c_ptr = .false.
        end if

    end function test_set_default_c_ptr

end module c_interface_unit_tests
