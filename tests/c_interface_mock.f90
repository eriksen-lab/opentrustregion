! Copyright (C) 2025- Jonas Greiner
!
! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.

module c_interface_mock

    use opentrustregion, only: stderr
    use c_interface, only: c_rp, c_ip, solver_c_wrapper_type, &
                           stability_check_c_wrapper_type
    use test_reference, only: ref_settings_c
    use, intrinsic :: iso_c_binding, only: c_bool, c_ptr, c_funptr, c_f_pointer, &
                                           c_f_procpointer, c_associated, c_null_char

    implicit none

    logical(c_bool) :: test_solver_interface = .true., &
                       test_stability_check_interface = .true.

    ! create function pointers to ensure that routines comply with interface
    procedure(solver_c_wrapper_type), pointer :: solver_c_wrapper_ptr => &
        mock_solver_c_wrapper
    procedure(stability_check_c_wrapper_type), pointer :: &
        stability_check_c_wrapper_ptr => mock_stability_check_c_wrapper

contains

    logical(c_bool) function test_solver_result() bind(C)
        !
        ! this function extracts whether the solver python interface test has passed
        !
        test_solver_result = test_solver_interface

    end function test_solver_result

    logical(c_bool) function test_stability_check_result() bind(C)
        !
        ! this function extracts whether the stability check python interface test has
        ! passed
        !
        test_stability_check_result = test_stability_check_interface

    end function test_stability_check_result

    function mock_solver_c_wrapper(update_orbs_c_funptr, obj_func_c_funptr, &
                                   n_param_c, settings_c) result(error_c) &
        bind(C, name="mock_solver")
        !
        ! this subroutine is a mock routine for the solver C wrapper subroutine
        !
        use c_interface, only: solver_settings_type_c, update_orbs_c_type, &
                               hess_x_c_type, obj_func_c_type, precond_c_type, &
                               conv_check_c_type, logger_c_type
        use test_reference, only: tol_c, operator(/=)

        type(c_funptr), intent(in), value :: update_orbs_c_funptr, obj_func_c_funptr
        integer(c_ip), intent(in), value :: n_param_c
        type(solver_settings_type_c), intent(in), value :: settings_c
        integer(c_ip) :: error_c

        type(c_funptr) :: hess_x_c_funptr
        procedure(update_orbs_c_type), pointer :: update_orbs_funptr
        procedure(hess_x_c_type), pointer :: hess_x_funptr
        procedure(obj_func_c_type), pointer :: obj_func_funptr
        real(c_rp) :: func
        real(c_rp), allocatable :: kappa(:), grad(:), h_diag(:), x(:), hess_x(:), &
                                   residual(:), precond_residual(:)
        procedure(precond_c_type), pointer :: precond_funptr
        procedure(conv_check_c_type), pointer :: conv_check_funptr
        procedure(logger_c_type), pointer :: logger_funptr
        character(:), allocatable, target :: message
        logical(c_bool) :: converged
        integer(c_ip) :: error

        ! get Fortran pointer to passed orbital update routine and call it
        call c_f_procpointer(cptr=update_orbs_c_funptr, fptr=update_orbs_funptr)
        allocate(kappa(n_param_c), grad(n_param_c), h_diag(n_param_c))
        kappa = 1.0_c_rp
        error = update_orbs_funptr(kappa, func, grad, h_diag, hess_x_c_funptr)

        ! check for error
        if (error /= 0) then
            write (stderr, *) "test_solver_py_interface failed: Passed orbital "// &
                "updating function returned error."
            test_solver_interface = .false.
        end if

        ! get Fortran pointer to Hessian linear transformation function and call it
        call c_f_procpointer(cptr=hess_x_c_funptr, fptr=hess_x_funptr)
        allocate(x(n_param_c), hess_x(n_param_c))
        x = 1.0_c_rp
        error = hess_x_funptr(x, hess_x)

        ! check for error
        if (error /= 0) then
            write (stderr, *) "test_solver_py_interface failed: Passed Hessian "// &
                "linear transformation function returned error."
            test_solver_interface = .false.
        end if

        ! get Fortran pointers to generated arrays and check whether these are filled
        ! correctly
        if (abs(func - 3.0_c_rp) > tol_c) then
            write (stderr, *) "test_solver_py_interface failed: Returned function "// &
                "value of passed orbital updating function wrong."
            test_solver_interface = .false.
        end if
        if (any(abs(grad - 2.0_c_rp) > tol_c)) then
            write (stderr, *) "test_solver_py_interface failed: Returned gradient "// &
                "of passed orbital updating function wrong."
            test_solver_interface = .false.
        end if
        if (any(abs(h_diag - 3.0_c_rp) > tol_c)) then
            write (stderr, *) "test_solver_py_interface failed: Returned Hessian "// &
                "diagonal of passed orbital updating function wrong."
            test_solver_interface = .false.
        end if
        if (any(abs(hess_x - 4.0_c_rp) > tol_c)) then
            write (stderr, *) "test_solver_py_interface failed: Returned Hessian "// &
                "linear transformation from Hessian linear transformation function "// &
                "of passed orbital updating function wrong."
            test_solver_interface = .false.
        end if
        deallocate(grad, h_diag, x, hess_x)

        ! get Fortran pointer to passed objective function and call it
        call c_f_procpointer(cptr=obj_func_c_funptr, fptr=obj_func_funptr)
        error = obj_func_funptr(kappa, func)

        ! check for error
        if (error /= 0) then
            write (stderr, *) "test_solver_py_interface failed: Passed objective "// &
                "function returned error."
            test_solver_interface = .false.
        end if

        ! check if generated function value is correct
        if (abs(func - 3.0_c_rp) > tol_c) then
            write (stderr, *) "test_solver_py_interface failed: Returned function "// &
                "value of passed objective function wrong."
            test_solver_interface = .false.
        end if
        deallocate(kappa)

        ! get Fortran pointer to passed preconditioner function, call it and check
        ! result
        call c_f_procpointer(cptr=settings_c%precond, fptr=precond_funptr)
        allocate(residual(n_param_c), precond_residual(n_param_c))
        residual = 1.0_c_rp
        error = precond_funptr(residual, 5.0_c_rp, precond_residual)
        if (error /= 0) then
            write (stderr, *) "test_solver_py_interface failed: Passed "// &
                "preconditioner function returned error."
            test_solver_interface = .false.
        end if
        if (any(abs(precond_residual - 5.0_c_rp) > tol_c)) then
            write (stderr, *) "test_solver_py_interface failed: Returned "// &
                "preconditioned residual from preconditioner function wrong."
            test_solver_interface = .false.
        end if
        deallocate(residual, precond_residual)

        ! get Fortran pointer to passed convergence check function, call it and check 
        ! result
        call c_f_procpointer(cptr=settings_c%conv_check, fptr=conv_check_funptr)
        error = conv_check_funptr(converged)
        if (error /= 0) then
            write (stderr, *) "test_solver_py_interface failed: Passed convergence "// &
                "check function returned error."
            test_solver_interface = .false.
        end if
        if (.not. converged) then
            write (stderr, *) "test_solver_py_interface failed: Returned "// &
                "convergence logical of convergence check function wrong."
            test_solver_interface = .false.
        end if

        ! get Fortran pointer to passed logging function and call it
        message = "test" // c_null_char
        call c_f_procpointer(cptr=settings_c%logger, fptr=logger_funptr)
        call logger_funptr(message)

        ! check optional settings against reference values
        if (settings_c /= ref_settings_c) then
            write (stderr, *) "test_solver_py_interface failed: Passed settings "// &
                "associated with wrong values."
            test_solver_interface = .false.
        end if

        ! set return arguments
        error_c = 0

    end function mock_solver_c_wrapper

    function mock_stability_check_c_wrapper(h_diag_c_ptr, hess_x_c_funptr, n_param_c, &
                                            stable_c, settings_c, kappa_c_ptr) &
        result(error_c) bind(C, name="mock_stability_check")
        !
        ! this subroutine is a mock routine for the stability check C wrapper
        ! subroutine
        !
        use c_interface, only: stability_settings_type_c, hess_x_c_type, &
                               precond_c_type, logger_c_type
        use test_reference, only: tol_c, operator(/=)

        type(c_ptr), intent(in), value :: h_diag_c_ptr, kappa_c_ptr
        type(c_funptr), intent(in), value :: hess_x_c_funptr
        integer(c_ip), intent(in), value :: n_param_c
        logical(c_bool), intent(out) :: stable_c
        type(stability_settings_type_c), intent(in), value :: settings_c
        integer(c_ip) :: error_c

        real(c_rp), pointer :: h_diag_ptr(:), kappa_ptr(:)
        procedure(hess_x_c_type), pointer :: hess_x_funptr
        procedure(precond_c_type), pointer :: precond_funptr
        real(c_rp), allocatable :: x(:), hess_x(:), residual(:), precond_residual(:)
        procedure(logger_c_type), pointer :: logger_funptr
        character(:), allocatable, target :: message
        integer(c_ip) :: error

        ! check if Hessian diagonal is passed correctly
        call c_f_pointer(h_diag_c_ptr, h_diag_ptr, [n_param_c])
        if (any(abs(h_diag_ptr - 3.0_c_rp) > tol_c)) then
            write (stderr, *) "test_stability_check_py_interface failed: Passed "// &
                "Hessian diagonal wrong."
            test_stability_check_interface = .false.
        end if

        ! get Fortran pointer to Hessian linear transformation function, call it and 
        ! check result
        call c_f_procpointer(cptr=hess_x_c_funptr, fptr=hess_x_funptr)
        allocate(x(n_param_c), hess_x(n_param_c))
        x = 1.0_c_rp
        error = hess_x_funptr(x, hess_x)
        if (error /= 0) then
            write (stderr, *) "test_stability_check_py_interface failed: Passed "// &
                "Hessian linear transformation function returned error."
            test_stability_check_interface = .false.
        end if
        if (any(abs(hess_x - 4.0_c_rp) > tol_c)) then
            write (stderr, *) "test_stability_check_py_interface failed: Returned "// &
                "Hessian linear transformation from passed Hessian linear "// &
                "transformation function wrong."
            test_stability_check_interface = .false.
        end if
        deallocate(x, hess_x)

        ! get Fortran pointer to passed precond function, call it and check result
        call c_f_procpointer(cptr=settings_c%precond, fptr=precond_funptr)
        allocate(residual(n_param_c), precond_residual(n_param_c))
        residual = 1.0_c_rp
        error = precond_funptr(residual, 5.0_c_rp, precond_residual)
        if (error /= 0) then
            write (stderr, *) "test_stability_check_py_interface failed: Passed "// &
                "preconditioner function returned error."
            test_stability_check_interface = .false.
        end if
        if (any(abs(precond_residual - 5.0_c_rp) > tol_c)) then
            write (stderr, *) "test_stability_check_py_interface failed: Returned "// &
                "preconditioned residual from preconditioner function wrong."
            test_stability_check_interface = .false.
        end if
        deallocate(residual, precond_residual)

        ! get Fortran pointer to passed logging function and call it
        message = "test" // c_null_char
        call c_f_procpointer(cptr=settings_c%logger, fptr=logger_funptr)
        call logger_funptr(message)

        ! check optional settings against reference values
        if (settings_c /= ref_settings_c) then
            write (stderr, *) "test_stability_check_py_interface failed: Passed "// &
                "settings associated with wrong values."
            test_stability_check_interface = .false.
        end if

        ! set return arguments
        stable_c = .false.
        if (c_associated(kappa_c_ptr)) then
            call c_f_pointer(kappa_c_ptr, kappa_ptr, [n_param_c])
            kappa_ptr = 1.0_c_rp
        end if
        error_c = 0

    end function mock_stability_check_c_wrapper

    subroutine mock_init_solver_settings_c(settings) &
        bind(C, name="mock_init_solver_settings")
        !
        ! this subroutine is a mock routine for the C solver setting initialization 
        ! subroutine
        !
        use c_interface, only: solver_settings_type_c
        use test_reference, only: assignment(=)

        type(solver_settings_type_c), intent(inout) :: settings

        ! set reference values
        settings = ref_settings_c

    end subroutine mock_init_solver_settings_c

    subroutine mock_init_stability_settings_c(settings) &
        bind(C, name="mock_init_stability_settings")
        !
        ! this subroutine is a mock routine for the C stability check setting 
        ! initialization subroutine
        !
        use c_interface, only: stability_settings_type_c
        use test_reference, only: assignment(=)

        type(stability_settings_type_c), intent(inout) :: settings

        ! set reference values
        settings = ref_settings_c

    end subroutine mock_init_stability_settings_c

end module c_interface_mock
