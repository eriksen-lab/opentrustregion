! Copyright (C) 2025- Jonas Greiner
!
! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.

module c_interface_mock

    use opentrustregion, only: stderr
    use c_interface, only: c_rp, c_ip, solver_c_wrapper, stability_check_c_wrapper, &
                           init_solver_settings_c, init_stability_settings_c
    use test_reference, only: ref_settings, n_param
    use, intrinsic :: iso_c_binding, only: c_bool, c_ptr, c_funptr, c_f_pointer, &
                                           c_f_procpointer, c_associated, c_null_char

    implicit none

    logical(c_bool), bind(C) :: test_solver_interface = .true., &
                                test_stability_check_interface = .true.

    ! create function pointers to ensure that routines comply with interface
    procedure(solver_c_wrapper), pointer :: mock_solver_c_wrapper_ptr => &
        mock_solver_c_wrapper
    procedure(stability_check_c_wrapper), pointer :: &
        mock_stability_check_c_wrapper_ptr => mock_stability_check_c_wrapper
    procedure(init_solver_settings_c), pointer :: mock_init_solver_settings_c_ptr &
        => mock_init_solver_settings_c
    procedure(init_stability_settings_c), pointer :: &
        mock_init_stability_settings_c_ptr => mock_init_stability_settings_c

contains

    function mock_solver_c_wrapper(update_orbs_c_funptr, obj_func_c_funptr, &
                                   n_param_c, settings_c) result(error_c) &
        bind(C, name="mock_solver")
        !
        ! this subroutine is a mock routine for the solver C wrapper subroutine
        !
        use c_interface, only: solver_settings_type_c, update_orbs_c_type, &
                               hess_x_c_type, obj_func_c_type, precond_c_type, &
                               project_c_type, conv_check_c_type, logger_c_type
        use test_reference, only: test_update_orbs_c_funptr, test_obj_func_c_funptr, &
                                  test_precond_c_funptr, test_project_c_funptr, &
                                  test_conv_check_c_funptr, operator(/=)

        type(c_funptr), intent(in), value :: update_orbs_c_funptr, obj_func_c_funptr
        integer(c_ip), intent(in), value :: n_param_c
        type(solver_settings_type_c), intent(in), value :: settings_c
        integer(c_ip) :: error_c

        procedure(logger_c_type), pointer :: logger_funptr
        character(:), allocatable, target :: message

        ! test passed orbital update function
        test_solver_interface = test_solver_interface .and. &
            test_update_orbs_c_funptr(update_orbs_c_funptr, "solver_py_interface", &
                                      " by given orbital updating function")

        ! test passed objective function
        test_solver_interface = test_solver_interface .and. &
            test_obj_func_c_funptr(obj_func_c_funptr, "solver_py_interface", &
                                   " by given objective function")

        ! check if passed number of parameters is correct
        if (n_param_c /= 3) then
            write (stderr, *) "test_solver_py_interface failed: Passed number of "// &
                "parameters wrong."
            test_solver_interface = .false.
        end if

        ! test passed preconditioner function
        test_solver_interface = test_solver_interface .and. &
            test_precond_c_funptr(settings_c%precond, "solver_py_interface", &
                                  " by given preconditioning function")

        ! test passed projection function
        test_solver_interface = test_solver_interface .and. &
            test_project_c_funptr(settings_c%project, "solver_py_interface", &
                                  " by given projection function")

        ! test passed convergence check function
        test_solver_interface = test_solver_interface .and. &
            test_conv_check_c_funptr(settings_c%conv_check, "solver_py_interface", &
                                     " by given convergence check function")

        ! get Fortran pointer to passed logging function and call it
        message = "test" // c_null_char
        call c_f_procpointer(cptr=settings_c%logger, fptr=logger_funptr)
        call logger_funptr(message)

        ! check optional settings against reference values
        if (settings_c /= ref_settings) then
            write (stderr, *) "test_solver_py_interface failed: Passed settings "// &
                "associated with wrong values."
            test_solver_interface = .false.
        end if

        ! set return arguments
        error_c = 0

    end function mock_solver_c_wrapper

    function mock_stability_check_c_wrapper(h_diag_c, hess_x_c_funptr, n_param_c, &
                                            stable_c, settings_c, kappa_c_ptr) &
        result(error_c) bind(C, name="mock_stability_check")
        !
        ! this subroutine is a mock routine for the stability check C wrapper
        ! subroutine
        !
        use c_interface, only: stability_settings_type_c, hess_x_c_type, &
                               precond_c_type, project_c_type, logger_c_type
        use test_reference, only: tol_c, test_hess_x_c_funptr, test_precond_c_funptr, &
                                  test_project_c_funptr, operator(/=)

        real(c_rp), intent(in), target :: h_diag_c(*)
        type(c_funptr), intent(in), value :: hess_x_c_funptr
        integer(c_ip), intent(in), value :: n_param_c
        logical(c_bool), intent(out) :: stable_c
        type(stability_settings_type_c), intent(in), value :: settings_c
        type(c_ptr), intent(in), value :: kappa_c_ptr
        integer(c_ip) :: error_c

        real(c_rp), pointer :: kappa_ptr(:)
        procedure(logger_c_type), pointer :: logger_funptr
        character(:), allocatable, target :: message

        ! check if Hessian diagonal is passed correctly
        if (any(abs(h_diag_c(:n_param_c) - 3.0_c_rp) > tol_c)) then
            write (stderr, *) "test_stability_check_py_interface failed: Passed "// &
                "Hessian diagonal wrong."
            test_stability_check_interface = .false.
        end if

        ! test passed Hessian linear transformation
        test_stability_check_interface = test_stability_check_interface .and. &
            test_hess_x_c_funptr(hess_x_c_funptr, "stability_check_py_interface", &
                                 " by given Hessian linear transformation function")

        ! check if passed number of parameters is correct
        if (n_param_c /= 3) then
            write (stderr, *) "test_stability_check_py_interface failed: Passed "// &
                "number of parameters wrong."
            test_stability_check_interface = .false.
        end if

        ! test passed preconditioner
        test_stability_check_interface = test_stability_check_interface .and. &
            test_precond_c_funptr(settings_c%precond, "stability_check_py_interface", &
                                  " by given preconditioning function")

        ! test passed projection function
        test_stability_check_interface = test_stability_check_interface .and. &
            test_project_c_funptr(settings_c%project, "stability_check_py_interface", &
                                  " by given projection function")

        ! get Fortran pointer to passed logging function and call it
        message = "test" // c_null_char
        call c_f_procpointer(cptr=settings_c%logger, fptr=logger_funptr)
        call logger_funptr(message)

        ! check optional settings against reference values
        if (settings_c /= ref_settings) then
            write (stderr, *) "test_stability_check_py_interface failed: Passed "// &
                "settings associated with wrong values."
            test_stability_check_interface = .false.
        end if

        ! set return arguments
        stable_c = .false.
        if (c_associated(kappa_c_ptr)) then
            call c_f_pointer(kappa_c_ptr, kappa_ptr, [n_param])
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
        settings = ref_settings

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
        settings = ref_settings

    end subroutine mock_init_stability_settings_c

end module c_interface_mock
