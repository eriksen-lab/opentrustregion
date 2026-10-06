! Copyright (C) 2025- Jonas Greiner
!
! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.

module opentrustregion_mock

    use opentrustregion, only: rp, ip, stderr, solver, stability_check, &
                               solver_settings_type
    use c_interface, only: c_callbacks_type
    use test_reference, only: tol, ref_stability_settings, operator(/=)

    implicit none

    logical :: test_passed

    ! output values the mock solver and stability check return or write into the
    ! settings object, so that a test can check the C wrappers hand them back
    integer(ip), parameter :: mock_error = 7, mock_n_update_orbs = 3, &
                              mock_n_hess_x = 5, mock_stability_n_hess_x = 2, &
                              mock_stability_check_n_hess_x = 4

    ! number of calls the mock solver makes to the optional callback functions of
    ! the internal stability check, the preconditioner, projection and logging
    ! function, with the context of the nested settings
    integer(ip), parameter :: mock_solver_n_nested_calls = 3

    ! settings and callback bundle of the nested settings the mock solver received,
    ! so that a test can check what the C wrapper handed over
    type(solver_settings_type) :: received_solver_settings
    type(c_callbacks_type) :: received_stability_callbacks

    ! create function pointers to ensure that routines comply with interface
    procedure(solver), pointer :: mock_solver_ptr => mock_solver
    procedure(stability_check), pointer :: mock_stability_check_ptr => &
        mock_stability_check

contains

    subroutine mock_solver(update_orbs_funptr, obj_func_funptr, n_param, error, &
                           settings)
        !
        ! this subroutine is a mock routine for solver to test the C interface
        !
        use opentrustregion, only: update_orbs_type, obj_func_type
        use test_reference, only: check_update_orbs_funptr, check_obj_func_funptr, &
                                  check_precond_funptr, check_project_funptr, &
                                  check_conv_check_funptr, n_param_ref => n_param

        procedure(update_orbs_type), intent(in), pointer :: update_orbs_funptr
        procedure(obj_func_type), intent(in), pointer :: obj_func_funptr
        integer(ip), intent(in) :: n_param
        integer(ip), intent(out) :: error
        type(solver_settings_type), intent(inout) :: settings

        ! initialize logical
        test_passed = .true.

        ! test passed orbital update subroutine
        test_passed = test_passed .and. check_update_orbs_funptr( &
            update_orbs_funptr, "solver_c_wrapper", &
            " by given orbital updating subroutine", settings%context)

        ! test passed objective function
        test_passed = test_passed .and. check_obj_func_funptr( &
            obj_func_funptr, "solver_c_wrapper", " by given objective function", &
            settings%context)

        ! check number of parameters
        if (n_param /= n_param_ref) then
            test_passed = .false.
            write(stderr, *) "test_solver_c_wrapper failed: Passed number of "// &
                "parameters wrong."
        end if

        ! check if optional preconditioner subroutine is correctly passed
        if (.not. associated(settings%precond)) then
            test_passed = .false.
            write(stderr, *) "test_solver_c_wrapper failed: Passed preconditioner "// &
                "function not associated with value."
        else
            test_passed = test_passed .and. check_precond_funptr( &
                settings%precond, "solver_c_wrapper", &
                " by given preconditioner subroutine", settings%context)
        end if

        ! check if optional projection subroutine is correctly passed
        if (.not. associated(settings%project)) then
            test_passed = .false.
            write(stderr, *) "test_solver_c_wrapper failed: Passed projection "// &
                "function not associated with value."
        else
            test_passed = test_passed .and. check_project_funptr( &
                settings%project, "solver_c_wrapper", &
                " by given projection subroutine", settings%context)
        end if

        ! check if optional convergence check function is correctly passed
        if (.not. associated(settings%conv_check)) then
            test_passed = .false.
            write(stderr, *) "test_solver_c_wrapper failed: Passed convergence "// &
                "check function not associated with value."
        else
            test_passed = test_passed .and. check_conv_check_funptr( &
                settings%conv_check, "solver_c_wrapper", &
                " by given convergence check function", settings%context)
        end if

        ! check if optional logging function is correctly passed
        if (.not. associated(settings%logger)) then
            test_passed = .false.
            write(stderr, *) "test_solver_c_wrapper failed: Passed logging "// &
                "function not associated with value."
        else
            call settings%logger("test", settings%context)
        end if

        ! call the optional callback functions as the internal stability check does
        ! with the nested settings' context, whose callback bundle decides whether the
        ! solver's C functions or the nested settings' own are called
        if (associated(settings%precond)) then
            test_passed = test_passed .and. check_precond_funptr( &
                settings%precond, "solver_c_wrapper", " by preconditioner of the "// &
                "internal stability check", settings%stability_settings%context)
        end if
        if (associated(settings%project)) then
            test_passed = test_passed .and. check_project_funptr( &
                settings%project, "solver_c_wrapper", " by projection of the "// &
                "internal stability check", settings%stability_settings%context)
        end if
        if (associated(settings%logger)) &
            call settings%logger("test", settings%stability_settings%context)

        ! record the settings and the callback bundle of the nested settings
        received_solver_settings = settings
        received_stability_callbacks = c_callbacks_type()
        select type (callbacks => settings%stability_settings%context)
        type is (c_callbacks_type)
            received_stability_callbacks = callbacks
        end select

        ! set output quantities and fields
        error = mock_error
        settings%max_precision_reached = .false.
        settings%n_update_orbs = mock_n_update_orbs
        settings%n_hess_x = mock_n_hess_x
        settings%stability_settings%n_hess_x = mock_stability_n_hess_x

    end subroutine mock_solver

    subroutine mock_stability_check(h_diag, hess_x_funptr, stable, error, settings, &
                                    kappa)
        !
        ! this subroutine is a mock routine for the stability check to test the C
        ! interface
        !
        use opentrustregion, only: stability_settings_type, hess_x_type
        use test_reference, only: check_hess_x_funptr, check_precond_funptr, &
                                  check_project_funptr, n_param

        real(rp), intent(in) :: h_diag(:)
        procedure(hess_x_type), intent(in), pointer :: hess_x_funptr
        logical, intent(out) :: stable
        integer(ip), intent(out) :: error
        type(stability_settings_type), intent(inout) :: settings
        real(rp), intent(out), optional :: kappa(:)

        ! initialize logical
        test_passed = .true.

        ! check Hessian diagonal
        if (size(h_diag) /= n_param) then
            test_passed = .false.
            write(stderr, *) "test_stability_check_c_wrapper failed: Passed "// &
                "Hessian diagonal has wrong size."
        else if (any(abs(h_diag - 3.0_rp) > tol)) then
            test_passed = .false.
            write(stderr, *) "test_stability_check_c_wrapper failed: Passed "// &
                "Hessian diagonal wrong."
        end if

        ! test passed Hessian linear transformation subroutine
        test_passed = test_passed .and. check_hess_x_funptr( &
            hess_x_funptr, "stability_check_c_wrapper", &
            " by given Hessian linear transformation subroutine", settings%context)

        ! check if optional preconditioner subroutine is correctly passed
        if (.not. associated(settings%precond)) then
            test_passed = .false.
            write(stderr, *) "test_stability_check_c_wrapper failed: Passed "// &
                "preconditioner function not associated with value."
        else
            test_passed = test_passed .and. check_precond_funptr( &
                settings%precond, "stability_check_c_wrapper", &
                " by given preconditioner subroutine", settings%context)
        end if

        ! check if optional projection subroutine is correctly passed
        if (.not. associated(settings%project)) then
            test_passed = .false.
            write(stderr, *) "test_stability_check_c_wrapper failed: Passed "// &
                "projection function not associated with value."
        else
            test_passed = test_passed .and. check_project_funptr( &
                settings%project, "stability_check_c_wrapper", &
                " by given projection subroutine", settings%context)
        end if

        ! check if optional logging function is correctly passed
        if (.not. associated(settings%logger)) then
            test_passed = .false.
            write(stderr, *) "test_stability_check_c_wrapper failed: Passed "// &
                "logging function not associated with value."
        else
            call settings%logger("test", settings%context)
        end if

        ! check if optional settings are correctly passed
        if (settings /= ref_stability_settings) then
            test_passed = .false.
            write(stderr, *) "test_stability_check_c_wrapper failed: Passed "// &
                "optional settings associated with wrong values."
        end if

        ! set output quantities and fields
        stable = .true.
        if (present(kappa)) kappa = 1.0_rp
        error = mock_error
        settings%n_hess_x = mock_stability_check_n_hess_x

    end subroutine mock_stability_check

end module opentrustregion_mock
