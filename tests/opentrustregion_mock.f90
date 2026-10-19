! Copyright (C) 2025- Jonas Greiner
!
! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.

module opentrustregion_mock

    use opentrustregion, only: rp, ip, stderr, solver, stability_check, &
                               solver_settings_type, stability_settings_type, &
                               solver_state_type, stability_state_type
    use c_interface, only: c_callbacks_type
    use test_reference, only: tol

    implicit none

    logical :: test_passed

    ! error code the mock solver and stability check return and states they write,
    ! whose fields hold values distinct from each other within their kind and from the
    ! same field of the reference states, so that a test can check the C wrappers hand
    ! them back
    integer(ip), parameter :: mock_error = 7
    type(stability_state_type), parameter :: mock_nested_state = stability_state_type( &
        iteration=7, min_eigval=-0.5_rp, residual_rms=0.25_rp, n_hess_x=2)
    type(solver_state_type), parameter :: mock_solver_state = solver_state_type( &
        macro_iteration=4, func=-1.5_rp, grad_rms=0.75_rp, trust_radius=0.5_rp, &
        micro_iteration=6, jacobi_davidson_iteration=1, relative_residual=0.375_rp, &
        level_shift=2.5_rp, step_size=0.125_rp, n_update_orbs=3, n_obj_func=8, &
        n_hess_x=5, max_precision_reached=.false., stability_state=mock_nested_state)
    type(stability_state_type), parameter :: mock_stability_state = &
        stability_state_type(iteration=9, min_eigval=-0.75_rp, residual_rms=0.0625_rp, &
                             n_hess_x=4)

    ! stability the mock stability check returns, which a test sets
    logical :: mock_stable = .true.

    ! settings and callback bundles the mock solver and stability check received,
    ! the bundle of the settings and that of the nested settings, so that a test can
    ! check what the C wrapper handed over
    type(solver_settings_type) :: received_solver_settings
    type(stability_settings_type) :: received_stability_settings
    type(c_callbacks_type) :: received_callbacks, received_stability_callbacks

    ! whether the callback bundles the mock solver and stability check received refer
    ! to the state record they were given, so that a test can check that the C
    ! wrappers copy this record into the host's
    logical :: solver_state_bundled, stability_state_bundled

    ! create function pointers to ensure that routines comply with interface
    procedure(solver), pointer :: mock_solver_ptr => mock_solver
    procedure(stability_check), pointer :: mock_stability_check_ptr => &
        mock_stability_check

contains

    subroutine mock_solver(update_orbs_funptr, obj_func_funptr, n_param, error, &
                           settings, state)
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
        type(solver_settings_type), intent(in), optional :: settings
        type(solver_state_type), intent(out), target, optional :: state

        ! assume tests pass
        test_passed = .true.

        ! set error code and check that the C wrapper passes the settings it converted
        error = mock_error
        if (.not. present(settings)) then
            write(stderr, *) "test_solver_c_wrapper failed: Settings not passed."
            test_passed = .false.
            return
        end if

        ! test passed orbital update subroutine
        if (.not. check_update_orbs_funptr(update_orbs_funptr, "solver_c_wrapper", &
                                           " by given orbital updating subroutine", &
                                           settings%context)) test_passed = .false.

        ! test passed objective function
        if (.not. check_obj_func_funptr(obj_func_funptr, "solver_c_wrapper", &
                                        " by given objective function", &
                                        settings%context)) test_passed = .false.

        ! check number of parameters
        if (n_param /= n_param_ref) then
            write(stderr, *) "test_solver_c_wrapper failed: Passed number of "// &
                "parameters wrong."
            test_passed = .false.
        end if

        ! check the optional callback functions that were passed, which of them are
        ! passed the test checks
        if (associated(settings%precond)) then
            if (.not. check_precond_funptr(settings%precond, "solver_c_wrapper", &
                                           " by given preconditioner subroutine", &
                                           settings%context)) test_passed = .false.
        end if
        if (associated(settings%project)) then
            if (.not. check_project_funptr(settings%project, "solver_c_wrapper", &
                                           " by given projection subroutine", &
                                           settings%context)) test_passed = .false.
        end if
        if (associated(settings%conv_check)) then
            if (.not. check_conv_check_funptr(settings%conv_check, "solver_c_wrapper", &
                                              " by given convergence check function", &
                                              settings%context)) test_passed = .false.
        end if
        if (associated(settings%logger)) call settings%logger("test", settings%context)

        ! record the settings and the callback bundles of the settings and the nested
        ! settings
        received_solver_settings = settings
        received_callbacks = c_callbacks_type()
        select type (callbacks => settings%context)
        type is (c_callbacks_type)
            received_callbacks = callbacks
        end select
        received_stability_callbacks = c_callbacks_type()
        select type (callbacks => settings%stability_settings%context)
        type is (c_callbacks_type)
            received_stability_callbacks = callbacks
        end select

        ! record whether the callback bundles of the settings and the nested settings
        ! both refer to the state record
        solver_state_bundled = .false.
        if (present(state)) then
            select type (callbacks => settings%context)
            type is (c_callbacks_type)
                solver_state_bundled = associated(callbacks%solver_state, state)
            end select
            select type (callbacks => settings%stability_settings%context)
            type is (c_callbacks_type)
                solver_state_bundled = solver_state_bundled .and. &
                                       associated(callbacks%solver_state, state)
            class default
                solver_state_bundled = .false.
            end select
        end if

        ! set state
        if (present(state)) state = mock_solver_state

    end subroutine mock_solver

    subroutine mock_stability_check(h_diag, hess_x_funptr, stable, error, settings, &
                                    state, kappa)
        !
        ! this subroutine is a mock routine for the stability check to test the C
        ! interface
        !
        use opentrustregion, only: hess_x_type
        use test_reference, only: check_hess_x_funptr, check_precond_funptr, &
                                  check_project_funptr, n_param

        real(rp), intent(in) :: h_diag(:)
        procedure(hess_x_type), intent(in), pointer :: hess_x_funptr
        logical, intent(out) :: stable
        integer(ip), intent(out) :: error
        type(stability_settings_type), intent(in), optional :: settings
        type(stability_state_type), intent(out), target, optional :: state
        real(rp), intent(out), optional :: kappa(:)

        ! assume tests pass
        test_passed = .true.

        ! set error code and stability and check that the C wrapper passes the settings
        ! it converted
        error = mock_error
        stable = mock_stable
        if (.not. present(settings)) then
            write(stderr, *) "test_stability_check_c_wrapper failed: Settings not "// &
                "passed."
            test_passed = .false.
            return
        end if

        ! check Hessian diagonal
        if (size(h_diag) /= n_param) then
            write(stderr, *) "test_stability_check_c_wrapper failed: Passed "// &
                "Hessian diagonal has wrong size."
            test_passed = .false.
        else if (any(abs(h_diag - 3.0_rp) > tol)) then
            write(stderr, *) "test_stability_check_c_wrapper failed: Passed "// &
                "Hessian diagonal wrong."
            test_passed = .false.
        end if

        ! test passed Hessian linear transformation subroutine
        if (.not. check_hess_x_funptr( &
            hess_x_funptr, "stability_check_c_wrapper", " by given Hessian linear "// &
            "transformation subroutine", settings%context)) test_passed = .false.

        ! check the optional callback functions that were passed, which of them are
        ! passed the test checks
        if (associated(settings%precond)) then
            if (.not. check_precond_funptr( &
                settings%precond, "stability_check_c_wrapper", " by given "// &
                "preconditioner subroutine", settings%context)) test_passed = .false.
        end if
        if (associated(settings%project)) then
            if (.not. check_project_funptr( &
                settings%project, "stability_check_c_wrapper", " by given "// &
                "projection subroutine", settings%context)) test_passed = .false.
        end if
        if (associated(settings%logger)) call settings%logger("test", settings%context)

        ! record the settings and their callback bundle
        received_stability_settings = settings
        received_callbacks = c_callbacks_type()
        select type (callbacks => settings%context)
        type is (c_callbacks_type)
            received_callbacks = callbacks
        end select

        ! record whether the callback bundle refers to the state record
        stability_state_bundled = .false.
        if (present(state)) then
            select type (callbacks => settings%context)
            type is (c_callbacks_type)
                stability_state_bundled = associated(callbacks%stability_state, state)
            end select
        end if

        ! set optional output quantities and state
        if (present(kappa)) kappa = 1.0_rp
        if (present(state)) state = mock_stability_state

    end subroutine mock_stability_check

end module opentrustregion_mock
