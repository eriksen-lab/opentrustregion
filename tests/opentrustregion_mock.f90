! Copyright (C) 2025- Jonas Greiner
!
! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.

module opentrustregion_mock

    use opentrustregion, only: rp, ip, stderr, solver_type, stability_check_type, &
                               update_orbs_type, hess_x_type, obj_func_type
    use test_reference, only: tol, ref_settings, operator(/=)

    implicit none

    logical :: test_passed

    ! create function pointers to ensure that routines comply with interface
    procedure(solver_type), pointer :: mock_solver_ptr => mock_solver
    procedure(stability_check_type), pointer :: mock_stability_check_ptr => &
        mock_stability_check

contains

    subroutine mock_solver(update_orbs_funptr, obj_func_funptr, n_param, error, &
                           settings)
        !
        ! this subroutine is a mock routine for solver to test the C interface
        !
        use opentrustregion, only: solver_settings_type

        procedure(update_orbs_type), intent(in), pointer :: update_orbs_funptr
        procedure(obj_func_type), intent(in), pointer :: obj_func_funptr
        integer(ip), intent(in) :: n_param
        integer(ip), intent(out) :: error
        type(solver_settings_type), intent(inout) :: settings

        real(rp), allocatable :: kappa(:), x(:), grad(:), h_diag(:), hess_x(:), &
                                 residual(:), precond_residual(:)
        real(rp) :: func
        procedure(hess_x_type), pointer :: hess_x_funptr
        logical :: converged

        ! initialize logical
        test_passed = .true.

        ! check if passed orbital updating subroutine produces correct quantities
        allocate(kappa(n_param), grad(n_param), h_diag(n_param))
        kappa = 1.0_rp
        call update_orbs_funptr(kappa, func, grad, h_diag, hess_x_funptr, error)
        if (error /= 0) then
            test_passed = .false.
            write (stderr, *) "test_solver_c_wrapper failed: Passed orbital "// &
                "updating function produced error."
        end if
        if (abs(func - 3.0_rp) > tol) then
            test_passed = .false.
            write (stderr, *) "test_solver_c_wrapper failed: Returned objective "// &
                "function value of passed orbital updating function wrong."
        end if
        if (any(abs(grad - 2.0_rp) > tol)) then
            test_passed = .false.
            write (stderr, *) "test_solver_c_wrapper failed: Returned gradient of "// &
                "passed orbital updating function wrong."
        end if
        if (any(abs(h_diag - 3.0_rp) > tol)) then
            test_passed = .false.
            write (stderr, *) "test_solver_c_wrapper failed: Returned Hessian "// &
                "diagonal of passed orbital updating function wrong."
        end if
        deallocate(grad, h_diag)
        allocate(x(n_param), hess_x(n_param))
        x = 1.0_rp
        call hess_x_funptr(x, hess_x, error)
        if (error /= 0) then
            test_passed = .false.
            write (stderr, *) "test_solver_c_wrapper failed: Hessian linear "// &
                "transformation function returned by passed orbital updating "// &
                "function produced error."
        end if
        if (any(abs(hess_x - 4.0_rp) > tol)) then
            test_passed = .false.
            write (stderr, *) "test_solver_c_wrapper failed: Returned Hessian "// &
                "linear transformation of Hessian linear transformation function "// &
                "returned by passed orbital updating function wrong."
        end if
        deallocate(x, hess_x)

        ! check if passed objective function produces correct quantities
        func = obj_func_funptr(kappa, error)
        if (error /= 0) then
            test_passed = .false.
            write (stderr, *) "test_solver_c_wrapper failed: Passed objective "// &
                "function produced error."
        end if
        if (abs(func - 3.0_rp) > tol) then
            test_passed = .false.
            write (stderr, *) "test_solver_c_wrapper failed: Returned objective "// &
                "function value of passed objective function wrong."
        end if
        deallocate(kappa)

        ! check number of parameters
        if (n_param /= 3) then
            test_passed = .false.
            write (stderr, *) "test_solver_c_wrapper failed: Passed number of "// &
                "parameters wrong."
        end if

        ! set output quantities
        error = 0

        ! check if optional preconditioner function is correctly passed
        if (.not. associated(settings%precond)) then
            test_passed = .false.
            write (stderr, *) "test_solver_c_wrapper failed: Passed preconditioner "// &
                "function not associated with value."
        else
            allocate(residual(n_param), precond_residual(n_param))
            residual = 1.0_rp
            call settings%precond(residual, 5.0_rp, precond_residual, error)
            if (error /= 0) then
                test_passed = .false.
                write (stderr, *) "test_solver_c_wrapper failed: Passed "// &
                    "preconditioner function produced error."
            end if
            if (any(abs(precond_residual - 5.0_rp) > tol)) then
                test_passed = .false.
                write (stderr, *) "test_solver_c_wrapper failed: Returned "// &
                    "preconditioner of passed preconditioner function wrong."
            end if
            deallocate(residual, precond_residual)
        end if

        ! check if optional convergence check function is correctly passed
        if (.not. associated(settings%conv_check)) then
            test_passed = .false.
            write (stderr, *) "test_solver_c_wrapper failed: Passed convergence "// &
                "check function not associated with value."
        else
            converged = settings%conv_check(error)
            if (error /= 0) then
                test_passed = .false.
                write (stderr, *) "test_solver_c_wrapper failed: Passed "// &
                    "convergence check function produced error."
            end if
            if (.not. converged) then
                test_passed = .false.
                write (stderr, *) "test_solver_c_wrapper failed: Returned "// &
                    "convergence logical of passed convergence check function wrong."
            end if
        end if

        ! check if optional logging function is correctly passed
        if (.not. associated(settings%logger)) then
            test_passed = .false.
            write (stderr, *) "test_solver_c_wrapper failed: Passed logging "// &
                "function not associated with value."
        else
            call settings%logger("test")
        end if

        ! check if optional settings are correctly passed
        if (settings /= ref_settings) then
            test_passed = .false.
            write (stderr, *) "test_solver_c_wrapper failed: Passed optional "// &
                "settings associated with wrong values."
        end if

    end subroutine mock_solver

    subroutine mock_stability_check(h_diag, hess_x_funptr, stable, error, settings, &
                                    kappa)
        !
        ! this subroutine is a mock routine for the stability check to test the C 
        ! interface
        !
        use opentrustregion, only: stability_settings_type

        real(rp), intent(in) :: h_diag(:)
        procedure(hess_x_type), intent(in), pointer :: hess_x_funptr
        logical, intent(out) :: stable
        integer(ip), intent(out) :: error
        type(stability_settings_type), intent(inout) :: settings
        real(rp), intent(out), optional :: kappa(:)

        real(rp), allocatable :: x(:), hess_x(:), residual(:), precond_residual(:)

        ! initialize logical
        test_passed = .true.

        ! check Hessian diagonal
        if (any(abs(h_diag - 3.0_rp) > tol)) then
            test_passed = .false.
            write (stderr, *) "test_stability_check_c_wrapper failed: Passed "// &
                "Hessian diagonal wrong."
        end if

        ! check if passed Hessian linear transformation function produces correct
        ! quantity
        allocate(x(size(h_diag)), hess_x(size(h_diag)))
        x = 1.0_rp
        call hess_x_funptr(x, hess_x, error)
        if (error /= 0) then
            test_passed = .false.
            write (stderr, *) "test_stability_check_c_wrapper failed: Passed "// &
                "Hessian linear transformation function produced error."
        end if
        if (any(abs(hess_x - 4.0_rp) > tol)) then
            test_passed = .false.
            write (stderr, *) "test_stability_check_c_wrapper failed: Hessian "// &
                "linear transformation returned by passed Hessian linear "// &
                "transformation function wrong."
        end if
        deallocate(x, hess_x)

        ! set output quantities
        stable = .false.
        if (present(kappa)) kappa = 1.0_rp
        error = 0

        ! check if optional preconditioner function is correctly passed
        if (.not. associated(settings%precond)) then
            test_passed = .false.
            write (stderr, *) "test_stability_check_c_wrapper failed: Passed "// &
                "preconditioner function not associated with value."
        else
            allocate(residual(size(h_diag)), precond_residual(size(h_diag)))
            residual = 1.0_rp
            call settings%precond(residual, 5.0_rp, precond_residual, error)
            if (error /= 0) then
                test_passed = .false.
                write (stderr, *) "test_stability_check_c_wrapper failed: Passed "// &
                    "preconditioner function produced error."
            end if
            if (any(abs(precond_residual - 5.0_rp) > tol)) then
                test_passed = .false.
                write (stderr, *) "test_stability_check_c_wrapper failed: Returned "// &
                    "preconditioner of passed preconditioner function wrong."
            end if
            deallocate(residual, precond_residual)
        end if

        ! check if optional logging function is correctly passed
        if (.not. associated(settings%logger)) then
            test_passed = .false.
            write (stderr, *) "test_stability_check_c_wrapper failed: Passed "// &
                "logging function not associated with value."
        else
            call settings%logger("test")
        end if

        ! check if optional settings are correctly passed
        if (settings /= ref_settings) then
            test_passed = .false.
            write (stderr, *) "test_stability_check_c_wrapper failed: Passed "// &
                "optional settings associated with wrong values."
        end if

    end subroutine mock_stability_check

end module opentrustregion_mock
