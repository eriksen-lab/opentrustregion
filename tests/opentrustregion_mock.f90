! Copyright (C) 2025- Jonas Greiner
!
! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.

module opentrustregion_mock

    use opentrustregion, only: rp, ip, stderr, update_orbs_type, hess_x_type, &
                               obj_func_type, precond_type, conv_check_type, logger_type

    implicit none

    real(rp), parameter :: tol = 1.d-10

    logical :: solver_default = .true., stability_check_default = .true., &
               test_passed

contains

    subroutine mock_solver(update_orbs_funptr, obj_func_funptr, n_param, error, &
                           precond_funptr, conv_check_funptr, stability, line_search, &
                           davidson, jacobi_davidson, prefer_jacobi_davidson, &
                           conv_tol, n_random_trial_vectors, start_trust_radius, &
                           n_macro, n_micro, global_red_factor, local_red_factor, &
                           seed, verbose, logger_funptr)
        !
        ! this subroutine is a mock routine for solver to test the C interface
        !
        use opentrustregion, only: solver_stability_default, &
                                   solver_line_search_default, &
                                   solver_davidson_default, &
                                   solver_jacobi_davidson_default, &
                                   solver_prefer_jacobi_davidson_default, &
                                   solver_conv_tol_default, &
                                   solver_n_random_trial_vectors_default, &
                                   solver_start_trust_radius_default, &
                                   solver_n_macro_default, solver_n_micro_default, &
                                   solver_global_red_factor_default, &
                                   solver_local_red_factor_default, &
                                   solver_seed_default, solver_verbose_default

        procedure(update_orbs_type), intent(in), pointer :: update_orbs_funptr
        procedure(obj_func_type), intent(in), pointer :: obj_func_funptr
        integer(ip), intent(in) :: n_param
        logical, intent(out) :: error
        procedure(precond_type), intent(in), pointer, optional :: precond_funptr
        procedure(conv_check_type), intent(in), pointer, optional :: conv_check_funptr
        logical, intent(in), optional :: stability, line_search, davidson, &
                                         jacobi_davidson, prefer_jacobi_davidson
        real(rp), intent(in), optional :: conv_tol, start_trust_radius, &
                                          global_red_factor, local_red_factor
        integer(ip), intent(in), optional :: n_random_trial_vectors, n_macro, n_micro, &
                                             seed, verbose
        procedure(logger_type), intent(in), pointer, optional :: logger_funptr

        real(rp), dimension(n_param) :: kappa, x, grad, h_diag, hess_x, residual
        real(rp) :: func
        procedure(hess_x_type), pointer :: hess_x_funptr
        logical :: converged

        ! initialize logical
        test_passed = .true.

        ! check if passed orbital updating subroutine produces correct quantities
        kappa = 1.d0
        call update_orbs_funptr(kappa, func, grad, h_diag, hess_x_funptr, error)
        if (error) then
            test_passed = .false.
            write (stderr, *) "test_solver_c_wrapper failed: Passed orbital "// &
                "updating function produced error."
        end if
        if (abs(func - 3.d0) > tol) then
            test_passed = .false.
            write (stderr, *) "test_solver_c_wrapper failed: Returned objective "// &
                "function value of passed orbital updating function wrong."
        end if
        if (any(abs(grad - 2.d0) > tol)) then
            test_passed = .false.
            write (stderr, *) "test_solver_c_wrapper failed: Returned gradient of "// &
                "passed orbital updating function wrong."
        end if
        if (any(abs(h_diag - 3.d0) > tol)) then
            test_passed = .false.
            write (stderr, *) "test_solver_c_wrapper failed: Returned Hessian "// &
                "diagonal of passed orbital updating function wrong."
        end if
        x = 1.d0
        hess_x = hess_x_funptr(x, error)
        if (error) then
            test_passed = .false.
            write (stderr, *) "test_solver_c_wrapper failed: Hessian linear "// &
                "transformation function returned by passed orbital updating "// &
                "function produced error."
        end if
        if (any(abs(hess_x - 4.d0) > tol)) then
            test_passed = .false.
            write (stderr, *) "test_solver_c_wrapper failed: Returned Hessian "// &
                "linear transformation of Hessian linear transformation function "// &
                "returned by passed orbital updating function wrong."
        end if

        ! check if passed objective function produces correct quantities
        func = obj_func_funptr(kappa, error)
        if (error) then
            test_passed = .false.
            write (stderr, *) "test_solver_c_wrapper failed: Passed objective "// &
                "function produced error."
        end if
        if (abs(func - 3.d0) > tol) then
            test_passed = .false.
            write (stderr, *) "test_solver_c_wrapper failed: Returned objective "// &
                "function value of passed objective function wrong."
        end if

        ! check number of parameters
        if (n_param /= 3) then
            test_passed = .false.
            write (stderr, *) "test_solver_c_wrapper failed: Passed number of "// &
                "parameters wrong."
        end if

        ! set output quantities
        error = .false.

        ! check if optional preconditioner function is correctly passed
        if (solver_default) then
            if (present(precond_funptr)) then
                if (associated(precond_funptr)) then
                    test_passed = .false.
                    write (stderr, *) "test_solver_c_wrapper failed: Passed "// &
                        "preconditioner function associated with value."
                end if
            end if
        else
            if (.not. present(precond_funptr)) then
                test_passed = .false.
                write (stderr, *) "test_solver_c_wrapper failed: Passed "// &
                    "preconditioner function not associated with value."
            else if (.not. associated(precond_funptr)) then
                test_passed = .false.
                write (stderr, *) "test_solver_c_wrapper failed: Passed "// &
                    "preconditioner function not associated with value."
            end if
            residual = 1.d0
            residual = precond_funptr(residual, 5.d0, error)
            if (error) then
                test_passed = .false.
                write (stderr, *) "test_solver_c_wrapper failed: Passed "// &
                    "preconditioner function produced error."
            end if
            if (any(abs(residual - 5.d0) > tol)) then
                test_passed = .false.
                write (stderr, *) "test_solver_c_wrapper failed: Returned "// &
                    "preconditioner of passed preconditioner function wrong."
            end if
        end if

        ! check if optional convergence check function is correctly passed
        if (solver_default) then
            if (present(conv_check_funptr)) then
                if (associated(conv_check_funptr)) then
                    test_passed = .false.
                    write (stderr, *) "test_solver_c_wrapper failed: Passed "// &
                        "convergence check function associated with value."
                end if
            end if
        else
            if (.not. present(conv_check_funptr)) then
                test_passed = .false.
                write (stderr, *) "test_solver_c_wrapper failed: Passed "// &
                    "convergence check function not associated with value."
            else if (.not. associated(conv_check_funptr)) then
                test_passed = .false.
                write (stderr, *) "test_solver_c_wrapper failed: Passed "// &
                    "convergence check function not associated with value."
            end if
            converged = conv_check_funptr(error)
            if (error) then
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

        ! check if optional arguments are associated with values
        if (.not. (present(stability) .and. present(line_search) .and. &
                   present(davidson) .and. present(jacobi_davidson) .and. &
                   present(prefer_jacobi_davidson) .and. present(conv_tol) .and. &
                   present(n_random_trial_vectors) .and. present(start_trust_radius) &
                   .and. present(n_macro) .and. present(n_micro) .and. &
                   present(global_red_factor) .and. present(local_red_factor) .and. &
                   present(seed) .and. present(verbose))) then
            test_passed = .false.
            write (stderr, *) "test_solver_c_wrapper failed: Passed optional "// &
                "arguments not associated with values."
        end if

        ! check if default arguments are set correctly
        if (solver_default) then
            if ((stability .neqv. solver_stability_default) .or. &
                (line_search .neqv. solver_line_search_default) .or. &
                (davidson .neqv. solver_davidson_default) .or. &
                (jacobi_davidson .neqv. solver_jacobi_davidson_default) .or. &
                (prefer_jacobi_davidson .neqv. solver_prefer_jacobi_davidson_default) &
                .or. abs(conv_tol - solver_conv_tol_default) > tol .or. &
                n_random_trial_vectors /= solver_n_random_trial_vectors_default .or. &
                abs(start_trust_radius - solver_start_trust_radius_default) > tol .or. &
                n_macro /= solver_n_macro_default .or. &
                n_micro /= solver_n_micro_default .or. &
                abs(global_red_factor - solver_global_red_factor_default) > tol .or. &
                abs(local_red_factor - solver_local_red_factor_default) > tol .or. &
                seed /= solver_seed_default .or. verbose /= solver_verbose_default) then
                test_passed = .false.
                write (stderr, *) "test_solver_c_wrapper failed: Passed optional "// &
                    "arguments associated with wrong values."
            end if
            ! check if optional arguments are correctly passed
        else
            if (stability .or. .not. line_search .or. davidson .or. jacobi_davidson &
                .or. .not. prefer_jacobi_davidson .or. abs(conv_tol - 1.d-3) > tol &
                .or. n_random_trial_vectors /= 5 .or. &
                abs(start_trust_radius - 0.2d0) > tol .or. n_macro /= 300 .or. &
                n_micro /= 200 .or. abs(global_red_factor - 1.d-2) > tol .or. &
                abs(local_red_factor - 1.d-3) > tol .or. seed /= 33 .or. verbose /= 3) &
                then
                test_passed = .false.
                write (stderr, *) "test_solver_c_wrapper failed: Passed optional "// &
                    "arguments associated with wrong values."
            end if
        end if

        ! check if optional logging function is correctly passed
        if (solver_default) then
            if (present(logger_funptr)) then
                if (associated(logger_funptr)) then
                    test_passed = .false.
                    write (stderr, *) "test_solver_c_wrapper failed: Passed "// &
                        "logging function associated with value."
                end if
            end if
        else
            if (.not. present(logger_funptr)) then
                test_passed = .false.
                write (stderr, *) "test_solver_c_wrapper failed: Passed "// &
                    "logging function not associated with value."
            else if (.not. associated(logger_funptr)) then
                test_passed = .false.
                write (stderr, *) "test_solver_c_wrapper failed: Passed "// &
                    "logging function not associated with value."
            end if
            call logger_funptr("test")
        end if

        ! move on to optional argument test for next test
        solver_default = .false.

    end subroutine mock_solver

    subroutine mock_stability_check(h_diag, hess_x_funptr, stable, kappa, error, &
                                    precond_funptr, jacobi_davidson, conv_tol, &
                                    n_random_trial_vectors, n_iter, verbose, &
                                    logger_funptr)
        !
        ! this subroutine performs a stability check
        !
        use opentrustregion, only: stability_jacobi_davidson_default, &
                                   stability_conv_tol_default, &
                                   stability_n_random_trial_vectors_default, &
                                   stability_n_iter_default, &
                                   stability_verbose_default

        real(rp), intent(in) :: h_diag(:)
        procedure(hess_x_type), intent(in), pointer :: hess_x_funptr
        logical, intent(out) :: stable, error
        real(rp), intent(out) :: kappa(:)
        procedure(precond_type), intent(in), pointer, optional :: precond_funptr
        logical, intent(in), optional :: jacobi_davidson
        real(rp), intent(in), optional :: conv_tol
        integer(ip), intent(in), optional :: n_random_trial_vectors, n_iter, verbose
        procedure(logger_type), intent(in), pointer, optional :: logger_funptr

        real(rp), dimension(size(h_diag)) :: x, hess_x, residual

        ! initialize logical
        test_passed = .true.

        ! check Hessian diagonal
        if (any(abs(h_diag - 3.d0) > tol)) then
            test_passed = .false.
            write (stderr, *) "test_stability_check_c_wrapper failed: Passed "// &
                "Hessian diagonal wrong."
        end if

        ! check if passed Hessian linear transformation function produces correct
        ! quantity
        x = 1.d0
        hess_x = hess_x_funptr(x, error)
        if (error) then
            test_passed = .false.
            write (stderr, *) "test_stability_check_c_wrapper failed: Passed "// &
                "Hessian linear transformation function produced error."
        end if
        if (any(abs(hess_x - 4.d0) > tol)) then
            test_passed = .false.
            write (stderr, *) "test_stability_check_c_wrapper failed: Hessian "// &
                "linear transformation returned by passed Hessian linear "// &
                "transformation function wrong."
        end if

        ! set output quantities
        stable = .false.
        kappa = 1.d0
        error = .false.

        ! check if optional preconditioner function is correctly passed
        if (stability_check_default) then
            if (present(precond_funptr)) then
                if (associated(precond_funptr)) then
                    test_passed = .false.
                    write (stderr, *) "test_stability_check_c_wrapper failed: "// &
                        "Passed preconditioner function associated with value."
                end if
            end if
        else
            if (.not. present(precond_funptr)) then
                test_passed = .false.
                write (stderr, *) "test_stability_check_c_wrapper failed: Passed "// &
                    "preconditioner function not associated with value."
            else if (.not. associated(precond_funptr)) then
                test_passed = .false.
                write (stderr, *) "test_stability_check_c_wrapper failed: Passed "// &
                    "preconditioner function not associated with value."
            end if
            residual = 1.d0
            residual = precond_funptr(residual, 5.d0, error)
            if (error) then
                test_passed = .false.
                write (stderr, *) "test_stability_check_c_wrapper failed: Passed "// &
                    "preconditioner function produced error."
            end if
            if (any(abs(residual - 5.d0) > tol)) then
                test_passed = .false.
                write (stderr, *) "test_stability_check_c_wrapper failed: Returned "// &
                    "preconditioner of passed preconditioner function wrong."
            end if
        end if

        ! check if optional arguments are associated with values
        if (.not. (present(jacobi_davidson) .and. present(conv_tol) .and. &
                   present(n_random_trial_vectors) .and. present(n_iter) .and. &
                   present(verbose))) then
            test_passed = .false.
            write (stderr, *) "test_stability_check_c_wrapper failed: Passed "// &
                "optional arguments not associated with values."
        end if

        ! check if default arguments are set correctly
        if (stability_check_default) then
            if ((jacobi_davidson .neqv. stability_jacobi_davidson_default) .or. &
                abs(conv_tol - stability_conv_tol_default) > tol .or. &
                n_random_trial_vectors /= stability_n_random_trial_vectors_default &
                .or. n_iter /= stability_n_iter_default .or. &
                verbose /= stability_verbose_default) then
                test_passed = .false.
                write (stderr, *) "test_stability_check_c_wrapper failed: Passed "// &
                    "optional arguments associated with wrong values."
            end if
            ! check if optional arguments are correctly passed
        else
            if (jacobi_davidson .or. abs(conv_tol - 1.d-3) > tol .or. &
                n_random_trial_vectors /= 3 .or. n_iter /= 50 .or. verbose /= 3) then
                test_passed = .false.
                write (stderr, *) "test_stability_check_c_wrapper failed: Passed "// &
                    "optional arguments associated with wrong values."
            end if
        end if

        ! check if optional logging function is correctly passed
        if (stability_check_default) then
            if (present(logger_funptr)) then
                if (associated(logger_funptr)) then
                    test_passed = .false.
                    write (stderr, *) "test_stability_check_c_wrapper failed: "// &
                        "Passed logging function associated with value."
                end if
            end if
        else
            if (.not. present(logger_funptr)) then
                test_passed = .false.
                write (stderr, *) "test_stability_check_c_wrapper failed: Passed "// &
                    "logging function not associated with value."
            else if (.not. associated(logger_funptr)) then
                test_passed = .false.
                write (stderr, *) "test_stability_check_c_wrapper failed: Passed "// &
                    "logging function not associated with value."
            end if
            call logger_funptr("test")
        end if

        ! move on to optional argument test for next test
        stability_check_default = .false.

    end subroutine mock_stability_check

end module opentrustregion_mock
