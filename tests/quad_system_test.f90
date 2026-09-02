! Copyright (C) 2025- Susi Lehtola
!
! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.
!
! Standalone end-to-end smoke test for the OpenTrustRegion solver. Runs
! Hartmann 6D in the library's working precision (whatever the library was
! built with) and asserts convergence to the known minimum, with tolerances
! scaled by epsilon(rp).
!
! Designed to validate the internal BLAS/LAPACK shim path: when the library
! is built with OPENTRUSTREGION_USE_INTERNAL_BLAS=ON, this test exercises
! every shim routine through the full solver pipeline. In quad-precision
! builds, this is the only test the testsuite has -- the Python testsuite
! cannot reach binary128 arrays through ctypes.
!
! Wired into CMake as `quad_system_test` (ctest target) when testing is
! enabled. Despite the name, it runs in whatever precision the library uses;
! "quad" reflects that the test was added together with quad-precision
! support.

module quad_system_test_mod
    use opentrustregion, only: rp, ip, hess_x_type
    implicit none

    real(rp), parameter :: alpha(4) = [1.0_rp, 1.2_rp, 3.0_rp, 3.2_rp]
    real(rp), parameter :: A(4, 6) = reshape([ &
        10.0_rp, 0.05_rp, 3.0_rp, 17.0_rp, &
        3.0_rp, 10.0_rp, 3.5_rp, 8.0_rp, &
        17.0_rp, 17.0_rp, 1.7_rp, 0.05_rp, &
        3.5_rp, 0.1_rp, 10.0_rp, 10.0_rp, &
        1.7_rp, 8.0_rp, 17.0_rp, 0.1_rp, &
        8.0_rp, 14.0_rp, 8.0_rp, 14.0_rp], [4, 6])
    real(rp), parameter :: P(4, 6) = reshape([ &
        0.1312_rp, 0.2329_rp, 0.2348_rp, 0.4047_rp, &
        0.1696_rp, 0.4135_rp, 0.1451_rp, 0.8828_rp, &
        0.5569_rp, 0.8307_rp, 0.3522_rp, 0.8732_rp, &
        0.0124_rp, 0.3736_rp, 0.2883_rp, 0.5743_rp, &
        0.8283_rp, 0.1004_rp, 0.3047_rp, 0.1091_rp, &
        0.5886_rp, 0.9991_rp, 0.6650_rp, 0.0381_rp], [4, 6])
    real(rp), parameter :: minimum1(6) = &
        [0.20168951_rp, 0.15001069_rp, 0.47687398_rp, &
         0.27533243_rp, 0.31165162_rp, 0.65730054_rp]

    real(rp) :: curr_vars(6), hess(6, 6)

contains

    function hartmann6d_func(vars) result(f)
        real(rp), intent(in) :: vars(:)
        real(rp) :: f, exp_term(4)
        integer(ip) :: i, j
        do i = 1, 4
            exp_term(i) = 0.0_rp
            do j = 1, 6
                exp_term(i) = exp_term(i) - A(i, j) * (vars(j) - P(i, j)) ** 2
            end do
            exp_term(i) = exp(exp_term(i))
        end do
        f = -sum(alpha * exp_term)
    end function hartmann6d_func

    subroutine hartmann6d_gradient(vars, grad)
        real(rp), intent(in) :: vars(:)
        real(rp), intent(out) :: grad(:)
        real(rp) :: exp_term(4)
        integer(ip) :: i, j
        do i = 1, 4
            exp_term(i) = 0.0_rp
            do j = 1, 6
                exp_term(i) = exp_term(i) - A(i, j) * (vars(j) - P(i, j)) ** 2
            end do
            exp_term(i) = exp(exp_term(i))
        end do
        do j = 1, 6
            grad(j) = sum(2.0_rp * alpha * A(:, j) * (vars(j) - P(:, j)) * exp_term)
        end do
    end subroutine hartmann6d_gradient

    subroutine hartmann6d_hessian(vars)
        real(rp), intent(in) :: vars(:)
        real(rp) :: exp_term(4)
        integer(ip) :: i, j
        do i = 1, 4
            exp_term(i) = 0.0_rp
            do j = 1, 6
                exp_term(i) = exp_term(i) - A(i, j) * (vars(j) - P(i, j)) ** 2
            end do
            exp_term(i) = exp(exp_term(i))
        end do
        do i = 1, 6
            hess(i, i) = 2.0_rp * sum(alpha * A(:, i) * exp_term * &
                (1.0_rp - 2.0_rp * A(:, i) * (vars(i) - P(:, i)) ** 2))
            do j = i + 1, 6
                hess(i, j) = -4.0_rp * sum(alpha * A(:, i) * A(:, j) * &
                    (vars(i) - P(:, i)) * (vars(j) - P(:, j)) * exp_term)
                hess(j, i) = hess(i, j)
            end do
        end do
    end subroutine hartmann6d_hessian

    subroutine hess_x_fun(x, hess_x_out, error)
        real(rp), intent(in), target :: x(:)
        real(rp), intent(out), target :: hess_x_out(:)
        integer(ip), intent(out) :: error
        error = 0
        hess_x_out = matmul(hess, x)
    end subroutine hess_x_fun

    subroutine update_orbs(delta_vars, func, grad, h_diag, hess_x_funptr, error)
        real(rp), intent(in), target :: delta_vars(:)
        real(rp), intent(out) :: func
        real(rp), intent(out), target :: grad(:), h_diag(:)
        procedure(hess_x_type), intent(out), pointer :: hess_x_funptr
        integer(ip), intent(out) :: error
        integer(ip) :: i

        error = 0
        curr_vars = curr_vars + delta_vars
        func = hartmann6d_func(curr_vars)
        call hartmann6d_gradient(curr_vars, grad)
        call hartmann6d_hessian(curr_vars)
        h_diag = [(hess(i, i), i=1, size(h_diag))]
        hess_x_funptr => hess_x_fun
    end subroutine update_orbs

    function obj_func(delta_vars, error) result(f)
        real(rp), intent(in), target :: delta_vars(:)
        integer(ip), intent(out) :: error
        real(rp) :: f
        error = 0
        f = hartmann6d_func(curr_vars + delta_vars)
    end function obj_func

end module quad_system_test_mod

program quad_system_test
    use opentrustregion, only: rp, ip, solver, solver_settings_type, &
                               update_orbs_type, obj_func_type
    use quad_system_test_mod, only: minimum1, curr_vars, &
                                    hartmann6d_func, hartmann6d_gradient, &
                                    update_orbs, obj_func
    use, intrinsic :: iso_fortran_env, only: error_unit, output_unit
    implicit none

    integer(ip), parameter :: n_param = 6
    integer(ip) :: error
    real(rp) :: f_initial, f_final, f_minimum1, dist_minimum1
    real(rp) :: final_grad(n_param)
    type(solver_settings_type) :: settings
    procedure(update_orbs_type), pointer :: update_orbs_fp
    procedure(obj_func_type), pointer :: obj_func_fp

    ! start in the quadratic region near minimum1
    curr_vars = [0.20_rp, 0.15_rp, 0.48_rp, 0.28_rp, 0.31_rp, 0.66_rp]
    f_initial = hartmann6d_func(curr_vars)
    update_orbs_fp => update_orbs
    obj_func_fp => obj_func

    call settings%init(error)
    if (error /= 0) then
        write(error_unit, '(A, I0)') &
            "quad_system_test: settings%init returned error ", error
        error stop 1
    end if
    ! tighten convergence so the test exercises tail accuracy at whatever
    ! precision the library was built with
    settings%conv_tol = sqrt(epsilon(0.0_rp)) * 10.0_rp

    call solver(update_orbs_fp, obj_func_fp, n_param, error, settings)
    if (error /= 0) then
        write(error_unit, '(A, I0)') &
            "quad_system_test: solver returned error ", error
        error stop 2
    end if

    ! gradient should be ~zero (this is the main accuracy claim)
    call hartmann6d_gradient(curr_vars, final_grad)
    if (norm2(final_grad) / sqrt(real(n_param, rp)) > settings%conv_tol) then
        write(error_unit, '(A, ES13.6, A, ES13.6)') &
            "quad_system_test: rms gradient ", &
            norm2(final_grad) / sqrt(real(n_param, rp)), &
            " exceeds conv_tol ", settings%conv_tol
        error stop 3
    end if

    ! function value should have decreased
    f_final = hartmann6d_func(curr_vars)
    if (.not. (f_final < f_initial)) then
        write(error_unit, '(A, ES23.15, A, ES23.15)') &
            "quad_system_test: function did not decrease: f_initial = ", &
            f_initial, ", f_final = ", f_final
        error stop 4
    end if

    ! sanity-check that the converged point is near the documented minimum1.
    ! Note that minimum1 literals are 8-digit values, so we can only check
    ! to ~1e-7 in any precision -- this is a basin-of-attraction check, not
    ! an accuracy check.
    dist_minimum1 = norm2(curr_vars - minimum1)
    if (dist_minimum1 > 1.0e-6_rp) then
        write(error_unit, '(A, ES13.6)') &
            "quad_system_test: converged outside minimum1 basin: dist = ", &
            dist_minimum1
        error stop 5
    end if
    f_minimum1 = hartmann6d_func(minimum1)
    if (.not. (f_final <= f_minimum1)) then
        ! the converged minimum should be at least as good as the 8-digit
        ! literal evaluation
        write(error_unit, '(A, ES23.15, A, ES23.15)') &
            "quad_system_test: f_final = ", f_final, &
            " is worse than f(minimum1 literal) = ", f_minimum1
        error stop 6
    end if

    write(output_unit, '(A, I0, A, ES23.15, A, ES13.6)') &
        "quad_system_test PASSED (rp kind=", rp, &
        ", f_final=", f_final, &
        ", |grad|=", norm2(final_grad)
end program quad_system_test
