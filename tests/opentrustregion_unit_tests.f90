! Copyright (C) 2025- Jonas Greiner
!
! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.

module opentrustregion_unit_tests

    use opentrustregion, only: rp, ip, stderr
    use test_reference, only: tol
    use, intrinsic :: iso_c_binding, only: c_bool

    implicit none

    ! parameters for 6D Hartmann function
    real(rp), parameter :: alpha(4) = [1.0_rp, 1.2_rp, 3.0_rp, 3.2_rp]
    real(rp), parameter :: A(4, 6) = reshape([10.0_rp, 0.05_rp, 3.0_rp, 17.0_rp, &
                                              3.0_rp, 10.0_rp, 3.5_rp, 8.0_rp, &
                                              17.0_rp, 17.0_rp, 1.7_rp, 0.05_rp, &
                                              3.5_rp, 0.1_rp, 10.0_rp, 10.0_rp, &
                                              1.7_rp, 8.0_rp, 17.0_rp, 0.1_rp, &
                                              8.0_rp, 14.0_rp, 8.0_rp, 14.0_rp], [4, 6])
    real(rp), parameter :: P(4, 6) = reshape([0.1312_rp, 0.2329_rp, 0.2348_rp, &
                                              0.4047_rp, 0.1696_rp, 0.4135_rp, &
                                              0.1451_rp, 0.8828_rp, 0.5569_rp, &
                                              0.8307_rp, 0.3522_rp, 0.8732_rp, &
                                              0.0124_rp, 0.3736_rp, 0.2883_rp, &
                                              0.5743_rp, 0.8283_rp, 0.1004_rp, &
                                              0.3047_rp, 0.1091_rp, 0.5886_rp, &
                                              0.9991_rp, 0.6650_rp, 0.0381_rp], [4, 6])

    ! stationary points of 6D Hartmann function
    real(rp), parameter :: minimum1(6) = [0.20168951_rp, 0.15001069_rp, 0.47687398_rp, &
                                          0.27533243_rp, 0.31165162_rp, 0.65730053_rp]
    real(rp), parameter :: minimum2(6) = [0.40465313_rp, 0.88244493_rp, 0.84610160_rp, &
                                          0.57398969_rp, 0.13892673_rp, 0.03849589_rp]
    real(rp), parameter :: saddle_point(6) = [0.35278250_rp, 0.59374767_rp, &
                                              0.47631257_rp, 0.40058250_rp, &
                                              0.31111531_rp, 0.32397158_rp]

    ! define global current variable and 6D Hartmann Hessian so that these can be 
    ! accessed by procedure pointers
    real(rp) :: curr_vars(6), hess(6, 6)

    ! global log message
    character(:), allocatable :: log_message

contains

    ! 6D Hartmann function definition

    function hartmann6d_func(vars) result(f)
        !
        ! this function defines the Hartmann 6D function
        !
        real(rp), intent(in) :: vars(:)
        real(rp) :: f, exp_term(4)
        integer(ip) :: i

        do i = 1, 4
            exp_term(i) = exp(-sum(A(i, :)*(vars - P(i, :))**2))
        end do

        f = -sum(alpha*exp_term)

    end function hartmann6d_func

    subroutine hartmann6d_gradient(vars, grad)
        !
        ! this subroutine defines the Hartmann 6D function's gradient
        !
        real(rp), intent(in) :: vars(:)
        real(rp), intent(out) :: grad(:)
        real(rp) :: exp_term(4)
        integer(ip) :: i, j

        do i = 1, 4
            exp_term(i) = exp(-sum(A(i, :)*(vars - P(i, :))**2))
        end do

        do j = 1, size(vars)
            grad(j) = sum(2.0_rp*alpha*A(:, j)*(vars(j) - P(:, j))*exp_term)
        end do

    end subroutine hartmann6d_gradient

    subroutine hartmann6d_hessian(vars)
        !
        ! this subroutine defines the Hartmann 6D function's Hessian
        !
        real(rp), intent(in) :: vars(:)
        real(rp) :: exp_term(4)
        integer(ip) :: i, j

        do i = 1, 4
            exp_term(i) = exp(-sum(A(i, :)*(vars - P(i, :))**2))
        end do

        do i = 1, size(vars)
            hess(i, i) = 2.0_rp*sum(alpha*A(:, i)*exp_term* &
                                    (1.0_rp - 2.0_rp*A(:, i)*(vars(i) - P(:, i))**2))
            do j = 1, i - 1
                hess(i, j) = -4.0_rp*sum(alpha*A(:, i)*A(:, j)*(vars(i) - P(:, i))* &
                                         (vars(j) - P(:, j))*exp_term)
                hess(j, i) = hess(i, j)
            end do
        end do

    end subroutine hartmann6d_hessian

    function hartmann6d_hess_x(x)
        !
        ! this function defines the Hessian linear transformation operation for the
        ! Hartmann 6D function
        !
        real(rp), intent(in) :: x(:)
        real(rp) :: hartmann6d_hess_x(size(x))

        hartmann6d_hess_x = matmul(hess, x)

    end function hartmann6d_hess_x

    subroutine hess_x_fun(x, hess_x, error)
        !
        ! this function describes the Hessian linear transformation operation for the
        ! Hartmann 6D function
        !
        real(rp), intent(in), target :: x(:)
        real(rp), intent(out), target :: hess_x(:)
        integer(ip), intent(out) :: error

        ! initialize error flag
        error = 0

        hess_x = hartmann6d_hess_x(x)

    end subroutine hess_x_fun

    function obj_func(delta_vars, error) result(func)
        !
        ! this function describes the objective function evaluation for the Hartmann
        ! 6D function
        !
        real(rp), intent(in), target :: delta_vars(:)
        integer(ip), intent(out) :: error
        real(rp) :: func

        ! initialize error flag
        error = 0

        func = hartmann6d_func(curr_vars + delta_vars)

    end function obj_func

    subroutine update_orbs(delta_vars, func, grad, h_diag, hess_x_funptr, error)
        !
        ! this function describes the orbital update equivalent for the Hartmann 6D
        ! function
        !
        use opentrustregion, only: hess_x_type

        real(rp), intent(in), target :: delta_vars(:)
        real(rp), intent(out) :: func
        real(rp), intent(out), target :: grad(:), h_diag(:)
        procedure(hess_x_type), intent(out), pointer :: hess_x_funptr
        integer(ip), intent(out) :: error

        integer(ip) :: i

        ! initialize error flag
        error = 0

        ! update variables
        curr_vars = curr_vars + delta_vars

        ! evaluate function, calculate gradient and Hessian diagonal and define
        ! Hessian linear transformation
        func = hartmann6d_func(curr_vars)
        call hartmann6d_gradient(curr_vars, grad)
        call hartmann6d_hessian(curr_vars)
        h_diag = [(hess(i, i), i=1, size(h_diag))]
        hess_x_funptr => hess_x_fun

    end subroutine update_orbs

    subroutine logger(message)
        !
        ! this subroutine is a mock logging subroutine
        !
        character(*), intent(in) :: message

        log_message = message

    end subroutine logger

    subroutine setup_settings(settings)
        !
        ! this subroutine sets up a settings object for tests
        !
        use opentrustregion, only: settings_type

        class(settings_type), intent(inout) :: settings
        integer(ip) :: error

        call settings%init(error)
        settings%verbose = 3
        settings%logger => logger

    end subroutine setup_settings

    logical(c_bool) function test_solver() bind(C)
        !
        ! this function tests the solver subroutine
        !
        use opentrustregion, only: update_orbs_type, obj_func_type, &
                                   solver_settings_type, solver, &
                                   default_settings => default_solver_settings

        integer(ip), parameter :: n_param = 6
        real(rp), parameter :: var_thres = 1e-6_rp
        integer(ip) :: error
        real(rp), allocatable :: final_grad(:)
        procedure(update_orbs_type), pointer :: update_orbs_funptr
        procedure(obj_func_type), pointer :: obj_func_funptr
        type(solver_settings_type) :: settings

        ! assume tests pass
        test_solver = .true.

        ! start in quadratic region near minimum
        curr_vars = [0.20_rp, 0.15_rp, 0.48_rp, 0.28_rp, 0.31_rp, 0.66_rp]
        update_orbs_funptr => update_orbs
        obj_func_funptr => obj_func

        ! initialize settings
        call settings%init(error)

        ! allocate space for the final gradient
        allocate(final_grad(n_param))

        ! run solver, check if error has occured and check whether gradient is zero and 
        ! agrees with correct minimum
        call solver(update_orbs_funptr, obj_func_funptr, n_param, error, settings)
        if (error /= 0) then
            write (stderr, *) "test_solver failed: Produced error."
            test_solver = .false.
        end if
        call hartmann6d_gradient(curr_vars, final_grad)
        if (norm2(final_grad)/sqrt(real(n_param, kind=rp)) > &
            default_settings%conv_tol) then
            write (stderr, *) "test_solver failed: Solver did not find stationary "// &
                "point."
            test_solver = .false.
        end if
        if (any(abs(curr_vars - minimum1) > var_thres)) then
            write (stderr, *) "test_solver failed: Solver did not find correct minimum."
            test_solver = .false.
        end if

        ! start near saddle point
        curr_vars = [0.35_rp, 0.59_rp, 0.48_rp, 0.40_rp, 0.31_rp, 0.32_rp]
        update_orbs_funptr => update_orbs
        obj_func_funptr => obj_func

        ! run solver, check if error has occured and check whether gradient is zero and 
        ! agrees with correct minimum
        call solver(update_orbs_funptr, obj_func_funptr, n_param, error, settings)
        if (error /= 0) then
            write (stderr, *) "test_solver failed: Produced error."
            test_solver = .false.
        end if
        call hartmann6d_gradient(curr_vars, final_grad)
        if (norm2(final_grad)/sqrt(real(n_param, kind=rp)) > &
            default_settings%conv_tol) then
            write (stderr, *) "test_solver failed: Solver did not find stationary "// &
                "point."
            test_solver = .false.
        end if
        if (any(abs(curr_vars - minimum1) > var_thres) .and. &
            any(abs(curr_vars - minimum2) > var_thres)) then
            write (stderr, *) "test_solver failed: Solver did not find correct minimum."
            test_solver = .false.
        end if

        ! start at saddle point
        curr_vars = saddle_point
        update_orbs_funptr => update_orbs
        obj_func_funptr => obj_func

        ! run solver, check if error has occured and check whether gradient is zero and 
        ! agrees with correct minimum
        call solver(update_orbs_funptr, obj_func_funptr, n_param, error, settings)
        if (error /= 0) then
            write (stderr, *) "test_solver failed: Produced error."
            test_solver = .false.
        end if
        call hartmann6d_gradient(curr_vars, final_grad)
        if (norm2(final_grad)/sqrt(real(n_param, kind=rp)) > &
            default_settings%conv_tol) then
            write (stderr, *) "test_solver failed: Solver did not find stationary "// &
                "point."
            test_solver = .false.
        end if
        if (any(abs(curr_vars - minimum1) > var_thres) .and. &
            any(abs(curr_vars - minimum2) > var_thres)) then
            write (stderr, *) "test_solver failed: Solver did not find minimum."
            test_solver = .false.
        end if

        ! deallocate space for the gradient
        deallocate(final_grad)

    end function test_solver

    logical(c_bool) function test_stability_check() bind(C)
        !
        ! this function tests the stability check subroutine
        !
        use opentrustregion, only: hess_x_type, stability_settings_type, stability_check

        real(rp) :: vars(6), h_diag(6), direction(6)
        procedure(hess_x_type), pointer :: hess_x_funptr
        logical :: stable
        integer(ip) :: error, i
        type(stability_settings_type) :: settings

        ! assume tests pass
        test_stability_check = .true.

        ! start at minimum and determine Hessian diagonal and define Hessian linear
        ! transformation
        vars = minimum1
        call hartmann6d_hessian(vars)
        h_diag = [(hess(i, i), i=1, size(h_diag))]
        hess_x_funptr => hess_x_fun

        ! initialize settings
        call settings%init(error)

        ! run stability, check if error has occured check and determine whether minimum 
        ! is stable and the returned direction vanishes
        call stability_check(h_diag, hess_x_funptr, stable, error, settings, direction)
        if (error /= 0) then
            write (stderr, *) "test_stability_check failed: Produced error."
            test_stability_check = .false.
        end if
        if (.not. stable) then
            write (stderr, *) "test_stability_check failed: Stability check "// &
                "incorrectly classifies stability of minimum."
            test_stability_check = .false.
        end if
        if (all(abs(direction) > tol)) then
            write (stderr, *) "test_stability_check failed: Stability check does "// &
                "not return zero vector for minimum"
            test_stability_check = .false.
        end if

        ! start at saddle point and determine Hessian diagonal and define linear
        ! linear transformation
        vars = saddle_point
        call hartmann6d_hessian(vars)
        h_diag = [(hess(i, i), i=1, size(h_diag))]
        hess_x_funptr => hess_x_fun

        ! run stability check, check if error has occured and determine whether saddle 
        ! point is unstable and the returned direction is correct
        call stability_check(h_diag, hess_x_funptr, stable, error, settings, direction)
        if (error /= 0) then
            write (stderr, *) "test_stability_check failed: Produced error."
            test_stability_check = .false.
        end if
        if (stable) then
            write (stderr, *) "test_stability_check failed: Stability check "// &
                "incorrectly classifies stability of saddle point."
            test_stability_check = .false.
        end if
        if (abs(abs(dot_product(direction, &
                                [-0.173375920238_rp, -0.518489821791_rp, &
                                 -6.432848975252e-3_rp, -0.340127852882_rp, &
                                 3.066460316955e-3_rp, 0.765095650196_rp])) - 1.0_rp) &
            > tol) then
            write (stderr, *) "test_stability_check failed: Stability check does "// &
                "not return correct direction for saddle point."
            test_stability_check = .false.
        end if

    end function test_stability_check

    logical(c_bool) function test_newton_step() bind(C)
        !
        ! this function tests the Newton step subroutine
        !
        use opentrustregion, only: solver_settings_type, newton_step

        type(solver_settings_type) :: settings
        real(rp) :: red_space_basis(6, 3), vars(6), grad(6), grad_norm, &
                    aug_hess(4, 4), solution(6), red_space_solution(3)
        integer(ip) :: i, j, error

        ! assume tests pass
        test_newton_step = .true.

        ! setup settings object
        call setup_settings(settings)

        ! defined a reduced space basis
        red_space_basis = &
            reshape([1.0_rp/sqrt(2.0_rp), -1.0_rp/sqrt(2.0_rp), 0.0_rp, 0.0_rp, &
                     0.0_rp, 0.0_rp, 1.0_rp/sqrt(6.0_rp), -1.0_rp/sqrt(6.0_rp), &
                     -2.0_rp/sqrt(6.0_rp), 0.0_rp, 0.0_rp, 0.0_rp, &
                     1.0_rp/sqrt(12.0_rp), -1.0_rp/sqrt(12.0_rp), &
                     1.0_rp/sqrt(12.0_rp), -3.0_rp/sqrt(12.0_rp), 0.0_rp, 0.0_rp], &
                    [6, 3])

        ! point in quadratic region near minimum
        vars = [0.20_rp, 0.15_rp, 0.48_rp, 0.28_rp, 0.31_rp, 0.66_rp]

        ! calculate gradient and Hessian to define augmented Hessian
        call hartmann6d_gradient(vars, grad)
        grad_norm = norm2(grad)
        call hartmann6d_hessian(vars)
        aug_hess = 0.0_rp
        do i = 1, 3
            do j = 1, 3
                aug_hess(i + 1, j + 1) = &
                    dot_product(red_space_basis(:, i), &
                                matmul(hess, red_space_basis(:, j)))
            end do
        end do

        ! perform Newton step, check if error has occured and determine whether 
        ! resulting solution is correct in reduced and full space
        call newton_step(aug_hess, grad_norm, red_space_basis, solution, &
                         red_space_solution, settings, error)
        if (error /= 0) then
            write (stderr, *) "test_newton_step failed: Produced error."
            test_newton_step = .false.
        end if
        if (any(abs(red_space_solution - [-2.555959788079e-2_rp, &
                                          1.565498761914e-2_rp, &
                                          4.727080080611e-3_rp]) > tol)) then
            write (stderr, *) "test_newton_step failed: Reduced space solution not "// &
                "correct."
            test_newton_step = .false.
        end if
        if (any(abs(solution - matmul(red_space_basis, red_space_solution)) > tol)) then
            write (stderr, *) "test_newton_step failed: Full space solution not "// &
                "correct."
            test_newton_step = .false.
        end if

    end function test_newton_step

    logical(c_bool) function test_bisection() bind(C)
        !
        ! this function tests the bisection subroutine
        !
        use opentrustregion, only: solver_settings_type, bisection

        type(solver_settings_type) :: settings
        real(rp) :: red_space_basis(6, 3), vars(6), grad(6), grad_norm, &
                    aug_hess(4, 4), solution(6), red_space_solution(3), trust_radius, mu
        integer(ip) :: i, j, error
        logical :: bracketed

        ! assume tests pass
        test_bisection = .true.

        ! setup settings object
        call setup_settings(settings)

        ! defined a reduced space basis
        red_space_basis = &
            reshape([1.0_rp/sqrt(2.0_rp), -1.0_rp/sqrt(2.0_rp), 0.0_rp, 0.0_rp, &
                     0.0_rp, 0.0_rp, 1.0_rp/sqrt(6.0_rp), -1.0_rp/sqrt(6.0_rp), &
                     -2.0_rp/sqrt(6.0_rp), 0.0_rp, 0.0_rp, 0.0_rp, &
                     1.0_rp/sqrt(12.0_rp), -1.0_rp/sqrt(12.0_rp), &
                     1.0_rp/sqrt(12.0_rp), -3.0_rp/sqrt(12.0_rp), 0.0_rp, 0.0_rp], &
                    [6, 3])

        ! choose target trust radius
        trust_radius = 0.4_rp

        ! point with strong negative curvature
        vars = [0.29_rp, 0.47_rp, 0.66_rp, 0.41_rp, 0.23_rp, 0.26_rp]

        ! calculate gradient and Hessian to define augmented Hessian
        call hartmann6d_gradient(vars, grad)
        grad_norm = norm2(grad)
        call hartmann6d_hessian(vars)
        aug_hess = 0.0_rp
        do i = 1, 3
            do j = 1, 3
                aug_hess(i + 1, j + 1) = dot_product(red_space_basis(:, i), &
                                                    matmul(hess, red_space_basis(:, j)))
            end do
        end do

        ! perform bisection, check whether error has occured, whether the correct trust 
        ! region was bracketed, determine whether resulting solution is correct in 
        ! reduced and full space and respects target trust radius
        call bisection(aug_hess, grad_norm, red_space_basis, trust_radius, solution, &
                       red_space_solution, mu, bracketed, settings, error)
        if (error /= 0) then
            write (stderr, *) "test_bisection failed: Produced error."
            test_bisection = .false.
        end if
        if (.not. bracketed) then
            write (stderr, *) "test_bisection failed: Unable to bracket trust region."
            test_bisection = .false.
        end if
        if (abs(norm2(solution) - trust_radius) > tol) then
            write (stderr, *) "test_bisection failed: Solution does not respect "// &
                "trust radius."
            test_bisection = .false.
        end if
        if (any(abs(red_space_solution - [-0.483593823965_rp, 0.482091645228_rp, &
                                          0.153783319727_rp]) > tol)) &
            then
            write (stderr, *) "test_bisection failed: Reduced space solution not "// &
                "correct."
            test_bisection = .false.
        end if
        if (any(abs(solution - matmul(red_space_basis, red_space_solution)) > tol)) then
            write (stderr, *) "test_bisection failed: Full space solution not correct."
            test_bisection = .false.
        end if

        ! point in quadratic region near minimum
        vars = [0.20_rp, 0.15_rp, 0.48_rp, 0.28_rp, 0.31_rp, 0.66_rp]

        ! calculate gradient and Hessian to define augmented Hessian
        call hartmann6d_gradient(vars, grad)
        grad_norm = norm2(grad)
        call hartmann6d_hessian(vars)
        aug_hess = 0.0_rp
        do i = 1, 3
            do j = 1, 3
                aug_hess(i + 1, j + 1) = &
                    dot_product(red_space_basis(:, i), &
                                matmul(hess, red_space_basis(:, j)))
            end do
        end do

        ! perform bisection and determine whether routine correctly throws error since
        ! minimum is closer than target trust radius and no level shift is necessary
        call bisection(aug_hess, grad_norm, red_space_basis, trust_radius, solution, &
                       red_space_solution, mu, bracketed, settings, error)
        if (error /= 0) then
            write (stderr, *) "test_bisection failed: Produced error."
            test_bisection = .false.
        end if
        if (bracketed) then
            write (stderr, *) "test_bisection failed: Bisection should fail if "// &
                "minimum is closer than trust radius."
            test_bisection = .false.
        end if

    end function test_bisection

    logical(c_bool) function test_bracket() bind(C)
        !
        ! this function tests the bisection subroutine
        !
        use opentrustregion, only: solver_settings_type, obj_func_type, bracket

        type(solver_settings_type) :: settings
        procedure(obj_func_type), pointer :: obj_func_funptr
        real(rp) :: vars(6), lower, upper, n
        integer(ip) :: error

        ! assume tests pass
        test_bracket = .true.

        ! setup settings object
        call setup_settings(settings)

        ! define procedure pointer
        obj_func_funptr => obj_func

        ! define direction
        vars = [0.20_rp, 0.15_rp, 0.48_rp, 0.28_rp, 0.31_rp, 0.66_rp]

        ! define lower and upper bound
        lower = 0.0_rp
        upper = 0.5_rp

        ! perform bracket and determine if new point decreases objective function in
        ! comparison to lower and upper bound
        n = bracket(obj_func_funptr, vars, lower, upper, settings, error)
        if (error /= 0) then
            write (stderr, *) "test_bracket failed: Produced error."
            test_bracket = .false.
        end if
        if (hartmann6d_func(n*vars) >= hartmann6d_func(lower*vars) .and. &
            hartmann6d_func(n*vars) >= hartmann6d_func(upper*vars)) then
            write (stderr, *) "test_bracket failed: Line search does not produce "// &
                "lower function value than starting points."
            test_bracket = .false.
        end if

        ! try different order
        n = bracket(obj_func_funptr, vars, upper, lower, settings, error)
        if (error /= 0) then
            write (stderr, *) "test_bracket failed: Produced error."
            test_bracket = .false.
        end if
        if (hartmann6d_func(n*vars) >= hartmann6d_func(lower*vars) .and. &
            hartmann6d_func(n*vars) >= hartmann6d_func(upper*vars)) then
            write (stderr, *) "test_bracket failed: Line search does not produce "// &
                "lower function value than starting points."
            test_bracket = .false.
        end if

    end function test_bracket

    logical(c_bool) function test_extend_symm_matrix() bind(C)
        !
        ! this function tests the subroutine for extending a symmetric matrix
        !
        use opentrustregion, only: extend_symm_matrix

        real(rp), allocatable :: matrix(:, :)
        real(rp) :: expected(3, 3), vector(3)

        ! assume tests pass
        test_extend_symm_matrix = .true.

        ! allocate and initialize symmetric matrix and vector to be added
        allocate(matrix(2, 2))
        matrix = reshape([1.0_rp, 2.0_rp, &
                          2.0_rp, 3.0_rp], [2, 2])
        vector = [4.0_rp, 5.0_rp, 6.0_rp]

        ! initialize expected matrix
        expected = reshape([1.0_rp, 2.0_rp, 4.0_rp, &
                            2.0_rp, 3.0_rp, 5.0_rp, &
                            4.0_rp, 5.0_rp, 6.0_rp], [3, 3])

        ! call routine and determine if dimensions and values of resulting matrix match
        call extend_symm_matrix(matrix, vector)
        if (size(matrix, 1) /= 3 .or. size(matrix, 2) /= 3) then
            write (stderr, *) "test_extend_symm_matrix failed: Incorrect matrix "// &
                "dimensions after extending."
            test_extend_symm_matrix = .false.
        end if
        if (norm2(matrix - expected) > tol) then
            write (stderr, *) "test_extend_symm_matrix failed: Incorrect matrix "// &
                "values after extending."
            test_extend_symm_matrix = .false.
        end if

        ! deallocate matrix
        deallocate(matrix)

    end function test_extend_symm_matrix

    logical(c_bool) function test_add_column() bind(C)
        !
        ! this function tests the subroutine for adding a column to a matrix
        !
        use opentrustregion, only: add_column

        real(rp), allocatable :: matrix(:, :)
        real(rp) :: expected(3, 3), new_col(3)

        ! assume tests pass
        test_add_column = .true.

        ! allocate and initialize matrix and vector to be added
        allocate(matrix(3, 2))
        matrix = reshape([1.0_rp, 2.0_rp, 3.0_rp, &
                          4.0_rp, 5.0_rp, 6.0_rp], [3, 2])
        new_col = [7.0_rp, 8.0_rp, 9.0_rp]

        ! initialize expected matrix
        expected = reshape([1.0_rp, 2.0_rp, 3.0_rp, &
                            4.0_rp, 5.0_rp, 6.0_rp, &
                            7.0_rp, 8.0_rp, 9.0_rp], [3, 3])

        ! call routine and determine if dimensions and values of resulting matrix match
        call add_column(matrix, new_col)
        if (size(matrix, 1) /= 3 .or. size(matrix, 2) /= 3) then
            write (stderr, *) "test_add_column failed: Incorrect matrix dimensions "// &
                "after adding column."
            test_add_column = .false.
        end if
        if (norm2(matrix - expected) > tol) then
            write (stderr, *) "test_add_column failed: Incorrect matrix values "// &
                "after adding column."
            test_add_column = .false.
        end if

        ! deallocate matrix
        deallocate(matrix)

    end function test_add_column

    logical(c_bool) function test_symm_mat_min_eig() bind(C)
        !
        ! this function tests the subroutine for determining the minimum eigenvalue and
        ! corresponding eigenvector for a symmetric matrix
        !
        use opentrustregion, only: solver_settings_type, symm_mat_min_eig

        type(solver_settings_type) :: settings
        real(rp) :: matrix(3, 3)
        real(rp) :: eigval, eigvec(3)
        integer(ip) :: error

        ! assume tests pass
        test_symm_mat_min_eig = .true.

        ! setup settings object
        call setup_settings(settings)

        ! initialize symmetric matrix
        matrix = reshape([3.0_rp, 1.0_rp, 1.0_rp, &
                          1.0_rp, 4.0_rp, 2.0_rp, &
                          1.0_rp, 2.0_rp, 5.0_rp], [3, 3])

        ! call routine and determine if lowest eigenvalue and corresponding eigenvector
        ! are found
        call symm_mat_min_eig(matrix, eigval, eigvec, settings, error)
        if (error /= 0) then
            write (stderr, *) "test_symm_mat_min_eig failed: Produced error."
            test_symm_mat_min_eig = .false.
        end if
        if (abs(eigval - 2.30797852837_rp) > tol) then
            write (stderr, *) "test_symm_mat_min_eig failed: Incorrect minimum "// &
                "eigenvalue for matrix."
            test_symm_mat_min_eig = .false.
        end if
        if (norm2(matmul(matrix, eigvec) - eigval*eigvec) > tol) then
            write (stderr, *) "test_symm_mat_min_eig failed: Incorrect eigenvector "// &
                "corresponding to minimum eigenvalue for matrix."
            test_symm_mat_min_eig = .false.
        end if

    end function test_symm_mat_min_eig

    logical(c_bool) function test_min_eigval() bind(C)
        !
        ! this function tests the function for determining the minimum eigenvalue for
        ! a symmetric matrix
        !
        use opentrustregion, only: solver_settings_type, min_eigval

        type(solver_settings_type) :: settings
        real(rp) :: matrix(3, 3), matrix_min_eigval
        integer(ip) :: error

        ! assume tests pass
        test_min_eigval = .true.

        ! setup settings object
        call setup_settings(settings)

        ! initialize symmetric matrix
        matrix = reshape([3.0_rp, 1.0_rp, 1.0_rp, &
                          1.0_rp, 4.0_rp, 2.0_rp, &
                          1.0_rp, 2.0_rp, 5.0_rp], [3, 3])

        ! call function and determine if lowest eigenvalue is found
        matrix_min_eigval = min_eigval(matrix, settings, error)
        if (error /= 0) then
            write (stderr, *) "test_min_eigval failed: Produced error."
            test_min_eigval = .false.
        end if
        if (abs(matrix_min_eigval - 2.30797852837_rp) > tol) then
            write (stderr, *) "test_min_eigval failed: Incorrect minimum "// &
                "eigenvalue for matrix."
            test_min_eigval = .false.
        end if

    end function test_min_eigval

    logical(c_bool) function test_init_rng() bind(C)
        !
        ! this function tests the initialization subroutine for the random number
        ! generator
        !
        use opentrustregion, only: init_rng

        integer(ip) :: seed1, seed2, i
        real(rp) :: rand_seq1(5), rand_seq2(5), rand_seq3(5)

        ! assume tests pass
        test_init_rng = .true.

        ! define seeds
        seed1 = 12345
        seed2 = 67890

        ! call rng with first seed
        call init_rng(seed1)
        do i = 1, 5
            call random_number(rand_seq1(i))
        end do

        ! call rng with first seed
        call init_rng(seed1)
        do i = 1, 5
            call random_number(rand_seq2(i))
        end do

        ! call rng with second seed
        call init_rng(seed2)
        do i = 1, 5
            call random_number(rand_seq3(i))
        end do

        ! check reproducibility
        if (any(abs(rand_seq1 - rand_seq2) > tol)) then
            write (stderr, *) "test_init_rng failed: RNG does not produce "// &
                "consistent sequences for the same seed."
            test_init_rng = .false.
        end if

        ! check variation
        if (all(abs(rand_seq1 - rand_seq3) < tol)) then
            write (stderr, *) "test_init_rng failed: RNG produces identical "// &
                "sequences for different seeds."
            test_init_rng = .false.
        end if

    end function test_init_rng

    logical(c_bool) function test_generate_trial_vectors() bind(C)
        !
        ! this function tests the function which generates trial vectors for the
        ! Davidson procedure
        !
        use opentrustregion, only: solver_settings_type, generate_trial_vectors

        type(solver_settings_type) :: settings
        real(rp), allocatable :: red_space_basis(:, :)
        real(rp) :: grad(4), h_diag(4), grad_norm
        integer(ip) :: error, i, j

        ! assume tests pass
        test_generate_trial_vectors = .true.

        ! setup settings object
        call setup_settings(settings)
        settings%n_random_trial_vectors = 2

        ! define gradient
        grad = [1.0_rp, 2.0_rp, 3.0_rp, 4.0_rp]
        grad_norm = norm2(grad)

        ! define all positive Hessian diagonal elements
        h_diag = [1.0_rp, 2.0_rp, 3.0_rp, 4.0_rp]

        ! generate trial vectors and determine whether function returns the correct
        ! number of orthonormal trial vectors
        red_space_basis = generate_trial_vectors(grad, grad_norm, h_diag, settings, &
                                                 error)
        if (error /= 0) then
            write (stderr, *) "test_generate_trial_vectors failed: Produced error."
            test_generate_trial_vectors = .false.
        end if
        if (.not. allocated(red_space_basis)) then
            write (stderr, *) "test_generate_trial_vectors failed: Reduced space "// &
                "basis not allocated."
            test_generate_trial_vectors = .false.
            return
        end if
        if (size(red_space_basis, 2) /= 1 + settings%n_random_trial_vectors) then
            write (stderr, *) "test_generate_trial_vectors failed: Incorrect "// &
                "number of vectors for Hessian with only positive diagonal elements."
            test_generate_trial_vectors = .false.
        end if
        do i = 1, size(red_space_basis, 2)
            do j = i + 1, size(red_space_basis, 2)
                if (abs(dot_product(red_space_basis(:, i), red_space_basis(:, j))) > &
                    tol) then
                    write (stderr, *) "test_generate_trial_vectors failed: "// &
                        "Generated vectors are not orthonormal for Hessian with "// &
                        "only positive diagonal elements."
                    test_generate_trial_vectors = .false.
                end if
            end do
        end do

        ! deallocate reduced space basis
        deallocate(red_space_basis)

        ! define Hessian diagonal with negative elements
        h_diag = [-1.0_rp, 2.0_rp, 3.0_rp, 4.0_rp]

        ! generate trial vectors and determine whether function returns the correct
        ! number of orthonormal trial vectors
        red_space_basis = generate_trial_vectors(grad, grad_norm, h_diag, settings, &
                                                 error)
        if (error /= 0) then
            write (stderr, *) "test_generate_trial_vectors failed: Produced error."
            test_generate_trial_vectors = .false.
        end if
        if (.not. allocated(red_space_basis)) then
            write (stderr, *) "test_generate_trial_vectors failed: Reduced space "// &
                "basis not allocated."
            test_generate_trial_vectors = .false.
            return
        end if
        if (size(red_space_basis, 2) /= 2 + settings%n_random_trial_vectors) then
            write (stderr, *) "test_generate_trial_vectors failed: Incorrect "// &
                "number of vectors for Hessian with diagonal elements."
            test_generate_trial_vectors = .false.
        end if
        do i = 1, size(red_space_basis, 2)
            do j = i + 1, size(red_space_basis, 2)
                if (abs(dot_product(red_space_basis(:, i), red_space_basis(:, j))) > &
                    tol) then
                    write (stderr, *) "test_generate_trial_vectors failed: "// &
                        "Generated vectors are not orthonormal for Hessian with "// &
                        "diagonal elements."
                    test_generate_trial_vectors = .false.
                end if
            end do
        end do

        ! deallocate reduced space basis
        deallocate(red_space_basis)

    end function test_generate_trial_vectors

    logical(c_bool) function test_generate_random_trial_vectors() bind(C)
        !
        ! this function tests the function which generates random trial vectors for the
        ! Davidson procedure
        !
        use opentrustregion, only: solver_settings_type, generate_random_trial_vectors

        type(solver_settings_type) :: settings
        real(rp), allocatable :: red_space_basis(:, :)
        integer(ip) :: error, i, j

        ! assume tests pass
        test_generate_random_trial_vectors = .true.

        ! setup settings object
        call setup_settings(settings)
        settings%n_random_trial_vectors = 2

        ! allocate reduced space basis and set first normalized basis vector
        allocate(red_space_basis(4, 3))
        red_space_basis(:, 1) = [1.0_rp, 2.0_rp, 3.0_rp, 4.0_rp]
        red_space_basis(:, 1) = red_space_basis(:, 1) / norm2(red_space_basis(:, 1))

        ! generate trial vectors and determine whether function returns orthonormal 
        ! trial vectors
        call generate_random_trial_vectors(red_space_basis, settings, error)
        if (error /= 0) then
            write (stderr, *) "test_generate_trial_vectors failed: Produced error."
            test_generate_random_trial_vectors = .false.
        end if
        do i = 1, size(red_space_basis, 2)
            do j = i + 1, size(red_space_basis, 2)
                if (abs(dot_product(red_space_basis(:, i), red_space_basis(:, j))) > &
                    tol) then
                    write (stderr, *) "test_generate_trial_vectors failed: "// &
                        "Generated vectors are not orthonormal for Hessian with "// &
                        "only positive diagonal elements."
                    test_generate_random_trial_vectors = .false.
                end if
            end do
        end do

        ! deallocate reduced space basis
        deallocate(red_space_basis)

    end function test_generate_random_trial_vectors

    logical(c_bool) function test_gram_schmidt() bind(C)
        !
        ! this function tests the Gram-Schmidt subroutine which orthonormalizes a 
        ! vector to a given basis
        !
        use opentrustregion, only: solver_settings_type, gram_schmidt

        type(solver_settings_type) :: settings
        real(rp) :: vector(4), lin_trans_vector(4), vector_small(2), space(4, 2), &
                    symm_matrix(4, 4), lin_trans_space(4, 2), space_small(2, 2)
        integer(ip) :: error

        ! assume tests pass
        test_gram_schmidt = .true.

        ! setup settings object
        call setup_settings(settings)

        ! define vector to be orthogonalized and space
        vector = [1.0_rp, 2.0_rp, 3.0_rp, 4.0_rp]
        space(:, 1) = [0.0_rp, 1.0_rp, 0.0_rp, 0.0_rp]
        space(:, 2) = [0.0_rp, 0.0_rp, 1.0_rp, 0.0_rp]

        ! perform Gram-Schmidt orthogonalization and determine whether added vector is
        ! orthonormalized
        call gram_schmidt(vector, space, settings, error)
        if (error /= 0) then
            write (stderr, *) "test_gram_schmidt failed: Produced error."
            test_gram_schmidt = .false.
        end if
        if (abs(dot_product(vector, space(:, 1))) > tol .or. &
            abs(dot_product(vector, space(:, 2))) > tol) then
            write (stderr, *) "test_gram_schmidt failed: Added vector not orthogonal."
            test_gram_schmidt = .false.
        end if
        if (abs(norm2(vector) - 1.0_rp) > tol) then
            write (stderr, *) "test_gram_schmidt failed: Added vector not normalized."
            test_gram_schmidt = .false.
        end if

        ! define vector to be orthogonalized and space
        vector = [1.0_rp, 2.0_rp, 3.0_rp, 4.0_rp]
        space(:, 1) = [0.0_rp, 1.0_rp, 0.0_rp, 0.0_rp]
        space(:, 2) = [0.0_rp, 0.0_rp, 1.0_rp, 0.0_rp]

        ! define symmetric linear transformation and corresponding vector and space
        symm_matrix = reshape([ 1.0_rp, -5.0_rp,  8.0_rp,  0.0_rp, &
                               -5.0_rp,  2.0_rp, -6.0_rp,  9.0_rp, &
                                8.0_rp, -6.0_rp,  3.0_rp, -7.0_rp, &
                                0.0_rp,  9.0_rp, -7.0_rp,  4.0_rp], &
                              shape(symm_matrix), order=[2,1])
        lin_trans_vector = matmul(symm_matrix, vector)
        lin_trans_space = matmul(symm_matrix, space)

        ! perform Gram-Schmidt orthogonalization and determine whether added vector is
        ! orthonormalized and linear transformation is correct
        call gram_schmidt(vector, space, settings, error, lin_trans_vector, &
                          lin_trans_space)
        if (error /= 0) then
            write (stderr, *) "test_gram_schmidt failed: Produced error."
            test_gram_schmidt = .false.
        end if
        if (abs(dot_product(vector, space(:, 1))) > tol .or. &
            abs(dot_product(vector, space(:, 2))) > tol) then
            write (stderr, *) "test_gram_schmidt failed: Added vector not orthogonal."
            test_gram_schmidt = .false.
        end if
        if (abs(norm2(vector) - 1.0_rp) > tol) then
            write (stderr, *) "test_gram_schmidt failed: Added vector not normalized."
            test_gram_schmidt = .false.
        end if
        if (abs(sum(abs(lin_trans_vector - matmul(symm_matrix, vector)))) > tol) then
            write (stderr, *) "test_gram_schmidt failed: Added linear "// &
                "transformation not correct."
            test_gram_schmidt = .false.
        end if

        ! define zero vector
        vector = [0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp]

        ! perform Gram-Schmidt orthogonalization and determine if function correctly
        ! throws error
        call gram_schmidt(vector, space, settings, error)
        if ((error == 0) .or. (log_message /= " Vector passed to Gram-Schmidt "// &
                                "procedure is numerically zero.")) then
            write (stderr, *) "test_gram_schmidt failed: No error returned during "// &
                "orthogonalization for zero vector."
            test_gram_schmidt = .false.
        end if

        ! define vector in space that is already complete
        vector_small = [1.0_rp, 2.0_rp]
        space_small(:, 1) = [1.0_rp, 0.0_rp]
        space_small(:, 2) = [0.0_rp, 1.0_rp]

        ! perform Gram-Schmidt orthogonalization and determine if function correctly
        ! throws error
        call gram_schmidt(vector_small, space_small, settings, error)
        if (error == 0 .or. (log_message /= " Number of vectors in Gram-Schmidt "// &
                              "procedure larger than dimension of vector space.")) then
            write (stderr, *) "test_gram_schmidt failed: No error returned during "// &
                "orthogonalization when number of vectors is larger than dimension "// &
                "of vector space."
            test_gram_schmidt = .false.
        end if

    end function test_gram_schmidt

    logical(c_bool) function test_init_solver_settings() bind(C)
        !
        ! this function tests the subroutine which initializes the solver settings
        !
        use opentrustregion, only: solver_settings_type, &
                                   default_settings => default_solver_settings
        use test_reference, only: operator(/=)

        type(solver_settings_type) :: settings
        integer(ip) :: error

        ! assume tests pass
        test_init_solver_settings = .true.

        ! initialize settings
        call settings%init(error)

        ! check function pointers
        if (error /= 0) then
            write (stderr, *) "test_init_solver_settings failed: Function raised "// &
                "error."
            test_init_solver_settings = .false.
        end if

        ! check function pointers
        if (associated(settings%precond) .or. associated(settings%conv_check) .or. &
            associated(settings%logger)) then
            write (stderr, *) "test_init_solver_settings failed: Function pointers "// &
                "should not be initialized."
            test_init_solver_settings = .false.
        end if

        ! check settings
        if (settings /= default_settings) then
            write (stderr, *) "test_init_solver_settings failed: Settings not "// &
                "initialized correctly."
            test_init_solver_settings = .false.
        end if

    end function test_init_solver_settings

    logical(c_bool) function test_init_stability_settings() bind(C)
        !
        ! this function tests the subroutine which initializes the stability check
        ! settings
        !
        use opentrustregion, only: stability_settings_type, &
                                   default_settings => default_stability_settings
        use test_reference, only: operator(/=)

        type(stability_settings_type) :: settings
        integer(ip) :: error

        ! assume tests pass
        test_init_stability_settings = .true.

        ! initialize settings
        call settings%init(error)

        ! check function pointers
        if (error /= 0) then
            write (stderr, *) "test_init_stability_settings failed: Function "// &
                "raised error."
            test_init_stability_settings = .false.
        end if

        ! check function pointers
        if (associated(settings%precond) .or. associated(settings%logger)) then
            write (stderr, *) "test_init_stability_settings failed: Function "// &
                "pointers should not be initialized."
            test_init_stability_settings = .false.
        end if

        ! check settings
        if (settings /= default_settings) then
            write (stderr, *) "test_init_stability_settings failed: Settings not "// &
                "initialized correctly."
            test_init_stability_settings = .false.
        end if

    end function test_init_stability_settings

    logical(c_bool) function test_level_shifted_diag_precond() bind(C)
        !
        ! this function tests the subroutine that constructs the level-shifted diagonal 
        ! preconditioner
        !
        use opentrustregion, only: level_shifted_diag_precond

        real(rp) :: residual(3), h_diag(3), precond_residual(3)

        ! assume tests pass
        test_level_shifted_diag_precond = .true.

        ! initialize quantities
        residual = [1.0_rp, 1.0_rp, 1.0_rp]
        h_diag = [-2.0_rp, 1.0_rp, 2.0_rp]

        ! call subroutine and check if results match
        call level_shifted_diag_precond(residual, -2.0_rp, h_diag, precond_residual)
        if (any(abs(precond_residual - [1e10_rp, 1.0_rp / 3, 0.25_rp]) > tol)) then
            write (stderr, *) "test_level_shifted_diag_precond failed: Returned "// &
                "preconditioned residual not correct."
            test_level_shifted_diag_precond = .false.
        end if

    end function test_level_shifted_diag_precond

    logical(c_bool) function test_abs_diag_precond() bind(C)
        !
        ! this function tests the subroutine that constructs the absolute diagonal 
        ! preconditioner
        !
        use opentrustregion, only: abs_diag_precond

        real(rp) :: residual(3), h_diag(3), precond_residual(3)

        ! assume tests pass
        test_abs_diag_precond = .true.

        ! initialize quantities
        residual = [1.0_rp, 1.0_rp, 1.0_rp]
        h_diag = [-2.0_rp, 0.0_rp, 2.0_rp]

        ! call subroutine and check if results match
        call abs_diag_precond(residual, h_diag, precond_residual)
        if (any(abs(precond_residual - [0.5_rp, 1e10_rp, 0.5_rp]) > tol)) then
            write (stderr, *) "test_abs_diag_precond failed: Returned "// &
                "preconditioned residual not correct."
            test_abs_diag_precond = .false.
        end if

    end function test_abs_diag_precond

    logical(c_bool) function test_orthogonal_projection() bind(C)
        !
        ! this function tests the orthogonal projection function which removes a 
        ! certain direction from a vector
        !
        use opentrustregion, only: orthogonal_projection

        real(rp), dimension(4) :: vector, direction

        ! assume tests pass
        test_orthogonal_projection = .true.

        ! define vector and direction to be projected out, the latter needs to be 
        ! normalized
        vector = [1.0_rp, 2.0_rp, 3.0_rp, 4.0_rp]
        direction = [0.0_rp, 1.0_rp, 2.0_rp, 0.0_rp] / sqrt(5.0_rp)

        ! perform orthogonal projection and determine whether vector contains direction
        vector = orthogonal_projection(vector, direction)
        if (abs(dot_product(vector, direction)) > tol) then
            write (stderr, *) "test_orthogonal_projection failed: Vector contains "// &
                "component from direction to be projected out."
            test_orthogonal_projection = .false.
        end if

    end function test_orthogonal_projection

    logical(c_bool) function test_jacobi_davidson_correction() bind(C)
        !
        ! this function tests the Jacobi-Davidson correction subroutine
        !
        use opentrustregion, only: solver_settings_type, hess_x_type, &
                                   jacobi_davidson_correction

        type(solver_settings_type) :: settings
        procedure(hess_x_type), pointer :: hess_x_funptr
        real(rp), dimension(6) :: vars, vector, solution, corr_vector, hess_vector
        integer(ip) :: error

        ! assume tests pass
        test_jacobi_davidson_correction = .true.

        ! setup settings object
        call setup_settings(settings)

        ! define point near saddle point, define trial vector, and solution to be 
        ! projected out
        vars = [0.35_rp, 0.59_rp, 0.48_rp, 0.40_rp, 0.31_rp, 0.32_rp]
        vector = [0.1_rp, 0.2_rp, 0.3_rp, 0.4_rp, 0.5_rp, 0.6_rp]
        solution = [1.0_rp, -2.0_rp, 2.0_rp, -1.0_rp, 1.0_rp, -2.0_rp]

        ! generate Hessian
        call hartmann6d_hessian(vars)

        ! define Hessian linear transformation
        hess_x_funptr => hess_x_fun

        ! calculate Jacobi-Davidson correction and compare values
        call jacobi_davidson_correction(hess_x_funptr, vector, solution, 0.5_rp, &
                                        corr_vector, hess_vector, settings, error)
        if (error /= 0) then
            write (stderr, *) "test_jacobi_davidson_correction failed: Returned error."
            test_jacobi_davidson_correction = .false.
        end if
        if (sum(abs(corr_vector - [-96.940677944_rp, 203.929698480_rp, &
                                   -216.199768920_rp, 100.656941418_rp, &
                                   -90.624469448_rp, 212.045768918_rp])) > 1e-8_rp) then
            write (stderr, *) "test_jacobi_davidson_correction failed: Returned "// &
                "correction vector wrong."
            test_jacobi_davidson_correction = .false.
        end if
        if (sum(abs(hess_vector - [14.407362159_rp, -18.566381727_rp, &
                                   6.546311286_rp, -10.441098685_rp, &
                                   20.923570656_rp, -10.250311288_rp])) > 1e-8_rp) then
            write (stderr, *) "test_jacobi_davidson_correction failed: Returned "// &
                "Hessian linear transformation wrong."
            test_jacobi_davidson_correction = .false.
        end if

    end function test_jacobi_davidson_correction

    logical(c_bool) function test_minres() bind(C)
        !
        ! this function tests the minimum residual method subroutine
        !
        use opentrustregion, only: solver_settings_type, hess_x_type, minres

        type(solver_settings_type) :: settings
        procedure(hess_x_type), pointer :: hess_x_funptr
        real(rp), dimension(6) :: vars, rhs, solution, vector, hess_vector, corr_vector
        real(rp) :: mu
        real(rp), parameter :: rtol = 1e-14_rp
        integer(ip) :: error

        ! assume tests pass
        test_minres = .true.

        ! setup settings object
        call setup_settings(settings)

        ! define point near saddle point
        vars = [0.35_rp, 0.59_rp, 0.48_rp, 0.40_rp, 0.31_rp, 0.32_rp]

        ! define solution to be projected out
        solution = [1.0_rp, 2.0_rp, 3.0_rp, 4.0_rp, 5.0_rp, 6.0_rp] / sqrt(91.0_rp)

        ! generate Hessian
        call hartmann6d_hessian(vars)

        ! define Hessian linear transformation
        hess_x_funptr => hess_x_fun

        ! define Rayleigh quotient
        mu = dot_product(solution, hartmann6d_hess_x(solution))

        ! define right hand side based on residual, this ensures rhs is orthogonal to 
        ! solution if mu describes the Rayleigh quotient and solution is normalized
        rhs = hartmann6d_hess_x(solution) - mu * solution

        ! run minimum residual method, check if Jacobi-Davidson correction equation is 
        ! solved and whether Hessian linear transformation is correct, if rhs and 
        ! solution are orthogonal (as in Jacobi-Davidson), the final vector will be 
        ! orthogonal to the solution vector and consequently the Hessian linear 
        ! transformation of the projected vector is equivalent to the Hessian linear 
        ! transformation of the vector itself
        call minres(-rhs, hess_x_funptr, solution, mu, rtol, vector, hess_vector, &
                    settings, error)
        if (error /= 0) then
            write (stderr, *) "test_minres failed: Returned error."
            test_minres = .false.
        end if
        corr_vector = vector - dot_product(vector, solution) * solution
        corr_vector = hartmann6d_hess_x(corr_vector) - mu * corr_vector
        corr_vector = corr_vector - dot_product(corr_vector, solution) * solution
        if (sum(abs(corr_vector + rhs)) > tol) then
            write (stderr, *) "test_minres failed: Returned solution does not "// & 
                "solve Jacobi-Davidson correction equation."
            test_minres = .false.
        end if
        if (sum(abs(hess_vector + dot_product(vector, solution) * &
                    hartmann6d_hess_x(solution) - hartmann6d_hess_x(vector))) > tol) &
            then
            write (stderr, *) "test_minres failed: Returned Hessian linear "// &
                "transformation wrong."
            test_minres = .false.
        end if

        ! run minimum residual method for vanishing right hand side
        rhs = 0.0_rp
        call minres(-rhs, hess_x_funptr, solution, mu, rtol, vector, hess_vector, &
                    settings, error)
        if (error /= 0) then
            write (stderr, *) "test_minres failed: Returned error."
            test_minres = .false.
        end if
        if (sum(abs(vector)) > tol) then
            write (stderr, *) "test_minres failed: Returned solution is not zero "// &
                "for a vanishing rhs."
            test_minres = .false.
        end if
        if (sum(abs(hess_vector)) > tol) then
            write (stderr, *) "test_minres failed: Returned Hessian linear "// &
                "transformation is not zero for a vanishing rhs."
            test_minres = .false.
        end if

    end function test_minres

    logical(c_bool) function test_print_results() bind(C)
        !
        ! this function tests the subroutine that prints the result table
        !
        use opentrustregion, only: solver_settings_type

        type(solver_settings_type) :: settings

        ! assume tests pass
        test_print_results = .true.

        ! setup settings object
        call setup_settings(settings)

        ! print row of results table without optional arguments and check if row is 
        ! correct
        call settings%print_results(1, 2.0_rp, 3.0_rp)
        if (log_message /= "        1   |     2.00000000000000E+00   "// &
            "|   3.00E+00   |      -      |        -   |        -     |      -   ") then
            write (stderr, *) "test_print_results failed: Printed row without "// &
                "optional arguments not correct."
            test_print_results = .false.
        end if

        ! print row of results table with optional arguments
        call settings%print_results(1, 2.0_rp, 3.0_rp, 4.0_rp, 5, 6, 7.0_rp, 8.0_rp)
        if (log_message /= "        1   |     2.00000000000000E+00   "// &
            "|   3.00E+00   |   4.00E+00  |   0 |   5  |    7.00E+00  |  8.00E+00") then
            write (stderr, *) "test_print_results failed: Printed row with "// &
                "optional arguments not correct."
            test_print_results = .false.
        end if

    end function test_print_results

    logical(c_bool) function test_print_message() bind(C)
        !
        ! this function tests the logging subroutine
        !
        use opentrustregion, only: solver_settings_type, print_message

        type(solver_settings_type) :: settings

        ! assume tests pass
        test_print_message = .true.

        ! setup settings object
        call setup_settings(settings)

        ! check if logging is correctly performed according to verbosity level when 
        ! logger is provided
        call print_message(settings, "This is a test message.", 1)
        if (trim(log_message) /= " This is a test message.") then
            write (stderr, *) "test_print_message failed: Log message is not "// &
                "printed correctly even though it should be according to verbosity "// &
                "level."
            test_print_message = .false.
        end if
        call print_message(settings, "This is another test message.", 4)
        if (log_message == " This is another test message.") then
            write (stderr, *) "test_print_message failed: Log message is printed "// &
                "even though it should not be according to verbosity level."
            test_print_message = .false.
        end if

    end function test_print_message

    logical(c_bool) function test_split_string_by_space() bind(C)
        !
        ! this function tests the subroutine which splits strings after a space if 
        ! they exceed a given maximum length
        !
        use opentrustregion, only: split_string_by_space

        character(23), parameter :: message = "This is a test message."
        character(:), allocatable :: substrings(:)

        ! assume tests pass
        test_split_string_by_space = .true.

        ! check if strings are split correctly on spaces
        call split_string_by_space(message, 8, substrings)
        if (size(substrings) == 3) then
            if (trim(substrings(1)) /= "This is" .or. trim(substrings(2)) /= "a test" &
                .or. trim(substrings(3)) /= "message.") then
                write (stderr, *) "test_split_string_by_space failed: Split "// &
                    "strings incorrect."
                test_split_string_by_space = .false.
            end if
        else
            write (stderr, *) "test_split_string_by_space failed: Number of "// &
                "substrings incorrect."
            test_split_string_by_space = .false.
        end if

        ! check if strings are split correctly if splitting on spaces is not possible
        call split_string_by_space(message, 5, substrings)
        if (size(substrings) == 5) then
            if (trim(substrings(1)) /= "This" .or. trim(substrings(2)) /= "is a" .or. &
                trim(substrings(3)) /= "test" .or. trim(substrings(4)) /= "messa" .or. &
                trim(substrings(5)) /= "ge.") then
                write (stderr, *) "test_split_string_by_space failed: Split "// &
                    "strings incorrect."
                test_split_string_by_space = .false.
            end if
        else
            write (stderr, *) "test_split_string_by_space failed: Number of "// &
                "substrings incorrect."
            test_split_string_by_space = .false.
        end if

    end function test_split_string_by_space

    logical(c_bool) function test_accept_trust_region_step() bind(C)
        !
        ! this function tests the subroutine which determines whether to accept a 
        ! trust-region step
        !
        use opentrustregion, only: solver_settings_type, accept_trust_region_step, &
                                   trust_radius_shrink_ratio, &
                                   trust_radius_expand_ratio, &
                                   trust_radius_shrink_factor, &
                                   trust_radius_expand_factor

        logical :: accept_step, max_precision_reached
        real(rp) :: solution(3), trust_radius
        type(solver_settings_type) :: settings

        ! assume tests pass
        test_accept_trust_region_step = .true.

        ! setup settings object
        call setup_settings(settings)

        solution = [0.3_rp, 0.3_rp, 0.3_rp]

        ! check if step is rejected and trust radius is correctly reduced if micro 
        ! iterations have not converged
        trust_radius = 1.0_rp
        accept_step = accept_trust_region_step(solution, 1.0_rp, .false., settings, &
                                               trust_radius, max_precision_reached)
        if (accept_step .or. abs(trust_radius - trust_radius_shrink_factor) > tol) then
            write(stderr, *) "test_accept_trust_region_step failed: Step accepted "// &
                "or trust radius not correctly reduced when micro iterations have "// &
                "not converged."
            test_accept_trust_region_step = .false.
        end if

        ! check if step is rejected and trust radius is correctly reduced if ratio is 
        ! negative
        trust_radius = 1.0_rp
        accept_step = accept_trust_region_step(solution, -1.0_rp, .true., settings, &
                                               trust_radius, max_precision_reached)
        if (accept_step .or. abs(trust_radius - trust_radius_shrink_factor) > tol) then
            write(stderr, *) "test_accept_trust_region_step failed: Step accepted "// & 
                "or trust radius not correctly reduced when ratio is negative."
            test_accept_trust_region_step = .false.
        end if

        ! check if step is rejected and trust radius is correctly reduced if 
        ! individual rotations are too large
        trust_radius = 1.0_rp
        solution(1) = 1.0_rp
        accept_step = accept_trust_region_step(solution, 1.0_rp, .true., settings, &
                                               trust_radius, max_precision_reached)
        if (accept_step .or. abs(trust_radius - trust_radius_shrink_factor) > tol) then
            write(stderr, *) "test_accept_trust_region_step failed: Step accepted "// &
                "or trust radius not correctly reduced when individual rotations "// &
                "are too large."
            test_accept_trust_region_step = .false.
        end if
        solution(1) = 0.3_rp

        ! check if step is accepted and trust radius is correctly reduced if ratio is 
        ! too small
        trust_radius = 1.0_rp
        accept_step = accept_trust_region_step(solution, &
                                               0.9_rp * trust_radius_shrink_ratio, &
                                               .true., settings, trust_radius, &
                                               max_precision_reached)
        if (.not. accept_step .or. abs(trust_radius - trust_radius_shrink_factor) > &
            tol) then
            write(stderr, *) "test_accept_trust_region_step failed: Step not "// &
                "accepted or trust radius not correctly reduced when ratio is too "// &
                "small."
            test_accept_trust_region_step = .false.
        end if

        ! check if step is accepted and trust radius is correctly reduced if ratio is 
        ! ok
        trust_radius = 1.0_rp
        accept_step = accept_trust_region_step(solution, &
                                               0.5_rp * (trust_radius_shrink_ratio + &
                                               trust_radius_expand_ratio), .true., &
                                               settings, trust_radius, &
                                               max_precision_reached)
        if (.not. accept_step .or. abs(trust_radius - 1.0_rp) > tol) then
            write(stderr, *) "test_accept_trust_region_step failed: Step not "// &
                "accepted or trust radius changed when ratio is acceptable."
            test_accept_trust_region_step = .false.
        end if

        ! check if step is accepted and trust radius is correctly expanded if ratio is 
        ! too large
        trust_radius = 1.0_rp
        accept_step = accept_trust_region_step(solution, &
                                               1.1_rp * trust_radius_expand_ratio, &
                                               .true., settings, trust_radius, &
                                               max_precision_reached)
        if (.not. accept_step .or. abs(trust_radius - trust_radius_expand_factor) > &
            tol) then
            write(stderr, *) "test_accept_trust_region_step failed: Step not "// &
                "accepted or trust radius not correctly expanded when ratio is too "// &
                "large."
            test_accept_trust_region_step = .false.
        end if

    end function test_accept_trust_region_step

    logical(c_bool) function test_solver_sanity_check() bind(C)
        !
        ! this function tests the subroutine which performs a sanity check for the 
        ! solver
        !
        use opentrustregion, only: solver_settings_type, solver_sanity_check

        type(solver_settings_type) :: settings
        real(rp) :: grad(3)
        integer(ip) :: error

        ! assume tests pass
        test_solver_sanity_check = .true.

        ! setup settings object
        call setup_settings(settings)

        ! check if error is incorrectly thrown for finite and non-negative number of 
        ! parameters
        settings%n_random_trial_vectors = 0
        call solver_sanity_check(settings, 3, grad, error)
        if (error /= 0) then
            write(stderr, *) "test_solver_sanity_check failed: Error thrown for "// &
                "non-negative and non-vanishing number of parameters."
            test_solver_sanity_check = .false.
        end if

        ! check if error is correctly thrown for vanishing number of parameters
        call solver_sanity_check(settings, 0, grad, error)
        if (error == 0) then
            write(stderr, *) "test_solver_sanity_check failed: Error not thrown "// &
                "for vanishing number of parameters."
            test_solver_sanity_check = .false.
        end if

        ! check if error is correctly thrown for negative number of parameters
        call solver_sanity_check(settings, -1, grad, error)
        if (error == 0) then
            write(stderr, *) "test_solver_sanity_check failed: Error not thrown "// &
                "for negative number of parameters."
            test_solver_sanity_check = .false.
        end if

        ! check if number of random trial vectors is reduced correctly
        settings%n_random_trial_vectors = 3
        call solver_sanity_check(settings, 3, grad, error)
        if (settings%n_random_trial_vectors /= 1) then
            write(stderr, *) "test_solver_sanity_check failed: Number of random "// &
                "trial not correctly set."
            test_solver_sanity_check = .false.
        end if

        ! check if gradient size is treated correctly
        call solver_sanity_check(settings, 3, grad, error)
        if (error /= 0) then
            write(stderr, *) "test_solver_sanity_check failed: Error thrown for "// &
                "gradient size."
            test_solver_sanity_check = .false.
        end if
        call solver_sanity_check(settings, 4, grad, error)
        if (error == 0) then
            write(stderr, *) "test_solver_sanity_check failed: Error not thrown "// &
                "for correct incorrect gradient size."
            test_solver_sanity_check = .false.
        end if

        ! check if subsystem solver is correctly checked
        settings%subsystem_solver = "davidson"
        call solver_sanity_check(settings, 3, grad, error)
        if (error /= 0) then
            write(stderr, *) "test_solver_sanity_check failed: Error thrown for "// &
                "davidson subsystem solver."
            test_solver_sanity_check = .false.
        end if
        settings%subsystem_solver = "jacobi-davidson"
        call solver_sanity_check(settings, 3, grad, error)
        if (error /= 0) then
            write(stderr, *) "test_solver_sanity_check failed: Error thrown for "// &
                "jacobi-davidson subsystem solver."
            test_solver_sanity_check = .false.
        end if
        settings%subsystem_solver = "tcg"
        call solver_sanity_check(settings, 3, grad, error)
        if (error /= 0) then
            write(stderr, *) "test_solver_sanity_check failed: Error thrown for "// &
                "tcg subsystem solver."
            test_solver_sanity_check = .false.
        end if
        settings%subsystem_solver = "unknown"
        call solver_sanity_check(settings, 3, grad, error)
        if (error == 0) then
            write(stderr, *) "test_solver_sanity_check failed: Error not thrown "// &
                "for unknown subsystem solver."
            test_solver_sanity_check = .false.
        end if

    end function test_solver_sanity_check

    logical(c_bool) function test_stability_sanity_check() bind(C)
        !
        ! this function tests the subroutine which performs a sanity check for the 
        ! stability check
        !
        use opentrustregion, only: stability_settings_type, stability_sanity_check

        type(stability_settings_type) :: settings
        integer(ip) :: error

        ! assume tests pass
        test_stability_sanity_check = .true.

        ! setup settings object
        call setup_settings(settings)

        ! check if number of random trial vectors is reduced correctly
        settings%n_random_trial_vectors = 3
        call stability_sanity_check(settings, 3, error)
        if (settings%n_random_trial_vectors /= 1) then
            write(stderr, *) "test_stability_sanity_check failed: Number of random "// &
                "trial not correctly set."
            test_stability_sanity_check = .false.
        end if

        ! check if subsystem solver is correctly checked
        settings%diag_solver = "davidson"
        call stability_sanity_check(settings, 3, error)
        if (error /= 0) then
            write(stderr, *) "test_stability_sanity_check failed: Error thrown for "// &
                "davidson diagonalization solver."
            test_stability_sanity_check = .false.
        end if
        settings%diag_solver = "jacobi-davidson"
        call stability_sanity_check(settings, 3, error)
        if (error /= 0) then
            write(stderr, *) "test_stability_sanity_check failed: Error thrown for "// &
                "jacobi-davidson diagonalization solver."
            test_stability_sanity_check = .false.
        end if
        settings%diag_solver = "unknown"
        call stability_sanity_check(settings, 3, error)
        if (error == 0) then
            write(stderr, *) "test_stability_sanity_check failed: Error not thrown "// &
                "for unknown diagonalization solver."
            test_stability_sanity_check = .false.
        end if

    end function test_stability_sanity_check

    logical(c_bool) function test_level_shifted_davidson() bind(C)
        !
        ! this function tests the level-shifted Davidson subroutine
        !
        use opentrustregion, only: obj_func_type, hess_x_type, solver_settings_type, &
                                   level_shifted_davidson, trust_radius_shrink_ratio, &
                                   trust_radius_expand_ratio, &
                                   trust_radius_shrink_factor, &
                                   trust_radius_expand_factor

        integer(ip), parameter :: n_param = 6
        real(rp) :: func, grad_norm, trust_radius, mu, ratio, solution_norm
        real(rp), dimension(n_param) :: grad, h_diag, solution
        integer(ip) :: i, imicro, imicro_jacobi_davidson, error
        procedure(obj_func_type), pointer :: obj_func_funptr
        procedure(hess_x_type), pointer :: hess_x_funptr
        type(solver_settings_type) :: settings
        logical :: jacobi_davidson_started, max_precision_reached

        ! assume tests pass
        test_level_shifted_davidson = .true.

        ! setup settings object
        call setup_settings(settings)

        ! initialize variables
        trust_radius = 0.4_rp
        obj_func_funptr => obj_func
        hess_x_funptr => hess_x_fun

        ! start in quadratic region near minimum
        curr_vars = [0.20_rp, 0.15_rp, 0.48_rp, 0.28_rp, 0.31_rp, 0.66_rp]
        func = hartmann6d_func(curr_vars)
        call hartmann6d_gradient(curr_vars, grad)
        grad_norm = norm2(grad)
        call hartmann6d_hessian(curr_vars)
        h_diag = [(hess(i, i), i=1, size(h_diag))]

        ! run level-shifted Davidson, check if error has occured, whether the level 
        ! shift vanishes and whether the solution stays within trust region and 
        ! describes the Newton step
        call level_shifted_davidson(func, grad, grad_norm, h_diag, n_param, &
                                    obj_func_funptr, hess_x_funptr, settings, &
                                    trust_radius, solution, mu, imicro, &
                                    imicro_jacobi_davidson, jacobi_davidson_started, &
                                    max_precision_reached, error)
        if (error /= 0) then
            write (stderr, *) "test_level_shifted_davidson failed: Produced error "// &
                "near minimum."
            test_level_shifted_davidson = .false.
        end if
        if (abs(mu) > tol) then
            write (stderr, *) "test_level_shifted_davidson failed: Level shift is "// &
                "not zero near minimum."
            test_level_shifted_davidson = .false.
        end if
        if (sum(abs(grad + hartmann6d_hess_x(solution))) > &
            settings%local_red_factor * grad_norm) then
            write (stderr, *) "test_level_shifted_davidson failed: Solution does "// &
                "not describe Newton step near minimum."
            test_level_shifted_davidson = .false.
        end if
        ratio = (hartmann6d_func(curr_vars + solution) - func) / &
                dot_product(solution, grad + 0.5_rp * hartmann6d_hess_x(solution))
        solution_norm = norm2(solution)
        if ((ratio < trust_radius_shrink_ratio .and. solution_norm > trust_radius &
             / trust_radius_shrink_factor) .or. &
            (trust_radius_shrink_ratio > ratio .and. ratio > trust_radius_expand_ratio &
             .and. solution_norm > trust_radius) .or. &
            (ratio > trust_radius_expand_ratio .and. solution_norm > trust_radius &
            / trust_radius_expand_factor)) then
            write (stderr, *) "test_level_shifted_davidson failed: Solution does "// &
                "not stay within trust region near minimum."
            test_level_shifted_davidson = .false.
        end if

        ! start near saddle point
        curr_vars = [0.35_rp, 0.59_rp, 0.48_rp, 0.40_rp, 0.31_rp, 0.32_rp]
        func = hartmann6d_func(curr_vars)
        call hartmann6d_gradient(curr_vars, grad)
        grad_norm = norm2(grad)
        call hartmann6d_hessian(curr_vars)
        h_diag = [(hess(i, i), i=1, size(h_diag))]
        trust_radius = 0.4_rp

        ! run level-shifted Davidson, check if error has occured, whether the level 
        ! shift is negative and whether the solution lies at the trust region boundary 
        ! and describes a level-shifted Newton step
        call level_shifted_davidson(func, grad, grad_norm, h_diag, n_param, &
                                    obj_func_funptr, hess_x_funptr, settings, &
                                    trust_radius, solution, mu, imicro, &
                                    imicro_jacobi_davidson, jacobi_davidson_started, &
                                    max_precision_reached, error)
        if (error /= 0) then
            write (stderr, *) "test_level_shifted_davidson failed: Produced error "// &
                "near saddle point."
            test_level_shifted_davidson = .false.
        end if
        if (mu >= 0.0_rp) then
            write (stderr, *) "test_level_shifted_davidson failed: Level shift is "// &
                "not negative near saddle point."
            test_level_shifted_davidson = .false.
        end if
        if (sum(abs(grad + hartmann6d_hess_x(solution) - mu * solution)) > &
            settings%global_red_factor * grad_norm) then
            write (stderr, *) "test_level_shifted_davidson failed: Solution does "// &
                "not describe level-shifted Newton step near saddle point."
            test_level_shifted_davidson = .false.
        end if
        ratio = (hartmann6d_func(curr_vars + solution) - func) / &
                dot_product(solution, grad + 0.5_rp * hartmann6d_hess_x(solution))
        solution_norm = norm2(solution)
        if ((trust_radius_shrink_ratio > ratio .and. &
             abs(solution_norm - (trust_radius / trust_radius_shrink_factor) ** 2) > &
             tol) .or. &
            (trust_radius_shrink_ratio > ratio .and. ratio > trust_radius_expand_ratio &
             .and. abs(solution_norm - trust_radius ** 2) < tol) .or. &
            (ratio > trust_radius_expand_ratio .and. &
             abs(solution_norm - (trust_radius / trust_radius_expand_factor) ** 2) < &
             tol)) then
            write (stderr, *) "test_level_shifted_davidson failed: Solution does "// &
                "not lie at trust region boundary near saddle point."
            test_level_shifted_davidson = .false.
        end if

        ! test Jacobi-Davidson near saddle point
        settings%subsystem_solver = "jacobi-davidson"
        trust_radius = 0.4_rp

        ! run level-shifted Jacobi-Davidson, check if error has occured, whether the 
        ! level shift is negative and whether the solution lies at the trust region 
        ! boundary and describes a level-shifted Newton step
        call level_shifted_davidson(func, grad, grad_norm, h_diag, n_param, &
                                    obj_func_funptr, hess_x_funptr, settings, &
                                    trust_radius, solution, mu, imicro, &
                                    imicro_jacobi_davidson, jacobi_davidson_started, &
                                    max_precision_reached, error)
        if (error /= 0) then
            write (stderr, *) "test_level_shifted_davidson failed: Produced error "// &
                "near saddle point with Jacobi-Davidson solver."
            test_level_shifted_davidson = .false.
        end if
        if (mu >= 0.0_rp) then
            write (stderr, *) "test_level_shifted_davidson failed: Level shift is "// &
                "not negative near saddle point with Jacobi-Davidson solver."
            test_level_shifted_davidson = .false.
        end if
        if (sum(abs(grad + hartmann6d_hess_x(solution) - mu * solution)) > &
            settings%global_red_factor * grad_norm) then
            write (stderr, *) "test_level_shifted_davidson failed: Solution does "// &
                "not describe level-shifted Newton step near saddle point with "// &
                "Jacobi-Davidson solver."
            test_level_shifted_davidson = .false.
        end if
        ratio = (hartmann6d_func(curr_vars + solution) - func) / &
                dot_product(solution, grad + 0.5_rp * hartmann6d_hess_x(solution))
        solution_norm = norm2(solution)
        if ((trust_radius_shrink_ratio > ratio .and. &
             abs(solution_norm - (trust_radius / trust_radius_shrink_factor) ** 2) > &
             tol) .or. &
            (trust_radius_shrink_ratio > ratio .and. ratio > trust_radius_expand_ratio &
             .and. abs(solution_norm - trust_radius ** 2) < tol) .or. &
            (ratio > trust_radius_expand_ratio .and. &
             abs(solution_norm - (trust_radius / trust_radius_expand_factor) ** 2) < &
             tol)) then
            write (stderr, *) "test_level_shifted_davidson failed: Solution does "// &
                "not lie at trust region boundary near saddle point with "// &
                "Jacobi-Davidson solver."
            test_level_shifted_davidson = .false.
        end if

    end function test_level_shifted_davidson

    logical(c_bool) function test_truncated_conjugate_gradient() bind(C)
        !
        ! this function tests the truncated conjugate gradient subroutine
        !
        use opentrustregion, only: obj_func_type, hess_x_type, solver_settings_type, &
                                   truncated_conjugate_gradient, &
                                   trust_radius_shrink_ratio, &
                                   trust_radius_expand_ratio, &
                                   trust_radius_shrink_factor, &
                                   trust_radius_expand_factor

        integer(ip), parameter :: n_param = 6
        real(rp) :: func, trust_radius, ratio, solution_norm
        real(rp), dimension(n_param) :: grad, h_diag, solution
        integer(ip) :: i, imicro, error
        procedure(obj_func_type), pointer :: obj_func_funptr
        procedure(hess_x_type), pointer :: hess_x_funptr
        type(solver_settings_type) :: settings
        logical :: max_precision_reached
        real(rp), parameter :: h_diag_floor = 1e-10_rp

        ! assume tests pass
        test_truncated_conjugate_gradient = .true.

        ! setup settings object
        call setup_settings(settings)
        settings%n_micro = 50

        ! initialize variables
        trust_radius = 0.4_rp
        obj_func_funptr => obj_func
        hess_x_funptr => hess_x_fun

        ! start in quadratic region near minimum
        curr_vars = [0.20_rp, 0.15_rp, 0.48_rp, 0.28_rp, 0.31_rp, 0.66_rp]
        func = hartmann6d_func(curr_vars)
        call hartmann6d_gradient(curr_vars, grad)
        call hartmann6d_hessian(curr_vars)
        h_diag = [(hess(i, i), i=1, size(h_diag))]

        ! run truncated conjugate gradient, check whether the solution lies at the 
        ! trust region boundary and reduces the function value
        call truncated_conjugate_gradient(func, grad, h_diag, n_param, &
                                          obj_func_funptr, hess_x_funptr, settings, &
                                          trust_radius, solution, imicro, &
                                          max_precision_reached, error)
        if (error /= 0) then
            write (stderr, *) "test_truncated_jacobi_davidson failed: Produced "// &
                "error near minimum."
            test_truncated_conjugate_gradient = .false.
        end if
        ratio = (hartmann6d_func(curr_vars + solution) - func) / &
                dot_product(solution, grad + 0.5_rp * hartmann6d_hess_x(solution))
        if (ratio <= 0.0_rp) then
            write (stderr, *) "test_truncated_conjugate_gradient failed: Solution "// &
                "does not reduce function value near minimum."
            test_truncated_conjugate_gradient = .false.
        end if
        solution_norm = dot_product(solution, solution / max(abs(h_diag), h_diag_floor))
        if ((ratio < trust_radius_shrink_ratio .and. solution_norm > (trust_radius / &
             trust_radius_shrink_factor) ** 2) .or. &
            (trust_radius_shrink_ratio > ratio .and. ratio > trust_radius_expand_ratio &
             .and. solution_norm > trust_radius ** 2) .or. &
            (ratio > trust_radius_expand_ratio .and. solution_norm > (trust_radius / &
             trust_radius_expand_factor) ** 2)) then
            write (stderr, *) "test_truncated_conjugate_gradient failed: Solution "// &
                "does not stay within trust region near minimum."
            test_truncated_conjugate_gradient = .false.
        end if

        ! start near saddle point
        curr_vars = [0.35_rp, 0.59_rp, 0.48_rp, 0.40_rp, 0.31_rp, 0.32_rp]
        func = hartmann6d_func(curr_vars)
        call hartmann6d_gradient(curr_vars, grad)
        call hartmann6d_hessian(curr_vars)
        h_diag = [(hess(i, i), i=1, size(h_diag))]
        trust_radius = 0.4_rp

        ! run truncated conjugate gradient, check whether the solution lies at the 
        ! trust region boundary and reduces the function value
        call truncated_conjugate_gradient(func, grad, h_diag, n_param, &
                                          obj_func_funptr, hess_x_funptr, settings, &
                                          trust_radius, solution, imicro, &
                                          max_precision_reached, error)
        if (error /= 0) then
            write (stderr, *) "test_truncated_jacobi_davidson failed: Produced "// &
                "error near saddle point."
            test_truncated_conjugate_gradient = .false.
        end if
        ratio = (hartmann6d_func(curr_vars + solution) - func) / &
                dot_product(solution, grad + 0.5_rp * hartmann6d_hess_x(solution))
        if (ratio <= 0.0_rp) then
            write (stderr, *) "test_truncated_conjugate_gradient failed: Solution "// &
                "does not reduce function value near saddle point."
            test_truncated_conjugate_gradient = .false.
        end if
        solution_norm = dot_product(solution, solution / max(abs(h_diag), h_diag_floor))
        if ((trust_radius_shrink_ratio > ratio .and. &
             abs(solution_norm - (trust_radius / trust_radius_shrink_factor) ** 2) > &
             tol) .or. &
            (trust_radius_shrink_ratio > ratio .and. ratio > trust_radius_expand_ratio &
             .and. abs(solution_norm - trust_radius ** 2) < tol) .or. &
            (ratio > trust_radius_expand_ratio .and. &
             abs(solution_norm - (trust_radius / trust_radius_expand_factor) ** 2) < &
             tol)) then
            write (stderr, *) "test_truncated_conjugate_gradient failed: Solution "// &
                "does not lie at trust region boundary near saddle point."
            test_truncated_conjugate_gradient = .false.
        end if

    end function test_truncated_conjugate_gradient

    logical(c_bool) function test_add_error_origin() bind(C)
        !
        ! this function tests the subroutine that adds the error origin to an error 
        ! code
        !
        use opentrustregion, only: solver_settings_type, add_error_origin

        type(solver_settings_type) :: settings
        integer(ip) :: error

        ! assume tests pass
        test_add_error_origin = .true.

        ! setup settings object
        call setup_settings(settings)

        ! check if subroutine adds error origin correctly if no origin is present
        error = 1
        call add_error_origin(error, 100, settings)
        if (error /= 101) then
            write (stderr, *) "test_add_error_origin failed: Error origin not "// &
                "correctly added."
            test_add_error_origin = .false.
        end if

        ! check if subroutine skips adding error origin if origin is already present
        call add_error_origin(error, 100, settings)
        if (error /= 101) then
            write (stderr, *) "test_add_error_origin failed: Error code modified "// &
                "even though error origin is already present."
            test_add_error_origin = .false.
        end if

        ! check if subroutine does not modify error code when no error is encountered
        error = 0
        call add_error_origin(error, 100, settings)
        if (error /= 0) then
            write (stderr, *) "test_add_error_origin failed: Error code modified "// &
                "even though error code of zero was passed."
            test_add_error_origin = .false.
        end if

        ! check if subroutine raises error for invalid error code
        error = -1
        call add_error_origin(error, 100, settings)
        if (error /= 101) then
            write (stderr, *) "test_add_error_origin failed: Error code not "// &
                "correctly returned for invalid (negative) error code."
            test_add_error_origin = .false.
        end if

    end function test_add_error_origin

    logical(c_bool) function test_string_to_lowercase() bind(C)
        !
        ! this function tests the function that transfers strings to lowercase
        !
        use opentrustregion, only: string_to_lowercase

        character(*), parameter :: input  = 'OpenTrustRegion123!', &
                                   expect = 'opentrustregion123!'

        ! assume tests pass
        test_string_to_lowercase = .true.

        ! test transfoer to lowercase
        if (string_to_lowercase(input) /= expect) then
            write (stderr, *) "test_string_to_lowercase failed: String not "// &
                "correctly transferred to lowercase."
            test_string_to_lowercase = .false.
        end if

    end function test_string_to_lowercase

end module opentrustregion_unit_tests
