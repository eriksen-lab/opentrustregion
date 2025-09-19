! Copyright (C) 2025- Jonas Greiner
!
! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.

module opentrustregion_unit_tests

    use opentrustregion, only: rp, ip, stderr, settings_type
    use iso_c_binding, only: c_bool

    implicit none

    real(rp), parameter :: tol = 1.d-10

    ! parameters for 6D Hartmann function
    real(rp), parameter :: alpha(4) = [1.0d0, 1.2d0, 3.0d0, 3.2d0]
    real(rp), parameter :: A(4, 6) = reshape([10.d0, 0.05d0, 3.0d0, 17.0d0, &
                                              3.d0, 10.d0, 3.5d0, 8.0d0, &
                                              17.d0, 17.d0, 1.7d0, 0.05d0, &
                                              3.5d0, 0.1d0, 10.d0, 10.0d0, &
                                              1.7d0, 8.d0, 17.d0, 0.1d0, &
                                              8.d0, 14.d0, 8.d0, 14.d0], [4, 6])
    real(rp), parameter :: P(4, 6) = reshape([0.1312d0, 0.2329d0, 0.2348d0, 0.4047d0, &
                                              0.1696d0, 0.4135d0, 0.1451d0, 0.8828d0, &
                                              0.5569d0, 0.8307d0, 0.3522d0, 0.8732d0, &
                                              0.0124d0, 0.3736d0, 0.2883d0, 0.5743d0, &
                                              0.8283d0, 0.1004d0, 0.3047d0, 0.1091d0, &
                                              0.5886d0, 0.9991d0, 0.6650d0, 0.0381d0], &
                                             [4, 6])

    ! stationary points of 6D Hartmann function
    real(rp), parameter :: minimum1(6) = [0.20168951d0, 0.15001069d0, 0.47687398d0, &
                                          0.27533243d0, 0.31165162d0, 0.65730053d0]
    real(rp), parameter :: minimum2(6) = [0.40465313d0, 0.88244493d0, 0.84610160d0, &
                                          0.57398969d0, 0.13892673d0, 0.03849589d0]
    real(rp), parameter :: saddle_point(6) = [0.35278250d0, 0.59374767d0, &
                                              0.47631257d0, 0.40058250d0, &
                                              0.31111531d0, 0.32397158d0]

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
        real(rp), intent(out) :: grad(size(vars))
        real(rp) :: exp_term(4)
        integer(ip) :: i, j

        do i = 1, 4
            exp_term(i) = exp(-sum(A(i, :)*(vars - P(i, :))**2))
        end do

        do j = 1, size(vars)
            grad(j) = sum(2.d0*alpha*A(:, j)*(vars(j) - P(:, j))*exp_term)
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
            hess(i, i) = 2.d0*sum(alpha*A(:, i)*exp_term* &
                                  (1.d0 - 2.d0*A(:, i)*(vars(i) - P(:, i))**2))
            do j = 1, i - 1
                hess(i, j) = -4.d0*sum(alpha*A(:, i)*A(:, j)*(vars(i) - P(:, i))* &
                                       (vars(j) - P(:, j))*exp_term)
                hess(j, i) = hess(i, j)
            end do
        end do

    end subroutine hartmann6d_hessian

    function hartmann6d_hess_x(x)
        !
        ! this function defiones the Hessian linear transformation operation for the
        ! Hartmann 6D function
        !
        real(rp), intent(in) :: x(:)

        real(rp) :: hartmann6d_hess_x(size(x))

        hartmann6d_hess_x = matmul(hess, x)

    end function hartmann6d_hess_x

    function hess_x(x, error)
        !
        ! this function describes the Hessian linear transformation operation for the
        ! Hartmann 6D function
        !
        real(rp), intent(in) :: x(:)
        logical, intent(out) :: error

        real(rp) :: hess_x(size(x))

        ! initialize error flag
        error = .false.

        hess_x = hartmann6d_hess_x(x)

    end function hess_x

    function obj_func(delta_vars, error) result(func)
        !
        ! this function describes the objective function evaluation for the Hartmann
        ! 6D function
        !
        real(rp), intent(in) :: delta_vars(:)
        logical, intent(out) :: error

        real(rp) :: func

        ! initialize error flag
        error = .false.

        func = hartmann6d_func(curr_vars + delta_vars)

    end function obj_func

    subroutine update_orbs(delta_vars, func, grad, h_diag, hess_x_funptr, error)
        !
        ! this function describes the orbital update equivalent for the Hartmann 6D
        ! function
        !
        use opentrustregion, only: hess_x_type

        real(rp), intent(in) :: delta_vars(:)
        real(rp), intent(out) :: func, grad(:), h_diag(:)
        procedure(hess_x_type), intent(out), pointer :: hess_x_funptr
        logical, intent(out) :: error

        integer(ip) :: i

        ! initialize error flag
        error = .false.

        ! update variables
        curr_vars = curr_vars + delta_vars

        ! evaluate function, calculate gradient and Hessian diagonal and define
        ! Hessian linear transformation
        func = hartmann6d_func(curr_vars)
        call hartmann6d_gradient(curr_vars, grad)
        call hartmann6d_hessian(curr_vars)
        h_diag = [(hess(i, i), i=1, size(h_diag))]
        hess_x_funptr => hess_x

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
        class(settings_type), intent(inout) :: settings

        settings%verbose = 3
        settings%logger => logger

    end subroutine setup_settings

    logical(c_bool) function test_solver() bind(C)
        !
        ! this function tests the solver subroutine
        !
        use opentrustregion, only: update_orbs_type, obj_func_type, solver, &
                                   solver_conv_tol_default

        integer(ip), parameter :: n_param = 6
        logical :: error
        real(rp) :: final_grad(n_param)
        procedure(update_orbs_type), pointer :: update_orbs_funptr
        procedure(obj_func_type), pointer :: obj_func_funptr

        ! assume tests pass
        test_solver = .true.

        ! start in quadratic region near minimum
        curr_vars = [0.20d0, 0.15d0, 0.48d0, 0.28d0, 0.31d0, 0.66d0]
        update_orbs_funptr => update_orbs
        obj_func_funptr => obj_func

        ! run solver, check if error has occured and check whether gradient is zero and 
        ! agrees with correct minimum
        call solver(update_orbs_funptr, obj_func_funptr, n_param, error)
        if (error) then
            write (stderr, *) "test_solver failed: Produced error."
            test_solver = .false.
        end if
        call hartmann6d_gradient(curr_vars, final_grad)
        if (norm2(final_grad)/sqrt(real(n_param, kind=rp)) > &
            solver_conv_tol_default) then
            write (stderr, *) "test_solver failed: Solver did not find stationary "// &
                "point."
            test_solver = .false.
        end if
        if (any(abs(curr_vars - minimum1) > 1d-8)) then
            write (stderr, *) "test_solver failed: Solver did not find correct minimum."
            test_solver = .false.
        end if

        ! start near saddle point
        curr_vars = [0.35d0, 0.59d0, 0.48d0, 0.40d0, 0.31d0, 0.32d0]
        update_orbs_funptr => update_orbs
        obj_func_funptr => obj_func

        ! run solver, check if error has occured and check whether gradient is zero and 
        ! agrees with correct minimum
        call solver(update_orbs_funptr, obj_func_funptr, n_param, error)
        if (error) then
            write (stderr, *) "test_solver failed: Produced error."
            test_solver = .false.
        end if
        call hartmann6d_gradient(curr_vars, final_grad)
        if (norm2(final_grad)/sqrt(real(n_param, kind=rp)) > &
            solver_conv_tol_default) then
            write (stderr, *) "test_solver failed: Solver did not find stationary "// &
                "point."
            test_solver = .false.
        end if
        if (any(abs(curr_vars - minimum1) > 1d-6) .and. &
            any(abs(curr_vars - minimum2) > 1d-6)) then
            write (stderr, *) "test_solver failed: Solver did not find correct minimum."
            test_solver = .false.
        end if

        ! start at saddle point
        curr_vars = saddle_point
        update_orbs_funptr => update_orbs
        obj_func_funptr => obj_func

        ! run solver, check if error has occured and check whether gradient is zero and 
        ! agrees with correct minimum
        call solver(update_orbs_funptr, obj_func_funptr, n_param, error)
        if (error) then
            write (stderr, *) "test_solver failed: Produced error."
            test_solver = .false.
        end if
        call hartmann6d_gradient(curr_vars, final_grad)
        if (norm2(final_grad)/sqrt(real(n_param, kind=rp)) > &
            solver_conv_tol_default) then
            write (stderr, *) "test_solver failed: Solver did not find stationary "// &
                "point."
            test_solver = .false.
        end if
        if (any(abs(curr_vars - minimum1) > 1d-6) .and. &
            any(abs(curr_vars - minimum2) > 1d-6)) then
            write (stderr, *) "test_solver failed: Solver did not find minimum."
            test_solver = .false.
        end if

    end function test_solver

    logical(c_bool) function test_stability_check() bind(C)
        !
        ! this function tests the stability check subroutine
        !
        use opentrustregion, only: hess_x_type, stability_check

        real(rp) :: vars(6), h_diag(6), direction(6)
        procedure(hess_x_type), pointer :: hess_x_funptr
        logical :: stable, error
        integer(ip) :: i

        ! assume tests pass
        test_stability_check = .true.

        ! start at minimum and determine Hessian diagonal and define Hessian linear
        ! transformation
        vars = minimum1
        call hartmann6d_hessian(vars)
        h_diag = [(hess(i, i), i=1, size(h_diag))]
        hess_x_funptr => hess_x

        ! run stability, check if error has occured check and determine whether minimum 
        ! is stable and the returned direction vanishes
        call stability_check(h_diag, hess_x_funptr, stable, direction, error)
        if (error) then
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
        hess_x_funptr => hess_x

        ! run stability check, check if error has occured and determine whether saddle 
        ! point is unstable and the returned direction is correct
        call stability_check(h_diag, hess_x_funptr, stable, direction, error)
        if (error) then
            write (stderr, *) "test_stability_check failed: Produced error."
            test_stability_check = .false.
        end if
        if (stable) then
            write (stderr, *) "test_stability_check failed: Stability check "// &
                "incorrectly classifies stability of saddle point."
            test_stability_check = .false.
        end if
        if (abs(abs(dot_product(direction, &
                                [-0.173375920238d0, -0.518489821791d0, &
                                 -6.432848975252d-3, -0.340127852882d0, &
                                 3.066460316955d-3, 0.765095650196d0])) - 1.d0) > tol) &
            then
            write (stderr, *) "test_stability_check failed: Stability check does "// &
                "not return correct direction for saddle point."
            test_stability_check = .false.
        end if

    end function test_stability_check

    logical(c_bool) function test_newton_step() bind(C)
        !
        ! this function tests the Newton step subroutine
        !
        use opentrustregion, only: newton_step

        type(settings_type) :: settings
        real(rp) :: red_space_basis(6, 3), vars(6), grad(6), grad_norm, &
                    aug_hess(4, 4), solution(6), red_space_solution(3)
        integer(ip) :: i, j
        logical :: error

        ! assume tests pass
        test_newton_step = .true.

        ! setup settings object
        call setup_settings(settings)

        ! defined a reduced space basis
        red_space_basis = &
            reshape([1.d0/sqrt(2.d0), -1.d0/sqrt(2.d0), 0.d0, 0.d0, 0.d0, 0.d0, &
                     1.d0/sqrt(6.d0), -1.d0/sqrt(6.d0), -2.d0/sqrt(6.d0), 0.d0, 0.d0, &
                     0.d0, 1.d0/sqrt(12.d0), -1.d0/sqrt(12.d0), 1.d0/sqrt(12.d0), &
                     -3.d0/sqrt(12.d0), 0.d0, 0.d0], [6, 3])

        ! point in quadratic region near minimum
        vars = [0.20d0, 0.15d0, 0.48d0, 0.28d0, 0.31d0, 0.66d0]

        ! calculate gradient and Hessian to define augmented Hessian
        call hartmann6d_gradient(vars, grad)
        grad_norm = norm2(grad)
        call hartmann6d_hessian(vars)
        aug_hess = 0.d0
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
        if (error) then
            write (stderr, *) "test_newton_step failed: Produced error."
            test_newton_step = .false.
        end if
        if (any(abs(red_space_solution - &
                    [-2.555959788079d-2, 1.565498761914d-2, 4.727080080611d-3]) > &
                tol)) then
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
        use opentrustregion, only: bisection

        type(settings_type) :: settings
        real(rp) :: red_space_basis(6, 3), vars(6), grad(6), grad_norm, &
                    aug_hess(4, 4), solution(6), red_space_solution(3), trust_radius, mu
        integer(ip) :: i, j
        logical :: bracketed, error

        ! assume tests pass
        test_bisection = .true.

        ! setup settings object
        call setup_settings(settings)

        ! defined a reduced space basis
        red_space_basis = &
            reshape([1.d0/sqrt(2.d0), -1.d0/sqrt(2.d0), 0.d0, 0.d0, 0.d0, 0.d0, &
                     1.d0/sqrt(6.d0), -1.d0/sqrt(6.d0), -2.d0/sqrt(6.d0), 0.d0, 0.d0, &
                     0.d0, 1.d0/sqrt(12.d0), -1.d0/sqrt(12.d0), 1.d0/sqrt(12.d0), &
                     -3.d0/sqrt(12.d0), 0.d0, 0.d0], [6, 3])

        ! choose target trust radius
        trust_radius = 0.4d0

        ! point with strong negative curvature
        vars = [0.29d0, 0.47d0, 0.66d0, 0.41d0, 0.23d0, 0.26d0]

        ! calculate gradient and Hessian to define augmented Hessian
        call hartmann6d_gradient(vars, grad)
        grad_norm = norm2(grad)
        call hartmann6d_hessian(vars)
        aug_hess = 0.d0
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
        if (error) then
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
        if (any(abs(red_space_solution - &
                    [-0.483593823965d0, 0.482091645228d0, 0.153783319727d0]) > tol)) &
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
        vars = [0.20d0, 0.15d0, 0.48d0, 0.28d0, 0.31d0, 0.66d0]

        ! calculate gradient and Hessian to define augmented Hessian
        call hartmann6d_gradient(vars, grad)
        grad_norm = norm2(grad)
        call hartmann6d_hessian(vars)
        aug_hess = 0.d0
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
        if (error) then
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
        use opentrustregion, only: obj_func_type, bracket

        type(settings_type) :: settings
        procedure(obj_func_type), pointer :: obj_func_funptr
        real(rp) :: vars(6), lower, upper, n
        logical :: error

        ! assume tests pass
        test_bracket = .true.

        ! setup settings object
        call setup_settings(settings)

        ! define procedure pointer
        obj_func_funptr => obj_func

        ! define direction
        vars = [0.20d0, 0.15d0, 0.48d0, 0.28d0, 0.31d0, 0.66d0]

        ! define lower and upper bound
        lower = 0.d0
        upper = 0.5d0

        ! perform bracket and determine if new point decreases objective function in
        ! comparison to lower and upper bound
        n = bracket(obj_func_funptr, vars, lower, upper, settings, error)
        if (error) then
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
        if (error) then
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
        allocate (matrix(2, 2))
        matrix = reshape([1.d0, 2.d0, &
                          2.d0, 3.d0], [2, 2])
        vector = [4.d0, 5.d0, 6.d0]

        ! initialize expected matrix
        expected = reshape([1.d0, 2.d0, 4.d0, &
                            2.d0, 3.d0, 5.d0, &
                            4.d0, 5.d0, 6.d0], [3, 3])

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
        deallocate (matrix)

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
        allocate (matrix(3, 2))
        matrix = reshape([1.d0, 2.d0, 3.d0, &
                          4.d0, 5.d0, 6.d0], [3, 2])
        new_col = [7.d0, 8.d0, 9.d0]

        ! initialize expected matrix
        expected = reshape([1.d0, 2.d0, 3.d0, &
                            4.d0, 5.d0, 6.d0, &
                            7.d0, 8.d0, 9.d0], [3, 3])

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
        deallocate (matrix)

    end function test_add_column

    logical(c_bool) function test_symm_mat_min_eig() bind(C)
        !
        ! this function tests the subroutine for determining the minimum eigenvalue and
        ! corresponding eigenvector for a symmetric matrix
        !
        use opentrustregion, only: symm_mat_min_eig

        type(settings_type) :: settings
        real(rp) :: matrix(3, 3)
        real(rp) :: eigval, eigvec(3)
        logical :: error

        ! assume tests pass
        test_symm_mat_min_eig = .true.

        ! setup settings object
        call setup_settings(settings)

        ! initialize symmetric matrix
        matrix = reshape([3.d0, 1.d0, 1.d0, &
                          1.d0, 4.d0, 2.d0, &
                          1.d0, 2.d0, 5.d0], [3, 3])

        ! call routine and determine if lowest eigenvalue and corresponding eigenvector
        ! are found
        call symm_mat_min_eig(matrix, eigval, eigvec, settings, error)
        if (error) then
            write (stderr, *) "test_symm_mat_min_eig failed: Produced error."
            test_symm_mat_min_eig = .false.
        end if
        if (abs(eigval - 2.30797852837d0) > tol) then
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
        use opentrustregion, only: min_eigval

        type(settings_type) :: settings
        real(rp) :: matrix(3, 3), matrix_min_eigval
        logical :: error

        ! assume tests pass
        test_min_eigval = .true.

        ! setup settings object
        call setup_settings(settings)

        ! initialize symmetric matrix
        matrix = reshape([3.d0, 1.d0, 1.d0, &
                          1.d0, 4.d0, 2.d0, &
                          1.d0, 2.d0, 5.d0], [3, 3])

        ! call function and determine if lowest eigenvalue is found
        matrix_min_eigval = min_eigval(matrix, settings, error)
        if (error) then
            write (stderr, *) "test_min_eigval failed: Produced error."
            test_min_eigval = .false.
        end if
        if (abs(matrix_min_eigval - 2.30797852837d0) > tol) then
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
        use opentrustregion, only: generate_trial_vectors

        type(settings_type) :: settings
        real(rp), allocatable :: red_space_basis(:, :)
        real(rp) :: grad(4), h_diag(4), grad_norm
        logical :: error
        integer(ip) :: i, j

        ! assume tests pass
        test_generate_trial_vectors = .true.

        ! setup settings object
        call setup_settings(settings)
        settings%n_random_trial_vectors = 2

        ! define gradient
        grad = [1.d0, 2.d0, 3.d0, 4.d0]
        grad_norm = norm2(grad)

        ! define all positive Hessian diagonal elements
        h_diag = [1.d0, 2.d0, 3.d0, 4.d0]

        ! generate trial vectors and determine whether function returns the correct
        ! number of orthonormal trial vectors
        red_space_basis = generate_trial_vectors(grad, grad_norm, h_diag, settings, &
                                                 error)
        if (error) then
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
        deallocate (red_space_basis)

        ! define Hessian diagonal with negative elements
        h_diag = [-1.d0, 2.d0, 3.d0, 4.d0]

        ! generate trial vectors and determine whether function returns the correct
        ! number of orthonormal trial vectors
        red_space_basis = generate_trial_vectors(grad, grad_norm, h_diag, settings, &
                                                 error)
        if (error) then
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
        deallocate (red_space_basis)

    end function test_generate_trial_vectors

    logical(c_bool) function test_generate_random_trial_vectors() bind(C)
        !
        ! this function tests the function which generates random trial vectors for the
        ! Davidson procedure
        !
        use opentrustregion, only: generate_random_trial_vectors

        type(settings_type) :: settings
        real(rp), allocatable :: red_space_basis(:, :)
        logical :: error
        integer(ip) :: i, j

        ! assume tests pass
        test_generate_random_trial_vectors = .true.

        ! setup settings object
        call setup_settings(settings)
        settings%n_random_trial_vectors = 2

        ! allocate reduced space basis and set first normalized basis vector
        allocate (red_space_basis(4, 3))
        red_space_basis(:, 1) = [1.d0, 2.d0, 3.d0, 4.d0]
        red_space_basis(:, 1) = red_space_basis(:, 1) / norm2(red_space_basis(:, 1))

        ! generate trial vectors and determine whether function returns orthonormal 
        ! trial vectors
        call generate_random_trial_vectors(red_space_basis, settings, error)
        if (error) then
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
        deallocate (red_space_basis)

    end function test_generate_random_trial_vectors

    logical(c_bool) function test_gram_schmidt() bind(C)
        !
        ! this function tests the Gram-Schmidt subroutine which orthonormalizes a 
        ! vector to a given basis
        !
        use opentrustregion, only: gram_schmidt

        type(settings_type) :: settings
        real(rp), dimension(4) :: vector, lin_trans_vector
        real(rp), dimension(2) :: vector_small
        real(rp) :: space(4, 2), symm_matrix(4, 4), lin_trans_space(4, 2), &
                    space_small(2, 2)
        logical :: error
        character(100) :: line

        ! assume tests pass
        test_gram_schmidt = .true.

        ! setup settings object
        call setup_settings(settings)

        ! define vector to be orthogonalized and space
        vector = [1.d0, 2.d0, 3.d0, 4.d0]
        space(:, 1) = [0.d0, 1.d0, 0.d0, 0.d0]
        space(:, 2) = [0.d0, 0.d0, 1.d0, 0.d0]

        ! perform Gram-Schmidt orthogonalization and determine whether added vector is
        ! orthonormalized
        call gram_schmidt(vector, space, settings, error)
        if (error) then
            write (stderr, *) "test_gram_schmidt failed: Produced error."
            test_gram_schmidt = .false.
        end if
        if (abs(dot_product(vector, space(:, 1))) > tol .or. &
            abs(dot_product(vector, space(:, 2))) > tol) then
            write (stderr, *) "test_gram_schmidt failed: Added vector not orthogonal."
            test_gram_schmidt = .false.
        end if
        if (abs(norm2(vector) - 1.d0) > tol) then
            write (stderr, *) "test_gram_schmidt failed: Added vector not normalized."
            test_gram_schmidt = .false.
        end if

        ! define vector to be orthogonalized and space
        vector = [1.d0, 2.d0, 3.d0, 4.d0]
        space(:, 1) = [0.d0, 1.d0, 0.d0, 0.d0]
        space(:, 2) = [0.d0, 0.d0, 1.d0, 0.d0]

        ! define symmetric linear transformation and corresponding vector and space
        symm_matrix = reshape([ 1.d0, -5.d0,  8.d0,  0.d0, &
                               -5.d0,  2.d0, -6.d0,  9.d0, &
                                8.d0, -6.d0,  3.d0, -7.d0, &
                                0.d0,  9.d0, -7.d0,  4.d0], &
                              shape(symm_matrix), order=[2,1])
        lin_trans_vector = matmul(symm_matrix, vector)
        lin_trans_space = matmul(symm_matrix, space)

        ! perform Gram-Schmidt orthogonalization and determine whether added vector is
        ! orthonormalized and linear transformation is correct
        call gram_schmidt(vector, space, settings, error, lin_trans_vector, &
                          lin_trans_space)
        if (error) then
            write (stderr, *) "test_gram_schmidt failed: Produced error."
            test_gram_schmidt = .false.
        end if
        if (abs(dot_product(vector, space(:, 1))) > tol .or. &
            abs(dot_product(vector, space(:, 2))) > tol) then
            write (stderr, *) "test_gram_schmidt failed: Added vector not orthogonal."
            test_gram_schmidt = .false.
        end if
        if (abs(norm2(vector) - 1.d0) > tol) then
            write (stderr, *) "test_gram_schmidt failed: Added vector not normalized."
            test_gram_schmidt = .false.
        end if
        if (abs(sum(abs(lin_trans_vector - matmul(symm_matrix, vector)))) > tol) then
            write (stderr, *) "test_gram_schmidt failed: Added linear "// &
                "transformation not correct."
            test_gram_schmidt = .false.
        end if

        ! define zero vector
        vector = [0.d0, 0.d0, 0.d0, 0.d0]

        ! perform Gram-Schmidt orthogonalization and determine if function correctly
        ! throws error
        call gram_schmidt(vector, space, settings, error)
        if ((.not. error) .or. (log_message /= " Vector passed to Gram-Schmidt "// &
                                "procedure is numerically zero.")) then
            write (stderr, *) "test_gram_schmidt failed: No error returned during "// &
                "orthogonalization for zero vector."
            test_gram_schmidt = .false.
        end if

        ! define vector in space that is already complete
        vector_small = [1.d0, 2.d0]
        space_small(:, 1) = [1.d0, 0.d0]
        space_small(:, 2) = [0.d0, 1.d0]

        ! perform Gram-Schmidt orthogonalization and determine if function correctly
        ! throws error
        call gram_schmidt(vector_small, space_small, settings, error)
        if (.not. error .or. (log_message /= " Number of vectors in Gram-Schmidt "// &
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
        use opentrustregion, only: logger_type, solver_settings_type, &
                                   init_solver_settings, solver_stability_default, &
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

        type(solver_settings_type) :: settings
        procedure(logger_type), pointer :: logger_funptr

        ! assume tests pass
        test_init_solver_settings = .true.

        ! call routine without optional arguments
        call init_solver_settings(settings)

        ! check if default values are correctly set
        if ((settings%stability .neqv. solver_stability_default) .or. &
            (settings%line_search .neqv. solver_line_search_default) .or. &
            (settings%davidson .neqv. solver_davidson_default) .or. &
            (settings%jacobi_davidson .neqv. solver_jacobi_davidson_default) .or. &
            (settings%prefer_jacobi_davidson .neqv. &
             solver_prefer_jacobi_davidson_default) .or. &
            abs(settings%conv_tol - solver_conv_tol_default) > tol .or. &
            settings%n_random_trial_vectors /= solver_n_random_trial_vectors_default &
            .or. abs(settings%start_trust_radius - solver_start_trust_radius_default) &
            > tol .or. settings%n_macro /= solver_n_macro_default .or. &
            settings%n_micro /= solver_n_micro_default .or. &
            abs(settings%global_red_factor - solver_global_red_factor_default) > tol &
            .or. abs(settings%local_red_factor - solver_local_red_factor_default) > &
            tol .or. settings%seed /= solver_seed_default .or. &
            settings%verbose /= solver_verbose_default .or. &
            associated(settings%logger)) then
            write (stderr, *) "test_init_solver_settings failed: Default arguments "// &
                "not set correctly."
            test_init_solver_settings = .false.
        end if

        ! set pointer for logging routine
        logger_funptr => logger

        ! call routine with optional arguments
        call init_solver_settings(settings, stability=.false., line_search=.true., &
                                  davidson=.false., jacobi_davidson=.false., &
                                  prefer_jacobi_davidson=.true., conv_tol=1.d-3, &
                                  n_random_trial_vectors=5, start_trust_radius=0.2d0, &
                                  n_macro=300, n_micro=200, global_red_factor=1.d-2, &
                                  local_red_factor=1.d-3, seed=33, verbose=3, &
                                  logger=logger_funptr)

        ! check if optional values are correctly set
        if ((settings%stability .neqv. .false.) .or. &
            (settings%line_search .neqv. .true.) .or. &
            (settings%davidson .neqv. .false.) .or. &
            (settings%jacobi_davidson .neqv. .false.) .or. &
            (settings%prefer_jacobi_davidson .neqv. .true.) .or. &
            abs(settings%conv_tol - 1.d-3) > tol .or. &
            settings%n_random_trial_vectors /= 5 .or. &
            abs(settings%start_trust_radius - 0.2d0) > tol .or. settings%n_macro &
            /= 300 .or. settings%n_micro /= 200 .or. &
            abs(settings%global_red_factor - 1.d-2) > tol .or. &
            abs(settings%local_red_factor - 1.d-3) > tol .or. &
            settings%seed /= 33 .or. settings%verbose /= 3 .or. &
            .not. associated(settings%logger)) then
            write (stderr, *) "test_init_solver_settings failed: Optional "// &
                "arguments not set correctly."
            test_init_solver_settings = .false.
        end if
        if (associated(settings%logger)) then
            call settings%logger("test")
            if (log_message /= "test") write (stderr, *) &
                "test_init_solver_settings failed: Optional logging subroutine "// &
                "not set correctly."
        end if

    end function test_init_solver_settings

    logical(c_bool) function test_init_stability_settings() bind(C)
        !
        ! this function tests the subroutine which initializes the stability check
        ! settings
        !
        use opentrustregion, only: logger_type, stability_settings_type, &
                                   init_stability_settings, &
                                   stability_jacobi_davidson_default, &
                                   stability_conv_tol_default, &
                                   stability_n_random_trial_vectors_default, &
                                   stability_n_iter_default, stability_verbose_default

        type(stability_settings_type) :: settings
        procedure(logger_type), pointer :: logger_funptr

        ! assume tests pass
        test_init_stability_settings = .true.

        ! call routine without optional arguments
        call init_stability_settings(settings)

        ! check if default values are correctly set
        if ((settings%jacobi_davidson .neqv. stability_jacobi_davidson_default) .or. &
            abs(settings%conv_tol - stability_conv_tol_default) > tol .or. &
            settings%n_random_trial_vectors /= &
            stability_n_random_trial_vectors_default .or. settings%n_iter /= &
            stability_n_iter_default .or. settings%verbose /= &
            stability_verbose_default .or. associated(settings%logger)) then
            write (stderr, *) "test_init_stability_settings failed: Default "// &
                "arguments not set correctly."
            test_init_stability_settings = .false.
        end if

        ! set pointer for logger routine
        logger_funptr => logger

        ! call routine with optional arguments
        call init_stability_settings(settings, jacobi_davidson=.false., &
                                     conv_tol=1.d-3, n_random_trial_vectors=3, &
                                     n_iter=50, verbose=3, logger=logger_funptr)

        ! check if optional values are correctly set
        if ((settings%jacobi_davidson .neqv. .false.) .or. &
            abs(settings%conv_tol - 1.d-3) > tol .or. settings%n_random_trial_vectors &
            /= 3 .or. settings%n_iter /= 50 .or. settings%verbose /= 3 .or. &
            .not. associated(settings%logger)) then
            write (stderr, *) "test_init_stability_settings failed: Optional "// &
                "arguments not set correctly."
            test_init_stability_settings = .false.
        end if
        if (associated(settings%logger)) then
            call settings%logger("test")
            if (log_message /= "test") write (stderr, *) &
                "test_init_stability_settings failed: Optional logging subroutine "// &
                "not set correctly."
        end if

    end function test_init_stability_settings

    logical(c_bool) function test_set_default() bind(C)
        !
        ! this function tests the subroutine which sets default values
        !
        use opentrustregion, only: set_default

        ! assume tests pass
        test_set_default = .true.

        ! check if optional arguments are correctly set for reals
        if (abs(set_default(2.d0, 1.d0) - 2.d0) > tol) then
            write (stderr, *) "test_set_default failed: Optional real argument not "// &
                "set correctly."
            test_set_default = .false.
        end if

        ! check if default arguments are correctly set for reals
        if (abs(set_default(default_value=1.d0) - 1.d0) > tol) then
            write (stderr, *) "test_set_default failed: Default real argument not "// &
                "set correctly."
            test_set_default = .false.
        end if

        ! check if optional arguments are correctly set for logicals
        if (set_default(.true., .false.) .neqv. .true.) then
            write (stderr, *) "test_set_default failed: Optional logical argument "// &
                "not set correctly."
            test_set_default = .false.
        end if

        ! check if default arguments are correctly set for logicals
        if (set_default(default_value=.false.) .neqv. .false.) then
            write (stderr, *) "test_set_default failed: Default logical argument "// &
                "not set correctly."
            test_set_default = .false.
        end if

        ! check if optional arguments are correctly set for integers
        if (set_default(2, 1) /= 2) then
            write (stderr, *) "test_set_default failed: Optional integer argument "// &
                "not set correctly."
            test_set_default = .false.
        end if

        ! check if default arguments are correctly set for integers
        if (set_default(default_value=1) /= 1) then
            write (stderr, *) "test_set_default failed: Default integer argument "// &
                "not set correctly."
            test_set_default = .false.
        end if

    end function test_set_default

    logical(c_bool) function test_level_shifted_diag_precond() bind(C)
        !
        ! this function tests the subroutine that constructs the level-shifted diagonal 
        ! preconditioner
        !
        use opentrustregion, only: level_shifted_diag_precond

        real(rp) :: residual(3), h_diag(3)

        ! assume tests pass
        test_level_shifted_diag_precond = .true.

        ! initialize quantities
        residual = [1.d0, 1.d0, 1.d0]
        h_diag = [-2.d0, 1.d0, 2.d0]

        ! call function and check if results match
        if (any(abs(level_shifted_diag_precond(residual, -2.d0, h_diag) - &
                    [1.d10, 1.d0 / 3, 0.25d0]) > tol)) then
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

        real(rp) :: residual(3), h_diag(3)

        ! assume tests pass
        test_abs_diag_precond = .true.

        ! initialize quantities
        residual = [1.d0, 1.d0, 1.d0]
        h_diag = [-2.d0, 0.d0, 2.d0]

        ! call function and check if results match
        if (any(abs(abs_diag_precond(residual, h_diag) - [0.5d0, 1.d10, 0.5d0]) > &
            tol)) then
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
        vector = [1.d0, 2.d0, 3.d0, 4.d0]
        direction = [0.d0, 1.d0, 2.d0, 0.d0] / sqrt(5.d0)

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
        use opentrustregion, only: hess_x_type, jacobi_davidson_correction

        procedure(hess_x_type), pointer :: hess_x_funptr
        real(rp), dimension(6) :: vars, vector, solution, corr_vector, hess_vector
        logical :: error

        ! assume tests pass
        test_jacobi_davidson_correction = .true.

        ! define point near saddle point, define trial vector, and solution to be 
        ! projected out
        vars = [0.35d0, 0.59d0, 0.48d0, 0.40d0, 0.31d0, 0.32d0]
        vector = [0.1d0, 0.2d0, 0.3d0, 0.4d0, 0.5d0, 0.6d0]
        solution = [1.d0, -2.d0, 2.d0, -1.d0, 1.d0, -2.d0]

        ! generate Hessian
        call hartmann6d_hessian(vars)

        ! define Hessian linear transformation
        hess_x_funptr => hess_x

        ! calculate Jacobi-Davidson correction and compare values
        call jacobi_davidson_correction(hess_x_funptr, vector, solution, 0.5d0, &
                                        corr_vector, hess_vector, error)
        if (error) then
            write (stderr, *) "test_jacobi_davidson_correction failed: Returned error."
            test_jacobi_davidson_correction = .false.
        end if
        if (sum(abs(corr_vector - [-96.940677944d0, 203.929698480d0, -216.199768920d0, &
                                   100.656941418d0, -90.624469448d0, 212.045768918d0]) &
                ) > 1d-8) then
            write (stderr, *) "test_jacobi_davidson_correction failed: Returned "// &
                "correction vector wrong."
            test_jacobi_davidson_correction = .false.
        end if
        if (sum(abs(hess_vector - [14.407362159d0, -18.566381727d0, 6.546311286d0, &
                                   -10.441098685d0, 20.923570656d0, -10.250311288d0]) &
                ) > 1d-8) then
            write (stderr, *) "test_jacobi_davidson_correction failed: Returned "// &
                "Hessian linear transformation wrong."
            test_jacobi_davidson_correction = .false.
        end if

    end function test_jacobi_davidson_correction

    logical(c_bool) function test_minres() bind(C)
        !
        ! this function tests the minimum residual method subroutine
        !
        use opentrustregion, only: hess_x_type, minres

        type(settings_type) :: settings
        procedure(hess_x_type), pointer :: hess_x_funptr
        real(rp), dimension(6) :: vars, rhs, solution, vector, hess_vector, corr_vector
        real(rp) :: mu
        logical :: error

        ! assume tests pass
        test_minres = .true.

        ! setup settings object
        call setup_settings(settings)

        ! define point near saddle point
        vars = [0.35d0, 0.59d0, 0.48d0, 0.40d0, 0.31d0, 0.32d0]

        ! define solution to be projected out
        solution = [1.d0, 2.d0, 3.d0, 4.d0, 5.d0, 6.d0] / sqrt(91.d0)

        ! generate Hessian
        call hartmann6d_hessian(vars)

        ! define Hessian linear transformation
        hess_x_funptr => hess_x

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
        call minres(-rhs, hess_x_funptr, solution, mu, 1d-14, vector, hess_vector, &
                    settings, error)
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
        rhs = 0.d0
        call minres(-rhs, hess_x_funptr, solution, mu, 1d-14, vector, hess_vector, &
                    settings, error)
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
        use opentrustregion, only: solver_settings_type, print_results

        type(solver_settings_type) :: settings

        ! assume tests pass
        test_print_results = .true.

        ! setup settings object
        call setup_settings(settings)

        ! print row of results table without optional arguments and check if row is 
        ! correct
        call print_results(settings, 1, 2.d0, 3.d0)
        if (log_message /= "        1   |     2.00000000000000E+00   "// &
            "|   3.00E+00   |      -      |        -   |        -     |      -   ") then
            write (stderr, *) "test_print_results failed: Printed row without "// &
                "optional arguments not correct."
            test_print_results = .false.
        end if

        ! print row of results table with optional arguments
        call print_results(settings, 1, 2.d0, 3.d0, 4.d0, 5, 6, 7.d0, 8.d0)
        if (log_message /= "        1   |     2.00000000000000E+00   "// &
            "|   3.00E+00   |   4.00E+00  |   0 |   5  |    7.00E+00  |  8.00E+00") then
            write (stderr, *) "test_print_results failed: Printed row with "// &
                "optional arguments not correct."
            test_print_results = .false.
        end if

    end function test_print_results

    logical(c_bool) function test_log() bind(C)
        !
        ! this function tests the logging subroutine
        !
        use opentrustregion, only: log

        type(settings_type) :: settings

        ! assume tests pass
        test_log = .true.

        ! setup settings object
        call setup_settings(settings)

        ! check if logging is correctly performed according to verbosity level when 
        ! logger is provided 
        call log(settings, "This is a test message.", 1)
        if (trim(log_message) /= " This is a test message.") then
            write (stderr, *) "test_log failed: Log message is not printed "// &
                "correctly even though it should be according to verbosity level."
            test_log = .false.
        end if
        call log(settings, "This is another test message.", 4)
        if (log_message == " This is another test message.") then
            write (stderr, *) "test_log failed: Log message is printed even "// &
                "though it should not be according to verbosity level."
            test_log = .false.
        end if

    end function test_log

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
        use opentrustregion, only: solver_settings_type, &
                                   accept_trust_region_step, &
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

        solution = [0.3d0, 0.3d0, 0.3d0]

        ! check if step is rejected and trust radius is correctly reduced if micro 
        ! iterations have not converged
        trust_radius = 1.d0
        accept_step = accept_trust_region_step(solution, 1.d0, .false., settings, &
                                               trust_radius, max_precision_reached)
        if (accept_step .or. abs(trust_radius - trust_radius_shrink_factor) > tol) then
            write(stderr, *) "test_accept_trust_region_step failed: Step accepted "// &
                "or trust radius not correctly reduced when micro iterations have "// &
                "not converged."
            test_accept_trust_region_step = .false.
        end if

        ! check if step is rejected and trust radius is correctly reduced if ratio is 
        ! negative
        trust_radius = 1.d0
        accept_step = accept_trust_region_step(solution, -1.d0, .true., settings, &
                                               trust_radius, max_precision_reached)
        if (accept_step .or. abs(trust_radius - trust_radius_shrink_factor) > tol) then
            write(stderr, *) "test_accept_trust_region_step failed: Step accepted "// & 
                "or trust radius not correctly reduced when ratio is negative."
            test_accept_trust_region_step = .false.
        end if

        ! check if step is rejected and trust radius is correctly reduced if 
        ! individual rotations are too large
        trust_radius = 1.d0
        solution(1) = 1.d0
        accept_step = accept_trust_region_step(solution, 1.d0, .true., settings, &
                                               trust_radius, max_precision_reached)
        if (accept_step .or. abs(trust_radius - trust_radius_shrink_factor) > tol) then
            write(stderr, *) "test_accept_trust_region_step failed: Step accepted "// &
                "or trust radius not correctly reduced when individual rotations "// &
                "are too large."
            test_accept_trust_region_step = .false.
        end if
        solution(1) = 0.3d0

        ! check if step is accepted and trust radius is correctly reduced if ratio is 
        ! too small
        trust_radius = 1.d0
        accept_step = accept_trust_region_step(solution, &
                                               0.9d0 * trust_radius_shrink_ratio, &
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
        trust_radius = 1.d0
        accept_step = accept_trust_region_step(solution, &
                                               0.5d0 * (trust_radius_shrink_ratio + &
                                               trust_radius_expand_ratio), .true., &
                                               settings, trust_radius, &
                                               max_precision_reached)
        if (.not. accept_step .or. abs(trust_radius - 1.d0) > tol) then
            write(stderr, *) "test_accept_trust_region_step failed: Step not "// &
                "accepted or trust radius changed when ratio is acceptable."
            test_accept_trust_region_step = .false.
        end if

        ! check if step is accepted and trust radius is correctly expanded if ratio is 
        ! too large
        trust_radius = 1.d0
        accept_step = accept_trust_region_step(solution, &
                                               1.1d0 * trust_radius_expand_ratio, &
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

    logical(c_bool) function test_sanity_check() bind(C)
        !
        ! this function tests the subroutine which performs a sanity check
        !
        use opentrustregion, only: solver_settings_type, sanity_check

        type(solver_settings_type) :: settings
        real(rp) :: grad(3)
        logical :: error

        ! assume tests pass
        test_sanity_check = .true.

        ! setup settings object
        call setup_settings(settings)

        ! check if error is incorrectly thrown for finite and non-negative number of 
        ! parameters
        settings%davidson = .true.
        settings%n_random_trial_vectors = 0
        call sanity_check(settings, 3, grad, error)
        if (error) then
            write(stderr, *) "test_sanity_check failed: Error thrown for "// &
                "non-negative and non-vanishing number of parameters."
            test_sanity_check = .false.
        end if

        ! check if error is correctly thrown for vanishing number of parameters
        call sanity_check(settings, 0, grad, error)
        if (.not. error) then
            write(stderr, *) "test_sanity_check failed: Error not thrown for "// &
                "vanishing number of parameters."
            test_sanity_check = .false.
        end if

        ! check if error is correctly thrown for negative number of parameters
        call sanity_check(settings, -1, grad, error)
        if (.not. error) then
            write(stderr, *) "test_sanity_check failed: Error not thrown for "// &
                "negative number of parameters."
            test_sanity_check = .false.
        end if

        ! check if number of random trial vectors is reduced correctly
        settings%n_random_trial_vectors = 3
        call sanity_check(settings, 3, grad, error)
        if (settings%n_random_trial_vectors /= 1) then
            write(stderr, *) "test_sanity_check failed: Number of random trial not "// &
                "correctly set."
            test_sanity_check = .false.
        end if

        ! check if gradient size is treated correctly
        call sanity_check(settings, 3, grad, error)
        if (error) then
            write(stderr, *) "test_sanity_check failed: Error thrown for correct "// &
                "gradient size."
            test_sanity_check = .false.
        end if
        call sanity_check(settings, 4, grad, error)
        if (.not. error) then
            write(stderr, *) "test_sanity_check failed: Error not thrown for "// &
                "incorrect gradient size."
            test_sanity_check = .false.
        end if

    end function test_sanity_check

    logical(c_bool) function test_level_shifted_davidson() bind(C)
        !
        ! this function tests the level-shifted Davidson subroutine
        !
        use opentrustregion, only: obj_func_type, hess_x_type, precond_type, &
                                   solver_settings_type, level_shifted_davidson, &
                                   trust_radius_shrink_ratio, &
                                   trust_radius_expand_ratio, &
                                   trust_radius_shrink_factor, &
                                   trust_radius_expand_factor

        integer(ip), parameter :: n_param = 6
        real(rp) :: func, grad_norm, trust_radius, mu, ratio, solution_norm
        real(rp), dimension(n_param) :: grad, h_diag, solution
        integer(ip) :: i, imicro, imicro_jacobi_davidson
        procedure(obj_func_type), pointer :: obj_func_funptr
        procedure(hess_x_type), pointer :: hess_x_funptr
        procedure(precond_type), pointer :: precond_funptr
        type(solver_settings_type) :: settings
        logical :: jacobi_davidson_started, max_precision_reached, error

        ! assume tests pass
        test_level_shifted_davidson = .true.

        ! setup settings object
        call setup_settings(settings)
        settings%n_random_trial_vectors = 1
        settings%n_micro = 50
        settings%global_red_factor = 1d-3
        settings%local_red_factor = 1d-4
        settings%jacobi_davidson = .false.
        settings%prefer_jacobi_davidson = .false.

        ! initialize variables
        trust_radius = 0.4d0
        obj_func_funptr => obj_func
        hess_x_funptr => hess_x

        ! start in quadratic region near minimum
        curr_vars = [0.20d0, 0.15d0, 0.48d0, 0.28d0, 0.31d0, 0.66d0]
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
                                    .false., precond_funptr, trust_radius, solution, &
                                    mu, imicro, imicro_jacobi_davidson, &
                                    jacobi_davidson_started, max_precision_reached, &
                                    error)
        if (error) then
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
                dot_product(solution, grad + 0.5d0*hartmann6d_hess_x(solution))
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
        curr_vars = [0.35d0, 0.59d0, 0.48d0, 0.40d0, 0.31d0, 0.32d0]
        func = hartmann6d_func(curr_vars)
        call hartmann6d_gradient(curr_vars, grad)
        grad_norm = norm2(grad)
        call hartmann6d_hessian(curr_vars)
        h_diag = [(hess(i, i), i=1, size(h_diag))]
        trust_radius = 0.4d0

        ! run level-shifted Davidson, check if error has occured, whether the level 
        ! shift is negative and whether the solution lies at the trust region boundary 
        ! and describes a level-shifted Newton step
        call level_shifted_davidson(func, grad, grad_norm, h_diag, n_param, &
                                    obj_func_funptr, hess_x_funptr, settings, &
                                    .false., precond_funptr, trust_radius, solution, &
                                    mu, imicro, imicro_jacobi_davidson, &
                                    jacobi_davidson_started, max_precision_reached, &
                                    error)
        if (error) then
            write (stderr, *) "test_level_shifted_davidson failed: Produced error "// &
                "near saddle point."
            test_level_shifted_davidson = .false.
        end if
        if (mu >= 0.d0) then
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
                dot_product(solution, grad + 0.5d0*hartmann6d_hess_x(solution))
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
        settings%jacobi_davidson = .true.
        trust_radius = 0.4d0

        ! run level-shifted Jacobi-Davidson, check if error has occured, whether the 
        ! level shift is negative and whether the solution lies at the trust region 
        ! boundary and describes a level-shifted Newton step
        call level_shifted_davidson(func, grad, grad_norm, h_diag, n_param, &
                                    obj_func_funptr, hess_x_funptr, settings, &
                                    .false., precond_funptr, trust_radius, solution, &
                                    mu, imicro, imicro_jacobi_davidson, &
                                    jacobi_davidson_started, max_precision_reached, &
                                    error)
        if (error) then
            write (stderr, *) "test_level_shifted_davidson failed: Produced error "// &
                "near saddle point with Jacobi-Davidson solver."
            test_level_shifted_davidson = .false.
        end if
        if (mu >= 0.d0) then
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
                dot_product(solution, grad + 0.5d0*hartmann6d_hess_x(solution))
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
        use opentrustregion, only: obj_func_type, hess_x_type, precond_type, &
                                   solver_settings_type, truncated_conjugate_gradient, &
                                   trust_radius_shrink_ratio, &
                                   trust_radius_expand_ratio, &
                                   trust_radius_shrink_factor, &
                                   trust_radius_expand_factor

        integer(ip), parameter :: n_param = 6
        real(rp) :: func, trust_radius, ratio, solution_norm
        real(rp), dimension(n_param) :: grad, h_diag, solution
        integer(ip) :: i, imicro
        procedure(obj_func_type), pointer :: obj_func_funptr
        procedure(hess_x_type), pointer :: hess_x_funptr
        procedure(precond_type), pointer :: precond_funptr
        type(solver_settings_type) :: settings
        logical :: max_precision_reached, error

        ! assume tests pass
        test_truncated_conjugate_gradient = .true.

        ! setup settings object
        call setup_settings(settings)
        settings%n_micro = 50

        ! initialize variables
        trust_radius = 0.4d0
        obj_func_funptr => obj_func
        hess_x_funptr => hess_x

        ! start in quadratic region near minimum
        curr_vars = [0.20d0, 0.15d0, 0.48d0, 0.28d0, 0.31d0, 0.66d0]
        func = hartmann6d_func(curr_vars)
        call hartmann6d_gradient(curr_vars, grad)
        call hartmann6d_hessian(curr_vars)
        h_diag = [(hess(i, i), i=1, size(h_diag))]

        ! run truncated conjugate gradient, check whether the solution lies at the 
        ! trust region boundary and reduces the function value
        call truncated_conjugate_gradient(func, grad, h_diag, n_param, &
                                          obj_func_funptr, hess_x_funptr, .false., &
                                          precond_funptr, settings, trust_radius, &
                                          solution, imicro, max_precision_reached, &
                                          error)
        if (error) then
            write (stderr, *) "test_truncated_jacobi_davidson failed: Produced "// &
                "error near minimum."
            test_truncated_conjugate_gradient = .false.
        end if
        ratio = (hartmann6d_func(curr_vars + solution) - func) / &
                dot_product(solution, grad + 0.5d0*hartmann6d_hess_x(solution))
        if (ratio <= 0.d0) then
            write (stderr, *) "test_truncated_conjugate_gradient failed: Solution "// &
                "does not reduce function value near minimum."
            test_truncated_conjugate_gradient = .false.
        end if
        solution_norm = dot_product(solution, solution / max(abs(h_diag), 1.d-10))
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
        curr_vars = [0.35d0, 0.59d0, 0.48d0, 0.40d0, 0.31d0, 0.32d0]
        func = hartmann6d_func(curr_vars)
        call hartmann6d_gradient(curr_vars, grad)
        call hartmann6d_hessian(curr_vars)
        h_diag = [(hess(i, i), i=1, size(h_diag))]
        trust_radius = 0.4d0

        ! run truncated conjugate gradient, check whether the solution lies at the 
        ! trust region boundary and reduces the function value
        call truncated_conjugate_gradient(func, grad, h_diag, n_param, &
                                          obj_func_funptr, hess_x_funptr, .false., &
                                          precond_funptr, settings, trust_radius, &
                                          solution, imicro, max_precision_reached, &
                                          error)
        if (error) then
            write (stderr, *) "test_truncated_jacobi_davidson failed: Produced "// &
                "error near saddle point."
            test_truncated_conjugate_gradient = .false.
        end if
        ratio = (hartmann6d_func(curr_vars + solution) - func) / &
                dot_product(solution, grad + 0.5d0*hartmann6d_hess_x(solution))
        if (ratio <= 0.d0) then
            write (stderr, *) "test_truncated_conjugate_gradient failed: Solution "// &
                "does not reduce function value near saddle point."
            test_truncated_conjugate_gradient = .false.
        end if
        solution_norm = dot_product(solution, solution / max(abs(h_diag), 1.d-10))
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

end module opentrustregion_unit_tests
