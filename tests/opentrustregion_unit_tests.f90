! Copyright (C) 2025- Jonas Greiner
!
! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.

module opentrustregion_unit_tests

    use opentrustregion, only: rp, ip, stderr
    use c_interface, only: c_rp, c_ip
    use test_reference, only: tol, host_context_type
    use, intrinsic :: iso_c_binding, only: c_bool

    implicit none

    ! parameters for 6D Hartmann function
    integer(ip), parameter :: n_param = 6, n_terms = 4
    real(rp), parameter :: alpha(n_terms) = [1.0_rp, 1.2_rp, 3.0_rp, 3.2_rp]
    real(rp), parameter :: A(n_terms, n_param) = &
        reshape([10.0_rp, 0.05_rp, 3.0_rp, 17.0_rp, &
                 3.0_rp, 10.0_rp, 3.5_rp, 8.0_rp, &
                 17.0_rp, 17.0_rp, 1.7_rp, 0.05_rp, &
                 3.5_rp, 0.1_rp, 10.0_rp, 10.0_rp, &
                 1.7_rp, 8.0_rp, 17.0_rp, 0.1_rp, &
                 8.0_rp, 14.0_rp, 8.0_rp, 14.0_rp], [n_terms, n_param])
    real(rp), parameter :: P(n_terms, n_param) = &
        reshape([0.1312_rp, 0.2329_rp, 0.2348_rp, 0.4047_rp, &
                 0.1696_rp, 0.4135_rp, 0.1451_rp, 0.8828_rp, &
                 0.5569_rp, 0.8307_rp, 0.3522_rp, 0.8732_rp, &
                 0.0124_rp, 0.3736_rp, 0.2883_rp, 0.5743_rp, &
                 0.8283_rp, 0.1004_rp, 0.3047_rp, 0.1091_rp, &
                 0.5886_rp, 0.9991_rp, 0.6650_rp, 0.0381_rp], [n_terms, n_param])
    integer(c_ip), bind(C, name="hartmann6d_n_param") :: n_param_c = n_param
    integer(c_ip), bind(C, name="hartmann6d_n_terms") :: n_terms_c = n_terms
    real(c_rp), bind(C, name="hartmann6d_alpha") :: alpha_c(n_terms) = alpha
    real(c_rp), bind(C, name="hartmann6d_A") :: A_c(n_terms, n_param) = A
    real(c_rp), bind(C, name="hartmann6d_P") :: P_c(n_terms, n_param) = P

    ! stationary points of 6D Hartmann function
    real(rp), parameter :: minimum1(n_param) = &
        [0.20168951_rp, 0.15001069_rp, 0.47687398_rp, 0.27533243_rp, 0.31165162_rp, &
         0.65730053_rp]
    real(rp), parameter :: minimum2(n_param) = &
        [0.40465313_rp, 0.88244493_rp, 0.84610160_rp, 0.57398969_rp, 0.13892673_rp, &
         0.03849589_rp]
    real(rp), parameter :: saddle_point(n_param) = &
        [0.35278250_rp, 0.59374767_rp, 0.47631257_rp, 0.40058250_rp, 0.31111531_rp, &
         0.32397158_rp]

    ! points of 6D Hartmann function in the quadratic regions near the first minimum
    ! and near the saddle point and a point far from both
    real(rp), parameter :: near_minimum(n_param) = &
        [0.20_rp, 0.15_rp, 0.48_rp, 0.28_rp, 0.31_rp, 0.66_rp]
    real(rp), parameter :: near_saddle_point(n_param) = &
        [0.35_rp, 0.59_rp, 0.48_rp, 0.40_rp, 0.31_rp, 0.32_rp]
    real(rp), parameter :: distant_point(n_param) = &
        [0.9_rp, 0.1_rp, 0.7_rp, 0.2_rp, 0.6_rp, 0.4_rp]

    real(c_rp), bind(C, name="hartmann6d_minimum1") :: minimum1_c(n_param) = minimum1
    real(c_rp), bind(C, name="hartmann6d_saddle_point") :: saddle_point_c(n_param) = &
        saddle_point
    real(c_rp), bind(C, name="hartmann6d_near_minimum") :: near_minimum_c(n_param) = &
        near_minimum

    ! define type for the host context handed to the mock callback functions, which
    ! holds the messages passed to the mock logging function
    type, extends(host_context_type) :: test_context_type
        character(len=:), allocatable :: log_message
    end type

    ! define type for the host context handed to the mock callback functions of the
    ! Hartmann 6D function, which holds its state and counts the callback calls
    type, extends(test_context_type) :: hartmann6d_context_type
        ! current variables and Hessian
        real(rp) :: vars(n_param) = 0.0_rp, hess(n_param, n_param) = 0.0_rp

        ! number of calls to the orbital update and the Hessian linear transformation,
        ! so that the counters reported by the library can be checked, and number of
        ! orbital updates at the previous convergence check
        integer(ip) :: n_update_orbs_calls = 0, n_hess_x_calls = 0, &
                       n_update_orbs_at_conv_check = -1
    end type

    ! define type for the host context of the Hartmann 6D function which additionally
    ! records whether the objective function was evaluated without displacement,
    ! whether the Hessian was applied to the normalized gradient, the last
    ! displacement passed to the orbital update and whether the preconditioner
    ! received this displacement
    type, extends(hartmann6d_context_type) :: hartmann6d_recording_context_type
        logical :: obj_func_zero_step = .false., hess_x_normalized_grad = .false., &
                   precond_received_step = .false.
        real(rp) :: last_delta_vars(n_param) = 0.0_rp
    end type

    ! callback functions of the Hartmann 6D function into which a fault can be
    ! injected, with their names for failure messages
    integer(ip), parameter :: fault_update_orbs = 1, fault_obj_func = 2, &
                              fault_hess_x = 3, fault_precond = 4, fault_project = 5, &
                              fault_conv_check = 6, n_fault_callbacks = 6
    character(len=*), parameter :: fault_callback_names(n_fault_callbacks) = &
        [character(len=29) :: "orbital update", "objective function", "Hessian "// &
         "linear transformation", "preconditioner", "projection", "convergence check"]

    ! define type for the host context of the Hartmann 6D function which counts the
    ! calls of every callback function and lets the call fault_call of the callback
    ! function fault_callback fail
    type, extends(hartmann6d_context_type) :: hartmann6d_fault_context_type
        integer(ip) :: fault_callback = 0, fault_call = 0, &
                       n_callback_calls(n_fault_callbacks) = 0
    end type

    ! define type for the host context handed to the mock callback functions of a
    ! quadratic model of any dimension, which holds its gradient and Hessian and the
    ! index of the next unit vector a preconditioner returns
    type, extends(test_context_type) :: quadratic_context_type
        real(rp), allocatable :: grad(:), hess(:, :)
        integer(ip) :: next_unit_vector = 0
    end type

    ! define type for the host context handed to the mock callback functions of the
    ! double well function x^2 - y^2 + y^4, which has a saddle point at the origin and
    ! minima at y = +-1/sqrt(2), and holds its current variables
    type, extends(test_context_type) :: double_well_context_type
        real(rp) :: vars(2) = 0.0_rp
    end type

    ! define type for the host context handed to the objective function used to test
    ! the bracketing, which selects one of the one-dimensional test functions
    type, extends(test_context_type) :: bracket_context_type
        integer(ip) :: func_case = 0
    end type

    ! error code returned by mock callback functions reached without their context,
    ! chosen so that it cannot be mistaken for an expected error code
    integer(ip), parameter :: missing_context_error = 42

contains

    ! 6D Hartmann function definition

    function hartmann6d_func(vars) result(f)
        !
        ! this function defines the Hartmann 6D function
        !
        real(rp), intent(in) :: vars(:)
        real(rp) :: f, exp_term(n_terms)
        integer(ip) :: i

        do i = 1, n_terms
            exp_term(i) = exp(-sum(A(i, :) * (vars - P(i, :))**2))
        end do

        f = -sum(alpha * exp_term)

    end function hartmann6d_func

    subroutine hartmann6d_gradient(vars, grad)
        !
        ! this subroutine defines the Hartmann 6D function's gradient
        !
        real(rp), intent(in) :: vars(:)
        real(rp), intent(out) :: grad(:)
        real(rp) :: exp_term(n_terms)
        integer(ip) :: i, j

        do i = 1, n_terms
            exp_term(i) = exp(-sum(A(i, :) * (vars - P(i, :))**2))
        end do

        do j = 1, n_param
            grad(j) = sum(2.0_rp * alpha * A(:, j) * (vars(j) - P(:, j)) * exp_term)
        end do

    end subroutine hartmann6d_gradient

    function hartmann6d_hessian(vars) result(hess)
        !
        ! this function defines the Hartmann 6D function's Hessian
        !
        real(rp), intent(in) :: vars(:)
        real(rp) :: hess(n_param, n_param)

        real(rp) :: exp_term(n_terms)
        integer(ip) :: i, j

        do i = 1, n_terms
            exp_term(i) = exp(-sum(A(i, :) * (vars - P(i, :))**2))
        end do

        do i = 1, n_param
            hess(i, i) = 2.0_rp * &
                         sum(alpha * A(:, i) * exp_term * &
                             (1.0_rp - 2.0_rp * A(:, i) * (vars(i) - P(:, i))**2))
            do j = 1, i - 1
                hess(i, j) = -4.0_rp * &
                             sum(alpha * A(:, i) * A(:, j) * (vars(i) - P(:, i)) * &
                                 (vars(j) - P(:, j)) * exp_term)
                hess(j, i) = hess(i, j)
            end do
        end do

    end function hartmann6d_hessian

    function resolve_hartmann6d_context(context, error) result(state)
        !
        ! this function returns the Hartmann 6D context handed to a mock callback
        ! function and an error if no such context was handed over
        !
        class(*), intent(in), pointer :: context
        integer(ip), intent(out) :: error
        class(hartmann6d_context_type), pointer :: state

        ! initialize error flag
        error = 0

        ! get context
        state => null()
        if (associated(context)) then
            select type (context)
            class is (hartmann6d_context_type)
                state => context
            end select
        end if

        ! report missing test context
        if (.not. associated(state)) error = missing_context_error

    end function resolve_hartmann6d_context

    subroutine hess_x_fun(x, hess_x, error, context)
        !
        ! this subroutine describes the Hessian linear transformation operation for the
        ! Hartmann 6D function
        !
        use test_reference, only: check_host_context

        real(rp), intent(in), target :: x(:)
        real(rp), intent(out), target :: hess_x(:)
        integer(ip), intent(out) :: error
        class(*), intent(in), pointer :: context

        class(hartmann6d_context_type), pointer :: state
        real(rp) :: grad(n_param)

        ! check host context
        call check_host_context(context)

        ! get Hartmann 6D state
        state => resolve_hartmann6d_context(context, error)
        if (error /= 0) return

        ! count call
        state%n_hess_x_calls = state%n_hess_x_calls + 1

        ! record whether the Hessian is applied to the normalized gradient at the
        ! current variables
        select type (state)
        class is (hartmann6d_recording_context_type)
            call hartmann6d_gradient(state%vars, grad)
            if (norm2(grad) > 0.0_rp) then
                if (norm2(x - grad / norm2(grad)) < tol) &
                    state%hess_x_normalized_grad = .true.
            end if
        end select

        hess_x = matmul(state%hess, x)

    end subroutine hess_x_fun

    subroutine hess_x_fun_asymmetric(x, hess_x, error, context)
        !
        ! this subroutine describes the Hessian linear transformation operation for the
        ! Hartmann 6D function with a small antisymmetric contribution, which mimics
        ! numerical noise large enough that the Jacobi-Davidson method has to
        ! recalculate linear transformations which no longer respect Hessian symmetry
        !
        use test_reference, only: check_host_context

        real(rp), intent(in), target :: x(:)
        real(rp), intent(out), target :: hess_x(:)
        integer(ip), intent(out) :: error
        class(*), intent(in), pointer :: context

        class(hartmann6d_context_type), pointer :: state
        real(rp), parameter :: asymmetry = 1e-8_rp

        ! check host context
        call check_host_context(context)

        ! get Hartmann 6D state
        state => resolve_hartmann6d_context(context, error)
        if (error /= 0) return

        ! count call
        state%n_hess_x_calls = state%n_hess_x_calls + 1

        hess_x = matmul(state%hess, x) + &
                 asymmetry * [x(2), -x(1), x(4), -x(3), x(6), -x(5)]

    end subroutine hess_x_fun_asymmetric

    subroutine hess_x_fun_failing(x, hess_x, error, context)
        !
        ! this subroutine describes a Hessian linear transformation which always fails
        !
        use test_reference, only: check_host_context

        real(rp), intent(in), target :: x(:)
        real(rp), intent(out), target :: hess_x(:)
        integer(ip), intent(out) :: error
        class(*), intent(in), pointer :: context

        class(hartmann6d_context_type), pointer :: state

        ! check host context
        call check_host_context(context)

        ! get Hartmann 6D state
        state => resolve_hartmann6d_context(context, error)
        if (error /= 0) return

        ! count call
        state%n_hess_x_calls = state%n_hess_x_calls + 1

        error = 1
        hess_x = x

    end subroutine hess_x_fun_failing

    function obj_func(delta_vars, error, context) result(func)
        !
        ! this function describes the objective function evaluation for the Hartmann 6D
        ! function
        !
        use test_reference, only: check_host_context

        real(rp), intent(in), target :: delta_vars(:)
        integer(ip), intent(out) :: error
        class(*), intent(in), pointer :: context

        class(hartmann6d_context_type), pointer :: state
        real(rp) :: func

        ! check host context
        call check_host_context(context)

        ! initialize result in case the test problem state is missing
        func = 0.0_rp

        ! get Hartmann 6D state
        state => resolve_hartmann6d_context(context, error)
        if (error /= 0) return

        ! record whether function is evaluated without displacement
        select type (state)
        class is (hartmann6d_recording_context_type)
            if (maxval(abs(delta_vars)) <= 0.0_rp) state%obj_func_zero_step = .true.
        end select

        func = hartmann6d_func(state%vars + delta_vars)

    end function obj_func

    function obj_func_failing(delta_vars, error, context) result(func)
        !
        ! this function describes an objective function evaluation which always fails
        !
        use test_reference, only: check_host_context

        real(rp), intent(in), target :: delta_vars(:)
        integer(ip), intent(out) :: error
        class(*), intent(in), pointer :: context
        real(rp) :: func

        ! check host context
        call check_host_context(context)

        func = sum(delta_vars)
        error = 1

    end function obj_func_failing

    subroutine update_orbs(delta_vars, func, grad, h_diag, hess_x_funptr, error, &
                           context)
        !
        ! this subroutine describes the orbital update equivalent for the Hartmann 6D
        ! function
        !
        use opentrustregion, only: hess_x_type
        use test_reference, only: check_host_context

        real(rp), intent(in), target :: delta_vars(:)
        real(rp), intent(out) :: func
        real(rp), intent(out), target :: grad(:), h_diag(:)
        procedure(hess_x_type), intent(inout), pointer :: hess_x_funptr
        integer(ip), intent(out) :: error
        class(*), intent(in), pointer :: context

        class(hartmann6d_context_type), pointer :: state
        integer(ip) :: i

        ! check host context
        call check_host_context(context)

        ! get Hartmann 6D state
        state => resolve_hartmann6d_context(context, error)
        if (error /= 0) return

        ! count call
        state%n_update_orbs_calls = state%n_update_orbs_calls + 1

        ! record and update variables
        select type (state)
        class is (hartmann6d_recording_context_type)
            state%last_delta_vars = delta_vars
        end select
        state%vars = state%vars + delta_vars

        ! evaluate function, calculate gradient and Hessian diagonal and define Hessian
        ! linear transformation
        func = hartmann6d_func(state%vars)
        call hartmann6d_gradient(state%vars, grad)
        state%hess = hartmann6d_hessian(state%vars)
        h_diag = [(state%hess(i, i), i=1, size(h_diag))]
        hess_x_funptr => hess_x_fun

    end subroutine update_orbs

    subroutine update_orbs_no_hess_x(delta_vars, func, grad, h_diag, hess_x_funptr, &
                                     error, context)
        !
        ! this subroutine describes an orbital update which reports success but does
        ! not provide a Hessian linear transformation
        !
        use opentrustregion, only: hess_x_type

        real(rp), intent(in), target :: delta_vars(:)
        real(rp), intent(out) :: func
        real(rp), intent(out), target :: grad(:), h_diag(:)
        procedure(hess_x_type), intent(inout), pointer :: hess_x_funptr
        integer(ip), intent(out) :: error
        class(*), intent(in), pointer :: context

        call update_orbs(delta_vars, func, grad, h_diag, hess_x_funptr, error, context)
        hess_x_funptr => null()

    end subroutine update_orbs_no_hess_x

    logical function inject_fault(context, callback)
        !
        ! this function counts a call of a callback function in a fault injection
        ! context and returns whether this call has to fail
        !
        class(*), intent(in), pointer :: context
        integer(ip), intent(in) :: callback

        ! count call and determine whether it fails
        inject_fault = .false.
        if (.not. associated(context)) return
        select type (context)
        class is (hartmann6d_fault_context_type)
            context%n_callback_calls(callback) = context%n_callback_calls(callback) + 1
            inject_fault = callback == context%fault_callback .and. &
                           context%n_callback_calls(callback) == context%fault_call
        end select

    end function inject_fault

    subroutine faulty_update_orbs(delta_vars, func, grad, h_diag, hess_x_funptr, &
                                  error, context)
        !
        ! this subroutine describes the orbital update of the Hartmann 6D function into
        ! which a fault can be injected, and which returns the Hessian linear
        ! transformation into which a fault can be injected
        !
        use opentrustregion, only: hess_x_type

        real(rp), intent(in), target :: delta_vars(:)
        real(rp), intent(out) :: func
        real(rp), intent(out), target :: grad(:), h_diag(:)
        procedure(hess_x_type), intent(inout), pointer :: hess_x_funptr
        integer(ip), intent(out) :: error
        class(*), intent(in), pointer :: context

        ! perform orbital update
        call update_orbs(delta_vars, func, grad, h_diag, hess_x_funptr, error, context)
        if (error /= 0) return
        hess_x_funptr => faulty_hess_x

        ! inject fault
        if (inject_fault(context, fault_update_orbs)) error = 1

    end subroutine faulty_update_orbs

    function faulty_obj_func(delta_vars, error, context) result(func)
        !
        ! this function describes the objective function evaluation of the Hartmann 6D
        ! function into which a fault can be injected
        !
        real(rp), intent(in), target :: delta_vars(:)
        integer(ip), intent(out) :: error
        class(*), intent(in), pointer :: context
        real(rp) :: func

        ! evaluate objective function
        func = obj_func(delta_vars, error, context)
        if (error /= 0) return

        ! inject fault
        if (inject_fault(context, fault_obj_func)) error = 1

    end function faulty_obj_func

    subroutine faulty_hess_x(x, hess_x, error, context)
        !
        ! this subroutine describes the Hessian linear transformation of the Hartmann 6D
        ! function into which a fault can be injected
        !
        real(rp), intent(in), target :: x(:)
        real(rp), intent(out), target :: hess_x(:)
        integer(ip), intent(out) :: error
        class(*), intent(in), pointer :: context

        ! perform Hessian linear transformation
        call hess_x_fun(x, hess_x, error, context)
        if (error /= 0) return

        ! inject fault
        if (inject_fault(context, fault_hess_x)) error = 1

    end subroutine faulty_hess_x

    subroutine faulty_precond(residual, mu, precond_residual, error, context)
        !
        ! this subroutine describes an identity preconditioner into which a fault can be
        ! injected
        !
        use test_reference, only: check_host_context

        real(rp), intent(in), target :: residual(:)
        real(rp), intent(in) :: mu
        real(rp), intent(out), target :: precond_residual(:)
        integer(ip), intent(out) :: error
        class(*), intent(in), pointer :: context

        ! check host context
        call check_host_context(context)

        ! apply preconditioner and inject fault
        precond_residual = residual + 0.0_rp * mu
        error = 0
        if (inject_fault(context, fault_precond)) error = 1

    end subroutine faulty_precond

    subroutine faulty_project(vector, error, context)
        !
        ! this subroutine describes an identity projection into which a fault can be
        ! injected
        !
        use test_reference, only: check_host_context

        real(rp), intent(inout), target :: vector(:)
        integer(ip), intent(out) :: error
        class(*), intent(in), pointer :: context

        ! check host context
        call check_host_context(context)

        ! keep vector and inject fault
        vector = vector + 0.0_rp
        error = 0
        if (inject_fault(context, fault_project)) error = 1

    end subroutine faulty_project

    function faulty_conv_check(error, context) result(converged)
        !
        ! this function describes a convergence check which never reports convergence
        ! and into which a fault can be injected
        !
        use test_reference, only: check_host_context

        integer(ip), intent(out) :: error
        class(*), intent(in), pointer :: context
        logical :: converged

        ! check host context
        call check_host_context(context)

        ! report no convergence and inject fault
        converged = .false.
        error = 0
        if (inject_fault(context, fault_conv_check)) error = 1

    end function faulty_conv_check

    subroutine step_recording_precond(residual, mu, precond_residual, error, context)
        !
        ! this subroutine is an identity preconditioner which records whether it
        ! received the last nonvanishing displacement passed to the Hartmann 6D orbital
        ! update
        !
        use test_reference, only: check_host_context

        real(rp), intent(in), target :: residual(:)
        real(rp), intent(in) :: mu
        real(rp), intent(out), target :: precond_residual(:)
        integer(ip), intent(out) :: error
        class(*), intent(in), pointer :: context

        class(hartmann6d_context_type), pointer :: state

        ! check host context
        call check_host_context(context)

        ! get Hartmann 6D state
        state => resolve_hartmann6d_context(context, error)
        if (error /= 0) return

        select type (state)
        class is (hartmann6d_recording_context_type)
            if (maxval(abs(state%last_delta_vars)) > 0.0_rp .and. &
                maxval(abs(residual - state%last_delta_vars)) <= 0.0_rp) &
                state%precond_received_step = .true.
        end select
        precond_residual = residual + 0.0_rp * mu

        error = 0

    end subroutine step_recording_precond

    function mock_conv_check_without_update(error, context) result(converged)
        !
        ! this function describes a convergence check which only passes in a macro
        ! iteration without an orbital update, which is the case once the maximum
        ! precision is reached
        !
        use test_reference, only: check_host_context

        integer(ip), intent(out) :: error
        class(*), intent(in), pointer :: context

        class(hartmann6d_context_type), pointer :: state
        logical :: converged

        ! check host context
        call check_host_context(context)

        ! initialize result in case the test problem state is missing
        converged = .false.

        ! get Hartmann 6D state
        state => resolve_hartmann6d_context(context, error)
        if (error /= 0) return

        converged = state%n_update_orbs_calls == state%n_update_orbs_at_conv_check
        state%n_update_orbs_at_conv_check = state%n_update_orbs_calls

    end function mock_conv_check_without_update

    function resolve_quadratic_context(context, error) result(state)
        !
        ! this function returns the quadratic model context handed to a mock callback
        ! function and an error if no such context was handed over
        !
        class(*), intent(in), pointer :: context
        integer(ip), intent(out) :: error
        class(quadratic_context_type), pointer :: state

        ! initialize error flag
        error = 0

        ! get context
        state => null()
        if (associated(context)) then
            select type (context)
            class is (quadratic_context_type)
                state => context
            end select
        end if

        ! report missing test context
        if (.not. associated(state)) error = missing_context_error

    end function resolve_quadratic_context

    subroutine quadratic_hess_x(x, hess_x, error, context)
        !
        ! this subroutine describes the Hessian linear transformation operation for the
        ! quadratic model
        !
        use test_reference, only: check_host_context

        real(rp), intent(in), target :: x(:)
        real(rp), intent(out), target :: hess_x(:)
        integer(ip), intent(out) :: error
        class(*), intent(in), pointer :: context

        class(quadratic_context_type), pointer :: state

        ! check host context
        call check_host_context(context)

        ! get quadratic model state
        state => resolve_quadratic_context(context, error)
        if (error /= 0) return

        hess_x = matmul(state%hess, x)

    end subroutine quadratic_hess_x

    function quadratic_obj_func(delta_vars, error, context) result(func)
        !
        ! this function describes the objective function evaluation for the quadratic
        ! model
        !
        use test_reference, only: check_host_context

        real(rp), intent(in), target :: delta_vars(:)
        integer(ip), intent(out) :: error
        class(*), intent(in), pointer :: context
        real(rp) :: func

        class(quadratic_context_type), pointer :: state

        ! initialize result in case the quadratic model state is missing
        func = 0.0_rp

        ! check host context
        call check_host_context(context)

        ! get quadratic model state
        state => resolve_quadratic_context(context, error)
        if (error /= 0) return

        func = dot_product(state%grad, delta_vars) + &
               0.5_rp * dot_product(delta_vars, matmul(state%hess, delta_vars))

    end function quadratic_obj_func

    subroutine identity_hess_x(x, hess_x, error, context)
        !
        ! this subroutine describes an identity Hessian linear transformation which
        ! needs no test context
        !
        use test_reference, only: check_host_context

        real(rp), intent(in), target :: x(:)
        real(rp), intent(out), target :: hess_x(:)
        integer(ip), intent(out) :: error
        class(*), intent(in), pointer :: context

        ! check host context
        call check_host_context(context)

        error = 0
        hess_x = x

    end subroutine identity_hess_x

    subroutine identity_precond(residual, mu, precond_residual, error, context)
        !
        ! this subroutine describes an identity preconditioner which needs no test
        ! context
        !
        use test_reference, only: check_host_context

        real(rp), intent(in), target :: residual(:)
        real(rp), intent(in) :: mu
        real(rp), intent(out), target :: precond_residual(:)
        integer(ip), intent(out) :: error
        class(*), intent(in), pointer :: context

        ! check host context
        call check_host_context(context)

        precond_residual = residual + 0.0_rp * mu

        error = 0

    end subroutine identity_precond

    subroutine identity_project(vector, error, context)
        !
        ! this subroutine describes an identity projection which needs no test context
        !
        use test_reference, only: check_host_context

        real(rp), intent(inout), target :: vector(:)
        integer(ip), intent(out) :: error
        class(*), intent(in), pointer :: context

        ! check host context
        call check_host_context(context)

        vector = vector + 0.0_rp

        error = 0

    end subroutine identity_project

    subroutine stationary_update_orbs(delta_vars, func, grad, h_diag, hess_x_funptr, &
                                      error, context)
        !
        ! this subroutine describes an orbital update at the minimum of a quadratic
        ! function with identity Hessian which needs no test context, since the
        ! gradient vanishes for every displacement passed to it here
        !
        use opentrustregion, only: hess_x_type
        use test_reference, only: check_host_context

        real(rp), intent(in), target :: delta_vars(:)
        real(rp), intent(out) :: func
        real(rp), intent(out), target :: grad(:), h_diag(:)
        procedure(hess_x_type), intent(inout), pointer :: hess_x_funptr
        integer(ip), intent(out) :: error
        class(*), intent(in), pointer :: context

        ! check host context
        call check_host_context(context)

        error = 0
        func = 0.5_rp * dot_product(delta_vars, delta_vars)
        grad = 0.0_rp
        h_diag = 1.0_rp
        hess_x_funptr => identity_hess_x

    end subroutine stationary_update_orbs

    function constant_obj_func(delta_vars, error, context) result(func)
        !
        ! this function describes an objective function evaluation that does not depend
        ! on the displacement
        !
        use test_reference, only: check_host_context

        real(rp), intent(in), target :: delta_vars(:)
        integer(ip), intent(out) :: error
        class(*), intent(in), pointer :: context
        real(rp) :: func

        ! check host context
        call check_host_context(context)

        error = 0
        func = 1.0_rp + 0.0_rp * sum(delta_vars)

    end function constant_obj_func

    ! double well function definition

    function double_well_func(vars) result(f)
        !
        ! this function defines the double well function
        !
        real(rp), intent(in) :: vars(:)
        real(rp) :: f

        f = vars(1)**2 - vars(2)**2 + vars(2)**4

    end function double_well_func

    function resolve_double_well_context(context, error) result(state)
        !
        ! this function returns the double well context handed to a mock callback
        ! function and an error if no such context was handed over
        !
        class(*), intent(in), pointer :: context
        integer(ip), intent(out) :: error
        class(double_well_context_type), pointer :: state

        ! initialize error flag
        error = 0

        ! get context
        state => null()
        if (associated(context)) then
            select type (context)
            class is (double_well_context_type)
                state => context
            end select
        end if

        ! report missing test context
        if (.not. associated(state)) error = missing_context_error

    end function resolve_double_well_context

    subroutine double_well_hess_x(x, hess_x, error, context)
        !
        ! this subroutine describes the Hessian linear transformation operation for the
        ! double well function
        !
        use test_reference, only: check_host_context

        real(rp), intent(in), target :: x(:)
        real(rp), intent(out), target :: hess_x(:)
        integer(ip), intent(out) :: error
        class(*), intent(in), pointer :: context

        class(double_well_context_type), pointer :: state

        ! check host context
        call check_host_context(context)

        ! get double well state
        state => resolve_double_well_context(context, error)
        if (error /= 0) return

        hess_x = [2.0_rp, -2.0_rp + 12.0_rp * state%vars(2)**2] * x

    end subroutine double_well_hess_x

    function double_well_obj_func(delta_vars, error, context) result(func)
        !
        ! this function describes the objective function evaluation for the double
        ! well function
        !
        use test_reference, only: check_host_context

        real(rp), intent(in), target :: delta_vars(:)
        integer(ip), intent(out) :: error
        class(*), intent(in), pointer :: context
        real(rp) :: func

        class(double_well_context_type), pointer :: state

        ! initialize result in case the double well state is missing
        func = 0.0_rp

        ! check host context
        call check_host_context(context)

        ! get double well state
        state => resolve_double_well_context(context, error)
        if (error /= 0) return

        func = double_well_func(state%vars + delta_vars)

    end function double_well_obj_func

    function raised_double_well_obj_func(delta_vars, error, context) result(func)
        !
        ! this function describes an objective function evaluation for the double
        ! well function that is inconsistent with the orbital update since it raises
        ! the function value, so that no displacement lowers the function value
        !
        use test_reference, only: check_host_context

        real(rp), intent(in), target :: delta_vars(:)
        integer(ip), intent(out) :: error
        class(*), intent(in), pointer :: context
        real(rp) :: func

        class(double_well_context_type), pointer :: state

        ! initialize result in case the double well state is missing
        func = 0.0_rp

        ! check host context
        call check_host_context(context)

        ! get double well state
        state => resolve_double_well_context(context, error)
        if (error /= 0) return

        func = double_well_func(state%vars + delta_vars) + 1.0_rp

    end function raised_double_well_obj_func

    subroutine double_well_update_orbs(delta_vars, func, grad, h_diag, hess_x_funptr, &
                                       error, context)
        !
        ! this subroutine describes the orbital update equivalent for the double well
        ! function
        !
        use opentrustregion, only: hess_x_type
        use test_reference, only: check_host_context

        real(rp), intent(in), target :: delta_vars(:)
        real(rp), intent(out) :: func
        real(rp), intent(out), target :: grad(:), h_diag(:)
        procedure(hess_x_type), intent(inout), pointer :: hess_x_funptr
        integer(ip), intent(out) :: error
        class(*), intent(in), pointer :: context

        class(double_well_context_type), pointer :: state

        ! check host context
        call check_host_context(context)

        ! get double well state
        state => resolve_double_well_context(context, error)
        if (error /= 0) return

        ! update variables
        state%vars = state%vars + delta_vars

        ! evaluate function, calculate gradient and Hessian diagonal and define Hessian
        ! linear transformation
        func = double_well_func(state%vars)
        grad = [2.0_rp * state%vars(1), &
                -2.0_rp * state%vars(2) + 4.0_rp * state%vars(2)**3]
        h_diag = [2.0_rp, -2.0_rp + 12.0_rp * state%vars(2)**2]
        hess_x_funptr => double_well_hess_x

    end subroutine double_well_update_orbs

    function resolve_bracket_context(context, error) result(state)
        !
        ! this function returns the bracketing context handed to a mock callback
        ! function and an error if no such context was handed over
        !
        class(*), intent(in), pointer :: context
        integer(ip), intent(out) :: error
        class(bracket_context_type), pointer :: state

        ! initialize error flag
        error = 0

        ! get context
        state => null()
        if (associated(context)) then
            select type (context)
            class is (bracket_context_type)
                state => context
            end select
        end if

        ! report missing test context
        if (.not. associated(state)) error = missing_context_error

    end function resolve_bracket_context

    function bracket_obj_func(delta_vars, error, context) result(func)
        !
        ! this function describes the objective function evaluation for the
        ! one-dimensional functions used to test the bracketing, selected by its
        ! context
        !
        use test_reference, only: check_host_context

        real(rp), intent(in), target :: delta_vars(:)
        integer(ip), intent(out) :: error
        class(*), intent(in), pointer :: context
        real(rp) :: func

        class(bracket_context_type), pointer :: state
        real(rp) :: x

        ! initialize result in case the bracketing state is missing
        func = 0.0_rp

        ! check host context
        call check_host_context(context)

        ! get bracketing state
        state => resolve_bracket_context(context, error)
        if (error /= 0) return

        x = delta_vars(1)
        select case (state%func_case)
        case (1)
            func = (x - 3.0_rp)**2
        case (2)
            func = (x - 2.0_rp)**4
        case (3)
            func = exp(-x) + 0.01_rp * x
        case (4)
            func = abs(x - 2.2_rp)
        case (5)
            func = sqrt(abs(x - 5.0_rp))
        case (6)
            func = (x - 1000.0_rp)**2
        case (7)
            func = (x - 2.0_rp)**2 + 10.0_rp * exp(-100.0_rp * (x - 2.0_rp)**2)
        case (8)
            func = -x
        case default
            func = 1.0_rp
        end select

    end function bracket_obj_func

    subroutine mock_precond(residual, mu, precond_residual, error, context)
        !
        ! this subroutine is a test subroutine for the preconditioner subroutine
        !
        use test_reference, only: check_host_context

        real(rp), intent(in), target :: residual(:)
        real(rp), intent(in) :: mu
        real(rp), intent(out), target :: precond_residual(:)
        integer(ip), intent(out) :: error
        class(*), intent(in), pointer :: context

        ! check host context
        call check_host_context(context)

        precond_residual = mu * residual

        error = 0

    end subroutine mock_precond

    subroutine mock_precond_first_unit_vector(residual, mu, precond_residual, error, &
                                              context)
        !
        ! this subroutine is a test subroutine for the preconditioner subroutine which
        ! always returns the first unit vector
        !
        use test_reference, only: check_host_context

        real(rp), intent(in), target :: residual(:)
        real(rp), intent(in) :: mu
        real(rp), intent(out), target :: precond_residual(:)
        integer(ip), intent(out) :: error
        class(*), intent(in), pointer :: context

        ! check host context
        call check_host_context(context)

        precond_residual = 0.0_rp * residual + 0.0_rp * mu
        precond_residual(1) = 1.0_rp

        error = 0

    end subroutine mock_precond_first_unit_vector

    subroutine mock_precond_next_unit_vector(residual, mu, precond_residual, error, &
                                             context)
        !
        ! this subroutine is a test subroutine for the preconditioner subroutine which
        ! returns the next unit vector of the quadratic model on every call and the
        ! last one once all have been returned
        !
        use test_reference, only: check_host_context

        real(rp), intent(in), target :: residual(:)
        real(rp), intent(in) :: mu
        real(rp), intent(out), target :: precond_residual(:)
        integer(ip), intent(out) :: error
        class(*), intent(in), pointer :: context

        class(quadratic_context_type), pointer :: state

        ! check host context
        call check_host_context(context)

        ! initialize result in case the quadratic model state is missing
        precond_residual = 0.0_rp * residual + 0.0_rp * mu

        ! get quadratic model state
        state => resolve_quadratic_context(context, error)
        if (error /= 0) return

        state%next_unit_vector = min(state%next_unit_vector + 1, size(residual))
        precond_residual(state%next_unit_vector) = 1.0_rp

    end subroutine mock_precond_next_unit_vector

    subroutine mock_precond_error(residual, mu, precond_residual, error, context)
        !
        ! this subroutine is a test subroutine for a preconditioner subroutine which
        ! produces an error
        !
        use test_reference, only: check_host_context

        real(rp), intent(in), target :: residual(:)
        real(rp), intent(in) :: mu
        real(rp), intent(out), target :: precond_residual(:)
        integer(ip), intent(out) :: error
        class(*), intent(in), pointer :: context

        ! check host context
        call check_host_context(context)

        precond_residual = mu * residual

        error = 1

    end subroutine mock_precond_error

    subroutine mock_project(vector, error, context)
        !
        ! this subroutine is a test subroutine for the projection subroutine
        !
        use test_reference, only: check_host_context

        real(rp), intent(inout), target :: vector(:)
        integer(ip), intent(out) :: error
        class(*), intent(in), pointer :: context

        ! check host context
        call check_host_context(context)

        vector = 2 * vector

        error = 0

    end subroutine mock_project

    subroutine mock_project_error(vector, error, context)
        !
        ! this subroutine is a test subroutine for a projection subroutine which
        ! produces an error
        !
        use test_reference, only: check_host_context

        real(rp), intent(inout), target :: vector(:)
        integer(ip), intent(out) :: error
        class(*), intent(in), pointer :: context

        ! check host context
        call check_host_context(context)

        vector = 2 * vector

        error = 1

    end subroutine mock_project_error

    subroutine mock_project_out_pair(vector, error, context)
        !
        ! this subroutine is a test subroutine for the projection subroutine which
        ! removes the component along the symmetric combination of the first two unit
        ! vectors
        !
        use test_reference, only: check_host_context

        real(rp), intent(inout), target :: vector(:)
        integer(ip), intent(out) :: error
        class(*), intent(in), pointer :: context

        ! check host context
        call check_host_context(context)

        vector(1:2) = vector(1:2) - 0.5_rp * sum(vector(1:2))

        error = 0

    end subroutine mock_project_out_pair

    subroutine mock_project_first_component(vector, error, context)
        !
        ! this subroutine is a test subroutine for the projection subroutine which
        ! keeps only the first component
        !
        use test_reference, only: check_host_context

        real(rp), intent(inout), target :: vector(:)
        integer(ip), intent(out) :: error
        class(*), intent(in), pointer :: context

        ! check host context
        call check_host_context(context)

        vector(2:) = 0.0_rp

        error = 0

    end subroutine mock_project_first_component

    function mock_conv_check(error, context) result(converged)
        !
        ! this function describes a convergence check which always passes
        !
        use test_reference, only: check_host_context

        integer(ip), intent(out) :: error
        class(*), intent(in), pointer :: context
        logical :: converged

        ! check host context
        call check_host_context(context)

        ! initialize error flag
        error = 0

        converged = .true.

    end function mock_conv_check

    subroutine logger(message, context)
        !
        ! this subroutine is a mock logging subroutine
        !
        use test_reference, only: check_host_context

        character(len=*), intent(in) :: message
        class(*), intent(in), pointer :: context

        ! check host context
        call check_host_context(context)

        ! append message to the log of the test context
        if (associated(context)) then
            select type (context)
            class is (test_context_type)
                if (.not. allocated(context%log_message)) context%log_message = ""
                context%log_message = context%log_message//trim(message)
            end select
        end if

    end subroutine logger

    function identity_matrix(n) result(matrix)
        !
        ! this function returns an identity matrix
        !
        integer(ip), intent(in) :: n
        real(rp) :: matrix(n, n)

        integer(ip) :: i

        matrix = 0.0_rp
        do i = 1, n
            matrix(i, i) = 1.0_rp
        end do

    end function identity_matrix

    function generate_random_symm_matrix(n) result(matrix)
        !
        ! this function generates a random symmetric matrix
        !
        integer(ip), intent(in) :: n
        real(rp) :: matrix(n, n)

        call random_number(matrix)
        matrix = matrix + transpose(matrix)

    end function generate_random_symm_matrix

    subroutine setup_settings(settings, context)
        !
        ! this subroutine sets up a settings object for tests, which hands the test
        ! context to the mock callback functions
        !
        use opentrustregion, only: settings_type

        class(settings_type), intent(inout) :: settings
        class(test_context_type), intent(inout), target :: context

        integer(ip) :: error

        call settings%init(error)
        settings%verbose = 3
        settings%logger => logger
        settings%context => context
        context%log_message = ""

    end subroutine setup_settings

    subroutine setup_error_logging(settings, context)
        !
        ! this subroutine sets up a settings object so that only error messages are
        ! passed to the mock logging function and clears the log of the test context,
        ! so that a test can check whether an error message was printed
        !
        use opentrustregion, only: settings_type, verbosity_error

        class(settings_type), intent(inout) :: settings
        class(test_context_type), intent(inout), target :: context

        settings%logger => logger
        settings%verbose = verbosity_error
        settings%context => context
        context%log_message = ""

    end subroutine setup_error_logging

    subroutine ref_symm_mat_diag(matrix, eigvals, eigvecs)
        !
        ! this subroutine reimplements the eigendecomposition of a symmetric matrix,
        ! returning the eigenvalues in ascending order
        !
        real(rp), intent(in) :: matrix(:, :)
        real(rp), intent(out) :: eigvals(:), eigvecs(:, :)

        integer(ip) :: n, lwork, info
        real(rp), allocatable :: work(:)
        external :: dsyev

        n = size(matrix, 1)
        eigvecs = matrix

        ! query optimal workspace size
        lwork = -1
        allocate(work(1))
        call dsyev("V", "U", n, eigvecs, n, eigvals, work, lwork, info)
        lwork = int(work(1), kind=ip)
        deallocate(work)
        allocate(work(lwork))

        ! perform eigendecomposition
        call dsyev("V", "U", n, eigvecs, n, eigvals, work, lwork, info)
        deallocate(work)

    end subroutine ref_symm_mat_diag

    function ref_step_trust_radius(trust_radius, ratio)
        !
        ! this function recovers the trust radius an accepted step was computed for
        ! from the trust radius and the ratio of actual to predicted function change of
        ! the step
        !
        use opentrustregion, only: &
            trust_radius_shrink_ratio, trust_radius_expand_ratio, &
            trust_radius_shrink_factor, trust_radius_expand_factor

        real(rp), intent(in) :: trust_radius, ratio
        real(rp) :: ref_step_trust_radius

        if (ratio < trust_radius_shrink_ratio) then
            ref_step_trust_radius = trust_radius / trust_radius_shrink_factor
        else if (ratio < trust_radius_expand_ratio) then
            ref_step_trust_radius = trust_radius
        else
            ref_step_trust_radius = trust_radius / trust_radius_expand_factor
        end if

    end function ref_step_trust_radius

    logical function check_call_counts(context, test_name, case_name, n_hess_x, &
                                       n_update_orbs)
        !
        ! this function checks that the counters reported by the library agree with
        ! the number of times the callback functions were actually called
        !
        class(hartmann6d_context_type), intent(in) :: context
        character(len=*), intent(in) :: test_name, case_name
        integer(ip), intent(in), optional :: n_hess_x, n_update_orbs

        ! assume test passes
        check_call_counts = .true.

        if (present(n_hess_x)) then
            if (n_hess_x /= context%n_hess_x_calls) then
                write(stderr, *) "test_"//test_name//" failed: Reported number of "// &
                    "Hessian linear transformations wrong "//case_name//"."
                check_call_counts = .false.
            end if
        end if
        if (present(n_update_orbs)) then
            if (n_update_orbs /= context%n_update_orbs_calls) then
                write(stderr, *) "test_"//test_name// &
                    " failed: Reported number of orbital updates wrong "//case_name//"."
                check_call_counts = .false.
            end if
        end if

    end function check_call_counts

    logical(c_bool) function test_solver() bind(C)
        !
        ! this function tests the solver subroutine
        !
        use opentrustregion, only: &
            update_orbs_type, obj_func_type, solver_settings_type, solver, &
            default_settings => default_solver_settings, error_solver_max_iter, &
            error_update_orbs, error_conv_check, error_obj_func, error_hess_x, &
            error_precond, error_project, verbosity_warning, subsystem_solver_options, &
            stability_settings_uninitialized_warning_msg
        use test_reference, only: arm_host_context, host_context_reached

        real(rp), parameter :: var_thres = 1e-6_rp
        integer(ip) :: error, i_solver, callback, i_call, n_calls(n_fault_callbacks)
        integer(ip), parameter :: fault_origins(n_fault_callbacks) = &
            [error_update_orbs, error_obj_func, error_hess_x, error_precond, &
             error_project, error_conv_check]
        character(len=20) :: call_number
        type(hartmann6d_fault_context_type), target :: fault_context
        real(rp) :: final_grad(n_param)
        procedure(update_orbs_type), pointer :: update_orbs_funptr
        procedure(obj_func_type), pointer :: obj_func_funptr
        type(solver_settings_type) :: settings, uninitialized_settings
        type(hartmann6d_recording_context_type), target :: context
        type(double_well_context_type), target :: double_well_context

        ! assume tests pass
        test_solver = .true.

        ! initialize settings
        call settings%init(error)

        ! start at saddle point
        context%vars = saddle_point
        update_orbs_funptr => update_orbs
        obj_func_funptr => obj_func

        ! check that the objective function is not evaluated without displacement,
        ! which only the line search does
        context%obj_func_zero_step = .false.

        ! hand the callback functions a host context, the internal stability check has
        ! to inherit it since it calls back into the same host
        call arm_host_context(settings, context)

        ! run solver, check if error has occured and check whether gradient is zero and
        ! agrees with correct minimum, that the reported number of orbital updates
        ! agrees with the orbital update calls and that convergence is not flagged as
        ! reaching maximum precision
        context%n_update_orbs_calls = 0
        call solver(update_orbs_funptr, obj_func_funptr, n_param, error, settings)
        if (error /= 0) then
            write(stderr, *) "test_solver failed: Produced error."
            test_solver = .false.
        end if
        call hartmann6d_gradient(context%vars, final_grad)
        if (norm2(final_grad) / sqrt(real(n_param, kind=rp)) > &
            default_settings%conv_tol) then
            write(stderr, *) "test_solver failed: Solver did not find stationary point."
            test_solver = .false.
        end if
        if (any(abs(context%vars - minimum1) > var_thres) .and. &
            any(abs(context%vars - minimum2) > var_thres)) then
            write(stderr, *) "test_solver failed: Solver did not find minimum."
            test_solver = .false.
        end if
        test_solver = test_solver .and. logical( &
            check_call_counts(context, "solver", "at saddle point", &
                              n_update_orbs=settings%n_update_orbs), kind=c_bool)
        if (settings%max_precision_reached) then
            write(stderr, *) "test_solver failed: Flagged that maximum precision "// &
                "was reached when convergence tolerance was met."
            test_solver = .false.
        end if
        if (context%obj_func_zero_step) then
            write(stderr, *) "test_solver failed: Line search performed although "// &
                "it was not requested."
            test_solver = .false.
        end if
        test_solver = test_solver .and. &
                      logical(host_context_reached("solver", context), kind=c_bool)

        ! start at saddle point again with nested stability check settings that were
        ! not initialized, the solver initializes them before handing down its context,
        ! so the internal stability check still calls the Hessian linear transformation
        ! with the solver's context, and prints a warning
        context%vars = saddle_point
        call arm_host_context(settings, context)
        settings%stability_settings%initialized = .false.
        settings%stability_settings%context => null()
        settings%logger => logger
        settings%verbose = verbosity_warning
        context%log_message = ""
        call solver(update_orbs_funptr, obj_func_funptr, n_param, error, settings)
        if (error /= 0) then
            write(stderr, *) "test_solver failed: Produced error when nested "// &
                "stability check settings were not initialized."
            test_solver = .false.
        end if
        test_solver = test_solver .and. &
                      logical(host_context_reached("solver", context), kind=c_bool)
        if (.not. settings%stability_settings%initialized) then
            write(stderr, *) "test_solver failed: Nested stability check settings "// &
                "were left uninitialized."
            test_solver = .false.
        end if
        if (index(context%log_message, &
                  " "//stability_settings_uninitialized_warning_msg) == 0) then
            write(stderr, *) "test_solver failed: Warning not printed when nested "// &
                "stability check settings were not initialized."
            test_solver = .false.
        end if

        ! start at saddle point but allow only a single macro iteration, the internal
        ! stability check finds the saddle point unstable and the solver stops before
        ! solving a trust region subproblem, so all Hessian linear transformations are
        ! the internal stability check's and have to be added to the solver's counter,
        ! a stale counter has to be reset first
        context%vars = saddle_point
        call settings%init(error)
        settings%n_macro = 1
        settings%n_hess_x = 1000
        call setup_error_logging(settings, context)
        call solver(update_orbs_funptr, obj_func_funptr, n_param, error, settings)
        if (error /= error_solver_max_iter) then
            write(stderr, *) "test_solver failed: Did not return maximum iteration "// &
                "error code when only the internal stability check runs."
            test_solver = .false.
        end if
        if (len_trim(context%log_message) == 0) then
            write(stderr, *) "test_solver failed: No error message printed when "// &
                "only the internal stability check runs."
            test_solver = .false.
        end if
        if (settings%stability_settings%n_hess_x <= 0) then
            write(stderr, *) "test_solver failed: Internal stability check "// &
                "performed no Hessian linear transformations."
            test_solver = .false.
        end if
        if (settings%n_hess_x /= settings%stability_settings%n_hess_x) then
            write(stderr, *) "test_solver failed: Hessian linear transformations "// &
                "of the internal stability check not added to the solver's counter."
            test_solver = .false.
        end if

        ! force non-convergence by allowing only a single macro iteration from a
        ! generic starting point and check that the specific maximum iteration error
        ! code is returned
        context%vars = distant_point
        call settings%init(error)
        settings%n_macro = 1
        context%n_update_orbs_calls = 0

        ! leave a stale internal stability check counter as a previous call would, no
        ! internal stability check is performed here so it has to be reset
        settings%stability_settings%n_hess_x = 1
        call setup_error_logging(settings, context)
        call solver(update_orbs_funptr, obj_func_funptr, n_param, error, settings)
        if (error /= error_solver_max_iter) then
            write(stderr, *) "test_solver failed: Did not return maximum iteration "// &
                "error code when exceeding n_macro."
            test_solver = .false.
        end if
        if (len_trim(context%log_message) == 0) then
            write(stderr, *) "test_solver failed: No error message printed when "// &
                "exceeding n_macro."
            test_solver = .false.
        end if
        if (settings%stability_settings%n_hess_x /= 0) then
            write(stderr, *) "test_solver failed: Internal stability check counter "// &
                "of a previous call was not reset."
            test_solver = .false.
        end if
        test_solver = test_solver .and. logical( &
            check_call_counts(context, "solver", "when exceeding n_macro", &
                              n_update_orbs=settings%n_update_orbs), kind=c_bool)

        ! force the maximum precision heuristic to trigger by requesting a convergence
        ! tolerance that floating-point noise in the gradient can never satisfy
        context%vars = near_minimum
        call settings%init(error)
        settings%context => context
        settings%conv_tol = 0.0_rp
        settings%subsystem_solver = "tcg"

        ! run solver, check that it still returns without error while flagging this in
        ! the settings object and that the reported number of orbital updates agrees
        ! with the calls
        context%n_update_orbs_calls = 0
        call solver(update_orbs_funptr, obj_func_funptr, n_param, error, settings)
        if (error /= 0) then
            write(stderr, *) "test_solver failed: Produced error when forcing "// &
                "maximum precision heuristic."
            test_solver = .false.
        end if
        if (.not. settings%max_precision_reached) then
            write(stderr, *) "test_solver failed: Did not flag that maximum "// &
                "precision was reached when convergence tolerance could not be met."
            test_solver = .false.
        end if
        test_solver = test_solver .and. logical( &
            check_call_counts(context, "solver", "when maximum precision is reached", &
                              n_update_orbs=settings%n_update_orbs), kind=c_bool)

        ! run solver again on the same settings object but stop it before convergence
        ! and check that the maximum precision flag is reset
        context%vars = distant_point
        settings%conv_tol = default_settings%conv_tol
        settings%n_macro = 1
        call setup_error_logging(settings, context)
        call solver(update_orbs_funptr, obj_func_funptr, n_param, error, settings)
        if (error /= error_solver_max_iter) then
            write(stderr, *) "test_solver failed: Did not return maximum iteration "// &
                "error code after maximum precision was reached in a previous call."
            test_solver = .false.
        end if
        if (len_trim(context%log_message) == 0) then
            write(stderr, *) "test_solver failed: No error message printed after "// &
                "maximum precision was reached in a previous call."
            test_solver = .false.
        end if
        if (settings%max_precision_reached) then
            write(stderr, *) "test_solver failed: Flag that maximum precision was "// &
                "reached was not reset by the next call."
            test_solver = .false.
        end if

        ! force the maximum precision heuristic to trigger again but let the
        ! convergence check pass in that same macro iteration, convergence then takes
        ! precedence and maximum precision is not flagged
        context%vars = near_minimum
        call settings%init(error)
        settings%context => context
        settings%conv_tol = 0.0_rp
        settings%subsystem_solver = "tcg"
        settings%conv_check => mock_conv_check_without_update
        context%n_update_orbs_at_conv_check = -1
        call solver(update_orbs_funptr, obj_func_funptr, n_param, error, settings)
        if (error /= 0) then
            write(stderr, *) "test_solver failed: Produced error when convergence "// &
                "check passes once maximum precision is reached."
            test_solver = .false.
        end if
        if (settings%max_precision_reached) then
            write(stderr, *) "test_solver failed: Flagged that maximum precision "// &
                "was reached when convergence check passed in the same iteration."
            test_solver = .false.
        end if

        ! run solver with a convergence check which always passes, the solver stops at
        ! the first check without taking a step
        context%vars = near_minimum
        call settings%init(error)
        settings%context => context
        settings%conv_check => mock_conv_check
        call solver(update_orbs_funptr, obj_func_funptr, n_param, error, settings)
        if (error /= 0) then
            write(stderr, *) "test_solver failed: Produced error when convergence "// &
                "check passes."
            test_solver = .false.
        end if
        if (any(abs(context%vars - near_minimum) > tol)) then
            write(stderr, *) "test_solver failed: Solver did not stop when "// &
                "convergence check passed."
            test_solver = .false.
        end if

        ! run solver, an orbital update which succeeds without providing a Hessian
        ! linear transformation is reported as an orbital update error, the orbital
        ! update is still counted
        context%vars = near_minimum
        update_orbs_funptr => update_orbs_no_hess_x
        call settings%init(error)
        context%n_update_orbs_calls = 0
        call setup_error_logging(settings, context)
        call solver(update_orbs_funptr, obj_func_funptr, n_param, error, settings)
        if (error /= error_update_orbs + 1) then
            write(stderr, *) "test_solver failed: Did not report a missing Hessian "// &
                "linear transformation."
            test_solver = .false.
        end if
        if (len_trim(context%log_message) == 0) then
            write(stderr, *) "test_solver failed: No error message printed for "// &
                "missing Hessian linear transformation."
            test_solver = .false.
        end if
        test_solver = test_solver .and. logical(check_call_counts( &
            context, "solver", "for missing Hessian linear transformation", &
            n_update_orbs=settings%n_update_orbs), kind=c_bool)

        ! run solver with settings that were not initialized and check that these are
        ! initialized, the initialization also resets their host context, so the
        ! callback functions are ones that need no test context
        update_orbs_funptr => stationary_update_orbs
        obj_func_funptr => constant_obj_func
        call solver(update_orbs_funptr, obj_func_funptr, n_param, error, &
                    uninitialized_settings)
        if (error /= 0) then
            write(stderr, *) "test_solver failed: Produced error for settings that "// &
                "were not initialized."
            test_solver = .false.
        end if
        if (.not. uninitialized_settings%initialized) then
            write(stderr, *) "test_solver failed: Settings were not initialized."
            test_solver = .false.
        end if
        update_orbs_funptr => update_orbs
        obj_func_funptr => obj_func

        ! check that the line search brackets the step starting from a vanishing step
        context%vars = near_minimum
        call settings%init(error)
        settings%context => context
        settings%line_search = .true.
        context%obj_func_zero_step = .false.
        call solver(update_orbs_funptr, obj_func_funptr, n_param, error, settings)
        if (error /= 0) then
            write(stderr, *) "test_solver failed: Produced error with line search."
            test_solver = .false.
        end if
        call hartmann6d_gradient(context%vars, final_grad)
        if (norm2(final_grad) / sqrt(real(n_param, kind=rp)) > settings%conv_tol) then
            write(stderr, *) "test_solver failed: Solver did not find stationary "// &
                "point with line search."
            test_solver = .false.
        end if
        if (.not. context%obj_func_zero_step) then
            write(stderr, *) "test_solver failed: Line search was not performed."
            test_solver = .false.
        end if

        ! switch to Jacobi-Davidson from the first micro iteration and check that the
        ! printed iterations split the micro iterations into Davidson and
        ! Jacobi-Davidson ones, which shows as zero Davidson micro iterations
        context%vars = near_minimum
        call setup_settings(settings, context)
        settings%subsystem_solver = "jacobi-davidson"
        settings%jacobi_davidson_start = 0
        call solver(update_orbs_funptr, obj_func_funptr, n_param, error, settings)
        if (error /= 0) then
            write(stderr, *) "test_solver failed: Produced error with Jacobi-Davidson."
            test_solver = .false.
        end if
        if (index(context%log_message, "|   0  |") == 0) then
            write(stderr, *) "test_solver failed: Printed iterations do not split "// &
                "micro iterations with Jacobi-Davidson."
            test_solver = .false.
        end if

        ! check that the truncated conjugate gradient solver is used instead of the
        ! Davidson solvers, whose first trial vector is the normalized gradient, and
        ! that the reported step size is measured with the preconditioner
        context%vars = near_minimum
        call settings%init(error)
        settings%context => context
        settings%subsystem_solver = "tcg"
        settings%precond => step_recording_precond
        context%hess_x_normalized_grad = .false.
        context%precond_received_step = .false.
        call solver(update_orbs_funptr, obj_func_funptr, n_param, error, settings)
        if (error /= 0) then
            write(stderr, *) "test_solver failed: Produced error with truncated "// &
                "conjugate gradient."
            test_solver = .false.
        end if
        call hartmann6d_gradient(context%vars, final_grad)
        if (norm2(final_grad) / sqrt(real(n_param, kind=rp)) > settings%conv_tol) then
            write(stderr, *) "test_solver failed: Solver did not find stationary "// &
                "point with truncated conjugate gradient."
            test_solver = .false.
        end if
        if (context%hess_x_normalized_grad) then
            write(stderr, *) "test_solver failed: Davidson solver used instead of "// &
                "truncated conjugate gradient."
            test_solver = .false.
        end if
        if (.not. context%precond_received_step) then
            write(stderr, *) "test_solver failed: Step size not measured with "// &
                "preconditioner for truncated conjugate gradient."
            test_solver = .false.
        end if

        ! converge to the saddle point of the double well function by excluding the
        ! direction of negative curvature from the reduced space and check that the
        ! requested stability check detects the saddle point and the optimization
        ! continues to a minimum
        call setup_settings(settings, double_well_context)
        settings%stability = .true.
        settings%n_random_trial_vectors = 0
        double_well_context%vars = [0.5_rp, 0.0_rp]
        update_orbs_funptr => double_well_update_orbs
        obj_func_funptr => double_well_obj_func
        call solver(update_orbs_funptr, obj_func_funptr, 2_ip, error, settings)
        if (error /= 0) then
            write(stderr, *) "test_solver failed: Produced error when saddle point "// &
                "is reached."
            test_solver = .false.
        end if
        if (abs(double_well_context%vars(2)) < 0.5_rp) then
            write(stderr, *) "test_solver failed: Did not continue to minimum "// &
                "after saddle point was reached."
            test_solver = .false.
        end if

        ! inject faults by running the solver from the saddle point, where it performs
        ! the internal stability check and the line search along the unstable mode,
        ! with a stability check at convergence, a line search and every optional
        ! callback function, once without a fault to count the calls of every callback
        ! function and then letting every one of these calls fail in turn, for every
        ! subsystem solver, and check that the error is reported with the origin of the
        ! failing callback function and that the reported counters still agree with the
        ! calls
        do i_solver = 1, size(subsystem_solver_options)
            call run_with_fault(0_ip, 0_ip)
            if (error /= 0) then
                write(stderr, *) "test_solver failed: Produced error without fault "// &
                    "with "//trim(subsystem_solver_options(i_solver))// &
                    " subsystem solver."
                test_solver = .false.
                cycle
            end if
            n_calls = fault_context%n_callback_calls
            do callback = 1, n_fault_callbacks
                do i_call = 1, n_calls(callback)
                    call run_with_fault(callback, i_call)
                    write(call_number, '(I0)') i_call
                    if (error /= fault_origins(callback) + 1) then
                        write(stderr, *) "test_solver failed: Error of failing "// &
                            trim(fault_callback_names(callback))//" at call "// &
                            trim(call_number)//" not reported with its origin with "// &
                            trim(subsystem_solver_options(i_solver))// &
                            " subsystem solver."
                        test_solver = .false.
                        exit
                    end if
                    if (.not. check_call_counts( &
                        fault_context, "solver", &
                        "for failing "//trim(fault_callback_names(callback))// &
                        " at call "//trim(call_number)//" with "// &
                        trim(subsystem_solver_options(i_solver))//" subsystem solver", &
                        n_hess_x=settings%n_hess_x, &
                        n_update_orbs=settings%n_update_orbs)) then
                        test_solver = .false.
                        exit
                    end if
                end do
            end do
        end do

    contains

        subroutine run_with_fault(fault_callback, fault_call)
            !
            ! this subroutine runs the solver from the saddle point with the subsystem
            ! solver of the current case and lets the given call of the given callback
            ! function fail
            !
            integer(ip), intent(in) :: fault_callback, fault_call

            procedure(update_orbs_type), pointer :: faulty_update_orbs_funptr
            procedure(obj_func_type), pointer :: faulty_obj_func_funptr

            ! set up fault injection context
            fault_context%vars = saddle_point
            fault_context%n_callback_calls = 0
            fault_context%n_update_orbs_calls = 0
            fault_context%n_hess_x_calls = 0
            fault_context%fault_callback = fault_callback
            fault_context%fault_call = fault_call

            ! set up settings with every optional callback function
            call settings%init(error)
            settings%context => fault_context
            settings%subsystem_solver = subsystem_solver_options(i_solver)
            settings%jacobi_davidson_start = 0
            settings%line_search = .true.
            settings%stability = .true.
            settings%precond => faulty_precond
            settings%project => faulty_project
            settings%conv_check => faulty_conv_check

            ! run solver
            faulty_update_orbs_funptr => faulty_update_orbs
            faulty_obj_func_funptr => faulty_obj_func
            call solver(faulty_update_orbs_funptr, faulty_obj_func_funptr, n_param, &
                        error, settings)

        end subroutine run_with_fault

    end function test_solver

    logical(c_bool) function test_stability_check() bind(C)
        !
        ! this function tests the stability check subroutine
        !
        use opentrustregion, only: hess_x_type, stability_settings_type, &
                                   stability_check, error_stability_check_max_iter, &
                                   verbosity_debug, verbosity_warning, error_hess_x, &
                                   error_precond, error_project, diag_solver_options, &
                                   unstable_warning_msg

        real(rp) :: vars(n_param), h_diag(n_param), direction(n_param), &
                    hess_eigvals(n_param), hess_eigvecs(n_param, n_param)
        procedure(hess_x_type), pointer :: hess_x_funptr
        logical :: stable
        integer(ip) :: error, i
        character(len=300) :: msg
        type(stability_settings_type) :: settings, uninitialized_settings
        integer(ip) :: i_solver, callback, i_call, n_calls(n_fault_callbacks)
        integer(ip), parameter :: &
            fault_callbacks(3) = [fault_hess_x, fault_precond, fault_project], &
            fault_origins(3) = [error_hess_x, error_precond, error_project]
        character(len=20) :: call_number
        type(hartmann6d_fault_context_type), target :: fault_context
        type(hartmann6d_context_type), target :: context

        ! assume tests pass
        test_stability_check = .true.

        ! start at minimum and determine Hessian diagonal and define Hessian linear
        ! transformation
        vars = minimum1
        context%hess = hartmann6d_hessian(vars)
        h_diag = [(context%hess(i, i), i=1, size(h_diag))]
        hess_x_funptr => hess_x_fun

        ! initialize settings
        call settings%init(error)
        settings%context => context

        ! run stability check, check if error has occured and determine whether minimum
        ! is stable and the returned direction vanishes
        context%n_hess_x_calls = 0
        call stability_check(h_diag, hess_x_funptr, stable, error, settings, direction)
        if (error /= 0) then
            write(stderr, *) "test_stability_check failed: Produced error for minimum."
            test_stability_check = .false.
        end if
        if (.not. stable) then
            write(stderr, *) "test_stability_check failed: Stability check "// &
                "incorrectly classifies stability of minimum."
            test_stability_check = .false.
        end if
        if (any(abs(direction) > tol)) then
            write(stderr, *) "test_stability_check failed: Stability check does "// &
                "not return zero vector for minimum."
            test_stability_check = .false.
        end if
        test_stability_check = test_stability_check .and. logical(check_call_counts( &
            context, "stability_check", "for minimum", n_hess_x=settings%n_hess_x), &
            kind=c_bool)

        ! start at saddle point and determine Hessian diagonal, define linear
        ! transformation and determine the eigenvector of the lowest Hessian eigenvalue
        ! independently
        vars = saddle_point
        context%hess = hartmann6d_hessian(vars)
        call ref_symm_mat_diag(context%hess, hess_eigvals, hess_eigvecs)
        h_diag = [(context%hess(i, i), i=1, size(h_diag))]
        hess_x_funptr => hess_x_fun

        ! run stability check, check if error has occured and determine whether saddle
        ! point is unstable and the returned direction is correct and a warning with
        ! the lowest eigenvalue is printed, the settings object is reused so the
        ! reported counter must not include the previous call
        context%n_hess_x_calls = 0
        settings%logger => logger
        settings%verbose = verbosity_warning
        context%log_message = ""
        call stability_check(h_diag, hess_x_funptr, stable, error, settings, direction)
        if (error /= 0) then
            write(stderr, *) "test_stability_check failed: Produced error for "// &
                "saddle point."
            test_stability_check = .false.
        end if
        if (stable) then
            write(stderr, *) "test_stability_check failed: Stability check "// &
                "incorrectly classifies stability of saddle point."
            test_stability_check = .false.
        end if
        if (abs(abs(dot_product(direction, hess_eigvecs(:, 1))) - 1.0_rp) > tol) then
            write(stderr, *) "test_stability_check failed: Stability check does "// &
                "not return correct direction for saddle point."
            test_stability_check = .false.
        end if
        write(msg, '(A, F0.4)') unstable_warning_msg, hess_eigvals(1)
        if (adjustl(context%log_message) /= trim(msg)) then
            write(stderr, *) "test_stability_check failed: Warning not printed for "// &
                "saddle point."
            test_stability_check = .false.
        end if
        test_stability_check = test_stability_check .and. logical( &
            check_call_counts(context, "stability_check", "for saddle point", &
                              n_hess_x=settings%n_hess_x), kind=c_bool)

        ! force non-convergence by allowing only a single iteration and check that the
        ! specific maximum iteration error code is returned
        call settings%init(error)
        settings%n_iter = 1
        call setup_error_logging(settings, context)
        call stability_check(h_diag, hess_x_funptr, stable, error, settings)
        if (error /= error_stability_check_max_iter) then
            write(stderr, *) "test_stability_check failed: Did not return maximum "// &
                "iteration error code when exceeding n_iter."
            test_stability_check = .false.
        end if
        if (len_trim(context%log_message) == 0) then
            write(stderr, *) "test_stability_check failed: No error message "// &
                "printed when exceeding n_iter."
            test_stability_check = .false.
        end if

        ! force the reduced space to grow until it spans the full parameter space by
        ! setting an unreachable convergence tolerance
        call settings%init(error)
        settings%context => context
        settings%conv_tol = 0.0_rp

        call stability_check(h_diag, hess_x_funptr, stable, error, settings, direction)
        if (error /= 0) then
            write(stderr, *) "test_stability_check failed: Produced error when "// &
                "reduced space grows to dimension of full parameter space."
            test_stability_check = .false.
        end if
        if (stable) then
            write(stderr, *) "test_stability_check failed: Stability check "// &
                "incorrectly classifies stability of saddle point when reduced "// &
                "space grows to dimension of full parameter space."
            test_stability_check = .false.
        end if
        if (abs(abs(dot_product(direction, hess_eigvecs(:, 1))) - 1.0_rp) > tol) then
            write(stderr, *) "test_stability_check failed: Stability check does "// &
                "not return correct direction for saddle point when reduced space "// &
                "grows to dimension of full parameter space."
            test_stability_check = .false.
        end if

        ! run stability check at saddle point with the Jacobi-Davidson method switched
        ! on after the first iteration, check that the correction equations are solved,
        ! which is logged at debug verbosity, so that the Hessian linear transformations
        ! performed in the Jacobi-Davidson correction are counted as well
        hess_x_funptr => hess_x_fun
        call setup_settings(settings, context)
        settings%verbose = verbosity_debug
        settings%diag_solver = "jacobi-davidson"
        settings%jacobi_davidson_start = 1
        context%n_hess_x_calls = 0
        call stability_check(h_diag, hess_x_funptr, stable, error, settings, direction)
        if (error /= 0) then
            write(stderr, *) "test_stability_check failed: Produced error with "// &
                "Jacobi-Davidson method."
            test_stability_check = .false.
        end if
        if (index(context%log_message, "MINRES") == 0) then
            write(stderr, *) "test_stability_check failed: Jacobi-Davidson "// &
                "correction equations not solved."
            test_stability_check = .false.
        end if
        if (stable) then
            write(stderr, *) "test_stability_check failed: Stability check "// &
                "incorrectly classifies stability of saddle point with "// &
                "Jacobi-Davidson method."
            test_stability_check = .false.
        end if
        test_stability_check = test_stability_check .and. logical(check_call_counts( &
            context, "stability_check", "with Jacobi-Davidson method", &
            n_hess_x=settings%n_hess_x), kind=c_bool)

        ! check that the Jacobi-Davidson method is used neither by the Davidson method
        ! nor before the iteration at which it is requested to start, which the
        ! stability check never reaches for a problem of this size
        call setup_settings(settings, context)
        settings%verbose = verbosity_debug
        settings%jacobi_davidson_start = 0
        call stability_check(h_diag, hess_x_funptr, stable, error, settings, direction)
        if (error /= 0 .or. index(context%log_message, "MINRES") /= 0) then
            write(stderr, *) "test_stability_check failed: Produced error or "// &
                "Jacobi-Davidson correction equations solved with Davidson method."
            test_stability_check = .false.
        end if
        call setup_settings(settings, context)
        settings%verbose = verbosity_debug
        settings%diag_solver = "jacobi-davidson"
        settings%jacobi_davidson_start = settings%n_iter
        call stability_check(h_diag, hess_x_funptr, stable, error, settings, direction)
        if (error /= 0 .or. index(context%log_message, "MINRES") /= 0) then
            write(stderr, *) "test_stability_check failed: Produced error or "// &
                "Jacobi-Davidson correction equations solved before "// &
                "Jacobi-Davidson method is started."
            test_stability_check = .false.
        end if

        ! run stability check with settings that were not initialized and check that
        ! these are initialized, the initialization also resets their host context, so
        ! the Hessian linear transformation is one that needs no test context
        hess_x_funptr => identity_hess_x
        call stability_check(h_diag, hess_x_funptr, stable, error, &
                             uninitialized_settings)
        if (error /= 0) then
            write(stderr, *) "test_stability_check failed: Produced error for "// &
                "settings that were not initialized."
            test_stability_check = .false.
        end if
        if (.not. uninitialized_settings%initialized) then
            write(stderr, *) "test_stability_check failed: Settings were not "// &
                "initialized."
            test_stability_check = .false.
        end if

        ! start at saddle point with a preconditioner which always returns the first
        ! unit vector and a Hessian diagonal whose minimum lies at the first element,
        ! so that the first new trial vector is linearly dependent on the first trial
        ! vector, and check that the stability check stops without error and without
        ! printing an error message
        context%hess = hartmann6d_hessian(saddle_point)
        h_diag = 1.0_rp
        h_diag(1) = 0.0_rp
        hess_x_funptr => hess_x_fun
        call settings%init(error)
        settings%precond => mock_precond_first_unit_vector
        call setup_error_logging(settings, context)
        call stability_check(h_diag, hess_x_funptr, stable, error, settings, direction)
        if (error /= 0) then
            write(stderr, *) "test_stability_check failed: Produced error when new "// &
                "trial vector is linearly dependent."
            test_stability_check = .false.
        end if
        if (len_trim(context%log_message) /= 0) then
            write(stderr, *) "test_stability_check failed: Error message printed "// &
                "when new trial vector is linearly dependent."
            test_stability_check = .false.
        end if

        ! inject faults by running the stability check at the saddle point with every
        ! optional callback function, once without a fault to count the calls of every
        ! callback function and then letting every one of these calls fail in turn, for
        ! every diagonalization solver, and check that the error is reported with the
        ! origin of the failing callback function and that the reported counter still
        ! agrees with the calls
        fault_context%hess = hartmann6d_hessian(saddle_point)
        h_diag = [(fault_context%hess(i, i), i=1, size(h_diag))]
        hess_x_funptr => faulty_hess_x
        do i_solver = 1, size(diag_solver_options)
            call run_with_fault(0_ip, 0_ip)
            if (error /= 0) then
                write(stderr, *) "test_stability_check failed: Produced error "// &
                    "without fault with "//trim(diag_solver_options(i_solver))// &
                    " diagonalization solver."
                test_stability_check = .false.
                cycle
            end if
            n_calls = fault_context%n_callback_calls
            do callback = 1, size(fault_callbacks)
                do i_call = 1, n_calls(fault_callbacks(callback))
                    call run_with_fault(fault_callbacks(callback), i_call)
                    write(call_number, '(I0)') i_call
                    if (error /= fault_origins(callback) + 1) then
                        write(stderr, *) "test_stability_check failed: Error of "// &
                            "failing "// &
                            trim(fault_callback_names(fault_callbacks(callback)))// &
                            " at call "//trim(call_number)//" not reported with "// &
                            "its origin with "//trim(diag_solver_options(i_solver))// &
                            " diagonalization solver."
                        test_stability_check = .false.
                        exit
                    end if
                    if (.not. check_call_counts( &
                        fault_context, "stability_check", "for failing "// &
                        trim(fault_callback_names(fault_callbacks(callback)))// &
                        " at call "//trim(call_number)//" with "// &
                        trim(diag_solver_options(i_solver))// &
                        " diagonalization solver", n_hess_x=settings%n_hess_x)) then
                        test_stability_check = .false.
                        exit
                    end if
                end do
            end do
        end do

    contains

        subroutine run_with_fault(fault_callback, fault_call)
            !
            ! this subroutine runs the stability check at the saddle point with the
            ! diagonalization solver of the current case and lets the given call of the
            ! given callback function fail
            !
            integer(ip), intent(in) :: fault_callback, fault_call

            ! set up fault injection context
            fault_context%n_callback_calls = 0
            fault_context%n_hess_x_calls = 0
            fault_context%fault_callback = fault_callback
            fault_context%fault_call = fault_call

            ! set up settings with every optional callback function
            call settings%init(error)
            settings%context => fault_context
            settings%diag_solver = diag_solver_options(i_solver)
            settings%jacobi_davidson_start = 0
            settings%precond => faulty_precond
            settings%project => faulty_project

            ! run stability check
            call stability_check(h_diag, hess_x_funptr, stable, error, settings, &
                                 direction)

        end subroutine run_with_fault

    end function test_stability_check

    logical(c_bool) function test_check_stationary_point() bind(C)
        !
        ! this function tests the subroutine which performs the stability check at a
        ! stationary point and the line search along an unstable mode
        !
        use opentrustregion, only: &
            solver_settings_type, hess_x_type, obj_func_type, check_stationary_point, &
            verbosity_debug, error_solver, error_stability_check, error_obj_func, &
            started_at_saddle_point_warning_msg, reached_saddle_point_warning_msg
        use test_reference, only: arm_host_context, host_context_reached

        type(solver_settings_type) :: settings
        type(double_well_context_type), target :: context, stability_context
        procedure(hess_x_type), pointer :: hess_x_funptr
        procedure(obj_func_type), pointer :: obj_func_funptr
        real(rp) :: kappa(2)
        logical :: stable
        integer(ip) :: error

        ! the double well function x^2 - y^2 + y^4 has the Hessian diagonal (2, -2) at
        ! its saddle point at the origin, where the logarithmic line search along the
        ! unstable mode first lowers the function at a tenth to the power of one half
        ! of the unit step, and the Hessian diagonal (2, 4) at its minima
        real(rp), parameter :: saddle_h_diag(2) = [2.0_rp, -2.0_rp], &
                               minimum_h_diag(2) = [2.0_rp, 4.0_rp], &
                               saddle_step = 10.0_rp**(-0.5_rp)

        ! assume tests pass
        test_check_stationary_point = .true.

        ! set function pointers
        hess_x_funptr => double_well_hess_x
        obj_func_funptr => double_well_obj_func

        ! check the saddle point in the first macro iteration with nested stability
        ! check settings without callback functions, verbosity or context of their own,
        ! which have to inherit the solver's for the duration of the check, and check
        ! that the Hessian linear transformations of the check are added to the
        ! solver's counter, that the step along the unstable mode is returned and that
        ! the warning for starting at a saddle point is printed
        call setup_settings(settings, context)
        settings%precond => identity_precond
        settings%project => identity_project
        settings%n_hess_x = 5
        context%vars = [0.0_rp, 0.0_rp]
        call arm_host_context(settings, context)
        call check_stationary_point(0.0_rp, saddle_h_diag, hess_x_funptr, &
                                    obj_func_funptr, 1_ip, settings, stable, kappa, &
                                    error)
        if (error /= 0 .or. stable) then
            write(stderr, *) "test_check_stationary_point failed: Produced error "// &
                "or saddle point not found to be unstable."
            test_check_stationary_point = .false.
        end if
        if (abs(kappa(1)) > tol .or. abs(abs(kappa(2)) - saddle_step) > tol) then
            write(stderr, *) "test_check_stationary_point failed: Step along "// &
                "unstable mode not correct."
            test_check_stationary_point = .false.
        end if
        if (.not. ( &
            associated(settings%stability_settings%precond, identity_precond) .and. &
            associated(settings%stability_settings%project, identity_project) .and. &
            associated(settings%stability_settings%logger, logger))) then
            write(stderr, *) "test_check_stationary_point failed: Solver's "// &
                "callback functions not inherited."
            test_check_stationary_point = .false.
        end if
        if (settings%stability_settings%verbose /= settings%verbose) then
            write(stderr, *) "test_check_stationary_point failed: Solver's "// &
                "verbosity not inherited."
            test_check_stationary_point = .false.
        end if
        test_check_stationary_point = test_check_stationary_point .and. logical( &
            host_context_reached("check_stationary_point", context), kind=c_bool)
        if (associated(settings%stability_settings%context)) then
            write(stderr, *) "test_check_stationary_point failed: Host context "// &
                "lent to the stability check left on the nested settings."
            test_check_stationary_point = .false.
        end if
        if (settings%stability_settings%n_hess_x <= 0 .or. &
            settings%n_hess_x /= 5 + settings%stability_settings%n_hess_x) then
            write(stderr, *) "test_check_stationary_point failed: Hessian linear "// &
                "transformations of the stability check not added to the solver's "// &
                "counter."
            test_check_stationary_point = .false.
        end if
        if (index(context%log_message, " "//started_at_saddle_point_warning_msg) == 0) &
            then
            write(stderr, *) "test_check_stationary_point failed: Warning not "// &
                "printed for starting at saddle point."
            test_check_stationary_point = .false.
        end if

        ! check the saddle point in a later macro iteration with nested stability check
        ! settings with a verbosity above the solver's and a context of their own,
        ! which they have to keep, and check that the warning for reaching a saddle
        ! point is printed
        call setup_settings(settings, context)
        context%vars = [0.0_rp, 0.0_rp]
        stability_context%vars = [0.0_rp, 0.0_rp]
        stability_context%n_calls = 0
        settings%stability_settings%verbose = verbosity_debug
        settings%stability_settings%context => stability_context
        call check_stationary_point(0.0_rp, saddle_h_diag, hess_x_funptr, &
                                    obj_func_funptr, 2_ip, settings, stable, kappa, &
                                    error)
        if (error /= 0 .or. stable) then
            write(stderr, *) "test_check_stationary_point failed: Produced error "// &
                "or saddle point not found to be unstable in later macro iteration."
            test_check_stationary_point = .false.
        end if
        if (settings%stability_settings%verbose /= verbosity_debug) then
            write(stderr, *) "test_check_stationary_point failed: Verbosity of "// &
                "nested settings above the solver's replaced."
            test_check_stationary_point = .false.
        end if
        if (.not. associated(settings%stability_settings%context, &
                             stability_context) .or. stability_context%n_calls == 0) &
            then
            write(stderr, *) "test_check_stationary_point failed: Context of "// &
                "nested settings replaced or not handed to the stability check."
            test_check_stationary_point = .false.
        end if
        if (index(context%log_message, " "//reached_saddle_point_warning_msg) == 0) then
            write(stderr, *) "test_check_stationary_point failed: Warning not "// &
                "printed for reaching saddle point."
            test_check_stationary_point = .false.
        end if

        ! check a minimum with nested stability check settings which provide their own
        ! logging function while the solver has none, which they have to keep, and
        ! check that the minimum is found to be stable without a step
        call setup_settings(settings, context)
        settings%logger => null()
        settings%stability_settings%logger => logger
        context%vars = [0.0_rp, 1.0_rp / sqrt(2.0_rp)]
        call check_stationary_point(double_well_func(context%vars), minimum_h_diag, &
                                    hess_x_funptr, obj_func_funptr, 2_ip, settings, &
                                    stable, kappa, error)
        if (error /= 0 .or. .not. stable) then
            write(stderr, *) "test_check_stationary_point failed: Produced error "// &
                "or minimum not found to be stable."
            test_check_stationary_point = .false.
        end if
        if (any(abs(kappa) > tol)) then
            write(stderr, *) "test_check_stationary_point failed: Step returned "// &
                "for minimum."
            test_check_stationary_point = .false.
        end if
        if (.not. associated(settings%stability_settings%logger, logger)) then
            write(stderr, *) "test_check_stationary_point failed: Logging function "// &
                "of nested settings replaced."
            test_check_stationary_point = .false.
        end if

        ! check that an error is returned and an error message is printed when no step
        ! along the unstable mode lowers the objective function
        call setup_settings(settings, context)
        call setup_error_logging(settings, context)
        context%vars = [0.0_rp, 0.0_rp]
        obj_func_funptr => raised_double_well_obj_func
        call check_stationary_point(0.0_rp, saddle_h_diag, hess_x_funptr, &
                                    obj_func_funptr, 1_ip, settings, stable, kappa, &
                                    error)
        if (error /= error_solver + 1) then
            write(stderr, *) "test_check_stationary_point failed: No error "// &
                "returned when line search along unstable mode fails."
            test_check_stationary_point = .false.
        end if
        if (len_trim(context%log_message) == 0) then
            write(stderr, *) "test_check_stationary_point failed: No error message "// &
                "printed when line search along unstable mode fails."
            test_check_stationary_point = .false.
        end if

        ! check that the errors of a failing objective function during the line search
        ! and of the stability check are reported with their origin
        call setup_settings(settings, context)
        obj_func_funptr => obj_func_failing
        call check_stationary_point(0.0_rp, saddle_h_diag, hess_x_funptr, &
                                    obj_func_funptr, 1_ip, settings, stable, kappa, &
                                    error)
        if (error /= error_obj_func + 1) then
            write(stderr, *) "test_check_stationary_point failed: Error of failing "// &
                "objective function not reported with its origin."
            test_check_stationary_point = .false.
        end if
        call setup_settings(settings, context)
        obj_func_funptr => double_well_obj_func
        settings%stability_settings%diag_solver = "unknown"
        call check_stationary_point(0.0_rp, saddle_h_diag, hess_x_funptr, &
                                    obj_func_funptr, 1_ip, settings, stable, kappa, &
                                    error)
        if (error /= error_stability_check + 1) then
            write(stderr, *) "test_check_stationary_point failed: Error of "// &
                "stability check not reported with its origin."
            test_check_stationary_point = .false.
        end if

    end function test_check_stationary_point

    logical(c_bool) function test_newton_step() bind(C)
        !
        ! this function tests the Newton step subroutine
        !
        use opentrustregion, only: solver_settings_type, newton_step

        type(solver_settings_type) :: settings
        integer(ip), parameter :: n_trial = 3
        real(rp) :: red_space_basis(n_param, n_trial), grad_norm, &
                    aug_hess(n_trial + 1, n_trial + 1), solution(n_param), &
                    red_space_solution(n_trial), red_space_grad(n_trial)
        integer(ip) :: error
        type(test_context_type), target :: context

        ! assume tests pass
        test_newton_step = .true.

        ! setup settings object
        call setup_settings(settings, context)

        ! generate reduced space basis, gradient norm and augmented Hessian whose
        ! reduced space Hessian is diagonally dominant and therefore nonsingular
        call random_number(red_space_basis)
        call random_number(grad_norm)
        aug_hess = 0.0_rp
        aug_hess(2:, 2:) = generate_random_symm_matrix(n_trial) + &
                           2 * n_trial * identity_matrix(n_trial)

        ! perform Newton step, check if error has occured and determine whether
        ! resulting solution solves the Newton equations in reduced space and is
        ! correctly transformed to the full space
        call newton_step(aug_hess, grad_norm, red_space_basis, solution, &
                         red_space_solution, settings, error)
        if (error /= 0) then
            write(stderr, *) "test_newton_step failed: Produced error."
            test_newton_step = .false.
        end if
        red_space_grad = 0.0_rp
        red_space_grad(1) = grad_norm
        if (any( &
            abs(matmul(aug_hess(2:, 2:), red_space_solution) + red_space_grad) > tol)) &
            then
            write(stderr, *) "test_newton_step failed: Reduced space solution does "// &
                "not solve Newton equations."
            test_newton_step = .false.
        end if
        if (any(abs(solution - matmul(red_space_basis, red_space_solution)) > tol)) then
            write(stderr, *) "test_newton_step failed: Full space solution not correct."
            test_newton_step = .false.
        end if

        ! use a vanishing reduced space Hessian for which the Newton equations cannot
        ! be solved and check that an error is returned and an error message is printed
        aug_hess = 0.0_rp
        call setup_error_logging(settings, context)
        call newton_step(aug_hess, grad_norm, red_space_basis, solution, &
                         red_space_solution, settings, error)
        if (error == 0) then
            write(stderr, *) "test_newton_step failed: No error returned for "// &
                "singular reduced space Hessian."
            test_newton_step = .false.
        end if
        if (len_trim(context%log_message) == 0) then
            write(stderr, *) "test_newton_step failed: No error message printed "// &
                "for singular reduced space Hessian."
            test_newton_step = .false.
        end if

    end function test_newton_step

    logical(c_bool) function test_bisection() bind(C)
        !
        ! this function tests the bisection subroutine
        !
        use opentrustregion, only: solver_settings_type, bisection

        type(solver_settings_type) :: settings
        integer(ip), parameter :: n_trial = 3
        real(rp) :: red_space_basis(n_param, n_trial), grad(n_param), grad_norm, &
                    hess(n_param, n_param), aug_hess(n_trial + 1, n_trial + 1), &
                    red_space_hess_eigvals(n_trial), &
                    red_space_hess_eigvecs(n_trial, n_trial), solution(n_param), &
                    red_space_solution(n_trial), trust_radius, mu, &
                    grad_coupled_component, newton_solution(n_param), &
                    newton_red_space_solution(n_trial)
        integer(ip) :: error
        type(test_context_type), target :: context

        ! assume tests pass
        test_bisection = .true.

        ! setup settings object
        call setup_settings(settings, context)

        ! define an orthonormal reduced space basis whose first vector is not the
        ! gradient direction
        red_space_basis = reshape( &
            [1.0_rp / sqrt(2.0_rp), -1.0_rp / sqrt(2.0_rp), 0.0_rp, 0.0_rp, 0.0_rp, &
             0.0_rp, 1.0_rp / sqrt(6.0_rp), 1.0_rp / sqrt(6.0_rp), &
             -2.0_rp / sqrt(6.0_rp), 0.0_rp, 0.0_rp, 0.0_rp, 1.0_rp / sqrt(12.0_rp), &
             1.0_rp / sqrt(12.0_rp), 1.0_rp / sqrt(12.0_rp), -3.0_rp / sqrt(12.0_rp), &
             0.0_rp, 0.0_rp], [n_param, n_trial])

        ! choose target trust radius
        trust_radius = 0.4_rp

        ! project the Hessian at a point with strong negative curvature, where the
        ! level shift has to be bisected
        call set_up_reduced_hessian( &
            [0.29_rp, 0.47_rp, 0.66_rp, 0.41_rp, 0.23_rp, 0.26_rp])
        test_bisection = &
            test_bisection .and. &
            bisected_step_correct("at point with strong negative curvature")

        ! project the Hessian in the quadratic region near minimum and determine
        ! whether routine correctly falls back to the Newton step since the minimum is
        ! closer than the target trust radius and no level shift is necessary
        call set_up_reduced_hessian(near_minimum)
        call bisection(aug_hess, grad_norm, red_space_basis, red_space_hess_eigvals, &
                       red_space_hess_eigvecs, trust_radius, solution, &
                       red_space_solution, mu, settings, error)
        if (error /= 0) then
            write(stderr, *) "test_bisection failed: Produced error instead of "// &
                "falling back to Newton step."
            test_bisection = .false.
        end if
        if (abs(mu) > tol) then
            write(stderr, *) "test_bisection failed: Level shift not zero for "// &
                "Newton step fallback."
            test_bisection = .false.
        end if
        newton_red_space_solution = &
            -grad_norm * matmul(red_space_hess_eigvecs, &
                                red_space_hess_eigvecs(1, :) / red_space_hess_eigvals)
        newton_solution = matmul(red_space_basis, newton_red_space_solution)
        if (any(abs(red_space_solution - newton_red_space_solution) > tol)) then
            write(stderr, *) "test_bisection failed: Newton step fallback reduced "// &
                "space solution not correct."
            test_bisection = .false.
        end if
        if (any(abs(solution - newton_solution) > tol)) then
            write(stderr, *) "test_bisection failed: Newton step fallback full "// &
                "space solution not correct."
            test_bisection = .false.
        end if

        ! set up a reduced space Hessian in an orthonormal basis with a negative
        ! eigenvalue whose eigenvector has a gradient component and which is large
        ! compared to the gradient, so that the step at the starting alpha is longer
        ! than the trust radius and alpha has to be increased repeatedly to bracket the
        ! trust radius
        red_space_basis = 0.0_rp
        red_space_basis(1, 1) = 1.0_rp
        red_space_basis(2, 2) = 1.0_rp
        red_space_basis(3, 3) = 1.0_rp
        grad_norm = 0.1_rp
        trust_radius = 0.5_rp
        call set_up_diagonal_reduced_hessian([-10.0_rp, 2.0_rp, 3.0_rp])
        test_bisection = test_bisection .and. &
                         bisected_step_correct("when alpha has to be increased")

        ! set up hard case, in which the lowest reduced space Hessian eigenvalue is
        ! negative and its eigenvector has no component along the gradient direction,
        ! so that no level shift can reproduce the trust region solution and it must be
        ! constructed directly from the eigendecomposition, the hard case step assumes
        ! an orthonormal basis
        grad_norm = 1.0_rp
        trust_radius = 0.6_rp
        call set_up_diagonal_reduced_hessian([5.0_rp, -2.0_rp, 3.0_rp])
        call bisection(aug_hess, grad_norm, red_space_basis, red_space_hess_eigvals, &
                       red_space_hess_eigvecs, trust_radius, solution, &
                       red_space_solution, mu, settings, error)
        if (error /= 0) then
            write(stderr, *) "test_bisection failed: Produced error for hard case."
            test_bisection = .false.
        end if
        if (abs(norm2(solution) - trust_radius) > tol) then
            write(stderr, *) "test_bisection failed: Hard case solution does not "// &
                "respect trust radius."
            test_bisection = .false.
        end if
        if (abs(mu - minval(red_space_hess_eigvals)) > tol) then
            write(stderr, *) "test_bisection failed: Hard case level shift not "// &
                "equal to lowest reduced space Hessian eigenvalue."
            test_bisection = .false.
        end if
        grad_coupled_component = grad_norm / (minval(red_space_hess_eigvals) - &
                                              maxval(red_space_hess_eigvals))
        if (any(abs(red_space_solution - [ &
            grad_coupled_component, sqrt(trust_radius**2 - grad_coupled_component**2), &
            0.0_rp]) > tol)) then
            write(stderr, *) "test_bisection failed: Hard case reduced space "// &
                "solution not correct."
            test_bisection = .false.
        end if
        if (any(abs(solution - matmul(red_space_basis, red_space_solution)) > tol)) then
            write(stderr, *) "test_bisection failed: Hard case full space solution "// &
                "not correct."
            test_bisection = .false.
        end if

        ! test hard case with a trust radius below the norm of the solution at the
        ! crossover point grad_norm / (5 - (-2)), where the solution cannot be filled
        ! up with the lowest eigenvector and the level shift is instead bisected
        ! starting from the crossover point
        trust_radius = 0.1_rp
        test_bisection = test_bisection .and. bisected_step_correct( &
            "for hard case with trust radius below crossover point")

        ! set up hard case with degenerate lowest eigenvalues, the eigenvector spanning
        ! the degenerate subspace used to fill the trust radius is not uniquely
        ! defined, so only invariant properties of the solution are checked rather than
        ! exact reduced space solution components
        trust_radius = 0.5_rp
        call set_up_diagonal_reduced_hessian([5.0_rp, -2.0_rp, -2.0_rp])
        call bisection(aug_hess, grad_norm, red_space_basis, red_space_hess_eigvals, &
                       red_space_hess_eigvecs, trust_radius, solution, &
                       red_space_solution, mu, settings, error)
        if (error /= 0) then
            write(stderr, *) "test_bisection failed: Produced error for degenerate "// &
                "hard case."
            test_bisection = .false.
        end if
        if (abs(norm2(solution) - trust_radius) > tol) then
            write(stderr, *) "test_bisection failed: Degenerate hard case solution "// &
                "does not respect trust radius."
            test_bisection = .false.
        end if
        if (abs(mu - minval(red_space_hess_eigvals)) > tol) then
            write(stderr, *) "test_bisection failed: Degenerate hard case level "// &
                "shift not equal to lowest reduced space Hessian eigenvalue."
            test_bisection = .false.
        end if
        grad_coupled_component = grad_norm / (minval(red_space_hess_eigvals) - &
                                              maxval(red_space_hess_eigvals))
        if (abs(red_space_solution(1) - grad_coupled_component) > tol) then
            write(stderr, *) "test_bisection failed: Degenerate hard case "// &
                "component along non-degenerate eigenvector not correct."
            test_bisection = .false.
        end if
        if (abs(norm2(red_space_solution(2:3)) - &
                sqrt(trust_radius**2 - grad_coupled_component**2)) > tol) then
            write(stderr, *) "test_bisection failed: Degenerate hard case fill "// &
                "magnitude within degenerate eigenspace not correct."
            test_bisection = .false.
        end if
        if (any(abs(solution - matmul(red_space_basis, red_space_solution)) > tol)) then
            write(stderr, *) "test_bisection failed: Degenerate hard case full "// &
                "space solution not correct."
            test_bisection = .false.
        end if

    contains

        subroutine set_up_reduced_hessian(vars)
            !
            ! this subroutine sets the gradient norm and the augmented Hessian with
            ! the Hartmann 6D Hessian at a point projected onto the reduced space and
            ! diagonalizes the reduced space Hessian
            !
            real(rp), intent(in) :: vars(n_param)

            call hartmann6d_gradient(vars, grad)
            grad_norm = norm2(grad)
            hess = hartmann6d_hessian(vars)
            aug_hess = 0.0_rp
            aug_hess(2:, 2:) = matmul(transpose(red_space_basis), &
                                      matmul(hess, red_space_basis))
            call ref_symm_mat_diag(aug_hess(2:, 2:), red_space_hess_eigvals, &
                                   red_space_hess_eigvecs)

        end subroutine set_up_reduced_hessian

        subroutine set_up_diagonal_reduced_hessian(diagonal)
            !
            ! this subroutine sets the augmented Hessian with a diagonal reduced space
            ! Hessian and diagonalizes the reduced space Hessian
            !
            real(rp), intent(in) :: diagonal(n_trial)

            integer(ip) :: i_diag

            aug_hess = 0.0_rp
            do i_diag = 1, n_trial
                aug_hess(i_diag + 1, i_diag + 1) = diagonal(i_diag)
            end do
            call ref_symm_mat_diag(aug_hess(2:, 2:), red_space_hess_eigvals, &
                                   red_space_hess_eigvecs)

        end subroutine set_up_diagonal_reduced_hessian

        logical function bisected_step_correct(case_name)
            !
            ! this function performs bisection and checks whether error has occured
            ! and whether the resulting solution respects the target trust radius,
            ! solves the level-shifted Newton equations in reduced space for a level
            ! shift below the lowest reduced space Hessian eigenvalue and is correctly
            ! transformed to the full space
            !
            character(len=*), intent(in) :: case_name

            real(rp) :: red_space_grad(n_trial)

            ! assume test passes
            bisected_step_correct = .true.

            call bisection(aug_hess, grad_norm, red_space_basis, &
                           red_space_hess_eigvals, red_space_hess_eigvecs, &
                           trust_radius, solution, red_space_solution, mu, settings, &
                           error)
            if (error /= 0) then
                write(stderr, *) "test_bisection failed: Produced error "//case_name// &
                    "."
                bisected_step_correct = .false.
            end if
            if (abs(norm2(solution) - trust_radius) > tol) then
                write(stderr, *) "test_bisection failed: Solution does not respect "// &
                    "trust radius "//case_name//"."
                bisected_step_correct = .false.
            end if
            red_space_grad = 0.0_rp
            red_space_grad(1) = grad_norm
            if (norm2(matmul(aug_hess(2:, 2:), red_space_solution) - &
                      mu * red_space_solution + red_space_grad) > tol) then
                write(stderr, *) "test_bisection failed: Reduced space solution "// &
                    "does not solve level-shifted Newton equations "//case_name//"."
                bisected_step_correct = .false.
            end if
            if (mu >= minval(red_space_hess_eigvals)) then
                write(stderr, *) "test_bisection failed: Level shift not below "// &
                    "lowest reduced space Hessian eigenvalue "//case_name//"."
                bisected_step_correct = .false.
            end if
            if (any(abs(solution - matmul(red_space_basis, red_space_solution)) > &
                    tol)) then
                write(stderr, *) "test_bisection failed: Full space solution not "// &
                    "correct "//case_name//"."
                bisected_step_correct = .false.
            end if

        end function bisected_step_correct

    end function test_bisection

    logical(c_bool) function test_bracket() bind(C)
        !
        ! this function tests the bracketing function
        !
        use opentrustregion, only: solver_settings_type, obj_func_type, bracket, &
                                   error_obj_func

        type(solver_settings_type) :: settings
        procedure(obj_func_type), pointer :: obj_func_funptr
        real(rp) :: n
        integer(ip) :: error, i
        integer(ip), parameter :: n_cases = 9
        real(rp), parameter :: direction(1) = [1.0_rp]
        character(len=*), parameter :: case_names(n_cases) = &
            [character(len=26) :: "(x - 3)^2", "(x - 2)^4", "exp(-x) + x / 100", &
             "|x - 2.2|", "sqrt(|x - 5|)", "(x - 1000)^2", "(x - 2)^2 with bump at 2", &
             "-x", "constant"]
        logical, parameter :: expect_error(n_cases) = &
            [.false., .false., .false., .false., .false., .false., .false., .true., &
             .true.]

        ! expected multipliers, the first one at the minimum where the parabolic fit
        ! is exact, the golden ratio steps 1 + phi and 3 + sqrt(5) for functions that
        ! reject the parabolic fits and the remaining ones reproduced with an
        ! independent implementation of the bracketing algorithm
        real(rp), parameter :: expected_n(n_cases) = &
            [3.0_rp, 1.8567627457812104_rp, 4.141367271592155_rp, &
             (3.0_rp + sqrt(5.0_rp)) / 2.0_rp, 3.0_rp + sqrt(5.0_rp), &
             1000.0000000003488_rp, 1.0_rp, 0.0_rp, 0.0_rp]
        type(bracket_context_type), target :: context

        ! assume tests pass
        test_bracket = .true.

        ! setup settings object
        call setup_settings(settings, context)

        ! define procedure pointer
        obj_func_funptr => bracket_obj_func

        ! bracket the minimum of one-dimensional functions which drive the bracketing
        ! through its parabolic extrapolation steps, its golden ratio steps and its
        ! error conditions and determine whether the expected multiplier or an error
        ! is returned
        do i = 1, n_cases
            context%func_case = i
            call setup_error_logging(settings, context)
            n = bracket(obj_func_funptr, direction, 0.0_rp, 1.0_rp, settings, error)
            if (expect_error(i)) then
                if (error == 0) then
                    write(stderr, *) "test_bracket failed: Did not return error "// &
                        "for the function "//trim(case_names(i))//"."
                    test_bracket = .false.
                end if
                if (len_trim(context%log_message) == 0) then
                    write(stderr, *) "test_bracket failed: No error message "// &
                        "printed for the function "//trim(case_names(i))//"."
                    test_bracket = .false.
                end if
            else
                if (error /= 0) then
                    write(stderr, *) "test_bracket failed: Produced error for the "// &
                        "function "//trim(case_names(i))//"."
                    test_bracket = .false.
                end if
                if (abs(n - expected_n(i)) > tol * abs(expected_n(i))) then
                    write(stderr, *) "test_bracket failed: Returned multiplier not "// &
                        "correct for the function "//trim(case_names(i))//"."
                    test_bracket = .false.
                end if
            end if
        end do

        ! swap lower and upper bound and determine whether the same minimum is found
        context%func_case = 1
        n = bracket(obj_func_funptr, direction, 1.0_rp, 0.0_rp, settings, error)
        if (error /= 0) then
            write(stderr, *) "test_bracket failed: Produced error for swapped bounds."
            test_bracket = .false.
        end if
        if (abs(n - expected_n(1)) > tol * abs(expected_n(1))) then
            write(stderr, *) "test_bracket failed: Returned multiplier not correct "// &
                "for swapped bounds."
            test_bracket = .false.
        end if

        ! bracket with an objective function which fails and check that its error is
        ! reported with its origin
        obj_func_funptr => obj_func_failing
        n = bracket(obj_func_funptr, direction, 0.0_rp, 1.0_rp, settings, error)
        if (error /= error_obj_func + 1) then
            write(stderr, *) "test_bracket failed: Did not report the error of a "// &
                "failing objective function with its origin."
            test_bracket = .false.
        end if

    end function test_bracket

    logical(c_bool) function test_extend_symm_matrix() bind(C)
        !
        ! this function tests the subroutine for extending a symmetric matrix
        !
        use opentrustregion, only: extend_symm_matrix

        real(rp), allocatable :: matrix(:, :)
        real(rp) :: initial_matrix(n_param, n_param), vector(n_param + 1)

        ! assume tests pass
        test_extend_symm_matrix = .true.

        ! generate symmetric matrix and vector to be added
        initial_matrix = generate_random_symm_matrix(n_param)
        matrix = initial_matrix
        call random_number(vector)

        ! call routine and determine if dimensions of resulting matrix match and whether
        ! the initial matrix is extended by the vector as last row and column
        call extend_symm_matrix(matrix, vector)
        if (size(matrix, 1) /= n_param + 1 .or. size(matrix, 2) /= n_param + 1) then
            write(stderr, *) "test_extend_symm_matrix failed: Incorrect matrix "// &
                "dimensions after extending."
            test_extend_symm_matrix = .false.
            return
        end if
        if (any(abs(matrix(:n_param, :n_param) - initial_matrix) > tol) .or. &
            any(abs(matrix(:, n_param + 1) - vector) > tol) .or. &
            any(abs(matrix(n_param + 1, :) - vector) > tol)) then
            write(stderr, *) "test_extend_symm_matrix failed: Incorrect matrix "// &
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
        real(rp) :: initial_matrix(n_param, 2), new_col(n_param)

        ! assume tests pass
        test_add_column = .true.

        ! generate matrix and column to be added
        call random_number(initial_matrix)
        matrix = initial_matrix
        call random_number(new_col)

        ! call routine and determine if dimensions of resulting matrix match and whether
        ! the initial matrix is extended by the column
        call add_column(matrix, new_col)
        if (size(matrix, 1) /= n_param .or. size(matrix, 2) /= 3) then
            write(stderr, *) "test_add_column failed: Incorrect matrix dimensions "// &
                "after adding column."
            test_add_column = .false.
            return
        end if
        if (any(abs(matrix(:, :2) - initial_matrix) > tol) .or. &
            any(abs(matrix(:, 3) - new_col) > tol)) then
            write(stderr, *) "test_add_column failed: Incorrect matrix values "// &
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
        real(rp) :: matrix(n_param, n_param), eigval, eigvec(n_param), &
                    eigvals(n_param), eigvecs(n_param, n_param)
        integer(ip) :: error
        type(test_context_type), target :: context

        ! assume tests pass
        test_symm_mat_min_eig = .true.

        ! setup settings object
        call setup_settings(settings, context)

        ! generate symmetric matrix and determine its eigenvalues independently
        matrix = generate_random_symm_matrix(n_param)
        call ref_symm_mat_diag(matrix, eigvals, eigvecs)

        ! call routine and determine if lowest eigenvalue and corresponding normalized
        ! eigenvector are found
        call symm_mat_min_eig(matrix, eigval, eigvec, settings, error)
        if (error /= 0) then
            write(stderr, *) "test_symm_mat_min_eig failed: Produced error."
            test_symm_mat_min_eig = .false.
        end if
        if (abs(eigval - minval(eigvals)) > tol) then
            write(stderr, *) "test_symm_mat_min_eig failed: Incorrect minimum "// &
                "eigenvalue for matrix."
            test_symm_mat_min_eig = .false.
        end if
        if (norm2(matmul(matrix, eigvec) - eigval * eigvec) > tol .or. &
            abs(norm2(eigvec) - 1.0_rp) > tol) then
            write(stderr, *) "test_symm_mat_min_eig failed: Incorrect eigenvector "// &
                "corresponding to minimum eigenvalue for matrix."
            test_symm_mat_min_eig = .false.
        end if

    end function test_symm_mat_min_eig

    logical(c_bool) function test_symm_mat_diag() bind(C)
        !
        ! this function tests the function for determining the eigenvalues and
        ! eigenvectors for a symmetric matrix
        !
        use opentrustregion, only: solver_settings_type, symm_mat_diag

        type(solver_settings_type) :: settings
        real(rp) :: matrix(n_param, n_param), eigvals(n_param), &
                    eigvecs(n_param, n_param)
        integer(ip) :: error
        type(test_context_type), target :: context

        ! assume tests pass
        test_symm_mat_diag = .true.

        ! setup settings object
        call setup_settings(settings, context)

        ! generate symmetric matrix
        matrix = generate_random_symm_matrix(n_param)

        ! call routine and determine if eigenvalues and orthonormal eigenvectors are
        ! found
        call symm_mat_diag(matrix, eigvals, eigvecs, settings, error)
        if (error /= 0) then
            write(stderr, *) "test_symm_mat_diag failed: Produced error."
            test_symm_mat_diag = .false.
        end if
        if (norm2(matmul(matrix, eigvecs) - eigvecs * &
                  spread(eigvals, dim=1, ncopies=size(eigvecs, 1))) > tol) then
            write(stderr, *) "test_symm_mat_diag failed: Incorrect eigenvectors "// &
                "and eigenvalues for matrix."
            test_symm_mat_diag = .false.
        end if
        if (norm2(matmul(transpose(eigvecs), eigvecs) - identity_matrix(n_param)) > &
            tol) then
            write(stderr, *) "test_symm_mat_diag failed: Eigenvectors are not "// &
                "orthonormal."
            test_symm_mat_diag = .false.
        end if

    end function test_symm_mat_diag

    logical(c_bool) function test_init_rng() bind(C)
        !
        ! this function tests the initialization subroutine for the random number
        ! generator
        !
        use opentrustregion, only: init_rng

        integer(ip) :: seed1, seed2
        real(rp) :: rand_seq1(5), rand_seq2(5), rand_seq3(5)

        ! assume tests pass
        test_init_rng = .true.

        ! define seeds
        seed1 = 12345
        seed2 = 67890

        ! call rng with first seed
        call init_rng(seed1)
        call random_number(rand_seq1)

        ! call rng with first seed
        call init_rng(seed1)
        call random_number(rand_seq2)

        ! call rng with second seed
        call init_rng(seed2)
        call random_number(rand_seq3)

        ! check reproducibility
        if (any(abs(rand_seq1 - rand_seq2) > tol)) then
            write(stderr, *) "test_init_rng failed: RNG does not produce "// &
                "consistent sequences for the same seed."
            test_init_rng = .false.
        end if

        ! check variation
        if (all(abs(rand_seq1 - rand_seq3) < tol)) then
            write(stderr, *) "test_init_rng failed: RNG produces identical "// &
                "sequences for different seeds."
            test_init_rng = .false.
        end if

    end function test_init_rng

    logical(c_bool) function test_generate_trial_vectors() bind(C)
        !
        ! this function tests the function which generates trial vectors for the
        ! Davidson procedure
        !
        use opentrustregion, only: solver_settings_type, generate_trial_vectors, &
                                   error_project

        type(solver_settings_type) :: settings
        real(rp), allocatable :: red_space_basis(:, :)
        real(rp) :: grad(n_param), h_diag(n_param), grad_norm, neg_curv_vec(n_param), &
                    grad_small(2), h_diag_small(2)
        integer(ip) :: error
        type(test_context_type), target :: context

        ! assume tests pass
        test_generate_trial_vectors = .true.

        ! setup settings object
        call setup_settings(settings, context)
        settings%n_random_trial_vectors = 2

        ! generate gradient
        call random_number(grad)
        grad_norm = norm2(grad)

        ! generate positive Hessian diagonal elements
        call random_number(h_diag)
        h_diag = h_diag + 1.0_rp

        ! generate trial vectors and determine whether function returns the normalized
        ! gradient followed by the requested number of random trial vectors
        red_space_basis = &
            generate_trial_vectors(grad, grad_norm, h_diag, settings, error)
        if (error /= 0) then
            write(stderr, *) "test_generate_trial_vectors failed: Produced error "// &
                "for Hessian with only positive diagonal elements."
            test_generate_trial_vectors = .false.
        end if
        if (.not. allocated(red_space_basis)) then
            write(stderr, *) "test_generate_trial_vectors failed: Reduced space "// &
                "basis not allocated for Hessian with only positive diagonal elements."
            test_generate_trial_vectors = .false.
            return
        end if
        if (norm2(matmul(transpose(red_space_basis), red_space_basis) - &
                  identity_matrix(size(red_space_basis, 2, kind=ip))) > tol) then
            write(stderr, *) "test_generate_trial_vectors failed: Reduced space "// &
                "basis not orthonormal for Hessian with only positive diagonal "// &
                "elements."
            test_generate_trial_vectors = .false.
        end if
        if (size(red_space_basis, 2) /= 1 + settings%n_random_trial_vectors) then
            write(stderr, *) "test_generate_trial_vectors failed: Incorrect number "// &
                "of vectors for Hessian with only positive diagonal elements."
            test_generate_trial_vectors = .false.
        end if
        if (any(abs(red_space_basis(:, 1) - grad / grad_norm) > tol)) then
            write(stderr, *) "test_generate_trial_vectors failed: First vector is "// &
                "not the normalized gradient for Hessian with only positive "// &
                "diagonal elements."
            test_generate_trial_vectors = .false.
        end if

        ! deallocate reduced space basis
        deallocate(red_space_basis)

        ! make one Hessian diagonal element negative
        h_diag(2) = -h_diag(2)

        ! generate trial vectors and determine whether function returns the normalized
        ! gradient, the unit vector along the most negative Hessian diagonal element
        ! orthonormalized against it and the requested number of random trial vectors
        red_space_basis = &
            generate_trial_vectors(grad, grad_norm, h_diag, settings, error)
        if (error /= 0) then
            write(stderr, *) "test_generate_trial_vectors failed: Produced error "// &
                "for Hessian with negative diagonal elements."
            test_generate_trial_vectors = .false.
        end if
        if (.not. allocated(red_space_basis)) then
            write(stderr, *) "test_generate_trial_vectors failed: Reduced space "// &
                "basis not allocated for Hessian with negative diagonal elements."
            test_generate_trial_vectors = .false.
            return
        end if
        if (norm2(matmul(transpose(red_space_basis), red_space_basis) - &
                  identity_matrix(size(red_space_basis, 2, kind=ip))) > tol) then
            write(stderr, *) "test_generate_trial_vectors failed: Reduced space "// &
                "basis not orthonormal for Hessian with negative diagonal elements."
            test_generate_trial_vectors = .false.
        end if
        if (size(red_space_basis, 2) /= 2 + settings%n_random_trial_vectors) then
            write(stderr, *) "test_generate_trial_vectors failed: Incorrect number "// &
                "of vectors for Hessian with negative diagonal elements."
            test_generate_trial_vectors = .false.
        end if
        if (any(abs(red_space_basis(:, 1) - grad / grad_norm) > tol)) then
            write(stderr, *) "test_generate_trial_vectors failed: First vector is "// &
                "not the normalized gradient for Hessian with negative diagonal "// &
                "elements."
            test_generate_trial_vectors = .false.
        end if
        neg_curv_vec = 0.0_rp
        neg_curv_vec(2) = 1.0_rp
        neg_curv_vec = neg_curv_vec - &
                       dot_product(neg_curv_vec, grad) / grad_norm**2 * grad
        neg_curv_vec = neg_curv_vec / norm2(neg_curv_vec)
        if (any(abs(red_space_basis(:, 2) - neg_curv_vec) > tol)) then
            write(stderr, *) "test_generate_trial_vectors failed: Second vector is "// &
                "not the orthonormalized direction of the negative Hessian "// &
                "diagonal element."
            test_generate_trial_vectors = .false.
        end if

        ! project the direction of the negative Hessian diagonal element with a
        ! projection which removes the component along the symmetric combination of the
        ! first two unit vectors and determine whether the second vector is the
        ! orthonormalized projected direction
        settings%project => mock_project_out_pair
        red_space_basis = &
            generate_trial_vectors(grad, grad_norm, h_diag, settings, error)
        if (error /= 0) then
            write(stderr, *) "test_generate_trial_vectors failed: Produced error "// &
                "with projection."
            test_generate_trial_vectors = .false.
        end if
        neg_curv_vec = 0.0_rp
        neg_curv_vec(2) = 1.0_rp
        neg_curv_vec(1:2) = neg_curv_vec(1:2) - 0.5_rp * sum(neg_curv_vec(1:2))
        neg_curv_vec = neg_curv_vec - &
                       dot_product(neg_curv_vec, grad) / grad_norm**2 * grad
        neg_curv_vec = neg_curv_vec / norm2(neg_curv_vec)
        if (any(abs(red_space_basis(:, 2) - neg_curv_vec) > tol)) then
            write(stderr, *) "test_generate_trial_vectors failed: Second vector is "// &
                "not the orthonormalized projected direction of the negative "// &
                "Hessian diagonal element."
            test_generate_trial_vectors = .false.
        end if
        settings%project => null()

        ! choose the gradient along the direction of the negative Hessian diagonal
        ! element, which is then linearly dependent on the gradient, and determine
        ! whether only the normalized gradient is followed by the random trial vectors
        grad = 0.0_rp
        call random_number(grad(2))
        grad(2) = grad(2) + 1.0_rp
        grad_norm = norm2(grad)
        call setup_error_logging(settings, context)
        red_space_basis = &
            generate_trial_vectors(grad, grad_norm, h_diag, settings, error)
        if (error /= 0) then
            write(stderr, *) "test_generate_trial_vectors failed: Produced error "// &
                "for gradient along negative Hessian diagonal element."
            test_generate_trial_vectors = .false.
        end if
        if (len_trim(context%log_message) /= 0) then
            write(stderr, *) "test_generate_trial_vectors failed: Error message "// &
                "printed for gradient along negative Hessian diagonal element."
            test_generate_trial_vectors = .false.
        end if
        if (size(red_space_basis, 2) /= 1 + settings%n_random_trial_vectors) then
            write(stderr, *) "test_generate_trial_vectors failed: Incorrect number "// &
                "of vectors for gradient along negative Hessian diagonal element."
            test_generate_trial_vectors = .false.
        end if
        if (any(abs(red_space_basis(:, 1) - grad / grad_norm) > tol)) then
            write(stderr, *) "test_generate_trial_vectors failed: First vector is "// &
                "not the normalized gradient for gradient along negative Hessian "// &
                "diagonal element."
            test_generate_trial_vectors = .false.
        end if

        ! use only two parameters with a negative Hessian diagonal element and
        ! determine whether no direction of the negative Hessian diagonal element is
        ! added
        settings%n_random_trial_vectors = 1
        call random_number(grad_small)
        call random_number(h_diag_small)
        h_diag_small(2) = -h_diag_small(2) - 1.0_rp
        red_space_basis = generate_trial_vectors(grad_small, norm2(grad_small), &
                                                 h_diag_small, settings, error)
        if (error /= 0) then
            write(stderr, *) "test_generate_trial_vectors failed: Produced error "// &
                "for two parameters."
            test_generate_trial_vectors = .false.
        end if
        if (size(red_space_basis, 2) /= 1 + settings%n_random_trial_vectors) then
            write(stderr, *) "test_generate_trial_vectors failed: Incorrect number "// &
                "of vectors for two parameters."
            test_generate_trial_vectors = .false.
        end if

        ! deallocate reduced space basis
        deallocate(red_space_basis)

        ! check that the function result is allocated even when the projection function
        ! produces an error before the reduced space basis is constructed
        settings%project => mock_project_error
        red_space_basis = &
            generate_trial_vectors(grad, grad_norm, h_diag, settings, error)
        if (error /= error_project + 1) then
            write(stderr, *) "test_generate_trial_vectors failed: Error of failing "// &
                "projection function not reported with its origin."
            test_generate_trial_vectors = .false.
        end if
        if (.not. allocated(red_space_basis)) then
            write(stderr, *) "test_generate_trial_vectors failed: Reduced space "// &
                "basis not allocated for failing projection function."
            test_generate_trial_vectors = .false.
            return
        end if

        ! deallocate reduced space basis
        deallocate(red_space_basis)

        ! project the direction of the negative Hessian diagonal element, which is the
        ! second unit vector, with a projection which keeps only the first component so
        ! that it vanishes, and check that the error of its orthonormalization is
        ! returned and the function result is still allocated
        settings%project => mock_project_first_component
        red_space_basis = &
            generate_trial_vectors(grad, grad_norm, h_diag, settings, error)
        if (error == 0) then
            write(stderr, *) "test_generate_trial_vectors failed: Error not "// &
                "produced for vanishing projected direction."
            test_generate_trial_vectors = .false.
        end if
        if (.not. allocated(red_space_basis)) then
            write(stderr, *) "test_generate_trial_vectors failed: Reduced space "// &
                "basis not allocated for vanishing projected direction."
            test_generate_trial_vectors = .false.
            return
        end if

        ! deallocate reduced space basis
        deallocate(red_space_basis)

    end function test_generate_trial_vectors

    logical(c_bool) function test_generate_random_trial_vectors() bind(C)
        !
        ! this function tests the function which generates random trial vectors for the
        ! Davidson procedure
        !
        use opentrustregion, only: solver_settings_type, &
                                   generate_random_trial_vectors, error_project

        type(solver_settings_type) :: settings
        real(rp), allocatable :: red_space_basis(:, :)
        real(rp) :: first_vector(n_param)
        integer(ip) :: error, i, j
        type(test_context_type), target :: context

        ! assume tests pass
        test_generate_random_trial_vectors = .true.

        ! setup settings object
        call setup_settings(settings, context)
        settings%n_random_trial_vectors = 2

        ! allocate reduced space basis and set first normalized basis vector
        allocate(red_space_basis(n_param, 3))
        call random_number(red_space_basis(:, 1))
        red_space_basis(:, 1) = red_space_basis(:, 1) / norm2(red_space_basis(:, 1))
        first_vector = red_space_basis(:, 1)

        ! generate trial vectors and determine whether function returns orthonormal
        ! trial vectors and leaves the first basis vector unchanged
        call generate_random_trial_vectors(red_space_basis, settings, error)
        if (error /= 0) then
            write(stderr, *) "test_generate_random_trial_vectors failed: Produced "// &
                "error."
            test_generate_random_trial_vectors = .false.
        end if
        if (any(abs(red_space_basis(:, 1) - first_vector) > tol)) then
            write(stderr, *) "test_generate_random_trial_vectors failed: First "// &
                "basis vector changed."
            test_generate_random_trial_vectors = .false.
        end if
        if (any(abs(norm2(red_space_basis(:, 2:), dim=1) - 1.0_rp) > tol)) then
            write(stderr, *) "test_generate_random_trial_vectors failed: Generated "// &
                "vectors are not normalized."
            test_generate_random_trial_vectors = .false.
        end if
        do i = 1, size(red_space_basis, 2)
            do j = i + 1, size(red_space_basis, 2)
                if (abs(dot_product(red_space_basis(:, i), red_space_basis(:, j))) > &
                    tol) then
                    write(stderr, *) "test_generate_random_trial_vectors failed: "// &
                        "Generated vectors are not orthogonal."
                    test_generate_random_trial_vectors = .false.
                end if
            end do
        end do

        ! project the random trial vectors with a projection which removes the
        ! component along the symmetric combination of the first two unit vectors,
        ! starting from a first basis vector without this component, and determine
        ! whether the generated vectors do not have this component
        settings%project => mock_project_out_pair
        call random_number(red_space_basis(:, 1))
        red_space_basis(1:2, 1) = red_space_basis(1:2, 1) - &
                                  0.5_rp * sum(red_space_basis(1:2, 1))
        red_space_basis(:, 1) = red_space_basis(:, 1) / norm2(red_space_basis(:, 1))
        call generate_random_trial_vectors(red_space_basis, settings, error)
        if (error /= 0) then
            write(stderr, *) "test_generate_random_trial_vectors failed: Produced "// &
                "error with projection."
            test_generate_random_trial_vectors = .false.
        end if
        if (any(abs(sum(red_space_basis(1:2, 2:), dim=1)) > tol)) then
            write(stderr, *) "test_generate_random_trial_vectors failed: Generated "// &
                "vectors are not projected."
            test_generate_random_trial_vectors = .false.
        end if

        ! check that the error of a failing projection function is reported with its
        ! origin
        settings%project => mock_project_error
        call generate_random_trial_vectors(red_space_basis, settings, error)
        if (error /= error_project + 1) then
            write(stderr, *) "test_generate_random_trial_vectors failed: Did not "// &
                "report the error of a failing projection function with its origin."
            test_generate_random_trial_vectors = .false.
        end if

        ! project the random trial vectors onto the first unit vector, which is the
        ! first basis vector, so that every attempt produces a linearly dependent
        ! vector, and determine whether an error is returned and an error message is
        ! printed once the maximum number of attempts is reached
        settings%project => mock_project_first_component
        red_space_basis(:, 1) = 0.0_rp
        red_space_basis(1, 1) = 1.0_rp
        call setup_error_logging(settings, context)
        call generate_random_trial_vectors(red_space_basis, settings, error)
        if (error == 0) then
            write(stderr, *) "test_generate_random_trial_vectors failed: No error "// &
                "returned when no linearly independent vector can be generated."
            test_generate_random_trial_vectors = .false.
        end if
        if (len_trim(context%log_message) == 0) then
            write(stderr, *) "test_generate_random_trial_vectors failed: No error "// &
                "message printed when no linearly independent vector can be generated."
            test_generate_random_trial_vectors = .false.
        end if

        ! deallocate reduced space basis
        deallocate(red_space_basis)

    end function test_generate_random_trial_vectors

    logical(c_bool) function test_gram_schmidt() bind(C)
        !
        ! this function tests the Gram-Schmidt subroutine which orthonormalizes a
        ! vector to a given basis
        !
        use opentrustregion, only: solver_settings_type, gram_schmidt, &
                                   error_gram_schmidt_lin_dep

        type(solver_settings_type) :: settings
        real(rp) :: vector(n_param), lin_trans_vector(n_param), vector_small(2), &
                    space(n_param, 2), symm_matrix(n_param, n_param), &
                    lin_trans_space(n_param, 2), space_small(2, 2), eigvals(n_param), &
                    eigvecs(n_param, n_param)
        integer(ip) :: error
        type(test_context_type), target :: context

        ! assume tests pass
        test_gram_schmidt = .true.

        ! setup settings object
        call setup_settings(settings, context)

        ! generate vector to be orthogonalized and orthonormal space spanned by
        ! eigenvectors of a symmetric matrix
        call random_number(vector)
        call ref_symm_mat_diag(generate_random_symm_matrix(n_param), eigvals, eigvecs)
        space = eigvecs(:, :2)

        ! perform Gram-Schmidt orthogonalization and determine whether added vector is
        ! orthonormalized
        call gram_schmidt(vector, space, settings, error)
        if (error /= 0) then
            write(stderr, *) "test_gram_schmidt failed: Produced error without "// &
                "linear transformation."
            test_gram_schmidt = .false.
        end if
        if (abs(dot_product(vector, space(:, 1))) > tol .or. &
            abs(dot_product(vector, space(:, 2))) > tol) then
            write(stderr, *) "test_gram_schmidt failed: Added vector not "// &
                "orthogonal without linear transformation."
            test_gram_schmidt = .false.
        end if
        if (abs(norm2(vector) - 1.0_rp) > tol) then
            write(stderr, *) "test_gram_schmidt failed: Added vector not "// &
                "normalized without linear transformation."
            test_gram_schmidt = .false.
        end if

        ! generate vector to be orthogonalized
        call random_number(vector)

        ! generate symmetric linear transformation and corresponding vector and space
        symm_matrix = generate_random_symm_matrix(n_param)
        lin_trans_vector = matmul(symm_matrix, vector)
        lin_trans_space = matmul(symm_matrix, space)

        ! perform Gram-Schmidt orthogonalization and determine whether added vector is
        ! orthonormalized and linear transformation is correct
        call gram_schmidt(vector, space, settings, error, lin_trans_vector, &
                          lin_trans_space)
        if (error /= 0) then
            write(stderr, *) "test_gram_schmidt failed: Produced error with linear "// &
                "transformation."
            test_gram_schmidt = .false.
        end if
        if (abs(dot_product(vector, space(:, 1))) > tol .or. &
            abs(dot_product(vector, space(:, 2))) > tol) then
            write(stderr, *) "test_gram_schmidt failed: Added vector not "// &
                "orthogonal with linear transformation."
            test_gram_schmidt = .false.
        end if
        if (abs(norm2(vector) - 1.0_rp) > tol) then
            write(stderr, *) "test_gram_schmidt failed: Added vector not "// &
                "normalized with linear transformation."
            test_gram_schmidt = .false.
        end if
        if (norm2(lin_trans_vector - matmul(symm_matrix, vector)) > tol) then
            write(stderr, *) "test_gram_schmidt failed: Added linear "// &
                "transformation not correct."
            test_gram_schmidt = .false.
        end if

        ! define zero vector
        vector = 0.0_rp

        ! perform Gram-Schmidt orthogonalization and determine if function correctly
        ! returns an error and prints an error message
        call setup_error_logging(settings, context)
        call gram_schmidt(vector, space, settings, error)
        if (error == 0) then
            write(stderr, *) "test_gram_schmidt failed: No error returned during "// &
                "orthogonalization for zero vector."
            test_gram_schmidt = .false.
        end if
        if (len_trim(context%log_message) == 0) then
            write(stderr, *) "test_gram_schmidt failed: No error message printed "// &
                "for zero vector."
            test_gram_schmidt = .false.
        end if

        ! define linearly dependent vector
        vector = space(:, 1)

        ! perform Gram-Schmidt orthogonalization and determine if function correctly
        ! returns an error and prints an error message
        call setup_error_logging(settings, context)
        call gram_schmidt(vector, space, settings, error)
        if (error /= error_gram_schmidt_lin_dep) then
            write(stderr, *) "test_gram_schmidt failed: No error returned during "// &
                "orthogonalization for linearly dependent vector."
            test_gram_schmidt = .false.
        end if
        if (len_trim(context%log_message) == 0) then
            write(stderr, *) "test_gram_schmidt failed: No error message printed "// &
                "for linearly dependent vector."
            test_gram_schmidt = .false.
        end if

        ! define linearly dependent vector
        vector = space(:, 1)

        ! perform Gram-Schmidt orthogonalization that stays silent on error and
        ! determine that the error is still returned but no error message is printed
        call setup_error_logging(settings, context)
        call gram_schmidt(vector, space, settings, error, silent_on_error=.true.)
        if (error /= error_gram_schmidt_lin_dep) then
            write(stderr, *) "test_gram_schmidt failed: No error returned during "// &
                "orthogonalization for linearly dependent vector when asked to "// &
                "stay silent on error."
            test_gram_schmidt = .false.
        end if
        if (len_trim(context%log_message) /= 0) then
            write(stderr, *) "test_gram_schmidt failed: Error message printed "// &
                "despite being asked to stay silent on error."
            test_gram_schmidt = .false.
        end if

        ! define vector in space that is already complete
        call random_number(vector_small)
        space_small = identity_matrix(2_ip)

        ! perform Gram-Schmidt orthogonalization and determine if function correctly
        ! returns an error and prints an error message
        call setup_error_logging(settings, context)
        call gram_schmidt(vector_small, space_small, settings, error)
        if (error == 0) then
            write(stderr, *) "test_gram_schmidt failed: No error returned during "// &
                "orthogonalization when number of vectors is larger than dimension "// &
                "of vector space."
            test_gram_schmidt = .false.
        end if
        if (len_trim(context%log_message) == 0) then
            write(stderr, *) "test_gram_schmidt failed: No error message printed "// &
                "when number of vectors is larger than dimension of vector space."
            test_gram_schmidt = .false.
        end if

    end function test_gram_schmidt

    logical(c_bool) function test_init_solver_settings() bind(C)
        !
        ! this function tests the subroutine which initializes the solver settings
        !
        use opentrustregion, only: solver_settings_type, &
                                   default_settings => default_solver_settings
        use test_reference, only: operator(/=), callbacks_unset

        type, extends(solver_settings_type) :: extended_solver_settings_type
        end type

        type(solver_settings_type) :: settings
        type(extended_solver_settings_type) :: extended_settings
        type(test_context_type), target :: context
        integer(ip) :: error

        ! assume tests pass
        test_init_solver_settings = .true.

        ! set callback functions, host context and a non-default value which the
        ! initialization has to discard
        settings%precond => mock_precond
        settings%project => mock_project
        settings%conv_check => mock_conv_check
        settings%logger => logger
        settings%context => context
        settings%conv_tol = 1.0_rp

        ! initialize settings
        call settings%init(error)

        ! check for error
        if (error /= 0) then
            write(stderr, *) "test_init_solver_settings failed: Function raised error."
            test_init_solver_settings = .false.
        end if

        ! check function pointers and host context
        if (.not. callbacks_unset(settings)) then
            write(stderr, *) "test_init_solver_settings failed: Function pointers "// &
                "not discarded."
            test_init_solver_settings = .false.
        end if
        if (associated(settings%context)) then
            write(stderr, *) "test_init_solver_settings failed: Host context not "// &
                "discarded."
            test_init_solver_settings = .false.
        end if

        ! check settings
        if (settings /= default_settings) then
            write(stderr, *) "test_init_solver_settings failed: Settings not "// &
                "initialized correctly."
            test_init_solver_settings = .false.
        end if

        ! initialize settings of a type extending the solver settings, which cannot be
        ! set to the default values, and check that an error is returned
        call extended_settings%init(error)
        if (error == 0) then
            write(stderr, *) "test_init_solver_settings failed: No error returned "// &
                "for extended settings type."
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
        use test_reference, only: operator(/=), callbacks_unset

        type, extends(stability_settings_type) :: extended_stability_settings_type
        end type

        type(stability_settings_type) :: settings
        type(extended_stability_settings_type) :: extended_settings
        type(test_context_type), target :: context
        integer(ip) :: error

        ! assume tests pass
        test_init_stability_settings = .true.

        ! set callback functions, host context and a non-default value which the
        ! initialization has to discard
        settings%precond => mock_precond
        settings%project => mock_project
        settings%logger => logger
        settings%context => context
        settings%conv_tol = 1.0_rp

        ! initialize settings
        call settings%init(error)

        ! check for error
        if (error /= 0) then
            write(stderr, *) "test_init_stability_settings failed: Function raised "// &
                "error."
            test_init_stability_settings = .false.
        end if

        ! check function pointers and host context
        if (.not. callbacks_unset(settings)) then
            write(stderr, *) "test_init_stability_settings failed: Function "// &
                "pointers not discarded."
            test_init_stability_settings = .false.
        end if
        if (associated(settings%context)) then
            write(stderr, *) "test_init_stability_settings failed: Host context "// &
                "not discarded."
            test_init_stability_settings = .false.
        end if

        ! check settings
        if (settings /= default_settings) then
            write(stderr, *) "test_init_stability_settings failed: Settings not "// &
                "initialized correctly."
            test_init_stability_settings = .false.
        end if

        ! initialize settings of a type extending the stability settings, which cannot
        ! be set to the default values, and check that an error is returned
        call extended_settings%init(error)
        if (error == 0) then
            write(stderr, *) "test_init_stability_settings failed: No error "// &
                "returned for extended settings type."
            test_init_stability_settings = .false.
        end if

    end function test_init_stability_settings

    logical(c_bool) function test_level_shifted_diag_precond() bind(C)
        !
        ! this function tests the subroutine that constructs the level-shifted diagonal
        ! preconditioner
        !
        use opentrustregion, only: solver_settings_type, level_shifted_diag_precond, &
                                   error_precond, error_project

        real(rp) :: vector(3), mu, h_diag(3), precond_vector(3)
        type(solver_settings_type) :: settings
        integer(ip) :: error
        type(test_context_type), target :: context

        ! assume tests pass
        test_level_shifted_diag_precond = .true.

        ! setup settings object
        call setup_settings(settings, context)

        ! initialize quantities
        vector = [1.0_rp, 1.0_rp, 1.0_rp]
        mu = -2.0_rp
        h_diag = [-2.0_rp, 1.0_rp, 2.0_rp]

        ! call subroutine and check if results match
        call level_shifted_diag_precond(vector, mu, h_diag, precond_vector, settings, &
                                        error)
        if (error /= 0) then
            write(stderr, *) "test_level_shifted_diag_precond failed: Returned "// &
                "error for default preconditioner."
            test_level_shifted_diag_precond = .false.
        end if
        if (any(abs(precond_vector - [1e10_rp, 1.0_rp / 3, 0.25_rp]) > tol)) then
            write(stderr, *) "test_level_shifted_diag_precond failed: Returned "// &
                "preconditioned vector not correct for default preconditioner."
            test_level_shifted_diag_precond = .false.
        end if

        ! test custom projection function
        settings%project => mock_project

        ! call subroutine and check if results match
        call level_shifted_diag_precond(vector, mu, h_diag, precond_vector, settings, &
                                        error)
        if (error /= 0) then
            write(stderr, *) "test_level_shifted_diag_precond failed: Returned "// &
                "error for custom projection function."
            test_level_shifted_diag_precond = .false.
        end if
        if (any(abs(precond_vector - [2e10_rp, 2.0_rp / 3, 0.5_rp]) > tol)) then
            write(stderr, *) "test_level_shifted_diag_precond failed: Returned "// &
                "preconditioned vector not correct for custom projection function."
            test_level_shifted_diag_precond = .false.
        end if

        ! test custom preconditioner
        settings%precond => mock_precond

        ! call subroutine and check if results match
        call level_shifted_diag_precond(vector, mu, h_diag, precond_vector, settings, &
                                        error)
        if (error /= 0) then
            write(stderr, *) "test_level_shifted_diag_precond failed: Returned "// &
                "error for custom preconditioner."
            test_level_shifted_diag_precond = .false.
        end if
        if (any(abs(precond_vector - [-2.0_rp, -2.0_rp, -2.0_rp]) > tol)) then
            write(stderr, *) "test_level_shifted_diag_precond failed: Returned "// &
                "preconditioned vector not correct for custom preconditioner."
            test_level_shifted_diag_precond = .false.
        end if

        ! check that the error of a failing custom preconditioner is reported with its
        ! origin
        settings%precond => mock_precond_error
        call level_shifted_diag_precond(vector, mu, h_diag, precond_vector, settings, &
                                        error)
        if (error /= error_precond + 1) then
            write(stderr, *) "test_level_shifted_diag_precond failed: Did not "// &
                "report the error of a failing preconditioner with its origin."
            test_level_shifted_diag_precond = .false.
        end if

        ! check that the error of a failing custom projection function is reported
        ! with its origin
        settings%precond => null()
        settings%project => mock_project_error
        call level_shifted_diag_precond(vector, mu, h_diag, precond_vector, settings, &
                                        error)
        if (error /= error_project + 1) then
            write(stderr, *) "test_level_shifted_diag_precond failed: Did not "// &
                "report the error of a failing projection function with its origin."
            test_level_shifted_diag_precond = .false.
        end if

    end function test_level_shifted_diag_precond

    logical(c_bool) function test_abs_diag_precond() bind(C)
        !
        ! this function tests the subroutine that constructs the absolute diagonal
        ! preconditioner
        !
        use opentrustregion, only: solver_settings_type, abs_diag_precond, &
                                   error_precond, error_project

        real(rp) :: vector(3), h_diag(3), precond_vector(3)
        type(solver_settings_type) :: settings
        integer(ip) :: error
        type(test_context_type), target :: context

        ! assume tests pass
        test_abs_diag_precond = .true.

        ! setup settings object
        call setup_settings(settings, context)

        ! initialize quantities
        vector = [1.0_rp, 1.0_rp, 1.0_rp]
        h_diag = [-2.0_rp, 0.0_rp, 2.0_rp]

        ! call subroutine and check if results match
        call abs_diag_precond(vector, h_diag, precond_vector, settings, error)
        if (error /= 0) then
            write(stderr, *) "test_abs_diag_precond failed: Returned error for "// &
                "default preconditioner."
            test_abs_diag_precond = .false.
        end if
        if (any(abs(precond_vector - [0.5_rp, 1e10_rp, 0.5_rp]) > tol)) then
            write(stderr, *) "test_abs_diag_precond failed: Returned "// &
                "preconditioned vector not correct for default preconditioner."
            test_abs_diag_precond = .false.
        end if

        ! test custom projection function
        settings%project => mock_project

        ! call subroutine and check if results match
        call abs_diag_precond(vector, h_diag, precond_vector, settings, error)
        if (error /= 0) then
            write(stderr, *) "test_abs_diag_precond failed: Returned error for "// &
                "custom projection function."
            test_abs_diag_precond = .false.
        end if
        if (any(abs(precond_vector - [1.0_rp, 2e10_rp, 1.0_rp]) > tol)) then
            write(stderr, *) "test_abs_diag_precond failed: Returned "// &
                "preconditioned vector not correct for custom projection function."
            test_abs_diag_precond = .false.
        end if

        ! test custom preconditioner
        settings%precond => mock_precond

        ! call subroutine and check if results match
        call abs_diag_precond(vector, h_diag, precond_vector, settings, error)
        if (error /= 0) then
            write(stderr, *) "test_abs_diag_precond failed: Returned error for "// &
                "custom preconditioner."
            test_abs_diag_precond = .false.
        end if
        if (any(abs(precond_vector - [0.0_rp, 0.0_rp, 0.0_rp]) > tol)) then
            write(stderr, *) "test_abs_diag_precond failed: Returned "// &
                "preconditioned vector not correct for custom preconditioner."
            test_abs_diag_precond = .false.
        end if

        ! check that the error of a failing custom preconditioner is reported with its
        ! origin
        settings%precond => mock_precond_error
        call abs_diag_precond(vector, h_diag, precond_vector, settings, error)
        if (error /= error_precond + 1) then
            write(stderr, *) "test_abs_diag_precond failed: Did not report the "// &
                "error of a failing preconditioner with its origin."
            test_abs_diag_precond = .false.
        end if

        ! check that the error of a failing custom projection function is reported
        ! with its origin
        settings%precond => null()
        settings%project => mock_project_error
        call abs_diag_precond(vector, h_diag, precond_vector, settings, error)
        if (error /= error_project + 1) then
            write(stderr, *) "test_abs_diag_precond failed: Did not report the "// &
                "error of a failing projection function with its origin."
            test_abs_diag_precond = .false.
        end if

    end function test_abs_diag_precond

    logical(c_bool) function test_orthogonal_projection() bind(C)
        !
        ! this function tests the orthogonal projection function which removes a
        ! certain direction from a vector
        !
        use opentrustregion, only: orthogonal_projection

        real(rp), dimension(n_param) :: vector, direction, complement

        ! assume tests pass
        test_orthogonal_projection = .true.

        ! generate vector and direction to be projected out, the latter needs to be
        ! normalized
        call random_number(vector)
        call random_number(direction)
        direction = direction / norm2(direction)

        ! perform orthogonal projection and determine whether the result contains the
        ! direction and whether only a component along the direction was removed
        complement = orthogonal_projection(vector, direction)
        if (abs(dot_product(complement, direction)) > tol) then
            write(stderr, *) "test_orthogonal_projection failed: Vector contains "// &
                "component from direction to be projected out."
            test_orthogonal_projection = .false.
        end if
        if (norm2(vector - complement - &
                  dot_product(vector - complement, direction) * direction) > tol) then
            write(stderr, *) "test_orthogonal_projection failed: Removed component "// &
                "not along direction to be projected out."
            test_orthogonal_projection = .false.
        end if

    end function test_orthogonal_projection

    logical(c_bool) function test_jacobi_davidson_correction() bind(C)
        !
        ! this function tests the Jacobi-Davidson correction subroutine
        !
        use opentrustregion, only: solver_settings_type, hess_x_type, &
                                   jacobi_davidson_correction, error_hess_x

        type(solver_settings_type) :: settings
        procedure(hess_x_type), pointer :: hess_x_funptr
        real(rp), dimension(n_param) :: vector, solution, corr_vector, hess_vector, &
                                        proj_vector, expected_corr_vector
        real(rp) :: eigval
        integer(ip) :: error
        type(hartmann6d_context_type), target :: context

        ! assume tests pass
        test_jacobi_davidson_correction = .true.

        ! setup settings object
        call setup_settings(settings, context)

        ! generate Hessian, trial vector, eigenvalue and normalized solution to be
        ! projected out
        context%hess = generate_random_symm_matrix(n_param)
        call random_number(vector)
        call random_number(solution)
        solution = solution / norm2(solution)
        call random_number(eigval)

        ! define Hessian linear transformation
        hess_x_funptr => hess_x_fun

        ! calculate Jacobi-Davidson correction and determine whether the Hessian linear
        ! transformation of the projected vector and the projected level-shifted
        ! Hessian linear transformation of the projected vector are returned and whether
        ! the reported number of Hessian linear transformations agrees with the calls
        context%n_hess_x_calls = 0
        call jacobi_davidson_correction(hess_x_funptr, vector, solution, eigval, &
                                        corr_vector, hess_vector, settings, error)
        if (error /= 0) then
            write(stderr, *) "test_jacobi_davidson_correction failed: Returned error."
            test_jacobi_davidson_correction = .false.
        end if
        proj_vector = vector - dot_product(vector, solution) * solution
        if (any(abs(hess_vector - matmul(context%hess, proj_vector)) > tol)) then
            write(stderr, *) "test_jacobi_davidson_correction failed: Returned "// &
                "Hessian linear transformation wrong."
            test_jacobi_davidson_correction = .false.
        end if
        expected_corr_vector = matmul(context%hess, proj_vector) - eigval * proj_vector
        expected_corr_vector = expected_corr_vector - &
                               dot_product(expected_corr_vector, solution) * solution
        if (any(abs(corr_vector - expected_corr_vector) > tol)) then
            write(stderr, *) "test_jacobi_davidson_correction failed: Returned "// &
                "correction vector wrong."
            test_jacobi_davidson_correction = .false.
        end if
        test_jacobi_davidson_correction = &
            test_jacobi_davidson_correction .and. logical(check_call_counts( &
                context, "jacobi_davidson_correction", "for valid input", &
                n_hess_x=settings%n_hess_x), kind=c_bool)

        ! calculate Jacobi-Davidson correction with a Hessian linear transformation
        ! which fails and check that the failing call is still counted
        hess_x_funptr => hess_x_fun_failing
        settings%n_hess_x = 0
        context%n_hess_x_calls = 0
        call jacobi_davidson_correction(hess_x_funptr, vector, solution, eigval, &
                                        corr_vector, hess_vector, settings, error)
        if (error /= error_hess_x + 1) then
            write(stderr, *) "test_jacobi_davidson_correction failed: Did not "// &
                "report the error of a failing Hessian linear transformation with "// &
                "its origin."
            test_jacobi_davidson_correction = .false.
        end if
        test_jacobi_davidson_correction = &
            test_jacobi_davidson_correction .and. logical(check_call_counts( &
                context, "jacobi_davidson_correction", "for failing Hessian linear "// &
                "transformation", n_hess_x=settings%n_hess_x), kind=c_bool)

    end function test_jacobi_davidson_correction

    logical(c_bool) function test_minres() bind(C)
        !
        ! this function tests the minimum residual method subroutine
        !
        use opentrustregion, only: solver_settings_type, hess_x_type, minres

        type(solver_settings_type) :: settings
        procedure(hess_x_type), pointer :: hess_x_funptr
        real(rp), dimension(n_param) :: rhs, solution, vector, hess_vector, &
                                        corr_vector, guess
        real(rp) :: mu
        real(rp), parameter :: minres_tol = 1e-14_rp
        integer(ip) :: error
        type(hartmann6d_context_type), target :: context

        ! assume tests pass
        test_minres = .true.

        ! setup settings object
        call setup_settings(settings, context)

        ! generate normalized solution to be projected out and symmetric Hessian
        call random_number(solution)
        solution = solution / norm2(solution)
        context%hess = generate_random_symm_matrix(n_param)

        ! define Hessian linear transformation
        hess_x_funptr => hess_x_fun

        ! define Rayleigh quotient
        mu = dot_product(solution, matmul(context%hess, solution))

        ! define right hand side based on residual, this ensures rhs is orthogonal to
        ! solution if mu describes the Rayleigh quotient and solution is normalized
        rhs = matmul(context%hess, solution) - mu * solution

        ! run minimum residual method, check if Jacobi-Davidson correction equation is
        ! solved and whether Hessian linear transformation is correct, if rhs and
        ! solution are orthogonal (as in Jacobi-Davidson), the final vector will be
        ! orthogonal to the solution vector and consequently the Hessian linear
        ! transformation of the projected vector is equivalent to the Hessian linear
        ! transformation of the vector itself
        call minres(-rhs, hess_x_funptr, solution, mu, minres_tol, vector, &
                    hess_vector, settings, error)
        if (error /= 0) then
            write(stderr, *) "test_minres failed: Returned error for correction "// &
                "equation."
            test_minres = .false.
        end if
        corr_vector = vector - dot_product(vector, solution) * solution
        corr_vector = matmul(context%hess, corr_vector) - mu * corr_vector
        corr_vector = corr_vector - dot_product(corr_vector, solution) * solution
        if (sum(abs(corr_vector + rhs)) > tol) then
            write(stderr, *) "test_minres failed: Returned solution does not solve "// &
                "Jacobi-Davidson correction equation."
            test_minres = .false.
        end if
        if (sum(abs(hess_vector + dot_product(vector, solution) * matmul( &
            context%hess, solution) - matmul(context%hess, vector))) > tol) then
            write(stderr, *) "test_minres failed: Returned Hessian linear "// &
                "transformation wrong."
            test_minres = .false.
        end if

        ! run minimum residual method for vanishing right hand side
        rhs = 0.0_rp
        call minres(-rhs, hess_x_funptr, solution, mu, minres_tol, vector, &
                    hess_vector, settings, error)
        if (error /= 0) then
            write(stderr, *) "test_minres failed: Returned error for vanishing rhs."
            test_minres = .false.
        end if
        if (sum(abs(vector)) > tol) then
            write(stderr, *) "test_minres failed: Returned solution is not zero "// &
                "for a vanishing rhs."
            test_minres = .false.
        end if
        if (sum(abs(hess_vector)) > tol) then
            write(stderr, *) "test_minres failed: Returned Hessian linear "// &
                "transformation is not zero for a vanishing rhs."
            test_minres = .false.
        end if

        ! run minimum residual method from an initial guess which only partially solves
        ! the Jacobi-Davidson correction equation and check if the equation is solved
        rhs = matmul(context%hess, solution) - mu * solution
        call minres(-rhs, hess_x_funptr, solution, mu, minres_tol, guess, hess_vector, &
                    settings, error)
        guess = 0.5_rp * guess
        call minres(-rhs, hess_x_funptr, solution, mu, minres_tol, vector, &
                    hess_vector, settings, error, guess=guess)
        if (error /= 0) then
            write(stderr, *) "test_minres failed: Returned error for initial guess."
            test_minres = .false.
        end if
        corr_vector = vector - dot_product(vector, solution) * solution
        corr_vector = matmul(context%hess, corr_vector) - mu * corr_vector
        corr_vector = corr_vector - dot_product(corr_vector, solution) * solution
        if (sum(abs(corr_vector + rhs)) > tol) then
            write(stderr, *) "test_minres failed: Returned solution does not solve "// &
                "Jacobi-Davidson correction equation for initial guess."
            test_minres = .false.
        end if

        ! run minimum residual method from an initial guess which solves the
        ! Jacobi-Davidson correction equation exactly and check that the guess is
        ! returned with only its own Hessian linear transformation
        call random_number(guess)
        guess = guess - dot_product(guess, solution) * solution
        corr_vector = matmul(context%hess, guess) - mu * guess
        rhs = -(corr_vector - dot_product(corr_vector, solution) * solution)
        context%n_hess_x_calls = 0
        call minres(-rhs, hess_x_funptr, solution, mu, minres_tol, vector, &
                    hess_vector, settings, error, guess=guess)
        if (error /= 0) then
            write(stderr, *) "test_minres failed: Returned error for initial guess "// &
                "which solves the equation."
            test_minres = .false.
        end if
        if (context%n_hess_x_calls /= 1 .or. any(abs(vector - guess) > tol)) then
            write(stderr, *) "test_minres failed: Initial guess which solves the "// &
                "equation not returned directly."
            test_minres = .false.
        end if

        ! allow only a single iteration and check that an error is returned and an
        ! error message is printed when the iteration limit is reached
        call setup_error_logging(settings, context)
        call minres(-rhs, hess_x_funptr, solution, mu, minres_tol, vector, &
                    hess_vector, settings, error, max_iter=1_ip)
        if (error == 0) then
            write(stderr, *) "test_minres failed: No error returned when iteration "// &
                "limit is reached."
            test_minres = .false.
        end if
        if (len_trim(context%log_message) == 0) then
            write(stderr, *) "test_minres failed: No error message printed when "// &
                "iteration limit is reached."
            test_minres = .false.
        end if

        ! run minimum residual method from an initial guess for a vanishing right hand
        ! side and check that the solution and its Hessian linear transformation vanish
        rhs = 0.0_rp
        call minres(-rhs, hess_x_funptr, solution, mu, minres_tol, vector, &
                    hess_vector, settings, error, guess=guess)
        if (error /= 0) then
            write(stderr, *) "test_minres failed: Returned error for initial guess "// &
                "and vanishing rhs."
            test_minres = .false.
        end if
        if (sum(abs(vector)) > tol .or. sum(abs(hess_vector)) > tol) then
            write(stderr, *) "test_minres failed: Returned solution or its Hessian "// &
                "linear transformation is not zero for initial guess and vanishing rhs."
            test_minres = .false.
        end if

    end function test_minres

    logical(c_bool) function test_add_trial_vector() bind(C)
        !
        ! this function tests the subroutine which adds a new trial vector to the
        ! reduced space basis
        !
        use opentrustregion, only: solver_settings_type, hess_x_type, &
                                   add_trial_vector, error_gram_schmidt_lin_dep, &
                                   error_hess_x

        type(solver_settings_type) :: settings
        procedure(hess_x_type), pointer :: hess_x_funptr
        real(rp), allocatable :: red_space_basis(:, :), h_basis(:, :)
        real(rp), dimension(n_param) :: residual, h_diag, solution, expected_vec, &
                                        eigvals
        real(rp) :: initial_basis(n_param, 2), eigvecs(n_param, n_param), &
                    proj_shifted_hess(n_param, n_param), level_shift, eigval
        real(rp), parameter :: minres_tol = 1e-14_rp
        integer(ip) :: error, i
        type(hartmann6d_context_type), target :: context

        ! assume tests pass
        test_add_trial_vector = .true.

        ! setup settings object
        call setup_settings(settings, context)

        ! generate Hessian, orthonormal basis spanned by eigenvectors of a symmetric
        ! matrix, residual, Hessian diagonal and level shift below the Hessian diagonal
        context%hess = generate_random_symm_matrix(n_param)
        call ref_symm_mat_diag(generate_random_symm_matrix(n_param), eigvals, eigvecs)
        initial_basis = eigvecs(:, :2)
        call random_number(residual)
        call random_number(h_diag)
        h_diag = h_diag + 1.0_rp
        call random_number(level_shift)
        level_shift = 0.5_rp * level_shift
        hess_x_funptr => hess_x_fun

        ! add trial vector from preconditioned residual and determine whether the
        ! reduced space basis is extended by the orthonormalized residual
        ! preconditioned with the level-shifted Hessian diagonal and its Hessian linear
        ! transformation
        red_space_basis = initial_basis
        h_basis = matmul(context%hess, initial_basis)
        context%n_hess_x_calls = 0
        settings%n_hess_x = 0
        call add_trial_vector(residual, level_shift, h_diag, .false., solution, &
                              0.0_rp, minres_tol, hess_x_funptr, red_space_basis, &
                              h_basis, settings, error)
        if (error /= 0) then
            write(stderr, *) "test_add_trial_vector failed: Produced error for "// &
                "Davidson."
            test_add_trial_vector = .false.
        end if
        if (size(red_space_basis, 2) /= 3 .or. size(h_basis, 2) /= 3) then
            write(stderr, *) "test_add_trial_vector failed: Reduced space basis "// &
                "not extended by one vector for Davidson."
            test_add_trial_vector = .false.
            return
        end if
        expected_vec = residual / (h_diag - level_shift)
        expected_vec = expected_vec - &
                       matmul(initial_basis, matmul(expected_vec, initial_basis))
        expected_vec = expected_vec / norm2(expected_vec)
        if (any(abs(red_space_basis(:, :2) - initial_basis) > tol) .or. &
            any(abs(red_space_basis(:, 3) - expected_vec) > tol)) then
            write(stderr, *) "test_add_trial_vector failed: New trial vector is "// &
                "not the orthonormalized preconditioned residual for Davidson."
            test_add_trial_vector = .false.
        end if
        if (any( &
            abs(h_basis(:, 3) - matmul(context%hess, red_space_basis(:, 3))) > tol)) &
            then
            write(stderr, *) "test_add_trial_vector failed: Hessian linear "// &
                "transformation of new trial vector wrong for Davidson."
            test_add_trial_vector = .false.
        end if
        test_add_trial_vector = test_add_trial_vector .and. logical(check_call_counts( &
            context, "add_trial_vector", "for Davidson", n_hess_x=settings%n_hess_x), &
            kind=c_bool)

        ! add trial vector from Jacobi-Davidson correction equations for the first
        ! basis vector and its Rayleigh quotient and determine whether the reduced
        ! space basis is extended by the orthonormalized solution of the correction
        ! equations, which solves the linear system with the projected shifted Hessian,
        ! made nonsingular along the solution, independently
        solution = initial_basis(:, 1)
        eigval = dot_product(solution, matmul(context%hess, solution))
        residual = matmul(context%hess, solution) - eigval * solution
        proj_shifted_hess = context%hess - eigval * identity_matrix(n_param)
        proj_shifted_hess = proj_shifted_hess - spread(solution, 2, n_param) * &
                            spread(matmul(solution, proj_shifted_hess), 1, n_param)
        proj_shifted_hess = proj_shifted_hess - spread(matmul( &
            proj_shifted_hess, solution), 2, n_param) * spread(solution, 1, n_param)
        proj_shifted_hess = proj_shifted_hess + &
                            spread(solution, 2, n_param) * spread(solution, 1, n_param)
        call ref_symm_mat_diag(proj_shifted_hess, eigvals, eigvecs)
        expected_vec = -matmul(eigvecs, matmul(residual, eigvecs) / eigvals)
        expected_vec = expected_vec - &
                       matmul(initial_basis, matmul(expected_vec, initial_basis))
        expected_vec = expected_vec / norm2(expected_vec)
        red_space_basis = initial_basis
        h_basis = matmul(context%hess, initial_basis)
        context%n_hess_x_calls = 0
        settings%n_hess_x = 0
        call add_trial_vector(residual, level_shift, h_diag, .true., solution, eigval, &
                              minres_tol, hess_x_funptr, red_space_basis, h_basis, &
                              settings, error)
        if (error /= 0) then
            write(stderr, *) "test_add_trial_vector failed: Produced error for "// &
                "Jacobi-Davidson."
            test_add_trial_vector = .false.
        end if
        if (size(red_space_basis, 2) /= 3 .or. size(h_basis, 2) /= 3) then
            write(stderr, *) "test_add_trial_vector failed: Reduced space basis "// &
                "not extended by one vector for Jacobi-Davidson."
            test_add_trial_vector = .false.
            return
        end if
        if (any(abs(red_space_basis(:, 3) - expected_vec) > tol)) then
            write(stderr, *) "test_add_trial_vector failed: New trial vector is "// &
                "not the orthonormalized solution of the correction equations for "// &
                "Jacobi-Davidson."
            test_add_trial_vector = .false.
        end if
        if (any( &
            abs(h_basis(:, 3) - matmul(context%hess, red_space_basis(:, 3))) > tol)) &
            then
            write(stderr, *) "test_add_trial_vector failed: Hessian linear "// &
                "transformation of new trial vector wrong for Jacobi-Davidson."
            test_add_trial_vector = .false.
        end if
        test_add_trial_vector = test_add_trial_vector .and. logical( &
            check_call_counts(context, "add_trial_vector", "for Jacobi-Davidson", &
                              n_hess_x=settings%n_hess_x), kind=c_bool)

        ! repeat with a slightly asymmetric Hessian linear transformation, whose new
        ! linear transformation then no longer respects Hessian symmetry with respect
        ! to the existing basis, and determine whether it is recalculated
        red_space_basis = initial_basis
        h_basis = matmul(context%hess, initial_basis)
        hess_x_funptr => hess_x_fun_asymmetric
        call add_trial_vector(residual, level_shift, h_diag, .true., solution, eigval, &
                              minres_tol, hess_x_funptr, red_space_basis, h_basis, &
                              settings, error)
        if (error /= 0) then
            write(stderr, *) "test_add_trial_vector failed: Produced error for "// &
                "asymmetric Hessian linear transformation."
            test_add_trial_vector = .false.
        end if
        if (size(h_basis, 2) /= 3) then
            write(stderr, *) "test_add_trial_vector failed: Reduced space basis "// &
                "not extended by one vector for asymmetric Hessian linear "// &
                "transformation."
            test_add_trial_vector = .false.
            return
        end if
        call hess_x_funptr(red_space_basis(:, 3), expected_vec, error, settings%context)
        if (any(abs(h_basis(:, 3) - expected_vec) > tol)) then
            write(stderr, *) "test_add_trial_vector failed: Linear transformation "// &
                "of new trial vector not recalculated for asymmetric Hessian "// &
                "linear transformation."
            test_add_trial_vector = .false.
        end if
        hess_x_funptr => hess_x_fun

        ! precondition a residual along the first basis vector with a constant Hessian
        ! diagonal so that the new vector is linearly dependent on the basis and
        ! determine whether this is returned without extending the basis or printing an
        ! error message
        red_space_basis = initial_basis
        h_basis = matmul(context%hess, initial_basis)
        residual = initial_basis(:, 1)
        h_diag = 1.0_rp
        call setup_error_logging(settings, context)
        call add_trial_vector(residual, 0.0_rp, h_diag, .false., solution, 0.0_rp, &
                              minres_tol, hess_x_funptr, red_space_basis, h_basis, &
                              settings, error)
        if (error /= error_gram_schmidt_lin_dep) then
            write(stderr, *) "test_add_trial_vector failed: Linear dependence not "// &
                "returned for Davidson."
            test_add_trial_vector = .false.
        end if
        if (size(red_space_basis, 2) /= 2 .or. size(h_basis, 2) /= 2) then
            write(stderr, *) "test_add_trial_vector failed: Reduced space basis "// &
                "extended despite linear dependence for Davidson."
            test_add_trial_vector = .false.
        end if
        if (len_trim(context%log_message) /= 0) then
            write(stderr, *) "test_add_trial_vector failed: Error message printed "// &
                "for linear dependence for Davidson."
            test_add_trial_vector = .false.
        end if

        ! choose a Hessian which couples the first unit vector, the solution, only to
        ! the second unit vector, which is an eigenvector of the Hessian projected onto
        ! the complement of the solution, so that the solution of the correction
        ! equations lies along the second unit vector, and determine whether its linear
        ! dependence on the basis of the first two unit vectors is returned without
        ! extending the basis or printing an error message
        context%hess = 0.0_rp
        context%hess(1, 1) = 1.0_rp
        context%hess(1, 2) = 0.5_rp
        context%hess(2, 1) = 0.5_rp
        do i = 2, n_param
            context%hess(i, i) = real(i, kind=rp)
        end do
        red_space_basis = identity_matrix(n_param)
        red_space_basis = red_space_basis(:, :2)
        h_basis = matmul(context%hess, red_space_basis)
        solution = red_space_basis(:, 1)
        eigval = context%hess(1, 1)
        residual = matmul(context%hess, solution) - eigval * solution
        call setup_error_logging(settings, context)
        call add_trial_vector(residual, 0.0_rp, h_diag, .true., solution, eigval, &
                              minres_tol, hess_x_funptr, red_space_basis, h_basis, &
                              settings, error)
        if (error /= error_gram_schmidt_lin_dep) then
            write(stderr, *) "test_add_trial_vector failed: Linear dependence not "// &
                "returned for Jacobi-Davidson."
            test_add_trial_vector = .false.
        end if
        if (size(red_space_basis, 2) /= 2 .or. size(h_basis, 2) /= 2) then
            write(stderr, *) "test_add_trial_vector failed: Reduced space basis "// &
                "extended despite linear dependence for Jacobi-Davidson."
            test_add_trial_vector = .false.
        end if
        if (len_trim(context%log_message) /= 0) then
            write(stderr, *) "test_add_trial_vector failed: Error message printed "// &
                "for linear dependence for Jacobi-Davidson."
            test_add_trial_vector = .false.
        end if

        ! add trial vector with a Hessian linear transformation which fails and check
        ! that its error is returned with its origin and the failing call is counted
        red_space_basis = initial_basis
        h_basis = matmul(context%hess, initial_basis)
        call random_number(residual)
        hess_x_funptr => hess_x_fun_failing
        context%n_hess_x_calls = 0
        settings%n_hess_x = 0
        call add_trial_vector(residual, 0.0_rp, h_diag, .false., solution, 0.0_rp, &
                              minres_tol, hess_x_funptr, red_space_basis, h_basis, &
                              settings, error)
        if (error /= error_hess_x + 1) then
            write(stderr, *) "test_add_trial_vector failed: Did not return error "// &
                "of failing Hessian linear transformation."
            test_add_trial_vector = .false.
        end if
        test_add_trial_vector = test_add_trial_vector .and. logical(check_call_counts( &
            context, "add_trial_vector", "for failing Hessian linear transformation", &
            n_hess_x=settings%n_hess_x), kind=c_bool)

    end function test_add_trial_vector

    logical(c_bool) function test_print_results() bind(C)
        !
        ! this function tests the subroutine that prints the result table
        !
        use opentrustregion, only: solver_settings_type

        type(solver_settings_type) :: settings
        type(test_context_type), target :: context

        ! assume tests pass
        test_print_results = .true.

        ! setup settings object
        call setup_settings(settings, context)

        ! print row of results table without optional arguments and check if row is
        ! correct
        call settings%print_results(1_ip, 2.0_rp, 3.0_rp)
        if (context%log_message /= "        1   |     2.00000000000000E+00   "// &
            "|   3.00E+00   |      -      |        -   |        -     |      -   ") then
            write(stderr, *) "test_print_results failed: Printed row without "// &
                "optional arguments not correct."
            test_print_results = .false.
        end if

        ! reset log message
        context%log_message = ""

        ! print row of results table with optional arguments
        call settings%print_results(1_ip, 2.0_rp, 3.0_rp, 4.0_rp, 5_ip, 6_ip, 7.0_rp, &
                                    8.0_rp)
        if (context%log_message /= "        1   |     2.00000000000000E+00   "// &
            "|   3.00E+00   |   4.00E+00  |   0 |   5  |    7.00E+00  |  8.00E+00") then
            write(stderr, *) "test_print_results failed: Printed row with optional "// &
                "arguments not correct."
            test_print_results = .false.
        end if

        ! reset log message
        context%log_message = ""

        ! print row of results table with micro iterations but without Jacobi-Davidson
        ! micro iterations
        call settings%print_results(1_ip, 2.0_rp, 3.0_rp, level_shift=4.0_rp, &
                                    n_micro=5_ip, trust_radius=7.0_rp, &
                                    kappa_norm=8.0_rp)
        if (context%log_message /= "        1   |     2.00000000000000E+00   "// &
            "|   3.00E+00   |   4.00E+00  |         5  |    7.00E+00  |  8.00E+00") then
            write(stderr, *) "test_print_results failed: Printed row without "// &
                "Jacobi-Davidson micro iterations not correct."
            test_print_results = .false.
        end if

    end function test_print_results

    logical(c_bool) function test_print_message() bind(C)
        !
        ! this function tests the logging subroutine
        !
        use opentrustregion, only: solver_settings_type, print_message, &
                                   verbosity_debug, verbosity_warning

        type(solver_settings_type) :: settings
        type(test_context_type), target :: context

        ! assume tests pass
        test_print_message = .true.

        ! setup settings object
        call setup_settings(settings, context)

        ! check if logging is correctly performed according to verbosity level when
        ! logger is provided
        call print_message(settings, "This is a test message.", verbosity_warning)
        if (trim(context%log_message) /= " This is a test message.") then
            write(stderr, *) "test_print_message failed: Log message is not "// &
                "printed correctly even though it should be according to verbosity "// &
                "level."
            test_print_message = .false.
        end if
        context%log_message = ""
        call print_message(settings, "This is another test message.", verbosity_debug)
        if (len_trim(context%log_message) /= 0) then
            write(stderr, *) "test_print_message failed: Log message is printed "// &
                "even though it should not be according to verbosity level."
            test_print_message = .false.
        end if

        ! check that a blank message does not abort
        context%log_message = ""
        call print_message(settings, "   ", verbosity_warning)
        if (trim(context%log_message) /= "") then
            write(stderr, *) "test_print_message failed: Blank message is not "// &
                "logged correctly."
            test_print_message = .false.
        end if

        ! check that a message longer than a line is split into several lines which
        ! together reproduce the message, each line is logged with a leading space
        context%log_message = ""
        call print_message(settings, repeat("This is a long test message. ", 8), &
                           verbosity_warning)
        if (context%log_message /= &
            " "//trim(repeat("This is a long test message. ", 8))) then
            write(stderr, *) "test_print_message failed: Long message is not "// &
                "logged correctly."
            test_print_message = .false.
        end if

    end function test_print_message

    logical(c_bool) function test_split_string_by_space() bind(C)
        !
        ! this function tests the subroutine which splits strings after a space if they
        ! exceed a given maximum length
        !
        use opentrustregion, only: split_string_by_space

        character(len=23), parameter :: message = "This is a test message."
        character(len=:), allocatable :: substrings(:)

        ! assume tests pass
        test_split_string_by_space = .true.

        ! check if strings are split correctly on spaces
        call split_string_by_space(message, 8_ip, substrings)
        if (size(substrings) == 3) then
            if (trim(substrings(1)) /= "This is" .or. &
                trim(substrings(2)) /= "a test" .or. &
                trim(substrings(3)) /= "message.") then
                write(stderr, *) "test_split_string_by_space failed: Split strings "// &
                    "incorrect when splitting on spaces."
                test_split_string_by_space = .false.
            end if
        else
            write(stderr, *) "test_split_string_by_space failed: Number of "// &
                "substrings incorrect when splitting on spaces."
            test_split_string_by_space = .false.
        end if

        ! check if strings are split correctly if splitting on spaces is not possible
        call split_string_by_space(message, 5_ip, substrings)
        if (size(substrings) == 5) then
            if (trim(substrings(1)) /= "This" .or. trim(substrings(2)) /= "is a" .or. &
                trim(substrings(3)) /= "test" .or. trim(substrings(4)) /= "messa" .or. &
                trim(substrings(5)) /= "ge.") then
                write(stderr, *) "test_split_string_by_space failed: Split strings "// &
                    "incorrect when splitting on spaces is not possible."
                test_split_string_by_space = .false.
            end if
        else
            write(stderr, *) "test_split_string_by_space failed: Number of "// &
                "substrings incorrect when splitting on spaces is not possible."
            test_split_string_by_space = .false.
        end if

        ! check that a blank string produces an allocated array of zero substrings
        call split_string_by_space("   ", 8_ip, substrings)
        if (.not. allocated(substrings)) then
            write(stderr, *) "test_split_string_by_space failed: Substrings not "// &
                "allocated for blank string."
            test_split_string_by_space = .false.
        else if (size(substrings) /= 0) then
            write(stderr, *) "test_split_string_by_space failed: Number of "// &
                "substrings incorrect for blank string."
            test_split_string_by_space = .false.
        end if

    end function test_split_string_by_space

    logical(c_bool) function test_accept_trust_region_step() bind(C)
        !
        ! this function tests the subroutine which determines whether to accept a
        ! trust-region step
        !
        use opentrustregion, only: &
            solver_settings_type, accept_trust_region_step, trust_radius_shrink_ratio, &
            trust_radius_too_small_warning_msg, trust_radius_expand_ratio, &
            trust_radius_shrink_factor, trust_radius_expand_factor

        logical :: accept_step, max_precision_reached
        real(rp) :: solution(3), trust_radius
        type(solver_settings_type) :: settings
        type(test_context_type), target :: context

        ! assume tests pass
        test_accept_trust_region_step = .true.

        ! setup settings object
        call setup_settings(settings, context)

        solution = [0.3_rp, 0.3_rp, 0.3_rp]

        ! check if step is rejected and trust radius is correctly reduced if micro
        ! iterations have not converged
        trust_radius = 1.0_rp
        accept_step = accept_trust_region_step(solution, 1.0_rp, .false., settings, &
                                               trust_radius, max_precision_reached)
        if (accept_step .or. max_precision_reached .or. &
            abs(trust_radius - trust_radius_shrink_factor) > tol) then
            write(stderr, *) "test_accept_trust_region_step failed: Step accepted, "// &
                "maximum precision reported as reached or trust radius not "// &
                "correctly reduced when micro iterations have not converged."
            test_accept_trust_region_step = .false.
        end if

        ! check if step is rejected and trust radius is correctly reduced if ratio is
        ! negative
        trust_radius = 1.0_rp
        accept_step = accept_trust_region_step(solution, -1.0_rp, .true., settings, &
                                               trust_radius, max_precision_reached)
        if (accept_step .or. max_precision_reached .or. &
            abs(trust_radius - trust_radius_shrink_factor) > tol) then
            write(stderr, *) "test_accept_trust_region_step failed: Step accepted, "// &
                "maximum precision reported as reached or trust radius not "// &
                "correctly reduced when ratio is negative."
            test_accept_trust_region_step = .false.
        end if

        ! check if step is rejected and trust radius is correctly reduced if individual
        ! rotations are too large
        trust_radius = 1.0_rp
        solution(1) = 1.0_rp
        accept_step = accept_trust_region_step(solution, 1.0_rp, .true., settings, &
                                               trust_radius, max_precision_reached)
        if (accept_step .or. max_precision_reached .or. &
            abs(trust_radius - trust_radius_shrink_factor) > tol) then
            write(stderr, *) "test_accept_trust_region_step failed: Step accepted, "// &
                "maximum precision reported as reached or trust radius not "// &
                "correctly reduced when individual rotations are too large."
            test_accept_trust_region_step = .false.
        end if
        solution(1) = 0.3_rp

        ! check if step is accepted and trust radius is correctly reduced if ratio is
        ! too small
        trust_radius = 1.0_rp
        accept_step = accept_trust_region_step( &
            solution, 0.9_rp * trust_radius_shrink_ratio, .true., settings, &
            trust_radius, max_precision_reached)
        if (.not. accept_step .or. max_precision_reached .or. &
            abs(trust_radius - trust_radius_shrink_factor) > tol) then
            write(stderr, *) "test_accept_trust_region_step failed: Step not "// &
                "accepted, maximum precision reported as reached or trust radius "// &
                "not correctly reduced when ratio is too small."
            test_accept_trust_region_step = .false.
        end if

        ! check if step is accepted and trust radius is correctly reduced if ratio is
        ! ok
        trust_radius = 1.0_rp
        accept_step = accept_trust_region_step( &
            solution, &
            0.5_rp * (trust_radius_shrink_ratio + trust_radius_expand_ratio), .true., &
            settings, trust_radius, max_precision_reached)
        if (.not. accept_step .or. max_precision_reached .or. &
            abs(trust_radius - 1.0_rp) > tol) then
            write(stderr, *) "test_accept_trust_region_step failed: Step not "// &
                "accepted, maximum precision reported as reached or trust radius "// &
                "changed when ratio is acceptable."
            test_accept_trust_region_step = .false.
        end if

        ! check if step is accepted and trust radius is correctly expanded if ratio is
        ! too large
        trust_radius = 1.0_rp
        accept_step = accept_trust_region_step( &
            solution, 1.1_rp * trust_radius_expand_ratio, .true., settings, &
            trust_radius, max_precision_reached)
        if (.not. accept_step .or. max_precision_reached .or. &
            abs(trust_radius - trust_radius_expand_factor) > tol) then
            write(stderr, *) "test_accept_trust_region_step failed: Step not "// &
                "accepted, maximum precision reported as reached or trust radius "// &
                "not correctly expanded when ratio is too large."
            test_accept_trust_region_step = .false.
        end if

        ! check if step is rejected and maximum precision is reported as reached if the
        ! reduced trust radius falls below numerical zero
        trust_radius = 1e-14_rp
        context%log_message = ""
        accept_step = accept_trust_region_step(solution, -1.0_rp, .true., settings, &
                                               trust_radius, max_precision_reached)
        if (accept_step .or. .not. max_precision_reached) then
            write(stderr, *) "test_accept_trust_region_step failed: Step accepted "// &
                "or maximum precision not reached when trust radius becomes too small."
            test_accept_trust_region_step = .false.
        end if
        if (adjustl(context%log_message) /= trust_radius_too_small_warning_msg) then
            write(stderr, *) "test_accept_trust_region_step failed: Warning not "// &
                "printed when trust radius becomes too small."
            test_accept_trust_region_step = .false.
        end if

    end function test_accept_trust_region_step

    logical(c_bool) function test_solver_sanity_check() bind(C)
        !
        ! this function tests the subroutine which performs a sanity check for the
        ! solver
        !
        use opentrustregion, only: &
            solver_settings_type, solver_sanity_check, project_warning_msg, &
            random_trial_vector_warning_msg, subsystem_solver_options

        type(solver_settings_type) :: settings
        real(rp) :: grad(3)
        integer(ip) :: error, i
        type(test_context_type), target :: context

        ! assume tests pass
        test_solver_sanity_check = .true.

        ! setup settings object
        call setup_settings(settings, context)

        ! check if error is incorrectly thrown for finite and non-negative number of
        ! parameters
        settings%n_random_trial_vectors = 0
        call solver_sanity_check(settings, 3_ip, grad, error)
        if (error /= 0) then
            write(stderr, *) "test_solver_sanity_check failed: Error thrown for "// &
                "non-negative and non-vanishing number of parameters."
            test_solver_sanity_check = .false.
        end if

        ! check if error is correctly thrown for vanishing number of parameters
        call setup_error_logging(settings, context)
        call solver_sanity_check(settings, 0_ip, grad, error)
        if (error == 0) then
            write(stderr, *) "test_solver_sanity_check failed: Error not thrown "// &
                "for vanishing number of parameters."
            test_solver_sanity_check = .false.
        end if
        if (len_trim(context%log_message) == 0) then
            write(stderr, *) "test_solver_sanity_check failed: No error message "// &
                "printed for vanishing number of parameters."
            test_solver_sanity_check = .false.
        end if

        ! check if error is correctly thrown for negative number of parameters
        call setup_error_logging(settings, context)
        call solver_sanity_check(settings, -1_ip, grad, error)
        if (error == 0) then
            write(stderr, *) "test_solver_sanity_check failed: Error not thrown "// &
                "for negative number of parameters."
            test_solver_sanity_check = .false.
        end if
        if (len_trim(context%log_message) == 0) then
            write(stderr, *) "test_solver_sanity_check failed: No error message "// &
                "printed for negative number of parameters."
            test_solver_sanity_check = .false.
        end if

        ! check if error is correctly thrown and an error message is printed for
        ! vanishing number of microiterations
        settings%n_micro = 0
        call setup_error_logging(settings, context)
        call solver_sanity_check(settings, 3_ip, grad, error)
        if (error == 0) then
            write(stderr, *) "test_solver_sanity_check failed: Error not thrown "// &
                "for vanishing number of microiterations."
            test_solver_sanity_check = .false.
        end if
        if (len_trim(context%log_message) == 0) then
            write(stderr, *) "test_solver_sanity_check failed: No error message "// &
                "printed for vanishing number of microiterations."
            test_solver_sanity_check = .false.
        end if
        settings%n_micro = 50

        ! check if number of random trial vectors is reduced correctly with a warning
        ! for the Davidson solvers
        settings%n_random_trial_vectors = 3
        settings%verbose = 3
        context%log_message = ""
        call solver_sanity_check(settings, 3_ip, grad, error)
        if (settings%n_random_trial_vectors /= 1) then
            write(stderr, *) "test_solver_sanity_check failed: Number of random "// &
                "trial vectors not correctly set."
            test_solver_sanity_check = .false.
        end if
        if (adjustl(context%log_message) /= &
            random_trial_vector_warning_msg//" Setting to 1.") then
            write(stderr, *) "test_solver_sanity_check failed: Warning message not "// &
                "correctly printed when number of random trial vectors is reduced."
            test_solver_sanity_check = .false.
        end if

        ! check that the number of random trial vectors is not reduced for the
        ! truncated conjugate gradient solver, which does not use them
        settings%subsystem_solver = "tcg"
        settings%n_random_trial_vectors = 3
        call solver_sanity_check(settings, 3_ip, grad, error)
        if (settings%n_random_trial_vectors /= 3) then
            write(stderr, *) "test_solver_sanity_check failed: Number of random "// &
                "trial vectors reduced for truncated conjugate gradient solver."
            test_solver_sanity_check = .false.
        end if
        settings%subsystem_solver = "davidson"
        settings%n_random_trial_vectors = 1

        ! check if gradient size is treated correctly
        call solver_sanity_check(settings, 3_ip, grad, error)
        if (error /= 0) then
            write(stderr, *) "test_solver_sanity_check failed: Error thrown for "// &
                "gradient size."
            test_solver_sanity_check = .false.
        end if
        call setup_error_logging(settings, context)
        call solver_sanity_check(settings, 4_ip, grad, error)
        if (error == 0) then
            write(stderr, *) "test_solver_sanity_check failed: Error not thrown "// &
                "for incorrect gradient size."
            test_solver_sanity_check = .false.
        end if
        if (len_trim(context%log_message) == 0) then
            write(stderr, *) "test_solver_sanity_check failed: No error message "// &
                "printed for incorrect gradient size."
            test_solver_sanity_check = .false.
        end if

        ! check if every subsystem solver option is accepted
        do i = 1, size(subsystem_solver_options)
            settings%subsystem_solver = subsystem_solver_options(i)
            call solver_sanity_check(settings, 3_ip, grad, error)
            if (error /= 0) then
                write(stderr, *) "test_solver_sanity_check failed: Error thrown "// &
                    "for "//trim(subsystem_solver_options(i))//" subsystem solver."
                test_solver_sanity_check = .false.
            end if
        end do
        settings%subsystem_solver = "Jacobi-Davidson"
        call solver_sanity_check(settings, 3_ip, grad, error)
        if (error /= 0 .or. settings%subsystem_solver /= "jacobi-davidson") then
            write(stderr, *) "test_solver_sanity_check failed: Error thrown or "// &
                "subsystem solver not converted to lowercase for mixed-case "// &
                "subsystem solver."
            test_solver_sanity_check = .false.
        end if
        settings%subsystem_solver = "unknown"
        call setup_error_logging(settings, context)
        call solver_sanity_check(settings, 3_ip, grad, error)
        if (error == 0) then
            write(stderr, *) "test_solver_sanity_check failed: Error not thrown "// &
                "for unknown subsystem solver."
            test_solver_sanity_check = .false.
        end if
        if (len_trim(context%log_message) == 0) then
            write(stderr, *) "test_solver_sanity_check failed: No error message "// &
                "printed for unknown subsystem solver."
            test_solver_sanity_check = .false.
        end if

        ! print warnings again and reset log message
        settings%verbose = 3
        context%log_message = ""

        ! check if a warning is printed when a custom projection function is set
        settings%subsystem_solver = "davidson"
        settings%project => mock_project
        call solver_sanity_check(settings, 3_ip, grad, error)
        if (adjustl(context%log_message) /= project_warning_msg) then
            write(stderr, *) "test_solver_sanity_check failed: Warning message not "// &
                "correctly printed when custom projecting function is set."
            test_solver_sanity_check = .false.
        end if

    end function test_solver_sanity_check

    logical(c_bool) function test_stability_sanity_check() bind(C)
        !
        ! this function tests the subroutine which performs a sanity check for the
        ! stability check
        !
        use opentrustregion, only: stability_settings_type, stability_sanity_check, &
                                   project_warning_msg, &
                                   random_trial_vector_warning_msg, diag_solver_options

        type(stability_settings_type) :: settings
        integer(ip) :: error, i
        type(test_context_type), target :: context

        ! assume tests pass
        test_stability_sanity_check = .true.

        ! setup settings object
        call setup_settings(settings, context)

        ! check if number of random trial vectors is reduced correctly with a warning
        settings%n_random_trial_vectors = 3
        call stability_sanity_check(settings, 3_ip, error)
        if (settings%n_random_trial_vectors /= 1) then
            write(stderr, *) "test_stability_sanity_check failed: Number of random "// &
                "trial vectors not correctly set."
            test_stability_sanity_check = .false.
        end if
        if (adjustl(context%log_message) /= &
            random_trial_vector_warning_msg//" Setting to 1.") then
            write(stderr, *) "test_stability_sanity_check failed: Warning message "// &
                "not correctly printed when number of random trial vectors is reduced."
            test_stability_sanity_check = .false.
        end if

        ! check if every diagonalization solver option is accepted
        do i = 1, size(diag_solver_options)
            settings%diag_solver = diag_solver_options(i)
            call stability_sanity_check(settings, 3_ip, error)
            if (error /= 0) then
                write(stderr, *) "test_stability_sanity_check failed: Error thrown "// &
                    "for "//trim(diag_solver_options(i))//" diagonalization solver."
                test_stability_sanity_check = .false.
            end if
        end do
        settings%diag_solver = "Jacobi-Davidson"
        call stability_sanity_check(settings, 3_ip, error)
        if (error /= 0 .or. settings%diag_solver /= "jacobi-davidson") then
            write(stderr, *) "test_stability_sanity_check failed: Error thrown or "// &
                "diagonalization solver not converted to lowercase for mixed-case "// &
                "diagonalization solver."
            test_stability_sanity_check = .false.
        end if
        settings%diag_solver = "unknown"
        call setup_error_logging(settings, context)
        call stability_sanity_check(settings, 3_ip, error)
        if (error == 0) then
            write(stderr, *) "test_stability_sanity_check failed: Error not thrown "// &
                "for unknown diagonalization solver."
            test_stability_sanity_check = .false.
        end if
        if (len_trim(context%log_message) == 0) then
            write(stderr, *) "test_stability_sanity_check failed: No error message "// &
                "printed for unknown diagonalization solver."
            test_stability_sanity_check = .false.
        end if

        ! print warnings again and reset log message
        settings%verbose = 3
        context%log_message = ""

        ! check if a warning is printed when a custom projection function is set
        settings%diag_solver = "davidson"
        settings%project => mock_project
        call stability_sanity_check(settings, 3_ip, error)
        if (adjustl(context%log_message) /= project_warning_msg) then
            write(stderr, *) "test_stability_sanity_check failed: Warning message "// &
                "not correctly printed when custom projecting function is set."
            test_stability_sanity_check = .false.
        end if

    end function test_stability_sanity_check

    logical(c_bool) function test_level_shifted_davidson() bind(C)
        !
        ! this function tests the level-shifted Davidson subroutine
        !
        use opentrustregion, only: obj_func_type, hess_x_type, solver_settings_type, &
                                   level_shifted_davidson, error_hess_x, &
                                   error_obj_func, error_precond

        real(rp) :: func, grad_norm, trust_radius, mu, ratio, input_trust_radius, &
                    overflow_residual_tol
        real(rp), dimension(n_param) :: grad, h_diag, solution
        integer(ip) :: i, i_case, imicro, imicro_jacobi_davidson, error
        real(rp) :: eigvals(n_param), rotation(n_param, n_param)
        integer(ip), parameter :: n_stagnation = 20
        real(rp) :: h_diag_stagnation(n_stagnation), solution_stagnation(n_stagnation)
        character(len=*), parameter :: &
            case_names(2) = [character(len=8) :: "diagonal", "rotated"], &
            stagnation_case_names(2) = &
                [character(len=15) :: "Davidson", "Jacobi-Davidson"]
        procedure(obj_func_type), pointer :: obj_func_funptr
        procedure(hess_x_type), pointer :: hess_x_funptr
        type(solver_settings_type) :: settings
        logical :: jacobi_davidson_started, max_precision_reached
        type(hartmann6d_context_type), target :: context
        type(quadratic_context_type), target :: quadratic_context

        ! assume tests pass
        test_level_shifted_davidson = .true.

        ! setup settings object
        call setup_settings(settings, context)

        ! initialize variables
        trust_radius = 0.4_rp
        obj_func_funptr => obj_func
        hess_x_funptr => hess_x_fun

        ! start in quadratic region near minimum
        context%vars = near_minimum
        func = hartmann6d_func(context%vars)
        call hartmann6d_gradient(context%vars, grad)
        grad_norm = norm2(grad)
        context%hess = hartmann6d_hessian(context%vars)
        h_diag = [(context%hess(i, i), i=1, size(h_diag))]

        ! run level-shifted Davidson, check if error has occured, whether the level
        ! shift vanishes, whether the solution stays within trust region and
        ! describes the Newton step and whether the reported number of Hessian linear
        ! transformations agrees with the calls
        settings%n_hess_x = 0
        context%n_hess_x_calls = 0
        call level_shifted_davidson( &
            func, grad, grad_norm, h_diag, n_param, obj_func_funptr, hess_x_funptr, &
            settings, trust_radius, solution, mu, imicro, imicro_jacobi_davidson, &
            jacobi_davidson_started, max_precision_reached, error)
        if (error /= 0) then
            write(stderr, *) "test_level_shifted_davidson failed: Produced error "// &
                "near minimum."
            test_level_shifted_davidson = .false.
        end if
        if (abs(mu) > tol) then
            write(stderr, *) "test_level_shifted_davidson failed: Level shift is "// &
                "not zero near minimum."
            test_level_shifted_davidson = .false.
        end if
        if (norm2(grad + matmul(context%hess, solution)) > &
            settings%local_red_factor * grad_norm) then
            write(stderr, *) "test_level_shifted_davidson failed: Solution does "// &
                "not describe Newton step near minimum."
            test_level_shifted_davidson = .false.
        end if
        ratio = (hartmann6d_func(context%vars + solution) - func) / &
                dot_product(solution, grad + 0.5_rp * matmul(context%hess, solution))
        if (norm2(solution) > ref_step_trust_radius(trust_radius, ratio) + tol) then
            write(stderr, *) "test_level_shifted_davidson failed: Solution does "// &
                "not stay within trust region near minimum."
            test_level_shifted_davidson = .false.
        end if
        test_level_shifted_davidson = test_level_shifted_davidson .and. logical( &
            check_call_counts(context, "level_shifted_davidson", "near minimum", &
                              n_hess_x=settings%n_hess_x), kind=c_bool)

        ! start near saddle point
        context%vars = near_saddle_point
        func = hartmann6d_func(context%vars)
        call hartmann6d_gradient(context%vars, grad)
        grad_norm = norm2(grad)
        context%hess = hartmann6d_hessian(context%vars)
        h_diag = [(context%hess(i, i), i=1, size(h_diag))]
        trust_radius = 0.4_rp

        ! run level-shifted Davidson, check if error has occured, whether the level
        ! shift is negative and whether the solution lies at the trust region boundary
        ! and describes a level-shifted Newton step
        call level_shifted_davidson( &
            func, grad, grad_norm, h_diag, n_param, obj_func_funptr, hess_x_funptr, &
            settings, trust_radius, solution, mu, imicro, imicro_jacobi_davidson, &
            jacobi_davidson_started, max_precision_reached, error)
        if (error /= 0) then
            write(stderr, *) "test_level_shifted_davidson failed: Produced error "// &
                "near saddle point."
            test_level_shifted_davidson = .false.
        end if
        if (mu >= 0.0_rp) then
            write(stderr, *) "test_level_shifted_davidson failed: Level shift is "// &
                "not negative near saddle point."
            test_level_shifted_davidson = .false.
        end if
        if (norm2(grad + matmul(context%hess, solution) - mu * solution) > &
            settings%global_red_factor * grad_norm) then
            write(stderr, *) "test_level_shifted_davidson failed: Solution does "// &
                "not describe level-shifted Newton step near saddle point."
            test_level_shifted_davidson = .false.
        end if
        ratio = (hartmann6d_func(context%vars + solution) - func) / &
                dot_product(solution, grad + 0.5_rp * matmul(context%hess, solution))
        if (abs(norm2(solution) - ref_step_trust_radius(trust_radius, ratio)) > tol) &
            then
            write(stderr, *) "test_level_shifted_davidson failed: Solution does "// &
                "not lie at trust region boundary near saddle point."
            test_level_shifted_davidson = .false.
        end if

        ! test Jacobi-Davidson near saddle point, switching from the first micro
        ! iteration since Davidson would otherwise converge before switching
        settings%subsystem_solver = "jacobi-davidson"
        settings%jacobi_davidson_start = 0
        trust_radius = 0.4_rp

        ! run level-shifted Jacobi-Davidson, check if error has occured, whether it
        ! switched to Jacobi-Davidson, whether the level shift is negative and whether
        ! the solution lies at the trust region boundary and describes a level-shifted
        ! Newton step
        call level_shifted_davidson( &
            func, grad, grad_norm, h_diag, n_param, obj_func_funptr, hess_x_funptr, &
            settings, trust_radius, solution, mu, imicro, imicro_jacobi_davidson, &
            jacobi_davidson_started, max_precision_reached, error)
        if (error /= 0) then
            write(stderr, *) "test_level_shifted_davidson failed: Produced error "// &
                "near saddle point with Jacobi-Davidson solver."
            test_level_shifted_davidson = .false.
        end if
        if (.not. jacobi_davidson_started) then
            write(stderr, *) "test_level_shifted_davidson failed: Did not switch "// &
                "to Jacobi-Davidson."
            test_level_shifted_davidson = .false.
        end if
        if (mu >= 0.0_rp) then
            write(stderr, *) "test_level_shifted_davidson failed: Level shift is "// &
                "not negative near saddle point with Jacobi-Davidson solver."
            test_level_shifted_davidson = .false.
        end if
        if (norm2(grad + matmul(context%hess, solution) - mu * solution) > &
            settings%global_red_factor * grad_norm) then
            write(stderr, *) "test_level_shifted_davidson failed: Solution does "// &
                "not describe level-shifted Newton step near saddle point with "// &
                "Jacobi-Davidson solver."
            test_level_shifted_davidson = .false.
        end if
        ratio = (hartmann6d_func(context%vars + solution) - func) / &
                dot_product(solution, grad + 0.5_rp * matmul(context%hess, solution))
        if (abs(norm2(solution) - ref_step_trust_radius(trust_radius, ratio)) > tol) &
            then
            write(stderr, *) "test_level_shifted_davidson failed: Solution does "// &
                "not lie at trust region boundary near saddle point with "// &
                "Jacobi-Davidson solver."
            test_level_shifted_davidson = .false.
        end if

        ! run level-shifted Davidson with a Hessian linear transformation which fails
        ! and check that the failing call is still counted
        hess_x_funptr => hess_x_fun_failing
        trust_radius = 0.4_rp
        settings%n_hess_x = 0
        context%n_hess_x_calls = 0
        call level_shifted_davidson( &
            func, grad, grad_norm, h_diag, n_param, obj_func_funptr, hess_x_funptr, &
            settings, trust_radius, solution, mu, imicro, imicro_jacobi_davidson, &
            jacobi_davidson_started, max_precision_reached, error)
        if (error /= error_hess_x + 1) then
            write(stderr, *) "test_level_shifted_davidson failed: Did not report "// &
                "the error of a failing Hessian linear transformation with its origin."
            test_level_shifted_davidson = .false.
        end if
        test_level_shifted_davidson = test_level_shifted_davidson .and. logical( &
            check_call_counts(context, "level_shifted_davidson", &
                              "for failing Hessian linear transformation", &
                              n_hess_x=settings%n_hess_x), kind=c_bool)

        ! run level-shifted Davidson with an objective function which fails and check
        ! that its error is reported with its origin
        call setup_settings(settings, context)
        hess_x_funptr => hess_x_fun
        obj_func_funptr => obj_func_failing
        trust_radius = 0.4_rp
        call level_shifted_davidson( &
            func, grad, grad_norm, h_diag, n_param, obj_func_funptr, hess_x_funptr, &
            settings, trust_radius, solution, mu, imicro, imicro_jacobi_davidson, &
            jacobi_davidson_started, max_precision_reached, error)
        if (error /= error_obj_func + 1) then
            write(stderr, *) "test_level_shifted_davidson failed: Did not report "// &
                "the error of a failing objective function with its origin."
            test_level_shifted_davidson = .false.
        end if
        obj_func_funptr => obj_func

        ! run level-shifted Davidson with a preconditioner which fails when the first
        ! trial vector is added and check that its error is returned
        settings%precond => mock_precond_error
        trust_radius = 0.4_rp
        call level_shifted_davidson( &
            func, grad, grad_norm, h_diag, n_param, obj_func_funptr, hess_x_funptr, &
            settings, trust_radius, solution, mu, imicro, imicro_jacobi_davidson, &
            jacobi_davidson_started, max_precision_reached, error)
        if (error /= error_precond + 1) then
            write(stderr, *) "test_level_shifted_davidson failed: Did not return "// &
                "the error of adding a trial vector."
            test_level_shifted_davidson = .false.
        end if
        settings%precond => null()

        ! let the residual stagnate with a quadratic model whose gradient only couples
        ! to the second unit vector and a preconditioner which returns the following
        ! unit vectors, which can never reduce the residual, and check that the
        ! Davidson method stops after ten micro iterations without sufficient residual
        ! reduction, so that its step is rejected and the trust radius reduced, while
        ! the Jacobi-Davidson method switches to the correction equations then, the
        ! first micro iteration sets the initial residual and the switch happens in the
        ! eleventh, long before the reduced space reaches the dimension of the full
        ! parameter space or the Jacobi-Davidson method would otherwise be started
        do i_case = 1, 2
            call setup_settings(settings, quadratic_context)
            settings%n_random_trial_vectors = 0
            settings%precond => mock_precond_next_unit_vector
            if (i_case == 2) settings%subsystem_solver = "jacobi-davidson"
            quadratic_context%hess = identity_matrix(n_stagnation)
            quadratic_context%hess(1:2, 1:2) = reshape([2.0_rp, 1.0_rp, &
                                                        1.0_rp, 3.0_rp], [2, 2])
            quadratic_context%grad = spread(0.0_rp, 1, n_stagnation)
            quadratic_context%grad(1) = 1.0_rp
            quadratic_context%next_unit_vector = 2
            h_diag_stagnation = [(quadratic_context%hess(i, i), i=1, n_stagnation)]
            obj_func_funptr => quadratic_obj_func
            hess_x_funptr => quadratic_hess_x
            trust_radius = 1.0_rp
            call level_shifted_davidson( &
                0.0_rp, quadratic_context%grad, 1.0_rp, h_diag_stagnation, &
                n_stagnation, obj_func_funptr, hess_x_funptr, settings, trust_radius, &
                solution_stagnation, mu, imicro, imicro_jacobi_davidson, &
                jacobi_davidson_started, max_precision_reached, error)
            if (error /= 0) then
                write(stderr, *) "test_level_shifted_davidson failed: Produced "// &
                    "error when residual stagnates with "// &
                    trim(stagnation_case_names(i_case))//" solver."
                test_level_shifted_davidson = .false.
            end if
            if (i_case == 1 .and. trust_radius >= 1.0_rp) then
                write(stderr, *) "test_level_shifted_davidson failed: Davidson "// &
                    "method did not stop and have its step rejected when residual "// &
                    "stagnates."
                test_level_shifted_davidson = .false.
            end if
            if (i_case == 2 .and. &
                (.not. jacobi_davidson_started .or. imicro_jacobi_davidson /= 11)) then
                write(stderr, *) "test_level_shifted_davidson failed: "// &
                    "Jacobi-Davidson method not switched to after ten micro "// &
                    "iterations when residual stagnates."
                test_level_shifted_davidson = .false.
            end if
        end do

        ! let the reduced space grow until it can no longer be expanded, a vanishing
        ! reduction factor prevents natural convergence, and one large Hessian
        ! eigenvalue ensures that even the exact full-rank solution's residual carries
        ! enough floating-point noise to stay above the solver's fixed convergence
        ! floor, for a diagonal Hessian with only two distinct eigenvalues the
        ! preconditioned residuals stay in a small subspace so that a new trial vector
        ! becomes linearly dependent, while for the same Hessian in a random
        ! orthonormal basis the reduced space grows to the dimension of the full
        ! parameter space
        do i_case = 1, 2
            call setup_settings(settings, quadratic_context)
            settings%local_red_factor = 0.0_rp
            settings%global_red_factor = 0.0_rp
            quadratic_context%hess = identity_matrix(n_param)
            quadratic_context%hess(1, 1) = 1e8_rp
            if (i_case == 2) then
                call ref_symm_mat_diag(generate_random_symm_matrix(n_param), eigvals, &
                                       rotation)
                quadratic_context%hess = matmul( &
                    rotation, matmul(quadratic_context%hess, transpose(rotation)))
            end if
            quadratic_context%grad = spread(1.0_rp, 1, n_param)
            overflow_residual_tol = 1e2_rp * epsilon(1.0_rp) * 1e8_rp
            h_diag = [(quadratic_context%hess(i, i), i=1, n_param)]
            grad_norm = norm2(quadratic_context%grad)
            func = 0.0_rp
            input_trust_radius = 1.0_rp
            trust_radius = input_trust_radius
            obj_func_funptr => quadratic_obj_func
            hess_x_funptr => quadratic_hess_x

            call level_shifted_davidson( &
                func, quadratic_context%grad, grad_norm, h_diag, n_param, &
                obj_func_funptr, hess_x_funptr, settings, trust_radius, solution, mu, &
                imicro, imicro_jacobi_davidson, jacobi_davidson_started, &
                max_precision_reached, error)
            if (error /= 0) then
                write(stderr, *) "test_level_shifted_davidson failed: Produced "// &
                    "error when reduced space stops growing for the "// &
                    trim(case_names(i_case))//" Hessian."
                test_level_shifted_davidson = .false.
            end if
            if (abs(norm2(solution) - input_trust_radius) > overflow_residual_tol) then
                write(stderr, *) "test_level_shifted_davidson failed: Solution "// &
                    "does not lie at trust region boundary when reduced space "// &
                    "stops growing for the "//trim(case_names(i_case))//" Hessian."
                test_level_shifted_davidson = .false.
            end if
            if (norm2(quadratic_context%grad + &
                      matmul(quadratic_context%hess, solution) - mu * solution) > &
                overflow_residual_tol) then
                write(stderr, *) "test_level_shifted_davidson failed: Solution "// &
                    "does not describe level-shifted Newton step when reduced "// &
                    "space stops growing for the "//trim(case_names(i_case))// &
                    " Hessian."
                test_level_shifted_davidson = .false.
            end if
        end do

    end function test_level_shifted_davidson

    logical(c_bool) function test_truncated_conjugate_gradient() bind(C)
        !
        ! this function tests the truncated conjugate gradient subroutine
        !
        use opentrustregion, only: obj_func_type, hess_x_type, solver_settings_type, &
                                   truncated_conjugate_gradient, error_hess_x, &
                                   error_obj_func, function_unchanged_warning_msg

        real(rp) :: func, trust_radius, ratio, solution_norm
        real(rp), dimension(n_param) :: grad, h_diag, solution
        integer(ip) :: i, imicro, error
        procedure(obj_func_type), pointer :: obj_func_funptr
        procedure(hess_x_type), pointer :: hess_x_funptr
        type(solver_settings_type) :: settings
        logical :: max_precision_reached
        real(rp), parameter :: h_diag_floor = 1e-10_rp
        type(hartmann6d_context_type), target :: context

        ! assume tests pass
        test_truncated_conjugate_gradient = .true.

        ! setup settings object
        call setup_settings(settings, context)

        ! initialize variables
        trust_radius = 0.4_rp
        obj_func_funptr => obj_func
        hess_x_funptr => hess_x_fun

        ! start in quadratic region near minimum
        context%vars = near_minimum
        func = hartmann6d_func(context%vars)
        call hartmann6d_gradient(context%vars, grad)
        context%hess = hartmann6d_hessian(context%vars)
        h_diag = [(context%hess(i, i), i=1, size(h_diag))]

        ! run truncated conjugate gradient, check whether the solution stays within the
        ! trust region and reduces the function value and whether the reported number
        ! of Hessian linear transformations agrees with the calls
        settings%n_hess_x = 0
        context%n_hess_x_calls = 0
        call truncated_conjugate_gradient( &
            func, grad, h_diag, n_param, obj_func_funptr, hess_x_funptr, settings, &
            trust_radius, solution, imicro, max_precision_reached, error)
        if (error /= 0) then
            write(stderr, *) "test_truncated_conjugate_gradient failed: Produced "// &
                "error near minimum."
            test_truncated_conjugate_gradient = .false.
        end if
        test_truncated_conjugate_gradient = &
            test_truncated_conjugate_gradient .and. logical(check_call_counts( &
                context, "truncated_conjugate_gradient", "near minimum", &
                n_hess_x=settings%n_hess_x), kind=c_bool)
        ratio = (hartmann6d_func(context%vars + solution) - func) / &
                dot_product(solution, grad + 0.5_rp * matmul(context%hess, solution))
        if (ratio <= 0.0_rp) then
            write(stderr, *) "test_truncated_conjugate_gradient failed: Solution "// &
                "does not reduce function value near minimum."
            test_truncated_conjugate_gradient = .false.
        end if
        solution_norm = sqrt(dot_product(solution, &
                                         solution / max(abs(h_diag), h_diag_floor)))
        if (solution_norm > ref_step_trust_radius(trust_radius, ratio) + tol) then
            write(stderr, *) "test_truncated_conjugate_gradient failed: Solution "// &
                "does not stay within trust region near minimum."
            test_truncated_conjugate_gradient = .false.
        end if

        ! start near saddle point
        context%vars = near_saddle_point
        func = hartmann6d_func(context%vars)
        call hartmann6d_gradient(context%vars, grad)
        context%hess = hartmann6d_hessian(context%vars)
        h_diag = [(context%hess(i, i), i=1, size(h_diag))]
        trust_radius = 0.4_rp

        ! run truncated conjugate gradient, check whether the solution lies at the
        ! trust region boundary and reduces the function value
        call truncated_conjugate_gradient( &
            func, grad, h_diag, n_param, obj_func_funptr, hess_x_funptr, settings, &
            trust_radius, solution, imicro, max_precision_reached, error)
        if (error /= 0) then
            write(stderr, *) "test_truncated_conjugate_gradient failed: Produced "// &
                "error near saddle point."
            test_truncated_conjugate_gradient = .false.
        end if
        ratio = (hartmann6d_func(context%vars + solution) - func) / &
                dot_product(solution, grad + 0.5_rp * matmul(context%hess, solution))
        if (ratio <= 0.0_rp) then
            write(stderr, *) "test_truncated_conjugate_gradient failed: Solution "// &
                "does not reduce function value near saddle point."
            test_truncated_conjugate_gradient = .false.
        end if
        solution_norm = sqrt(dot_product(solution, &
                                         solution / max(abs(h_diag), h_diag_floor)))
        if (abs(solution_norm - ref_step_trust_radius(trust_radius, ratio)) > tol) then
            write(stderr, *) "test_truncated_conjugate_gradient failed: Solution "// &
                "does not lie at trust region boundary near saddle point."
            test_truncated_conjugate_gradient = .false.
        end if

        ! run truncated conjugate gradient with a Hessian linear transformation which
        ! fails and check that the failing call is still counted
        hess_x_funptr => hess_x_fun_failing
        trust_radius = 0.4_rp
        settings%n_hess_x = 0
        context%n_hess_x_calls = 0
        call truncated_conjugate_gradient( &
            func, grad, h_diag, n_param, obj_func_funptr, hess_x_funptr, settings, &
            trust_radius, solution, imicro, max_precision_reached, error)
        if (error /= error_hess_x + 1) then
            write(stderr, *) "test_truncated_conjugate_gradient failed: Did not "// &
                "report the error of a failing Hessian linear transformation with "// &
                "its origin."
            test_truncated_conjugate_gradient = .false.
        end if
        test_truncated_conjugate_gradient = &
            test_truncated_conjugate_gradient .and. logical(check_call_counts( &
                context, "truncated_conjugate_gradient", "for failing Hessian "// &
                "linear transformation", n_hess_x=settings%n_hess_x), kind=c_bool)

        ! run truncated conjugate gradient with an objective function which fails and
        ! check that its error is reported with its origin
        trust_radius = 0.4_rp
        obj_func_funptr => obj_func_failing
        hess_x_funptr => hess_x_fun
        call truncated_conjugate_gradient( &
            func, grad, h_diag, n_param, obj_func_funptr, hess_x_funptr, settings, &
            trust_radius, solution, imicro, max_precision_reached, error)
        if (error /= error_obj_func + 1) then
            write(stderr, *) "test_truncated_conjugate_gradient failed: Did not "// &
                "report the error of a failing objective function with its origin."
            test_truncated_conjugate_gradient = .false.
        end if

        ! evaluate an objective function which does not change with the step and
        ! check that the calculation is reported as converged up to floating point
        ! precision and that a warning is printed
        func = 1.0_rp
        trust_radius = 0.4_rp
        obj_func_funptr => constant_obj_func
        context%log_message = ""
        call truncated_conjugate_gradient( &
            func, grad, h_diag, n_param, obj_func_funptr, hess_x_funptr, settings, &
            trust_radius, solution, imicro, max_precision_reached, error)
        if (error /= 0) then
            write(stderr, *) "test_truncated_conjugate_gradient failed: Produced "// &
                "error when function value does not change."
            test_truncated_conjugate_gradient = .false.
        end if
        if (.not. max_precision_reached) then
            write(stderr, *) "test_truncated_conjugate_gradient failed: Maximum "// &
                "precision not reached when function value does not change."
            test_truncated_conjugate_gradient = .false.
        end if
        if (index(context%log_message, " "//function_unchanged_warning_msg) == 0) then
            write(stderr, *) "test_truncated_conjugate_gradient failed: Warning "// &
                "not printed when function value does not change."
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
        type(test_context_type), target :: context

        ! assume tests pass
        test_add_error_origin = .true.

        ! setup settings object
        call setup_settings(settings, context)

        ! check if subroutine adds error origin correctly if no origin is present
        error = 1
        call add_error_origin(error, 100_ip, settings)
        if (error /= 101) then
            write(stderr, *) "test_add_error_origin failed: Error origin not "// &
                "correctly added."
            test_add_error_origin = .false.
        end if

        ! check if subroutine skips adding error origin if origin is already present
        call add_error_origin(error, 100_ip, settings)
        if (error /= 101) then
            write(stderr, *) "test_add_error_origin failed: Error code modified "// &
                "even though error origin is already present."
            test_add_error_origin = .false.
        end if

        ! check if subroutine does not modify error code when no error is encountered
        error = 0
        call add_error_origin(error, 100_ip, settings)
        if (error /= 0) then
            write(stderr, *) "test_add_error_origin failed: Error code modified "// &
                "even though error code of zero was passed."
            test_add_error_origin = .false.
        end if

        ! check if subroutine raises error for invalid error code
        error = -1
        call setup_error_logging(settings, context)
        call add_error_origin(error, 100_ip, settings)
        if (error /= 101) then
            write(stderr, *) "test_add_error_origin failed: Error code not "// &
                "correctly returned for invalid (negative) error code."
            test_add_error_origin = .false.
        end if
        if (len_trim(context%log_message) == 0) then
            write(stderr, *) "test_add_error_origin failed: No error message "// &
                "printed for invalid (negative) error code."
            test_add_error_origin = .false.
        end if

    end function test_add_error_origin

    logical(c_bool) function test_string_to_lowercase() bind(C)
        !
        ! this function tests the function that transfers strings to lowercase
        !
        use opentrustregion, only: string_to_lowercase

        character(len=*), parameter :: input = "OpenTrustRegion123!", &
                                       expect = "opentrustregion123!"

        ! assume tests pass
        test_string_to_lowercase = .true.

        ! test transfer to lowercase
        if (string_to_lowercase(input) /= expect) then
            write(stderr, *) "test_string_to_lowercase failed: String not "// &
                "correctly transferred to lowercase."
            test_string_to_lowercase = .false.
        end if

    end function test_string_to_lowercase

    logical(c_bool) function test_option_list() bind(C)
        !
        ! this function tests the function which lists options
        !
        use opentrustregion, only: option_list

        ! assume tests pass
        test_option_list = .true.

        ! check that options are quoted, stripped of trailing blanks and separated by
        ! commas
        if (option_list([character(len=5) :: "ab", "c", "d e"]) /= &
            """ab"", ""c"", ""d e""") then
            write(stderr, *) "test_option_list failed: Several options not listed "// &
                "correctly."
            test_option_list = .false.
        end if
        if (option_list(["ab"]) /= """ab""") then
            write(stderr, *) "test_option_list failed: Single option not listed "// &
                "correctly."
            test_option_list = .false.
        end if

    end function test_option_list

end module opentrustregion_unit_tests
