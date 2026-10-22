! Copyright (C) 2025- Jonas Greiner
!
! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.

module opentrustregion_unit_tests

    use opentrustregion, only: rp, ip, stderr, solver_state_type, stability_state_type
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
    ! and near the saddle point, a point far from both and a point with strong
    ! negative curvature
    real(rp), parameter :: near_minimum(n_param) = &
        [0.20_rp, 0.15_rp, 0.48_rp, 0.28_rp, 0.31_rp, 0.66_rp]
    real(rp), parameter :: near_saddle_point(n_param) = &
        [0.35_rp, 0.59_rp, 0.48_rp, 0.40_rp, 0.31_rp, 0.32_rp]
    real(rp), parameter :: distant_point(n_param) = &
        [0.9_rp, 0.1_rp, 0.7_rp, 0.2_rp, 0.6_rp, 0.4_rp]
    real(rp), parameter :: negative_curvature_point(n_param) = &
        [0.29_rp, 0.47_rp, 0.66_rp, 0.41_rp, 0.23_rp, 0.26_rp]

    ! dimension of the reduced space of subproblems
    integer(ip), parameter :: n_trial = 3

    real(c_rp), bind(C, name="hartmann6d_minimum1") :: minimum1_c(n_param) = minimum1
    real(c_rp), bind(C, name="hartmann6d_saddle_point") :: saddle_point_c(n_param) = &
        saddle_point
    real(c_rp), bind(C, name="hartmann6d_near_minimum") :: near_minimum_c(n_param) = &
        near_minimum

    ! define type for the host context handed to the mock callback functions, which
    ! holds the messages passed to the mock logging function and points at the state of
    ! the stability check that the mock convergence check of the stability check reads,
    ! with the RMS residual this convergence check saw last
    type, extends(host_context_type) :: test_context_type
        character(len=:), allocatable :: log_message
        type(stability_state_type), pointer :: stability_state => null()
        real(rp) :: seen_residual_rms = -1.0_rp
    end type

    ! define type for the host context handed to the mock callback functions of the
    ! Hartmann 6D function, which holds its state and counts the callback calls
    type, extends(test_context_type) :: hartmann6d_context_type
        ! current variables and Hessian
        real(rp) :: vars(n_param) = 0.0_rp, hess(n_param, n_param) = 0.0_rp

        ! number of calls to the orbital update, the objective function and the
        ! Hessian linear transformation, so that the counters reported by the library
        ! can be checked, and number of orbital updates at the previous convergence
        ! check
        integer(ip) :: n_update_orbs_calls = 0, n_obj_func_calls = 0, &
                       n_hess_x_calls = 0, n_update_orbs_at_conv_check = -1
    end type

    ! define type for the host context of the Hartmann 6D function which additionally
    ! records whether the objective function was evaluated without displacement,
    ! whether the Hessian was applied to the normalized gradient and the last
    ! displacement passed to the orbital update
    type, extends(hartmann6d_context_type) :: hartmann6d_recording_context_type
        logical :: obj_func_zero_step = .false., hess_x_normalized_grad = .false.
        real(rp) :: last_delta_vars(n_param) = 0.0_rp
    end type

    ! callback functions of the Hartmann 6D function into which a fault can be
    ! injected, with their names for failure messages
    integer(ip), parameter :: fault_update_orbs = 1, fault_obj_func = 2, &
                              fault_hess_x = 3, fault_precond = 4, fault_project = 5, &
                              fault_conv_check = 6, fault_get_extra_trial_vectors = 7, &
                              fault_stability_conv_check = 8, n_fault_callbacks = 8
    character(len=*), parameter :: fault_callback_names(n_fault_callbacks) = &
        [character(len=33) :: "orbital update", "objective function", "Hessian "// &
         "linear transformation", "preconditioner", "projection", "convergence check", &
         "extra trial vector function", "stability check convergence check"]

    ! define type for the host context of the Hartmann 6D function which counts the
    ! calls of every callback function and lets the call fault_call of the callback
    ! function fault_callback fail
    type, extends(hartmann6d_context_type) :: hartmann6d_fault_context_type
        integer(ip) :: fault_callback = 0, fault_call = 0, &
                       n_callback_calls(n_fault_callbacks) = 0
    end type

    ! define type for the host context handed to the mock callback functions of a
    ! quadratic model of any dimension, which holds its gradient and Hessian, the
    ! index of the next unit vector a preconditioner returns and the level shift a
    ! preconditioner was last called with
    type, extends(test_context_type) :: quadratic_context_type
        real(rp), allocatable :: grad(:), hess(:, :)
        integer(ip) :: next_unit_vector = 0
        real(rp) :: recorded_mu = 0.0_rp
    end type

    ! define type for the host context handed to the mock callback functions of the
    ! double well function x^2 - y^2 + y^4, which has a saddle point at the origin and
    ! minima at y = +-1/sqrt(2), and holds its current variables, records whether the
    ! Hessian linear transformation received it, counts the objective function
    ! evaluations and keeps the objective function at the previous convergence check
    type, extends(test_context_type) :: double_well_context_type
        real(rp) :: vars(2) = 0.0_rp, last_conv_check_func = huge(1.0_rp)
        logical :: hess_x_called = .false.
        integer(ip) :: n_obj_func_calls = 0
    end type

    ! define type for the host context of the double well function whose orbital update
    ! reads the solver state through a pointer during the optimization and records it
    ! with the size of the step it receives, in the Euclidean norm and in the norm of
    ! the absolute diagonal preconditioner of the previous point, the metric of the
    ! trust region of truncated conjugate gradient, for which it keeps the last Hessian
    ! diagonal
    type, extends(double_well_context_type) :: double_well_state_context_type
        type(solver_state_type), pointer :: state => null()
        type(solver_state_type), allocatable :: seen_states(:)
        real(rp), allocatable :: step_norms(:), precond_step_norms(:)
        real(rp) :: h_diag(2) = 1.0_rp
    end type

    ! define type for the host context of the quadratic model 1/2 |x|^2 whose
    ! objective function is evaluated at half the displacement, which holds its
    ! current variables and records the first nonvanishing displacement passed to the
    ! orbital update
    type, extends(test_context_type) :: half_step_context_type
        real(rp) :: vars(2) = 0.0_rp, first_step(2) = 0.0_rp
    end type

    ! define type for the host context handed to the objective function used to test
    ! the bracketing, which selects one of the one-dimensional test functions and the
    ! call of the objective function that fails, none if zero, the calls are counted
    ! by the host context
    type, extends(test_context_type) :: bracket_context_type
        integer(ip) :: func_case = 0, fail_call = 0
    end type

    ! Hartmann 6D context of the mock callback functions of a solver called without
    ! settings, which hands them no host context, so that they reach it here instead
    type(hartmann6d_context_type), target :: settingless_context

    ! error codes returned by mock callback functions reached without their context or
    ! with a context they should not receive, chosen so that they cannot be mistaken
    ! for an expected error code
    integer(ip), parameter :: missing_context_error = 42, unexpected_context_error = 43

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
        ! recalculate linear transformations which no longer respect Hessian symmetry,
        ! and into which a fault can be injected
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
        if (inject_fault(context, fault_hess_x)) error = 1

    end subroutine hess_x_fun_asymmetric

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

        ! count call
        state%n_obj_func_calls = state%n_obj_func_calls + 1

        ! record whether function is evaluated without displacement
        select type (state)
        class is (hartmann6d_recording_context_type)
            if (maxval(abs(delta_vars)) <= 0.0_rp) state%obj_func_zero_step = .true.
        end select

        func = hartmann6d_func(state%vars + delta_vars)

    end function obj_func

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

    subroutine update_orbs_hess_x_once(delta_vars, func, grad, h_diag, hess_x_funptr, &
                                       error, context)
        !
        ! this subroutine describes an orbital update which reports success but only
        ! provides a Hessian linear transformation on its first call and leaves the
        ! argument as it received it on later calls
        !
        use opentrustregion, only: hess_x_type

        real(rp), intent(in), target :: delta_vars(:)
        real(rp), intent(out) :: func
        real(rp), intent(out), target :: grad(:), h_diag(:)
        procedure(hess_x_type), intent(inout), pointer :: hess_x_funptr
        integer(ip), intent(out) :: error
        class(*), intent(in), pointer :: context

        procedure(hess_x_type), pointer :: received_hess_x_funptr
        class(hartmann6d_context_type), pointer :: state

        ! keep the received argument
        received_hess_x_funptr => hess_x_funptr

        ! perform orbital update, which counts the call
        call update_orbs(delta_vars, func, grad, h_diag, hess_x_funptr, error, context)
        if (error /= 0) return

        ! restore the received argument on every call but the first
        state => resolve_hartmann6d_context(context, error)
        if (error /= 0) return
        if (state%n_update_orbs_calls > 1) hess_x_funptr => received_hess_x_funptr

    end subroutine update_orbs_hess_x_once

    subroutine settingless_update_orbs(delta_vars, func, grad, h_diag, hess_x_funptr, &
                                       error, context)
        !
        ! this subroutine describes the orbital update for the Hartmann 6D function of
        ! a solver called without settings, which reaches the variables through the
        ! module context since it receives no host context, and returns a Hessian
        ! linear transformation that does the same
        !
        use opentrustregion, only: hess_x_type

        real(rp), intent(in), target :: delta_vars(:)
        real(rp), intent(out) :: func
        real(rp), intent(out), target :: grad(:), h_diag(:)
        procedure(hess_x_type), intent(inout), pointer :: hess_x_funptr
        integer(ip), intent(out) :: error
        class(*), intent(in), pointer :: context

        class(*), pointer :: module_context

        ! check that the solver handed no host context
        if (associated(context)) then
            error = unexpected_context_error
            return
        end if

        ! update the orbitals with the module context
        module_context => settingless_context
        call update_orbs(delta_vars, func, grad, h_diag, hess_x_funptr, error, &
                         module_context)
        hess_x_funptr => settingless_hess_x

    end subroutine settingless_update_orbs

    subroutine settingless_hess_x(x, hess_x, error, context)
        !
        ! this subroutine describes the Hessian linear transformation for the Hartmann
        ! 6D function of a solver called without settings, which reaches the variables
        ! through the module context since it receives no host context
        !
        real(rp), intent(in), target :: x(:)
        real(rp), intent(out), target :: hess_x(:)
        integer(ip), intent(out) :: error
        class(*), intent(in), pointer :: context

        class(*), pointer :: module_context

        ! check that the solver or stability check handed no host context
        if (associated(context)) then
            error = unexpected_context_error
            return
        end if

        ! apply the Hessian with the module context
        module_context => settingless_context
        call hess_x_fun(x, hess_x, error, module_context)

    end subroutine settingless_hess_x

    function settingless_obj_func(delta_vars, error, context) result(func)
        !
        ! this function describes the objective function evaluation for the Hartmann
        ! 6D function of a solver called without settings, which reaches the variables
        ! through the module context since it receives no host context
        !
        real(rp), intent(in), target :: delta_vars(:)
        integer(ip), intent(out) :: error
        class(*), intent(in), pointer :: context
        real(rp) :: func

        class(*), pointer :: module_context

        ! check that the solver handed no host context
        if (associated(context)) then
            func = 0.0_rp
            error = unexpected_context_error
            return
        end if

        ! evaluate the function with the module context
        module_context => settingless_context
        func = obj_func(delta_vars, error, module_context)

    end function settingless_obj_func

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

    function first_step_rejecting_obj_func(delta_vars, error, context) result(func)
        !
        ! this function describes the objective function evaluation of the Hartmann 6D
        ! function into which a fault can be injected and which raises the value of its
        ! first evaluation, so that the first step is rejected
        !
        real(rp), intent(in), target :: delta_vars(:)
        integer(ip), intent(out) :: error
        class(*), intent(in), pointer :: context
        real(rp) :: func

        ! evaluate objective function and inject fault
        func = faulty_obj_func(delta_vars, error, context)
        if (error /= 0) return

        ! raise the value of the first evaluation
        select type (context)
        class is (hartmann6d_fault_context_type)
            if (context%n_callback_calls(fault_obj_func) == 1) func = func + 1.0_rp
        end select

    end function first_step_rejecting_obj_func

    subroutine faulty_hess_x(x, hess_x, error, context)
        !
        ! this subroutine describes the Hessian linear transformation of the
        ! Hartmann 6D function into which a fault can be injected
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
        ! this subroutine describes the identity preconditioner into which a fault can
        ! be injected
        !
        real(rp), intent(in), target :: residual(:)
        real(rp), intent(in) :: mu
        real(rp), intent(out), target :: precond_residual(:)
        integer(ip), intent(out) :: error
        class(*), intent(in), pointer :: context

        ! apply preconditioner
        call identity_precond(residual, mu, precond_residual, error, context)
        if (error /= 0) return

        ! inject fault
        if (inject_fault(context, fault_precond)) error = 1

    end subroutine faulty_precond

    subroutine faulty_project(vector, error, context)
        !
        ! this subroutine describes the identity projection into which a fault can be
        ! injected
        !
        real(rp), intent(inout), target :: vector(:)
        integer(ip), intent(out) :: error
        class(*), intent(in), pointer :: context

        ! apply projection
        call identity_project(vector, error, context)
        if (error /= 0) return

        ! inject fault
        if (inject_fault(context, fault_project)) error = 1

    end subroutine faulty_project

    function faulty_conv_check(error, context) result(converged)
        !
        ! this function describes a convergence check of the Hartmann 6D function which
        ! reports convergence once the RMS gradient at the current variables drops
        ! below the default convergence threshold, since it replaces the gradient
        ! criterion, and into which a fault can be injected
        !
        use opentrustregion, only: default_solver_settings
        use test_reference, only: check_host_context

        integer(ip), intent(out) :: error
        class(*), intent(in), pointer :: context
        logical :: converged

        class(hartmann6d_context_type), pointer :: state
        real(rp) :: grad(n_param)

        ! check host context
        call check_host_context(context)

        ! initialize result in case the test problem state is missing
        converged = .false.

        ! get Hartmann 6D state
        state => resolve_hartmann6d_context(context, error)
        if (error /= 0) return

        ! check the gradient and inject fault
        call hartmann6d_gradient(state%vars, grad)
        converged = norm2(grad) / sqrt(real(n_param, kind=rp)) < &
                    default_solver_settings%conv_tol
        if (inject_fault(context, fault_conv_check)) error = 1

    end function faulty_conv_check

    subroutine faulty_get_extra_trial_vectors(trial_vectors, error, context)
        !
        ! this subroutine describes the extra trial vector function which returns
        ! scaled unit vectors into which a fault can be injected
        !
        real(rp), intent(out), target :: trial_vectors(:, :)
        integer(ip), intent(out) :: error
        class(*), intent(in), pointer :: context

        ! return extra trial vectors
        call mock_get_extra_trial_vectors(trial_vectors, error, context)
        if (error /= 0) return

        ! inject fault
        if (inject_fault(context, fault_get_extra_trial_vectors)) error = 1

    end subroutine faulty_get_extra_trial_vectors

    function faulty_stability_conv_check(error, context) result(converged)
        !
        ! this function describes the convergence check of the stability check which
        ! accepts a residual in the state of the stability check that vanishes relative
        ! to the lowest eigenvalue into which a fault can be injected
        !
        integer(ip), intent(out) :: error
        class(*), intent(in), pointer :: context
        logical :: converged

        ! check the state
        converged = mock_stability_conv_check(error, context)
        if (error /= 0) return

        ! inject fault
        if (inject_fault(context, fault_stability_conv_check)) error = 1

    end function faulty_stability_conv_check

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

    subroutine mu_recording_precond(residual, mu, precond_residual, error, context)
        !
        ! this subroutine describes a preconditioner for the quadratic model which
        ! scales the residual by the level shift it is called with and records it
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

        ! get quadratic model state
        precond_residual = 0.0_rp
        state => resolve_quadratic_context(context, error)
        if (error /= 0) return

        ! record the level shift and scale the residual by it
        state%recorded_mu = mu
        precond_residual = mu * residual

    end subroutine mu_recording_precond

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

        state%hess_x_called = .true.
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

        state%n_obj_func_calls = state%n_obj_func_calls + 1
        func = double_well_func(state%vars + delta_vars)

    end function double_well_obj_func

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

    subroutine state_reading_hess_x(x, hess_x, error, context)
        !
        ! this subroutine describes the Hessian linear transformation of the double
        ! well function, which also records the solver state it reads through the
        ! pointer in its context
        !
        real(rp), intent(in), target :: x(:)
        real(rp), intent(out), target :: hess_x(:)
        integer(ip), intent(out) :: error
        class(*), intent(in), pointer :: context

        ! apply the Hessian of the double well function
        call double_well_hess_x(x, hess_x, error, context)
        if (error /= 0) return

        ! record the state
        select type (context)
        class is (double_well_state_context_type)
            if (.not. allocated(context%seen_states)) &
                allocate(context%seen_states(0), context%step_norms(0), &
                         context%precond_step_norms(0))
            context%seen_states = [context%seen_states, context%state]
        end select

    end subroutine state_reading_hess_x

    subroutine state_reading_update_orbs(delta_vars, func, grad, h_diag, &
                                         hess_x_funptr, error, context)
        !
        ! this subroutine describes the orbital update of the double well function,
        ! which also records the solver state it reads through the pointer in its
        ! context and the size of the step it receives in the Euclidean norm and in
        ! the norm of the absolute diagonal preconditioner of the previous point
        !
        use opentrustregion, only: hess_x_type, precond_floor

        real(rp), intent(in), target :: delta_vars(:)
        real(rp), intent(out) :: func
        real(rp), intent(out), target :: grad(:), h_diag(:)
        procedure(hess_x_type), intent(inout), pointer :: hess_x_funptr
        integer(ip), intent(out) :: error
        class(*), intent(in), pointer :: context

        ! update the orbitals of the double well function
        call double_well_update_orbs(delta_vars, func, grad, h_diag, hess_x_funptr, &
                                     error, context)
        if (error /= 0) return

        ! record the state and the size of the step in both norms and keep the
        ! Hessian diagonal for the norm of the next step
        select type (context)
        class is (double_well_state_context_type)
            if (.not. allocated(context%seen_states)) &
                allocate(context%seen_states(0), context%step_norms(0), &
                         context%precond_step_norms(0))
            context%seen_states = [context%seen_states, context%state]
            context%step_norms = [context%step_norms, norm2(delta_vars)]
            context%precond_step_norms = [ &
                context%precond_step_norms, &
                sqrt(sum(delta_vars**2 / max(abs(context%h_diag), precond_floor)))]
            context%h_diag = h_diag
        end select

    end subroutine state_reading_update_orbs

    function double_well_conv_check(error, context) result(converged)
        !
        ! this function describes a convergence check of the double well function
        ! which, like a host's check of the change of the objective function, reports
        ! convergence once the objective function has changed by less than 1e-8 since
        ! the previous check, so that it can not pass at a stationary point the solver
        ! has just reached
        !
        use test_reference, only: check_host_context

        integer(ip), intent(out) :: error
        class(*), intent(in), pointer :: context
        logical :: converged

        class(double_well_context_type), pointer :: state
        real(rp) :: func

        ! check host context
        call check_host_context(context)

        ! initialize result in case the double well state is missing
        converged = .false.

        ! get double well state
        state => resolve_double_well_context(context, error)
        if (error /= 0) return

        ! compare the objective function against that of the previous check and keep
        ! it for the next check
        func = double_well_func(state%vars)
        converged = abs(func - state%last_conv_check_func) < 1e-8_rp
        state%last_conv_check_func = func

    end function double_well_conv_check

    function resolve_half_step_context(context, error) result(state)
        !
        ! this function returns the context of the half step model handed to a mock
        ! callback function and an error if no such context was handed over
        !
        class(*), intent(in), pointer :: context
        integer(ip), intent(out) :: error
        class(half_step_context_type), pointer :: state

        ! initialize error flag
        error = 0

        ! get context
        state => null()
        if (associated(context)) then
            select type (context)
            class is (half_step_context_type)
                state => context
            end select
        end if

        ! report missing test context
        if (.not. associated(state)) error = missing_context_error

    end function resolve_half_step_context

    subroutine half_step_update_orbs(delta_vars, func, grad, h_diag, hess_x_funptr, &
                                     error, context)
        !
        ! this subroutine describes the orbital update of the quadratic model
        ! 1/2 |x|^2, whose Newton step leads to the origin, and records the first
        ! nonvanishing displacement
        !
        use opentrustregion, only: hess_x_type
        use test_reference, only: check_host_context

        real(rp), intent(in), target :: delta_vars(:)
        real(rp), intent(out) :: func
        real(rp), intent(out), target :: grad(:), h_diag(:)
        procedure(hess_x_type), intent(inout), pointer :: hess_x_funptr
        integer(ip), intent(out) :: error
        class(*), intent(in), pointer :: context

        class(half_step_context_type), pointer :: state

        ! check host context
        call check_host_context(context)

        ! get half step model state
        state => resolve_half_step_context(context, error)
        if (error /= 0) return

        ! record first nonvanishing displacement and update variables
        if (maxval(abs(state%first_step)) <= 0.0_rp) state%first_step = delta_vars
        state%vars = state%vars + delta_vars

        ! evaluate model, gradient and Hessian diagonal and define Hessian linear
        ! transformation
        func = 0.5_rp * dot_product(state%vars, state%vars)
        grad = state%vars
        h_diag = 1.0_rp
        hess_x_funptr => identity_hess_x

    end subroutine half_step_update_orbs

    function half_step_obj_func(delta_vars, error, context) result(func)
        !
        ! this function describes an objective function which evaluates the quadratic
        ! model 1/2 |x|^2 at half the displacement, so that its minimum along a step
        ! lies at twice the Newton step of the model
        !
        use test_reference, only: check_host_context

        real(rp), intent(in), target :: delta_vars(:)
        integer(ip), intent(out) :: error
        class(*), intent(in), pointer :: context
        real(rp) :: func

        class(half_step_context_type), pointer :: state

        ! initialize result in case the half step model state is missing
        func = 0.0_rp

        ! check host context
        call check_host_context(context)

        ! get half step model state
        state => resolve_half_step_context(context, error)
        if (error /= 0) return

        func = 0.5_rp * sum((state%vars + 0.5_rp * delta_vars)**2)

    end function half_step_obj_func

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

        ! get bracketing state and report missing test context
        error = 0
        state => null()
        if (associated(context)) then
            select type (context)
            class is (bracket_context_type)
                state => context
            end select
        end if
        if (.not. associated(state)) then
            error = missing_context_error
            return
        end if

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

        ! let the requested call fail
        if (state%n_calls == state%fail_call) error = 1

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

    subroutine vanishing_project(vector, error, context)
        !
        ! this subroutine describes a projection which removes every vector and needs
        ! no test context
        !
        use test_reference, only: check_host_context

        real(rp), intent(inout), target :: vector(:)
        integer(ip), intent(out) :: error
        class(*), intent(in), pointer :: context

        ! check host context
        call check_host_context(context)

        vector = 0.0_rp

        error = 0

    end subroutine vanishing_project

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

    subroutine mock_get_extra_trial_vectors(trial_vectors, error, context)
        !
        ! this subroutine is a test subroutine for the extra trial vector subroutine
        !
        use test_reference, only: check_host_context

        real(rp), intent(out), target :: trial_vectors(:, :)
        integer(ip), intent(out) :: error
        class(*), intent(in), pointer :: context

        integer(ip) :: i

        ! check host context
        call check_host_context(context)

        trial_vectors = 0.0_rp
        do i = 1, size(trial_vectors, 2)
            trial_vectors(i, i) = real(i, kind=rp)
        end do

        error = 0

    end subroutine mock_get_extra_trial_vectors

    function mock_stability_conv_check(error, context) result(converged)
        !
        ! this function describes a convergence check of the stability check which
        ! reads the state of the stability check the test context points at, records
        ! the RMS residual it sees there and accepts a residual that vanishes relative
        ! to the lowest eigenvalue
        !
        use test_reference, only: check_host_context

        integer(ip), intent(out) :: error
        class(*), intent(in), pointer :: context
        logical :: converged

        ! check host context
        call check_host_context(context)

        ! report a missing test context or state
        converged = .false.
        error = missing_context_error
        if (.not. associated(context)) return

        ! check the lowest eigenpair in the state
        select type (context)
        class is (test_context_type)
            if (.not. associated(context%stability_state)) return
            context%seen_residual_rms = context%stability_state%residual_rms
            converged = context%stability_state%residual_rms < &
                        tol * abs(context%stability_state%min_eigval)
            error = 0
        end select

    end function mock_stability_conv_check

    subroutine mock_vanishing_extra_trial_vectors(trial_vectors, error, context)
        !
        ! this subroutine is a test subroutine for an extra trial vector subroutine
        ! which has no direction to contribute
        !
        use test_reference, only: check_host_context

        real(rp), intent(out), target :: trial_vectors(:, :)
        integer(ip), intent(out) :: error
        class(*), intent(in), pointer :: context

        ! check host context
        call check_host_context(context)

        trial_vectors = 0.0_rp

        error = 0

    end subroutine mock_vanishing_extra_trial_vectors

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

    function generate_random_upper_triangular_matrix(n) result(matrix)
        !
        ! this function generates a random upper triangular matrix, whose eigenvalues
        ! are its diagonal elements
        !
        integer(ip), intent(in) :: n
        real(rp) :: matrix(n, n)

        integer(ip) :: i

        call random_number(matrix)
        do i = 1, n - 1
            matrix(i + 1:, i) = 0.0_rp
        end do

    end function generate_random_upper_triangular_matrix

    subroutine setup_settings(settings, context)
        !
        ! this subroutine sets up a settings object for tests, which hands the test
        ! context to the mock callback functions
        !
        use opentrustregion, only: settings_type, verbosity_info

        class(settings_type), intent(out) :: settings
        class(test_context_type), intent(inout), target :: context

        settings%verbose = verbosity_info
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

    subroutine setup_hartmann6d_point(context, vars, func, grad, h_diag)
        !
        ! this subroutine moves a Hartmann 6D context to a point and evaluates the
        ! function, its gradient and Hessian and the Hessian diagonal there
        !
        class(hartmann6d_context_type), intent(inout) :: context
        real(rp), intent(in) :: vars(:)
        real(rp), intent(out) :: func, grad(:), h_diag(:)

        integer(ip) :: i

        context%vars = vars
        func = hartmann6d_func(vars)
        call hartmann6d_gradient(vars, grad)
        context%hess = hartmann6d_hessian(vars)
        h_diag = [(context%hess(i, i), i=1, n_param)]

    end subroutine setup_hartmann6d_point

    subroutine setup_double_well_state_context(context, state, func, grad, h_diag)
        !
        ! this subroutine moves a double well context, whose callback functions read
        ! the solver state, to a point of negative curvature, evaluates the function,
        ! its gradient and its Hessian diagonal there and points the context to a
        ! state holding the values of an earlier subproblem, which a subsystem solver
        ! has to reset before it calls the callback functions
        !
        type(double_well_state_context_type), intent(inout) :: context
        type(solver_state_type), intent(inout), target :: state
        real(rp), intent(out) :: func, grad(2), h_diag(2)

        ! move to the point and evaluate the function there
        context%vars = [0.5_rp, 0.1_rp]
        func = double_well_func(context%vars)
        grad = [2.0_rp * context%vars(1), &
                -2.0_rp * context%vars(2) + 4.0_rp * context%vars(2)**3]
        h_diag = [2.0_rp, -2.0_rp + 12.0_rp * context%vars(2)**2]

        ! point to a state of an earlier subproblem
        context%state => state
        state%trust_radius = 0.4_rp
        state%micro_iteration = 7
        state%jacobi_davidson_iteration = 3
        state%relative_residual = 0.3_rp

    end subroutine setup_double_well_state_context

    subroutine setup_reduced_space(vars, grad, grad_norm, red_space_basis, &
                                   red_space_hess)
        !
        ! this subroutine sets up the reduced space of the subproblem at a point of the
        ! Hartmann 6D function, the gradient, an orthonormal basis which starts with
        ! the gradient direction and continues with the leading unit vectors and the
        ! Hessian in this basis
        !
        real(rp), intent(in) :: vars(:)
        real(rp), intent(out) :: grad(:), grad_norm, red_space_basis(:, :), &
                                 red_space_hess(:, :)

        integer(ip) :: i, j

        ! calculate gradient
        call hartmann6d_gradient(vars, grad)
        grad_norm = norm2(grad)

        ! orthonormalize the leading unit vectors against the gradient direction
        red_space_basis = 0.0_rp
        red_space_basis(:, 1) = grad / grad_norm
        do i = 2, size(red_space_basis, 2)
            red_space_basis(i - 1, i) = 1.0_rp
            do j = 1, i - 1
                red_space_basis(:, i) = red_space_basis(:, i) - dot_product( &
                    red_space_basis(:, i), red_space_basis(:, j)) * &
                                        red_space_basis(:, j)
            end do
            red_space_basis(:, i) = red_space_basis(:, i) / norm2(red_space_basis(:, i))
        end do

        ! project Hessian onto the basis
        red_space_hess = matmul(transpose(red_space_basis), &
                                matmul(hartmann6d_hessian(vars), red_space_basis))

    end subroutine setup_reduced_space

    subroutine setup_hard_case(hess_diag, red_space_basis, grad_norm, red_space_hess)
        !
        ! this subroutine sets up the reduced space of a hard case of the subproblem,
        ! a basis of the leading unit vectors, since the hard case step assumes an
        ! orthonormal basis, a unit gradient norm and a diagonal Hessian, whose lowest
        ! eigenvalue is negative and has no eigenvector along the gradient direction,
        ! the first basis vector
        !
        real(rp), intent(in) :: hess_diag(:)
        real(rp), intent(out) :: red_space_basis(:, :), grad_norm, red_space_hess(:, :)

        integer(ip) :: i

        red_space_basis = 0.0_rp
        red_space_hess = 0.0_rp
        do i = 1, size(hess_diag)
            red_space_basis(i, i) = 1.0_rp
            red_space_hess(i, i) = hess_diag(i)
        end do
        grad_norm = 1.0_rp

    end subroutine setup_hard_case

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

    subroutine ref_general_mat_diag(matrix, eigvals_re, eigvals_im, right_eigvecs, &
                                    left_eigvecs)
        !
        ! this subroutine reimplements the eigendecomposition of a general matrix,
        ! returning the real and imaginary parts of the eigenvalues and left
        ! eigenvectors scaled to be biorthonormal to the right ones
        !
        real(rp), intent(in) :: matrix(:, :)
        real(rp), intent(out) :: eigvals_re(:), eigvals_im(:), right_eigvecs(:, :), &
                                 left_eigvecs(:, :)

        integer(ip) :: n, lwork, info, i
        real(rp), allocatable :: work(:), temp(:, :)
        complex(rp), allocatable :: left(:)
        external :: dgeev

        n = size(matrix, 1)
        temp = matrix

        ! query optimal workspace size
        lwork = -1
        allocate(work(1))
        call dgeev("V", "V", n, temp, n, eigvals_re, eigvals_im, left_eigvecs, n, &
                   right_eigvecs, n, work, lwork, info)
        lwork = int(work(1), kind=ip)
        deallocate(work)
        allocate(work(lwork))

        ! perform eigendecomposition
        call dgeev("V", "V", n, temp, n, eigvals_re, eigvals_im, left_eigvecs, n, &
                   right_eigvecs, n, work, lwork, info)
        deallocate(work, temp)

        ! scale left eigenvectors to biorthonormality, the two columns holding the real
        ! and imaginary parts of a complex eigenvector together
        i = 1
        do while (i <= n)
            if (abs(eigvals_im(i)) > 0.0_rp) then
                left = cmplx(left_eigvecs(:, i), left_eigvecs(:, i + 1), kind=rp)
                left = left / conjg(dot_product( &
                    left, cmplx(right_eigvecs(:, i), right_eigvecs(:, i + 1), kind=rp)))
                left_eigvecs(:, i) = real(left)
                left_eigvecs(:, i + 1) = aimag(left)
                i = i + 2
            else
                left_eigvecs(:, i) = left_eigvecs(:, i) / dot_product( &
                    left_eigvecs(:, i), right_eigvecs(:, i))
                i = i + 1
            end if
        end do

    end subroutine ref_general_mat_diag

    function ref_hartmann6d_ratio(context, step) result(ratio)
        !
        ! this function reimplements the ratio of the actual to the predicted change of
        ! the Hartmann 6D function for a step from the point of a Hartmann 6D context,
        ! which the trust-region solvers hand the step acceptance
        !
        class(hartmann6d_context_type), intent(in) :: context
        real(rp), intent(in) :: step(:)
        real(rp) :: ratio

        real(rp) :: grad(n_param)

        call hartmann6d_gradient(context%vars, grad)
        ratio = &
            (hartmann6d_func(context%vars + step) - hartmann6d_func(context%vars)) / &
            dot_product(step, grad + 0.5_rp * matmul(context%hess, step))

    end function ref_hartmann6d_ratio

    function ref_step_trust_radius(trust_radius, ratio, step_norm)
        !
        ! this function recovers the trust radius an accepted step was computed for
        ! from the trust radius, the ratio of actual to predicted function change of
        ! the step and its norm, since a high ratio only expands the trust radius for a
        ! step which nearly reaches it, taking it as expanded where the norm of the
        ! step lies within both the expanded and the unchanged trust radius
        !
        use opentrustregion, only: &
            trust_radius_shrink_ratio, trust_radius_expand_ratio, &
            trust_radius_shrink_factor, trust_radius_expand_factor

        real(rp), intent(in) :: trust_radius, ratio, step_norm
        real(rp) :: ref_step_trust_radius

        if (ratio < trust_radius_shrink_ratio) then
            ref_step_trust_radius = trust_radius / trust_radius_shrink_factor
        else if ( &
            ratio >= trust_radius_expand_ratio .and. &
            step_norm >= 0.99_rp * trust_radius / trust_radius_expand_factor .and. &
            step_norm <= trust_radius / trust_radius_expand_factor + tol) then
            ref_step_trust_radius = trust_radius / trust_radius_expand_factor
        else
            ref_step_trust_radius = trust_radius
        end if

    end function ref_step_trust_radius

    logical function check_call_counts(context, test_name, case_name, n_hess_x, &
                                       n_update_orbs, n_obj_func)
        !
        ! this function checks that the counters reported by the library agree with
        ! the number of times the callback functions were actually called
        !
        class(hartmann6d_context_type), intent(in) :: context
        character(len=*), intent(in) :: test_name, case_name
        integer(ip), intent(in), optional :: n_hess_x, n_update_orbs, n_obj_func

        ! assume tests pass
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
        if (present(n_obj_func)) then
            if (n_obj_func /= context%n_obj_func_calls) then
                write(stderr, *) "test_"//test_name//" failed: Reported number of "// &
                    "objective function evaluations wrong "//case_name//"."
                check_call_counts = .false.
            end if
        end if

    end function check_call_counts

    logical function check_reset_state_seen(context, test_name, callback_name)
        !
        ! this function checks that the first call of a callback function which reads
        ! the solver state through the pointer of a double well context saw the state
        ! reset by the subsystem solver, without micro iterations and with the
        ! relative residual of the vanishing solution, which is the gradient
        !
        type(double_well_state_context_type), intent(in) :: context
        character(len=*), intent(in) :: test_name, callback_name

        ! assume tests pass
        check_reset_state_seen = .true.

        if (.not. allocated(context%seen_states)) then
            write(stderr, *) "test_"//test_name//" failed: State not read by the "// &
                callback_name//"."
            check_reset_state_seen = .false.
        else if (context%seen_states(1)%micro_iteration /= 0 .or. &
                 context%seen_states(1)%jacobi_davidson_iteration /= 0 .or. &
                 abs(context%seen_states(1)%relative_residual - 1.0_rp) > tol) then
            write(stderr, *) "test_"//test_name//" failed: State of an earlier "// &
                "subproblem seen by the first call of the "//callback_name//"."
            check_reset_state_seen = .false.
        end if

    end function check_reset_state_seen

    logical(c_bool) function test_solver() bind(C)
        !
        ! this function tests the solver subroutine
        !
        use opentrustregion, only: &
            update_orbs_type, obj_func_type, solver_settings_type, solver, &
            default_settings => default_solver_settings, error_solver_max_iter, &
            error_update_orbs, error_conv_check, error_obj_func, error_hess_x, &
            error_solver, error_precond, error_project, error_get_extra_trial_vectors, &
            subsystem_solver_options

        real(rp), parameter :: var_thres = 1e-4_rp
        integer(ip) :: error, i_solver, i_case, callback, i_call, &
                       n_calls(n_fault_callbacks)
        logical :: grad_noise
        integer(ip), parameter :: fault_origins(n_fault_callbacks) = &
            [error_update_orbs, error_obj_func, error_hess_x, error_precond, &
             error_project, error_conv_check, error_get_extra_trial_vectors, &
             error_conv_check]
        character(len=20) :: call_number
        type(hartmann6d_fault_context_type), target :: fault_context
        real(rp) :: final_grad(n_param), start_vars(2), final_vars(n_param)
        procedure(update_orbs_type), pointer :: update_orbs_funptr
        procedure(obj_func_type), pointer :: obj_func_funptr
        type(solver_settings_type) :: settings
        type(hartmann6d_recording_context_type), target :: context
        type(double_well_context_type), target :: double_well_context
        type(half_step_context_type), target :: half_step_context
        type(double_well_state_context_type), target :: state_context
        type(solver_state_type), target :: state
        type(solver_state_type) :: explicit_state
        integer(ip) :: n_updates, i_update
        real(rp), allocatable :: step_norms(:)
        character(len=60) :: case_name

        ! assume tests pass
        test_solver = .true.

        ! start at saddle point
        context%vars = saddle_point
        update_orbs_funptr => update_orbs
        obj_func_funptr => obj_func

        ! hand the callback functions the host context, which they need to evaluate
        ! the function, so that the internal stability check has to inherit it
        settings%context => context

        ! run solver, check if error has occurred and check whether gradient is zero
        ! and agrees with correct minimum and that convergence is not flagged as
        ! reaching maximum precision
        call solver(update_orbs_funptr, obj_func_funptr, n_param, error, settings, &
                    state)
        if (error /= 0) then
            write(stderr, *) "test_solver failed: Produced error from saddle point."
            test_solver = .false.
        end if
        call hartmann6d_gradient(context%vars, final_grad)
        if (norm2(final_grad) / sqrt(real(n_param, kind=rp)) > &
            default_settings%conv_tol) then
            write(stderr, *) "test_solver failed: Solver did not find stationary "// &
                "point from saddle point."
            test_solver = .false.
        end if
        if (any(abs(context%vars - minimum1) > var_thres) .and. &
            any(abs(context%vars - minimum2) > var_thres)) then
            write(stderr, *) "test_solver failed: Solver did not find minimum from "// &
                "saddle point."
            test_solver = .false.
        end if
        if (state%max_precision_reached) then
            write(stderr, *) "test_solver failed: Flagged that maximum precision "// &
                "was reached when convergence tolerance was met from saddle point."
            test_solver = .false.
        end if
        if (context%obj_func_zero_step) then
            write(stderr, *) "test_solver failed: Line search performed although "// &
                "it was not requested from saddle point."
            test_solver = .false.
        end if

        ! start at saddle point but allow only a single macro iteration, the internal
        ! stability check finds the saddle point unstable and the solver stops before
        ! solving a trust region subproblem, so that all Hessian linear transformations
        ! are the internal stability check's, whose addition to the solver's counter
        ! test_check_stationary_point covers, and check that a stale counter in the
        ! state record is reset so that it equals the internal stability check's
        context%vars = saddle_point
        settings = solver_settings_type()
        settings%n_macro = 1
        state%n_hess_x = 1000
        call setup_error_logging(settings, context)
        call solver(update_orbs_funptr, obj_func_funptr, n_param, error, settings, &
                    state)
        if (error /= error_solver_max_iter) then
            write(stderr, *) "test_solver failed: Did not return maximum iteration "// &
                "error code when only the internal stability check runs."
            test_solver = .false.
        end if
        if (state%n_hess_x /= state%stability_state%n_hess_x) then
            write(stderr, *) "test_solver failed: Stale counter of Hessian linear "// &
                "transformations not reset."
            test_solver = .false.
        end if

        ! run solver with settings which the sanity check rejects, nested settings for
        ! the Jacobi-Davidson method, which are only rejected once they inherit the
        ! solver's non-symmetric Hessian, which the level-shifted Davidson subsystem
        ! solver supports, and check that its error is returned with the solver's
        ! origin before the orbitals are updated
        context%vars = near_minimum
        settings = solver_settings_type()
        settings%subsystem_solver = "davidson_ls"
        settings%hess_symm = .false.
        settings%stability_settings%diag_solver = "jacobi-davidson"
        call setup_error_logging(settings, context)
        call solver(update_orbs_funptr, obj_func_funptr, n_param, error, settings, &
                    state)
        if (error /= error_solver + 1) then
            write(stderr, *) "test_solver failed: Did not return error of sanity check."
            test_solver = .false.
        end if
        if (state%n_update_orbs /= 0) then
            write(stderr, *) "test_solver failed: Orbitals updated before settings "// &
                "the sanity check rejects were refused."
            test_solver = .false.
        end if

        ! force non-convergence by allowing only two macro iterations from a generic
        ! starting point with a start trust radius below the default and check that the
        ! specific maximum iteration error code is returned and that the step taken
        ! after the first macro iteration stays within the start trust radius
        context%vars = distant_point
        settings = solver_settings_type()
        settings%n_macro = 2
        settings%start_trust_radius = 0.1_rp

        ! leave a stale internal stability check counter in the state record as a
        ! previous call would, no internal stability check is performed here so it has
        ! to be reset
        state%stability_state%n_hess_x = 1
        call setup_error_logging(settings, context)
        call solver(update_orbs_funptr, obj_func_funptr, n_param, error, settings, &
                    state)
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
        if (norm2(context%last_delta_vars) > settings%start_trust_radius + tol) then
            write(stderr, *) "test_solver failed: First step exceeds start trust "// &
                "radius."
            test_solver = .false.
        end if
        if (state%stability_state%n_hess_x /= 0) then
            write(stderr, *) "test_solver failed: Internal stability check counter "// &
                "of a previous call was not reset."
            test_solver = .false.
        end if

        ! force the maximum precision heuristic to trigger by requesting a convergence
        ! tolerance that floating-point noise in the gradient can never satisfy
        context%vars = near_minimum
        settings = solver_settings_type()
        settings%context => context
        settings%conv_tol = 0.0_rp
        settings%subsystem_solver = "tcg"

        ! run solver, check that it still returns without error while flagging this in
        ! the state and that the reported number of orbital updates agrees with the
        ! calls
        context%n_update_orbs_calls = 0
        call solver(update_orbs_funptr, obj_func_funptr, n_param, error, settings, &
                    state)
        if (error /= 0) then
            write(stderr, *) "test_solver failed: Produced error when forcing "// &
                "maximum precision heuristic."
            test_solver = .false.
        end if
        if (.not. state%max_precision_reached) then
            write(stderr, *) "test_solver failed: Did not flag that maximum "// &
                "precision was reached when convergence tolerance could not be met."
            test_solver = .false.
        end if
        if (.not. check_call_counts( &
            context, "solver", "when maximum precision is reached", &
            n_update_orbs=state%n_update_orbs)) test_solver = .false.

        ! run solver again on the same state record but stop it before convergence and
        ! check that the maximum precision flag is reset, that the reported number of
        ! orbital updates does not include those of the previous call and that no
        ! trust region subproblem is solved in the only macro iteration, whose step
        ! would never be evaluated
        context%vars = distant_point
        settings%conv_tol = default_settings%conv_tol
        settings%n_macro = 1
        context%n_update_orbs_calls = 0
        context%n_hess_x_calls = 0
        call setup_error_logging(settings, context)
        call solver(update_orbs_funptr, obj_func_funptr, n_param, error, settings, &
                    state)
        if (error /= error_solver_max_iter) then
            write(stderr, *) "test_solver failed: Did not return maximum iteration "// &
                "error code after maximum precision was reached in a previous call."
            test_solver = .false.
        end if
        if (state%max_precision_reached) then
            write(stderr, *) "test_solver failed: Flag that maximum precision was "// &
                "reached was not reset by the next call."
            test_solver = .false.
        end if
        if (.not. check_call_counts( &
            context, "solver", "on reused state record", n_hess_x=state%n_hess_x, &
            n_update_orbs=state%n_update_orbs)) test_solver = .false.
        if (state%n_hess_x /= 0) then
            write(stderr, *) "test_solver failed: Trust region subproblem solved "// &
                "in the last macro iteration."
            test_solver = .false.
        end if

        ! force the maximum precision heuristic to trigger again but let the
        ! convergence check pass in that same macro iteration, convergence then takes
        ! precedence and maximum precision is not flagged
        context%vars = near_minimum
        settings = solver_settings_type()
        settings%context => context
        settings%subsystem_solver = "tcg"
        settings%conv_check => mock_conv_check_without_update
        call solver(update_orbs_funptr, obj_func_funptr, n_param, error, settings, &
                    state)
        if (error /= 0) then
            write(stderr, *) "test_solver failed: Produced error when convergence "// &
                "check passes once maximum precision is reached."
            test_solver = .false.
        end if
        if (state%max_precision_reached) then
            write(stderr, *) "test_solver failed: Flagged that maximum precision "// &
                "was reached when convergence check passed in the same iteration."
            test_solver = .false.
        end if

        ! run solver with a convergence check which always passes, the solver stops at
        ! the first check without taking a step
        context%vars = near_minimum
        settings = solver_settings_type()
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

        ! run solver with a convergence threshold which the RMS gradient already
        ! satisfies and a convergence check of the gradient with the default threshold,
        ! which replaces the criterion on the RMS gradient, so that the solver does not
        ! stop at the first check but continues to the minimum
        context%vars = near_minimum
        settings = solver_settings_type()
        settings%context => context
        settings%conv_tol = 1.0_rp
        settings%conv_check => faulty_conv_check
        call solver(update_orbs_funptr, obj_func_funptr, n_param, error, settings)
        if (error /= 0) then
            write(stderr, *) "test_solver failed: Produced error when convergence "// &
                "check replaces gradient criterion."
            test_solver = .false.
        end if
        if (any(abs(context%vars - minimum1) > var_thres)) then
            write(stderr, *) "test_solver failed: Solver stopped by gradient "// &
                "criterion which convergence check replaces."
            test_solver = .false.
        end if

        ! run solver with an orbital update which succeeds but only provides a Hessian
        ! linear transformation in the first macro iteration and leaves the argument as
        ! it received it afterwards, which violates the interface, and check that the
        ! solver, which disassociates the argument before every orbital update, reports
        ! this as its own error in the second macro iteration
        context%vars = distant_point
        context%n_update_orbs_calls = 0
        update_orbs_funptr => update_orbs_hess_x_once
        settings = solver_settings_type()
        call setup_error_logging(settings, context)
        call solver(update_orbs_funptr, obj_func_funptr, n_param, error, settings)
        if (error /= error_solver + 1) then
            write(stderr, *) "test_solver failed: Did not report a missing Hessian "// &
                "linear transformation."
            test_solver = .false.
        end if
        if (len_trim(context%log_message) == 0) then
            write(stderr, *) "test_solver failed: No error message printed for "// &
                "missing Hessian linear transformation."
            test_solver = .false.
        end if
        if (context%n_update_orbs_calls /= 2) then
            write(stderr, *) "test_solver failed: Missing Hessian linear "// &
                "transformation not reported in the second macro iteration."
            test_solver = .false.
        end if

        ! run solver from the saddle point of the Hartmann 6D function, whose path
        ! depends on the settings, with the default settings and without settings and
        ! check that both take exactly the same steps to the same point, the callback
        ! functions reach the variables through the module context since the solver
        ! hands them no host context without settings
        update_orbs_funptr => settingless_update_orbs
        obj_func_funptr => settingless_obj_func
        settingless_context%vars = saddle_point
        call solver(update_orbs_funptr, obj_func_funptr, n_param, error, &
                    default_settings, explicit_state)
        if (error /= 0) then
            write(stderr, *) "test_solver failed: Produced error with default settings."
            test_solver = .false.
        end if
        final_vars = settingless_context%vars
        settingless_context%vars = saddle_point
        call solver(update_orbs_funptr, obj_func_funptr, n_param, error, state=state)
        if (error /= 0) then
            write(stderr, *) "test_solver failed: Produced error without settings."
            test_solver = .false.
        end if
        if (any(abs(settingless_context%vars - final_vars) > 0.0_rp) .or. &
            state%n_update_orbs /= explicit_state%n_update_orbs .or. &
            state%n_obj_func /= explicit_state%n_obj_func .or. &
            state%n_hess_x /= explicit_state%n_hess_x .or. &
            state%stability_state%n_hess_x /= explicit_state%stability_state%n_hess_x) &
            then
            write(stderr, *) "test_solver failed: Did not use the default settings "// &
                "without settings."
            test_solver = .false.
        end if
        update_orbs_funptr => update_orbs
        obj_func_funptr => obj_func

        ! check that the line search brackets the step starting from a vanishing step
        context%vars = near_minimum
        settings = solver_settings_type()
        settings%context => context
        settings%line_search = .true.
        context%obj_func_zero_step = .false.
        call solver(update_orbs_funptr, obj_func_funptr, n_param, error, settings)
        if (error /= 0) then
            write(stderr, *) "test_solver failed: Produced error with line search."
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
        settings%subsystem_solver = "jacobi-davidson_ls"
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
        ! Davidson solvers, whose first trial vector is the normalized gradient, which
        ! the gradient perturbation would hide
        context%vars = near_minimum
        settings = solver_settings_type()
        settings%context => context
        settings%subsystem_solver = "tcg"
        settings%grad_noise = .false.
        context%hess_x_normalized_grad = .false.
        call solver(update_orbs_funptr, obj_func_funptr, n_param, error, settings)
        if (error /= 0) then
            write(stderr, *) "test_solver failed: Produced error with truncated "// &
                "conjugate gradient."
            test_solver = .false.
        end if
        if (context%hess_x_normalized_grad) then
            write(stderr, *) "test_solver failed: Davidson solver used instead of "// &
                "truncated conjugate gradient."
            test_solver = .false.
        end if

        ! perturb the gradient with reduction factors of the local and the global
        ! region orders of magnitude apart and start near the minimum, where the Newton
        ! step lies in the local region, and stop after the first subproblem, which
        ! counts as one in the global region and is therefore solved only to its
        ! gradient perturbation, a tenth of the gradient norm, so that its relative
        ! residual lies above the reduction factor of the local region
        context%vars = near_minimum
        settings = solver_settings_type()
        call setup_error_logging(settings, context)
        settings%grad_noise = .true.
        settings%local_red_factor = 1e-8_rp
        settings%global_red_factor = 1e-1_rp
        settings%n_macro = 2
        call solver(update_orbs_funptr, obj_func_funptr, n_param, error, settings, &
                    state)
        if (error /= error_solver_max_iter) then
            write(stderr, *) "test_solver failed: Did not stop after the first "// &
                "subproblem with gradient perturbation."
            test_solver = .false.
        end if
        if (state%relative_residual <= settings%local_red_factor) then
            write(stderr, *) "test_solver failed: First subproblem with gradient "// &
                "perturbation not solved for the global region."
            test_solver = .false.
        end if

        ! let the same optimization converge, the subproblems after the first one,
        ! which follow a Newton step of the local region, are perturbed and solved for
        ! it, so that the relative residual of the last one lies below the reduction
        ! factor of the local region
        context%vars = near_minimum
        settings%n_macro = default_settings%n_macro
        call solver(update_orbs_funptr, obj_func_funptr, n_param, error, settings, &
                    state)
        if (error /= 0) then
            write(stderr, *) "test_solver failed: Produced error with gradient "// &
                "perturbation."
            test_solver = .false.
        end if
        if (state%n_update_orbs < 3) then
            write(stderr, *) "test_solver failed: Converged with gradient "// &
                "perturbation before a second subproblem was solved."
            test_solver = .false.
        end if
        if (state%relative_residual > settings%local_red_factor) then
            write(stderr, *) "test_solver failed: Subproblem with gradient "// &
                "perturbation not solved for the local region of the previous one."
            test_solver = .false.
        end if

        ! perform a line search for an objective function whose minimum along the step
        ! lies at twice the Newton step of the model, which the gradient perturbation
        ! would move, and check that the step and its size in the state are scaled by
        ! the multiplier the bracketing returns, the solver is stopped by the maximum
        ! number of macro iterations once the first step is taken
        call setup_settings(settings, half_step_context)
        settings%line_search = .true.
        settings%grad_noise = .false.
        settings%n_macro = 2
        call random_number(start_vars)
        start_vars = 2.0_rp * start_vars - 1.0_rp
        half_step_context%vars = start_vars
        update_orbs_funptr => half_step_update_orbs
        obj_func_funptr => half_step_obj_func
        call solver(update_orbs_funptr, obj_func_funptr, 2_ip, error, settings, state)
        if (error /= error_solver_max_iter) then
            write(stderr, *) "test_solver failed: Did not stop at maximum number "// &
                "of macro iterations for line search."
            test_solver = .false.
        end if
        if (any(abs(half_step_context%first_step + 2.0_rp * start_vars) > tol)) then
            write(stderr, *) "test_solver failed: Step not scaled by multiplier of "// &
                "line search."
            test_solver = .false.
        end if
        if (abs(state%step_size - norm2(half_step_context%first_step)) > tol) then
            write(stderr, *) "test_solver failed: Step size in the state not "// &
                "scaled by multiplier of line search."
            test_solver = .false.
        end if

        ! minimize the double well function x^2 - y^2 + y^4 from (0.5, 0) with every
        ! subsystem solver and with a state record, which the orbital update reads
        ! during the optimization through the pointer in its context, and without
        ! random or extra trial vectors and gradient perturbation, so that the solver
        ! first converges to the saddle point at the origin, which the requested
        ! stability check finds unstable, and then continues to a minimum
        update_orbs_funptr => state_reading_update_orbs
        obj_func_funptr => double_well_obj_func
        state_context%state => state
        do i_solver = 1, size(subsystem_solver_options)
            settings = solver_settings_type()
            settings%subsystem_solver = subsystem_solver_options(i_solver)
            settings%stability = .true.
            settings%grad_noise = .false.
            settings%n_random_trial_vectors = 0
            settings%n_extra_trial_vectors = 0
            settings%start_trust_radius = 0.3_rp
            settings%context => state_context
            case_name = trim(subsystem_solver_options(i_solver))//" subsystem solver"
            state_context%vars = [0.5_rp, 0.0_rp]
            state_context%n_obj_func_calls = 0
            if (allocated(state_context%seen_states)) &
                deallocate(state_context%seen_states, state_context%step_norms, &
                           state_context%precond_step_norms)
            call solver(update_orbs_funptr, obj_func_funptr, 2_ip, error, settings, &
                        state)
            if (error /= 0) then
                write(stderr, *) "test_solver failed: Produced error with state "// &
                    "record for "//trim(case_name)//"."
                test_solver = .false.
            end if
            if (abs(abs(state_context%vars(2)) - 1.0_rp / sqrt(2.0_rp)) > var_thres) &
                then
                write(stderr, *) "test_solver failed: Did not continue from the "// &
                    "saddle point of the double well function to a minimum for "// &
                    trim(case_name)//"."
                test_solver = .false.
            end if

            ! check the state the solver returns, at the converged point, against the
            ! orbital updates the orbital update recorded
            n_updates = 0
            if (allocated(state_context%seen_states)) &
                n_updates = size(state_context%seen_states)
            if (state%macro_iteration /= n_updates - 1 .or. &
                abs(state%func - double_well_func(state_context%vars)) > tol .or. &
                abs(state%grad_rms - norm2([ &
                    2.0_rp * state_context%vars(1), -2.0_rp * state_context%vars(2) + &
                    4.0_rp * state_context%vars(2)**3]) / sqrt(2.0_rp)) > tol .or. &
                state%trust_radius <= 0.0_rp) then
                write(stderr, *) "test_solver failed: Returned state not that of "// &
                    "the converged point for "//trim(case_name)//"."
                test_solver = .false.
            end if
            if (state%n_update_orbs /= n_updates .or. &
                state%n_obj_func /= state_context%n_obj_func_calls) then
                write(stderr, *) "test_solver failed: Counters in the returned "// &
                    "state wrong for "//trim(case_name)//"."
                test_solver = .false.
            end if
            if (state%max_precision_reached) then
                write(stderr, *) "test_solver failed: Maximum precision flagged in "// &
                    "the returned state for "//trim(case_name)//"."
                test_solver = .false.
            end if

            ! check the state the orbital update read during the optimization: the
            ! first update sees the starting trust radius before anything was
            ! evaluated, every later one the previous macro iteration and the step it
            ! is given, whose size is measured in the metric of the trust region, the
            ! Euclidean norm for (Jacobi-)Davidson and the norm of the absolute
            ! diagonal preconditioner of the previous point for truncated conjugate
            ! gradient, while the step along an unstable mode, which needs no micro
            ! iterations, is measured in the Euclidean norm, the second one the level
            ! shift of the first step, which the trust radius truncates for
            ! (Jacobi-)Davidson and which vanishes for truncated conjugate gradient,
            ! the level shift never being negative or a negative zero, and one update
            ! the step along the unstable mode, without micro iterations, level shift
            ! and residual, after steps of the subsystem solver
            if (n_updates < 3) then
                write(stderr, *) "test_solver failed: Too few orbital updates to "// &
                    "check the state during the optimization for "//trim(case_name)//"."
                test_solver = .false.
                cycle
            end if
            associate (seen => state_context%seen_states)
                step_norms = state_context%step_norms
                if (settings%subsystem_solver == "tcg") &
                    where (seen%micro_iteration > 0) step_norms = &
                    state_context%precond_step_norms
                if (seen(1)%macro_iteration /= 0 .or. seen(1)%n_update_orbs /= 0 .or. &
                    abs(seen(1)%trust_radius - settings%start_trust_radius) > tol .or. &
                    seen(1)%micro_iteration /= 0 .or. abs(seen(1)%step_size) > tol) then
                    write(stderr, *) "test_solver failed: State wrong before the "// &
                        "starting point is evaluated for "//trim(case_name)//"."
                    test_solver = .false.
                end if
                if (any([(seen(i_update)%macro_iteration /= i_update - 2 .or. &
                          seen(i_update)%n_update_orbs /= i_update - 1, i_update=2, &
                          n_updates)])) then
                    write(stderr, *) "test_solver failed: Macro iteration or "// &
                        "number of orbital updates in the state wrong during the "// &
                        "optimization for "//trim(case_name)//"."
                    test_solver = .false.
                end if
                if (any([(abs(seen(i_update)%step_size - step_norms(i_update)) > tol, &
                          i_update=2, n_updates)])) then
                    write(stderr, *) "test_solver failed: Step size in the state "// &
                        "wrong during the optimization for "//trim(case_name)//"."
                    test_solver = .false.
                end if
                if (settings%subsystem_solver == "tcg") then
                    if (any([(abs(seen(i_update)%level_shift) > 0.0_rp, i_update=2, &
                              n_updates)])) then
                        write(stderr, *) "test_solver failed: Level shift in the "// &
                            "state does not vanish for "//trim(case_name)//"."
                        test_solver = .false.
                    end if
                else if (seen(2)%level_shift <= 0.0_rp) then
                    write(stderr, *) "test_solver failed: Level shift of the step "// &
                        "the trust radius truncates not in the state for "// &
                        trim(case_name)//"."
                    test_solver = .false.
                end if
                if (any([(sign(1.0_rp, seen(i_update)%level_shift) < 0.0_rp, &
                          i_update=2, n_updates)])) then
                    write(stderr, *) "test_solver failed: Negative level shift or "// &
                        "negative zero in the state for "//trim(case_name)//"."
                    test_solver = .false.
                end if
                if (.not. any([(seen(i_update)%micro_iteration == 0 .and. &
                                seen(i_update)%jacobi_davidson_iteration == 0 .and. &
                                abs(seen(i_update)%level_shift) < tol .and. &
                                abs(seen(i_update)%relative_residual) < tol .and. &
                                seen(i_update)%step_size > 0.0_rp .and. &
                                seen(i_update - 1)%micro_iteration > 0, i_update=3, &
                                n_updates)])) then
                    write(stderr, *) "test_solver failed: Step along the unstable "// &
                        "mode after steps of the subsystem solver not in the state "// &
                        "for "//trim(case_name)//"."
                    test_solver = .false.
                end if
            end associate
        end do

        ! minimize the double well function from (0.5, 0) as above, but without the
        ! stability check, so that the solver reaches the saddle point at an exactly
        ! vanishing gradient, once with a convergence check of the change of the
        ! objective function, which can not pass at a stationary point the solver has
        ! just reached, so that the solver takes the saddle point as a point at maximum
        ! precision and resets the flag once it continues, and which would pass at the
        ! point the step along the unstable mode reaches, since this step changes the
        ! objective function by far less than the tolerance of the check for a small
        ! convergence threshold, and once with the criterion on the RMS gradient, which
        ! passes there, and check that the solver checks the stability of the saddle
        ! point nevertheless, continues along the unstable mode and finally stops at a
        ! minimum
        update_orbs_funptr => double_well_update_orbs
        obj_func_funptr => double_well_obj_func
        do i_case = 1, 2
            settings = solver_settings_type()
            settings%grad_noise = .false.
            settings%n_random_trial_vectors = 0
            settings%n_extra_trial_vectors = 0
            settings%start_trust_radius = 0.3_rp
            settings%context => double_well_context
            if (i_case == 1) then
                settings%conv_check => double_well_conv_check
                settings%conv_tol = 1e-7_rp
                case_name = "with convergence check of the function change"
            else
                case_name = "with gradient criterion"
            end if
            double_well_context%vars = [0.5_rp, 0.0_rp]
            double_well_context%last_conv_check_func = huge(1.0_rp)
            call solver(update_orbs_funptr, obj_func_funptr, 2_ip, error, settings)
            if (error /= 0) then
                write(stderr, *) "test_solver failed: Produced error at exactly "// &
                    "vanishing gradient "//trim(case_name)//"."
                test_solver = .false.
            end if
            if (abs(abs(double_well_context%vars(2)) - 1.0_rp / sqrt(2.0_rp)) > &
                var_thres) then
                write(stderr, *) "test_solver failed: Did not continue from the "// &
                    "saddle point of the double well function at an exactly "// &
                    "vanishing gradient to a minimum "//trim(case_name)//"."
                test_solver = .false.
            end if
        end do

        ! start at the minimum of the quadratic model 1/2 |x|^2 at the origin, where
        ! the gradient vanishes exactly, with a vanishing convergence threshold, which
        ! the gradient can therefore not pass, and check that the solver stops there
        ! without an error and flags that maximum precision was reached, since the
        ! stability check finds the point stable
        call setup_settings(settings, half_step_context)
        settings%conv_tol = 0.0_rp
        half_step_context%vars = 0.0_rp
        update_orbs_funptr => half_step_update_orbs
        obj_func_funptr => half_step_obj_func
        call solver(update_orbs_funptr, obj_func_funptr, 2_ip, error, settings, state)
        if (error /= 0) then
            write(stderr, *) "test_solver failed: Produced error at a stable point "// &
                "with an exactly vanishing gradient."
            test_solver = .false.
        end if
        if (.not. state%max_precision_reached) then
            write(stderr, *) "test_solver failed: Maximum precision not flagged at "// &
                "a stable point with an exactly vanishing gradient which can not "// &
                "pass the convergence check."
            test_solver = .false.
        end if

        ! inject faults by running the solver from the saddle point, where it performs
        ! the internal stability check, with a stability check at convergence, a line
        ! search and every optional callback function, also of the nested settings,
        ! once without a fault to count the calls of every callback function and then
        ! letting every one of these calls fail in turn, for every subsystem solver
        ! without and with a perturbed gradient, whose perturbation is projected, and
        ! check that the unfaulted run reaches every callback function, that the error
        ! is reported with the origin of the failing callback function and that the
        ! reported counters still agree with the calls
        do i_solver = 1, size(subsystem_solver_options)
            do i_case = 1, 2
                grad_noise = i_case == 2
                case_name = trim(subsystem_solver_options(i_solver))// &
                            " subsystem solver"
                if (grad_noise) case_name = trim(case_name)//" and perturbed gradient"
                call run_with_fault(0_ip, 0_ip)
                if (error /= 0) then
                    write(stderr, *) "test_solver failed: Produced error without "// &
                        "fault with "//trim(case_name)//"."
                    test_solver = .false.
                    cycle
                end if
                n_calls = fault_context%n_callback_calls
                do callback = 1, n_fault_callbacks
                    if (n_calls(callback) == 0) then
                        write(stderr, *) "test_solver failed: Unfaulted run did "// &
                            "not call "//trim(fault_callback_names(callback))// &
                            " with "//trim(case_name)//"."
                        test_solver = .false.
                    end if
                    do i_call = 1, n_calls(callback)
                        call run_with_fault(callback, i_call)
                        write(call_number, '(I0)') i_call
                        if (error /= fault_origins(callback) + 1) then
                            write(stderr, *) "test_solver failed: Error of failing "// &
                                trim(fault_callback_names(callback))//" at call "// &
                                trim(call_number)//" not reported with its origin "// &
                                "with "//trim(case_name)//"."
                            test_solver = .false.
                        end if
                        if (.not. check_call_counts( &
                            fault_context, "solver", &
                            "for failing "//trim(fault_callback_names(callback))// &
                            " at call "//trim(call_number)//" with "//trim(case_name), &
                            n_hess_x=state%n_hess_x, &
                            n_update_orbs=state%n_update_orbs, &
                            n_obj_func=state%n_obj_func)) then
                            test_solver = .false.
                        end if
                    end do
                end do
            end do
        end do

    contains

        subroutine run_with_fault(fault_callback, fault_call)
            !
            ! this subroutine runs the solver from the saddle point with the subsystem
            ! solver and the gradient perturbation of the current case and lets the
            ! given call of the given callback function fail
            !
            integer(ip), intent(in) :: fault_callback, fault_call

            procedure(update_orbs_type), pointer :: faulty_update_orbs_funptr
            procedure(obj_func_type), pointer :: faulty_obj_func_funptr

            ! set up fault injection context
            fault_context%vars = saddle_point
            fault_context%n_callback_calls = 0
            fault_context%n_update_orbs_calls = 0
            fault_context%n_obj_func_calls = 0
            fault_context%n_hess_x_calls = 0
            fault_context%fault_callback = fault_callback
            fault_context%fault_call = fault_call
            fault_context%stability_state => state%stability_state

            ! set up settings with every optional callback function, also of the nested
            ! settings
            settings = solver_settings_type()
            settings%context => fault_context
            settings%subsystem_solver = subsystem_solver_options(i_solver)
            settings%jacobi_davidson_start = 0
            settings%line_search = .true.
            settings%stability = .true.
            settings%grad_noise = grad_noise
            settings%precond => faulty_precond
            settings%project => faulty_project
            settings%conv_check => faulty_conv_check
            settings%get_extra_trial_vectors => faulty_get_extra_trial_vectors
            settings%stability_settings%conv_check => faulty_stability_conv_check
            settings%stability_settings%get_extra_trial_vectors => &
                faulty_get_extra_trial_vectors

            ! run solver
            faulty_update_orbs_funptr => faulty_update_orbs
            faulty_obj_func_funptr => faulty_obj_func
            call solver(faulty_update_orbs_funptr, faulty_obj_func_funptr, n_param, &
                        error, settings, state)

        end subroutine run_with_fault

    end function test_solver

    logical(c_bool) function test_stability_check() bind(C)
        !
        ! this function tests the stability check subroutine
        !
        use opentrustregion, only: &
            hess_x_type, stability_settings_type, stability_check, &
            error_stability_check, error_stability_check_max_iter, verbosity_warning, &
            error_hess_x, error_precond, error_project, error_conv_check, &
            error_get_extra_trial_vectors, diag_solver_options, unstable_warning_msg, &
            default_stability_settings
        use test_reference, only: arm_host_context, host_context_reached, &
                                  n_extra_trial_vectors

        real(rp), parameter :: points(n_param, 2) = &
            reshape([minimum1, saddle_point], [n_param, 2_ip])
        character(len=*), parameter :: point_names(2) = &
            [character(len=12) :: "minimum", "saddle point"]
        integer(ip), parameter :: &
            fault_callbacks(5) = [fault_hess_x, fault_precond, fault_project, &
                                  fault_get_extra_trial_vectors, &
                                  fault_stability_conv_check], &
            fault_origins(5) = [error_hess_x, error_precond, error_project, &
                                error_get_extra_trial_vectors, error_conv_check]
        real(rp) :: h_diag(n_param), direction(n_param), hess_eigvals(n_param), &
                    hess_eigvecs(n_param, n_param)
        procedure(hess_x_type), pointer :: hess_x_funptr
        logical :: stable, expected_calls(n_fault_callbacks)
        integer(ip) :: error, i, i_point, i_solver, callback, i_call, &
                       n_calls(n_fault_callbacks)
        character(len=300) :: msg
        character(len=20) :: call_number
        character(len=100) :: case_name
        type(stability_settings_type) :: settings
        type(stability_state_type), target :: state
        type(stability_state_type) :: explicit_state
        type(hartmann6d_fault_context_type), target :: fault_context
        type(hartmann6d_context_type), target :: context

        ! assume tests pass
        test_stability_check = .true.

        ! run the stability check at the minimum and at the saddle point and check that
        ! it classifies the point correctly and returns the lowest Hessian eigenvalue
        ! and its eigenvector, which are determined independently, at the minimum as
        ! well, that the state holds the iteration, this eigenvalue and its converged
        ! residual and that the reported number of Hessian linear transformations
        ! agrees with the calls, the state record is reused so that the counter must
        ! not include the previous calls
        hess_x_funptr => hess_x_fun
        settings = stability_settings_type()
        settings%context => context
        do i_point = 1, size(points, 2)
            context%hess = hartmann6d_hessian(points(:, i_point))
            call ref_symm_mat_diag(context%hess, hess_eigvals, hess_eigvecs)
            h_diag = [(context%hess(i, i), i=1, size(h_diag))]
            case_name = point_names(i_point)
            context%n_hess_x_calls = 0
            call stability_check(h_diag, hess_x_funptr, stable, error, settings, &
                                 state, direction)
            if (error /= 0) then
                write(stderr, *) "test_stability_check failed: Produced error for "// &
                    trim(case_name)//"."
                test_stability_check = .false.
            end if
            if (stable .neqv. i_point == 1) then
                write(stderr, *) "test_stability_check failed: Stability "// &
                    "incorrectly classified for "//trim(case_name)//"."
                test_stability_check = .false.
            end if
            if (abs(state%min_eigval - hess_eigvals(1)) > tol) then
                write(stderr, *) "test_stability_check failed: Lowest eigenvalue "// &
                    "not held by the state for "//trim(case_name)//"."
                test_stability_check = .false.
            end if
            if (abs(abs(dot_product(direction, hess_eigvecs(:, 1))) - 1.0_rp) > tol) &
                then
                write(stderr, *) "test_stability_check failed: Eigenvector of the "// &
                    "lowest eigenvalue not returned for "//trim(case_name)//"."
                test_stability_check = .false.
            end if
            if (state%iteration < 1) then
                write(stderr, *) "test_stability_check failed: Iteration not in "// &
                    "the state for "//trim(case_name)//"."
                test_stability_check = .false.
            end if
            if (state%residual_rms >= settings%conv_tol) then
                write(stderr, *) "test_stability_check failed: Converged residual "// &
                    "of the lowest eigenvalue not in the state for "// &
                    trim(case_name)//"."
                test_stability_check = .false.
            end if
            if (.not. check_call_counts( &
                context, "stability_check", "for "//trim(case_name), &
                n_hess_x=state%n_hess_x)) test_stability_check = .false.
        end do

        ! check that a warning with the lowest eigenvalue is printed at the saddle point
        settings = stability_settings_type()
        settings%context => context
        settings%logger => logger
        settings%verbose = verbosity_warning
        context%log_message = ""
        call stability_check(h_diag, hess_x_funptr, stable, error, settings)
        if (error /= 0) then
            write(stderr, *) "test_stability_check failed: Produced error for "// &
                "saddle point with logging function."
            test_stability_check = .false.
        end if
        write(msg, '(A, F0.4)') unstable_warning_msg, hess_eigvals(1)
        if (index(context%log_message, " "//trim(msg)) == 0) then
            write(stderr, *) "test_stability_check failed: Warning not printed for "// &
                "saddle point."
            test_stability_check = .false.
        end if

        ! check that the Hessian linear transformation, the preconditioner, the
        ! projection, the extra trial vector function, the convergence check and the
        ! logging function receive the host context of the settings and that the
        ! convergence check reads the final state in the caller's record
        context%stability_state => state
        call setup_settings(settings, context)
        settings%precond => identity_precond
        settings%project => identity_project
        settings%conv_check => mock_stability_conv_check
        settings%n_extra_trial_vectors = n_extra_trial_vectors
        settings%get_extra_trial_vectors => mock_get_extra_trial_vectors
        call arm_host_context(settings, context)
        context%seen_residual_rms = -1.0_rp
        call stability_check(h_diag, hess_x_funptr, stable, error, settings, state, &
                             direction)
        if (error /= 0) then
            write(stderr, *) "test_stability_check failed: Produced error for "// &
                "saddle point with every callback function."
            test_stability_check = .false.
        end if
        if (.not. host_context_reached("stability_check", context)) &
            test_stability_check = .false.
        if (abs(context%seen_residual_rms - state%residual_rms) > tol * &
            state%residual_rms .or. state%residual_rms <= 0.0_rp) then
            write(stderr, *) "test_stability_check failed: Convergence check did "// &
                "not read the final state in the caller's record."
            test_stability_check = .false.
        end if
        context%stability_state => null()

        ! check that the saddle point is found to be unstable without a returned
        ! direction and eigenvalue
        call setup_settings(settings, context)
        call stability_check(h_diag, hess_x_funptr, stable, error, settings)
        if (error /= 0) then
            write(stderr, *) "test_stability_check failed: Produced error for "// &
                "saddle point without returned direction."
            test_stability_check = .false.
        end if
        if (stable) then
            write(stderr, *) "test_stability_check failed: Stability check "// &
                "incorrectly classifies stability of saddle point without returned "// &
                "direction."
            test_stability_check = .false.
        end if

        ! check that a stability threshold below the lowest eigenvalue counts the
        ! saddle point as stable
        settings%stability_thresh = 2.0_rp * hess_eigvals(1)
        call stability_check(h_diag, hess_x_funptr, stable, error, settings)
        if (error /= 0) then
            write(stderr, *) "test_stability_check failed: Produced error for "// &
                "saddle point with stability threshold below its lowest eigenvalue."
            test_stability_check = .false.
        end if
        if (.not. stable) then
            write(stderr, *) "test_stability_check failed: Saddle point not stable "// &
                "with stability threshold below its lowest eigenvalue."
            test_stability_check = .false.
        end if

        ! run stability check with settings which the sanity check rejects and check
        ! that its error is returned with the stability check's origin
        settings = stability_settings_type()
        settings%diag_solver = "unknown"
        call setup_error_logging(settings, context)
        call stability_check(h_diag, hess_x_funptr, stable, error, settings, &
                             kappa=direction)
        if (error /= error_stability_check + 1) then
            write(stderr, *) "test_stability_check failed: Did not return error of "// &
                "sanity check."
            test_stability_check = .false.
        end if

        ! force non-convergence by allowing only a single iteration and check that the
        ! specific maximum iteration error code is returned and that the state holds
        ! this iteration and its unconverged residual after the error
        settings = stability_settings_type()
        settings%n_iter = 1
        call setup_error_logging(settings, context)
        call stability_check(h_diag, hess_x_funptr, stable, error, settings, state)
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
        if (state%iteration /= 1 .or. state%residual_rms < settings%conv_tol) then
            write(stderr, *) "test_stability_check failed: Iteration or "// &
                "unconverged residual not in the state when exceeding n_iter."
            test_stability_check = .false.
        end if

        ! force the reduced space to grow until it spans the full parameter space by
        ! setting an unreachable convergence tolerance and check that its exact
        ! eigenpair is returned without error
        settings = stability_settings_type()
        settings%context => context
        settings%conv_tol = 0.0_rp
        call stability_check(h_diag, hess_x_funptr, stable, error, settings, &
                             kappa=direction)
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

        ! run stability check at saddle point with the Jacobi-Davidson method and a
        ! slightly asymmetric Hessian linear transformation, whose products are
        ! recalculated once the symmetry is found to be violated, and check that these
        ! are counted as well, the antisymmetric contribution leaves a residual of its
        ! own size, so the convergence threshold has to lie above it
        hess_x_funptr => hess_x_fun_asymmetric
        settings = stability_settings_type()
        settings%context => context
        settings%diag_solver = "jacobi-davidson"
        settings%jacobi_davidson_start = 1
        settings%conv_tol = 1e-6_rp
        context%n_hess_x_calls = 0
        call stability_check(h_diag, hess_x_funptr, stable, error, settings, state, &
                             direction)
        if (error /= 0) then
            write(stderr, *) "test_stability_check failed: Produced error with "// &
                "asymmetric Hessian linear transformation."
            test_stability_check = .false.
        end if
        if (stable) then
            write(stderr, *) "test_stability_check failed: Stability check "// &
                "incorrectly classifies stability of saddle point with asymmetric "// &
                "Hessian linear transformation."
            test_stability_check = .false.
        end if
        if (.not. check_call_counts( &
            context, "stability_check", "with asymmetric Hessian linear "// &
            "transformation", n_hess_x=state%n_hess_x)) test_stability_check = .false.
        hess_x_funptr => hess_x_fun

        ! run stability check at the saddle point, whose iterations depend on the
        ! settings, with the default settings and without settings and check that both
        ! reach exactly the same eigenvalue and residual with the same Hessian linear
        ! transformations, the Hessian linear transformation reaches the Hessian
        ! through the module context since the stability check hands it no host context
        ! without settings
        settingless_context%hess = context%hess
        hess_x_funptr => settingless_hess_x
        call stability_check(h_diag, hess_x_funptr, stable, error, &
                             default_stability_settings, explicit_state)
        if (error /= 0) then
            write(stderr, *) "test_stability_check failed: Produced error with "// &
                "default settings."
            test_stability_check = .false.
        end if
        call stability_check(h_diag, hess_x_funptr, stable, error, state=state)
        if (error /= 0) then
            write(stderr, *) "test_stability_check failed: Produced error without "// &
                "settings."
            test_stability_check = .false.
        end if
        if (state%n_hess_x /= explicit_state%n_hess_x .or. &
            state%iteration /= explicit_state%iteration .or. &
            abs(state%min_eigval - explicit_state%min_eigval) > 0.0_rp .or. &
            abs(state%residual_rms - explicit_state%residual_rms) > 0.0_rp) then
            write(stderr, *) "test_stability_check failed: Did not use the default "// &
                "settings without settings."
            test_stability_check = .false.
        end if

        ! start at saddle point with a preconditioner which always returns the first
        ! unit vector and a Hessian diagonal whose minimum lies at the first element,
        ! so that no new trial vector can be added before the eigenpair converges, and
        ! check that the maximum iteration error code is returned before the
        ! iterations are used up and that an error message is printed
        h_diag = 1.0_rp
        h_diag(1) = 0.0_rp
        hess_x_funptr => hess_x_fun
        settings = stability_settings_type()
        settings%precond => mock_precond_first_unit_vector
        call setup_error_logging(settings, context)
        call stability_check(h_diag, hess_x_funptr, stable, error, settings, state)
        if (error /= error_stability_check_max_iter) then
            write(stderr, *) "test_stability_check failed: Did not return maximum "// &
                "iteration error code when trial space can not be expanded."
            test_stability_check = .false.
        end if
        if (len_trim(context%log_message) == 0) then
            write(stderr, *) "test_stability_check failed: No error message "// &
                "printed when trial space can not be expanded."
            test_stability_check = .false.
        end if
        if (state%iteration >= settings%n_iter) then
            write(stderr, *) "test_stability_check failed: Iterations used up when "// &
                "trial space can not be expanded."
            test_stability_check = .false.
        end if

        ! run the stability check with the Hessian treated as non-symmetric for an
        ! upper triangular Hessian, whose eigenvalues are its diagonal elements, one
        ! apart from each other and lowest at the last element, so that the
        ! corresponding eigenvector has no zero component, and check that the lowest
        ! eigenvalue and a right eigenvector, whose residual lies below the convergence
        ! threshold, are returned, which only a diagonalization of the non-symmetric
        ! reduced space yields
        call random_number(context%hess)
        do i = 1, n_param
            context%hess(i + 1:, i) = 0.0_rp
            context%hess(i, i) = real(n_param - i, kind=rp) - 1.5_rp
        end do
        h_diag = [(context%hess(i, i), i=1, size(h_diag))]
        hess_x_funptr => hess_x_fun
        settings = stability_settings_type()
        settings%context => context
        settings%hess_symm = .false.
        call stability_check(h_diag, hess_x_funptr, stable, error, settings, state, &
                             direction)
        if (error /= 0) then
            write(stderr, *) "test_stability_check failed: Produced error for "// &
                "non-symmetric Hessian."
            test_stability_check = .false.
        end if
        if (abs(state%min_eigval - minval(h_diag)) > 0.5_rp .or. &
            norm2(matmul(context%hess, direction) - state%min_eigval * direction) / &
            sqrt(real(n_param, kind=rp)) >= settings%conv_tol) then
            write(stderr, *) "test_stability_check failed: Lowest eigenpair of "// &
                "non-symmetric Hessian not returned."
            test_stability_check = .false.
        end if

        ! inject faults by running the stability check at the saddle point with every
        ! optional callback function, once without a fault to count the calls of every
        ! callback function and then letting every one of these calls fail in turn, for
        ! every diagonalization solver, and check that the unfaulted run reaches
        ! exactly the callback functions the case uses, that the error is reported with
        ! the origin of the failing callback function and that the reported counter
        ! still agrees with the calls
        fault_context%hess = hartmann6d_hessian(saddle_point)
        h_diag = [(fault_context%hess(i, i), i=1, size(h_diag))]
        hess_x_funptr => faulty_hess_x
        do i_solver = 1, size(diag_solver_options)
            case_name = trim(diag_solver_options(i_solver))//" diagonalization solver"
            call run_with_fault(0_ip, 0_ip)
            if (error /= 0) then
                write(stderr, *) "test_stability_check failed: Produced error "// &
                    "without fault for "//trim(case_name)//"."
                test_stability_check = .false.
                cycle
            end if
            n_calls = fault_context%n_callback_calls

            ! the unfaulted run reaches every callback function of the stability check
            ! and the preconditioner, which the Jacobi-Davidson method started at the
            ! first iteration replaces, only with the Davidson method
            expected_calls = .false.
            expected_calls(fault_callbacks) = .true.
            expected_calls(fault_precond) = diag_solver_options(i_solver) == "davidson"
            do callback = 1, n_fault_callbacks
                if ((n_calls(callback) > 0) .neqv. expected_calls(callback)) then
                    write(stderr, *) "test_stability_check failed: Unfaulted run "// &
                        trim(merge("did not call", "called      ", &
                                   expected_calls(callback)))//" "// &
                        trim(fault_callback_names(callback))//" for "// &
                        trim(case_name)//"."
                    test_stability_check = .false.
                end if
            end do

            ! let every call fail in turn
            do callback = 1, size(fault_callbacks)
                do i_call = 1, n_calls(fault_callbacks(callback))
                    call run_with_fault(fault_callbacks(callback), i_call)
                    write(call_number, '(I0)') i_call
                    if (error /= fault_origins(callback) + 1) then
                        write(stderr, *) "test_stability_check failed: Error of "// &
                            "failing "// &
                            trim(fault_callback_names(fault_callbacks(callback)))// &
                            " at call "//trim(call_number)// &
                            " not reported with its origin for "//trim(case_name)//"."
                        test_stability_check = .false.
                    end if
                    if (.not. check_call_counts( &
                        fault_context, "stability_check", "for failing "// &
                        trim(fault_callback_names(fault_callbacks(callback)))// &
                        " at call "//trim(call_number)//" for "//trim(case_name), &
                        n_hess_x=state%n_hess_x)) test_stability_check = .false.
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
            fault_context%stability_state => state

            ! set up settings with every optional callback function
            settings = stability_settings_type()
            settings%context => fault_context
            settings%diag_solver = diag_solver_options(i_solver)
            settings%jacobi_davidson_start = 0
            settings%precond => faulty_precond
            settings%project => faulty_project
            settings%conv_check => faulty_stability_conv_check
            settings%n_extra_trial_vectors = n_extra_trial_vectors
            settings%get_extra_trial_vectors => faulty_get_extra_trial_vectors

            ! run stability check
            call stability_check(h_diag, hess_x_funptr, stable, error, settings, &
                                 state, direction)

        end subroutine run_with_fault

    end function test_stability_check

    logical(c_bool) function test_check_stationary_point() bind(C)
        !
        ! this function tests the subroutine which performs the stability check at a
        ! stationary point and determines the step along an unstable mode
        !
        use opentrustregion, only: solver_settings_type, hess_x_type, &
                                   check_stationary_point, error_stability_check, &
                                   started_at_saddle_point_warning_msg, &
                                   reached_saddle_point_warning_msg, &
                                   reached_saddle_point_no_grad_noise_warning_msg
        use test_reference, only: arm_host_context, host_context_reached

        type(solver_settings_type) :: settings
        type(solver_state_type) :: state
        type(double_well_context_type), target :: context, stability_context
        procedure(hess_x_type), pointer :: hess_x_funptr
        real(rp) :: kappa(2), kappa_norm, grad(2), grad_rms
        logical :: stable
        integer(ip) :: error

        ! the double well function x^2 - y^2 + y^4 has a vanishing gradient and the
        ! Hessian diagonal (2, -2) at its saddle point at the origin, whose unstable
        ! mode is the second unit vector, and the Hessian diagonal (2, 4) at its minima
        real(rp), parameter :: stationary_grad(2) = [0.0_rp, 0.0_rp], &
                               saddle_h_diag(2) = [2.0_rp, -2.0_rp], &
                               minimum_h_diag(2) = [2.0_rp, 4.0_rp]

        ! assume tests pass
        test_check_stationary_point = .true.

        ! set function pointer
        hess_x_funptr => double_well_hess_x

        ! check the saddle point in the first macro iteration with a converged gradient
        ! along the unstable mode and a context in the nested stability check settings
        ! and check that this context reaches the callback functions of the stability
        ! check, that its Hessian linear transformations are added to the solver's
        ! counter, that a step along the unstable mode is returned with its norm and
        ! the length the gradient and the eigenvalue set and that the warning for
        ! starting at a saddle point is printed
        call setup_settings(settings, context)
        state%n_hess_x = 5
        state%trust_radius = 1.0_rp
        context%vars = [0.0_rp, 0.0_rp]
        grad = [0.0_rp, 0.5_rp * settings%conv_tol * sqrt(2.0_rp)]
        grad_rms = norm2(grad) / sqrt(2.0_rp)
        call arm_host_context(settings%stability_settings, context)
        call check_stationary_point(grad, grad_rms, saddle_h_diag, hess_x_funptr, &
                                    1_ip, settings, state, stable, kappa, kappa_norm, &
                                    error)
        if (error /= 0) then
            write(stderr, *) "test_check_stationary_point failed: Produced error "// &
                "at saddle point."
            test_check_stationary_point = .false.
        end if
        if (stable) then
            write(stderr, *) "test_check_stationary_point failed: Saddle point not "// &
                "found to be unstable."
            test_check_stationary_point = .false.
        end if
        if (abs(kappa(1)) > tol .or. abs(kappa(2)) < tol) then
            write(stderr, *) "test_check_stationary_point failed: Step not along "// &
                "unstable mode."
            test_check_stationary_point = .false.
        end if
        if (abs(kappa_norm - norm2(kappa)) > tol) then
            write(stderr, *) "test_check_stationary_point failed: Norm of step "// &
                "along unstable mode not returned."
            test_check_stationary_point = .false.
        end if

        ! the step along the unstable mode, whose eigenvalue -2 the stability check
        ! leaves in its state, has to raise the RMS of the linearized gradient from the
        ! given converged one to ten times the convergence threshold, pointing against
        ! the gradient, comparing such that a NaN fails
        if (.not. abs(norm2(grad - 2.0_rp * kappa) / sqrt(2.0_rp) - &
                      10.0_rp * settings%conv_tol) <= tol) then
            write(stderr, *) "test_check_stationary_point failed: Step along "// &
                "unstable mode does not raise the linearized gradient RMS to ten "// &
                "times the convergence threshold."
            test_check_stationary_point = .false.
        end if
        if (.not. dot_product(grad, kappa) < 0.0_rp) then
            write(stderr, *) "test_check_stationary_point failed: Step along "// &
                "unstable mode does not point against the gradient."
            test_check_stationary_point = .false.
        end if
        if (.not. host_context_reached("check_stationary_point", context)) &
            test_check_stationary_point = .false.
        if (state%stability_state%n_hess_x <= 0) then
            write(stderr, *) "test_check_stationary_point failed: Hessian linear "// &
                "transformations of the stability check not counted."
            test_check_stationary_point = .false.
        end if
        if (state%n_hess_x /= 5 + state%stability_state%n_hess_x) then
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

        ! check the saddle point in a later macro iteration with the gradient
        ! perturbation and with nested stability check settings with a context of their
        ! own, which the stability check hands to every callback function it calls,
        ! which its own test covers, so that the Hessian linear transformation has to
        ! receive it instead of the solver's, and check that the warning for reaching a
        ! saddle point is printed
        call setup_settings(settings, context)
        settings%grad_noise = .true.
        context%vars = [0.0_rp, 0.0_rp]
        stability_context%vars = [0.0_rp, 0.0_rp]
        stability_context%hess_x_called = .false.
        settings%stability_settings%context => stability_context
        call check_stationary_point(stationary_grad, 0.0_rp, saddle_h_diag, &
                                    hess_x_funptr, 2_ip, settings, state, stable, &
                                    kappa, kappa_norm, error)
        if (error /= 0) then
            write(stderr, *) "test_check_stationary_point failed: Produced error "// &
                "at saddle point in later macro iteration."
            test_check_stationary_point = .false.
        end if
        if (.not. stability_context%hess_x_called) then
            write(stderr, *) "test_check_stationary_point failed: Context of "// &
                "nested settings not handed to the Hessian linear transformation."
            test_check_stationary_point = .false.
        end if
        if (index(context%log_message, " "//reached_saddle_point_warning_msg) == 0) then
            write(stderr, *) "test_check_stationary_point failed: Warning not "// &
                "printed for reaching saddle point."
            test_check_stationary_point = .false.
        end if

        ! check the saddle point in a later macro iteration without the gradient
        ! perturbation, with truncated conjugate gradient, which also only perturbs the
        ! gradient with it, and with a vanishing convergence threshold, as the solver
        ! reaches it at an exactly vanishing gradient, and check that the step along
        ! the unstable mode reaches the trust radius of the state and that the warning
        ! for reaching a saddle point advises enabling the perturbation
        call setup_settings(settings, context)
        context%vars = [0.0_rp, 0.0_rp]
        settings%grad_noise = .false.
        settings%subsystem_solver = "tcg"
        settings%conv_tol = 0.0_rp
        settings%stability_settings%context => context
        state%trust_radius = 0.7_rp
        call check_stationary_point(stationary_grad, 0.0_rp, saddle_h_diag, &
                                    hess_x_funptr, 2_ip, settings, state, stable, &
                                    kappa, kappa_norm, error)
        if (error /= 0) then
            write(stderr, *) "test_check_stationary_point failed: Produced error "// &
                "at saddle point without gradient perturbation."
            test_check_stationary_point = .false.
        end if
        if (abs(kappa_norm - state%trust_radius) > tol) then
            write(stderr, *) "test_check_stationary_point failed: Step along "// &
                "unstable mode does not reach trust radius of the state at "// &
                "vanishing gradient and convergence threshold."
            test_check_stationary_point = .false.
        end if
        if (index(context%log_message, &
                  " "//reached_saddle_point_no_grad_noise_warning_msg) == 0) then
            write(stderr, *) "test_check_stationary_point failed: Warning not "// &
                "printed for reaching saddle point without gradient perturbation."
            test_check_stationary_point = .false.
        end if

        ! check that a minimum is found to be stable
        call setup_settings(settings, context)
        settings%stability_settings%context => context
        context%vars = [0.0_rp, 1.0_rp / sqrt(2.0_rp)]
        call check_stationary_point(stationary_grad, 0.0_rp, minimum_h_diag, &
                                    hess_x_funptr, 2_ip, settings, state, stable, &
                                    kappa, kappa_norm, error)
        if (error /= 0) then
            write(stderr, *) "test_check_stationary_point failed: Produced error "// &
                "at minimum."
            test_check_stationary_point = .false.
        end if
        if (.not. stable) then
            write(stderr, *) "test_check_stationary_point failed: Minimum not "// &
                "found to be stable."
            test_check_stationary_point = .false.
        end if

        ! check that nested stability check settings which its sanity check rejects,
        ! here an unknown diagonalization solver, are reported with the origin of the
        ! stability check
        call setup_settings(settings, context)
        settings%stability_settings%diag_solver = "unknown"
        call check_stationary_point(stationary_grad, 0.0_rp, saddle_h_diag, &
                                    hess_x_funptr, 2_ip, settings, state, stable, &
                                    kappa, kappa_norm, error)
        if (error /= error_stability_check + 1) then
            write(stderr, *) "test_check_stationary_point failed: Error of "// &
                "rejected nested stability check settings not reported with the "// &
                "origin of the stability check."
            test_check_stationary_point = .false.
        end if

    end function test_check_stationary_point

    logical(c_bool) function test_saddle_point_step() bind(C)
        !
        ! this function tests the saddle-point step subroutine
        !
        use opentrustregion, only: solver_settings_type, saddle_point_step
        use, intrinsic :: ieee_arithmetic, only: ieee_is_nan

        ! cases with a gradient RMS below and well above the convergence threshold and
        ! a direction of negative curvature orthogonal to the gradient or overlapping
        ! it positively or negatively
        character(len=*), parameter :: case_names(4) = [ &
            character(len=53) :: "below the convergence threshold", "above the "// &
            "convergence threshold", "above the convergence threshold with "// &
            "positive overlap", "above the convergence threshold with negative overlap"]
        real(rp) :: grad_dir(n_param), kappa_dir(n_param), unstable_dir(n_param), &
                    grad(n_param), kappa(n_param), kappa_norm, min_eigval, grad_rms, &
                    target_grad_rms, predicted_grad_rms, trust_radius
        integer(ip) :: i
        type(solver_settings_type) :: settings

        ! assume tests pass
        test_saddle_point_step = .true.

        ! generate a trust radius, which only a vanishing target uses
        call random_number(trust_radius)
        trust_radius = trust_radius + 0.1_rp

        ! set convergence threshold
        settings%conv_tol = 1e-6_rp

        ! generate a random gradient direction and a random direction of negative
        ! curvature orthogonal to it, mimicking a (near-)saddle point where the
        ! negative-curvature eigenvector has no overlap with the gradient
        call random_number(grad_dir)
        grad_dir = grad_dir / norm2(grad_dir)
        call random_number(kappa_dir)
        kappa_dir = kappa_dir - dot_product(kappa_dir, grad_dir) * grad_dir
        kappa_dir = kappa_dir / norm2(kappa_dir)
        min_eigval = -0.5_rp

        ! test a gradient RMS below the convergence threshold, as when starting exactly
        ! at a stationary point, and well above it, as when maximum precision is
        ! reached at an unstable point before the gradient has converged, with the
        ! direction of negative curvature orthogonal to the gradient, and well above it
        ! with a direction overlapping the gradient, along which the step has to point
        ! against the gradient
        do i = 1, size(case_names)
            if (i == 1) then
                grad = 0.1_rp * settings%conv_tol * sqrt(real(n_param, kind=rp)) * &
                       grad_dir
            else
                grad = 50.0_rp * settings%conv_tol * sqrt(real(n_param, kind=rp)) * &
                       grad_dir
            end if
            grad_rms = norm2(grad) / sqrt(real(n_param, kind=rp))
            target_grad_rms = 10.0_rp * max(settings%conv_tol, grad_rms)
            if (i <= 2) then
                unstable_dir = kappa_dir
            else
                unstable_dir = kappa_dir + merge(0.5_rp, -0.5_rp, i == 3) * grad_dir
                unstable_dir = unstable_dir / norm2(unstable_dir)
            end if

            kappa = unstable_dir
            call saddle_point_step(grad, grad_rms, min_eigval, trust_radius, settings, &
                                   kappa, kappa_norm)

            ! the linearized gradient after taking the step should have exactly the
            ! targeted gradient RMS
            predicted_grad_rms = norm2(grad + min_eigval * kappa) / &
                                 sqrt(real(n_param, kind=rp))
            if (ieee_is_nan(predicted_grad_rms) .or. &
                abs(predicted_grad_rms - target_grad_rms) > tol) then
                write(stderr, *) "test_saddle_point_step failed: Step does not "// &
                    "reach targeted gradient RMS "//trim(case_names(i))//"."
                test_saddle_point_step = .false.
            end if
            if (abs(kappa_norm - norm2(kappa)) > tol) then
                write(stderr, *) "test_saddle_point_step failed: Returned norm "// &
                    "does not match step "//trim(case_names(i))//"."
                test_saddle_point_step = .false.
            end if
            if (i > 2 .and. dot_product(grad, kappa) >= 0.0_rp) then
                write(stderr, *) "test_saddle_point_step failed: Step does not "// &
                    "point against the gradient "//trim(case_names(i))//"."
                test_saddle_point_step = .false.
            end if
        end do

        ! with a vanishing gradient and convergence threshold no gradient RMS can be
        ! targeted, and the step along the direction of negative curvature has to reach
        ! the trust radius instead
        settings%conv_tol = 0.0_rp
        grad = 0.0_rp
        kappa = kappa_dir
        call saddle_point_step(grad, 0.0_rp, min_eigval, trust_radius, settings, &
                               kappa, kappa_norm)
        if (norm2(kappa - trust_radius * kappa_dir) > tol .or. &
            abs(kappa_norm - trust_radius) > tol) then
            write(stderr, *) "test_saddle_point_step failed: Step does not reach "// &
                "trust radius at vanishing gradient and convergence threshold."
            test_saddle_point_step = .false.
        end if

        ! with a trust radius below the length of the step that reaches the targeted
        ! gradient RMS, the step, which is taken without evaluating the objective
        ! function, is limited to the trust radius and still points against the
        ! gradient
        settings%conv_tol = 1e-6_rp
        grad = 50.0_rp * settings%conv_tol * sqrt(real(n_param, kind=rp)) * grad_dir
        grad_rms = norm2(grad) / sqrt(real(n_param, kind=rp))
        unstable_dir = kappa_dir + 0.5_rp * grad_dir
        kappa = unstable_dir / norm2(unstable_dir)
        trust_radius = 1e-5_rp
        call saddle_point_step(grad, grad_rms, min_eigval, trust_radius, settings, &
                               kappa, kappa_norm)
        if (.not. abs(norm2(kappa) - trust_radius) <= tol) then
            write(stderr, *) "test_saddle_point_step failed: Step not limited to "// &
                "the trust radius."
            test_saddle_point_step = .false.
        end if
        if (.not. dot_product(grad, kappa) < 0.0_rp) then
            write(stderr, *) "test_saddle_point_step failed: Step limited to the "// &
                "trust radius does not point against the gradient."
            test_saddle_point_step = .false.
        end if

    end function test_saddle_point_step

    logical(c_bool) function test_newton_step() bind(C)
        !
        ! this function tests the Newton step subroutine
        !
        use opentrustregion, only: newton_step

        character(len=*), parameter :: case_names(3) = &
            [character(len=10) :: "symmetric", "asymmetric", "singular"]
        real(rp) :: red_space_basis(n_param, n_trial), grad(n_param), grad_norm, &
                    red_space_hess(n_trial, n_trial), &
                    symm_red_space_hess(n_trial, n_trial), eigvals(n_trial), &
                    eigvals_im(n_trial), right_eigvecs(n_trial, n_trial), &
                    left_eigvecs(n_trial, n_trial), solution(n_param), &
                    red_space_solution(n_trial)
        integer(ip) :: i_case

        ! assume tests pass
        test_newton_step = .true.

        ! set up the reduced space near a minimum, where its Hessian is positive
        ! definite
        call setup_reduced_space(near_minimum, grad, grad_norm, red_space_basis, &
                                 symm_red_space_hess)

        ! perform the Newton step with the eigendecomposition of the reduced space
        ! Hessian, of this Hessian made asymmetric and of this Hessian with the last
        ! basis vector, which has no gradient component, decoupled with a vanishing
        ! eigenvalue, whose eigenvector is given a gradient component at the level of
        ! rounding errors, along which the step then has no component, and determine
        ! whether it solves the Newton equation in the reduced space and is returned in
        ! the full space
        do i_case = 1, size(case_names)
            red_space_hess = symm_red_space_hess
            if (i_case == 1) then
                call ref_symm_mat_diag(red_space_hess, eigvals, right_eigvecs)
                left_eigvecs = right_eigvecs
            else if (i_case == 2) then
                red_space_hess(1, 2) = red_space_hess(1, 2) + 0.5_rp
                call ref_general_mat_diag(red_space_hess, eigvals, eigvals_im, &
                                          right_eigvecs, left_eigvecs)
            else
                red_space_hess(:, n_trial) = 0.0_rp
                red_space_hess(n_trial, :) = 0.0_rp
                eigvals = 0.0_rp
                right_eigvecs = 0.0_rp
                call ref_symm_mat_diag(red_space_hess(:n_trial - 1, :n_trial - 1), &
                                       eigvals(:n_trial - 1), &
                                       right_eigvecs(:n_trial - 1, :n_trial - 1))
                right_eigvecs(n_trial, n_trial) = 1.0_rp
                left_eigvecs = right_eigvecs
                left_eigvecs(1, n_trial) = 1e-16_rp
            end if
            call newton_step(grad_norm, red_space_basis, eigvals, right_eigvecs, &
                             left_eigvecs, solution, red_space_solution)
            if (.not. norm2(matmul(red_space_basis, matmul( &
                red_space_hess, red_space_solution)) + grad) <= tol) then
                write(stderr, *) "test_newton_step failed: Newton step does not "// &
                    "satisfy Newton equation for "//trim(case_names(i_case))// &
                    " Hessian."
                test_newton_step = .false.
            end if
            if (any(abs(solution - matmul(red_space_basis, red_space_solution)) > &
                    tol)) then
                write(stderr, *) "test_newton_step failed: Full space solution not "// &
                    "correct for "//trim(case_names(i_case))//" Hessian."
                test_newton_step = .false.
            end if
        end do

    end function test_newton_step

    logical(c_bool) function test_bisection_ah() bind(C)
        !
        ! this function tests the augmented Hessian bisection subroutine
        !
        use opentrustregion, only: solver_settings_type, bisection_ah
        use, intrinsic :: ieee_arithmetic, only: ieee_is_nan

        type(solver_settings_type) :: settings
        real(rp) :: red_space_basis(n_param, n_trial), grad(n_param), grad_norm, &
                    red_space_hess(n_trial, n_trial), red_space_hess_eigvals(n_trial), &
                    red_space_hess_left_eigvecs(n_trial, n_trial), &
                    red_space_hess_right_eigvecs(n_trial, n_trial), solution(n_param), &
                    red_space_solution(n_trial), trust_radius, mu, &
                    grad_coupled_component, newton_solution(n_param), &
                    newton_red_space_solution(n_trial)
        integer(ip) :: error, i
        type(test_context_type), target :: context

        ! assume tests pass
        test_bisection_ah = .true.

        ! setup settings object
        call setup_settings(settings, context)

        ! set up the reduced space at a point with strong negative curvature
        trust_radius = 1.0_rp
        call setup_reduced_space(negative_curvature_point, grad, grad_norm, &
                                 red_space_basis, red_space_hess)
        call ref_symm_mat_diag(red_space_hess, red_space_hess_eigvals, &
                               red_space_hess_right_eigvecs)
        red_space_hess_left_eigvecs = red_space_hess_right_eigvecs

        ! perform bisection, check whether error has occured and determine whether
        ! resulting solution is correct in reduced and full space and respects target
        ! trust radius
        call bisection_ah(red_space_hess, grad_norm, red_space_basis, &
                          red_space_hess_eigvals, red_space_hess_right_eigvecs, &
                          red_space_hess_left_eigvecs, trust_radius, solution, &
                          red_space_solution, mu, settings, error)
        if (error /= 0) then
            write(stderr, *) "test_bisection_ah failed: Produced error."
            test_bisection_ah = .false.
        end if
        if (abs(norm2(solution) - trust_radius) > tol) then
            write(stderr, *) "test_bisection_ah failed: Solution does not respect "// &
                "trust radius."
            test_bisection_ah = .false.
        end if
        if (norm2(matmul(red_space_basis, matmul(red_space_hess, red_space_solution) - &
                         mu * red_space_solution) + grad) > tol) then
            write(stderr, *) "test_bisection_ah failed: Step does not satisfy "// &
                "level-shifted Newton equation."
            test_bisection_ah = .false.
        end if
        if (any(abs(solution - matmul(red_space_basis, red_space_solution)) > tol)) then
            write(stderr, *) "test_bisection_ah failed: Full space solution not "// &
                "correct for symmetric Hessian."
            test_bisection_ah = .false.
        end if
        if (mu >= minval(red_space_hess_eigvals)) then
            write(stderr, *) "test_bisection_ah failed: Level shift not below "// &
                "lowest reduced space Hessian eigenvalue."
            test_bisection_ah = .false.
        end if

        ! set up a reduced space Hessian in an orthonormal basis with a negative
        ! eigenvalue whose eigenvector lies along the gradient and which is large
        ! compared to the gradient, so that the step at the starting alpha is longer
        ! than the trust radius and alpha has to be increased repeatedly to bracket the
        ! trust radius
        red_space_basis = 0.0_rp
        red_space_hess = 0.0_rp
        do i = 1, n_trial
            red_space_basis(i, i) = 1.0_rp
        end do
        red_space_hess(1, 1) = -10.0_rp
        red_space_hess(2, 2) = 2.0_rp
        red_space_hess(3, 3) = 3.0_rp
        grad_norm = 0.1_rp
        grad = 0.0_rp
        grad(1) = grad_norm
        trust_radius = 0.5_rp
        call ref_symm_mat_diag(red_space_hess, red_space_hess_eigvals, &
                               red_space_hess_right_eigvecs)
        red_space_hess_left_eigvecs = red_space_hess_right_eigvecs

        ! perform bisection and check whether the solution respects the trust radius
        ! and satisfies the level-shifted Newton equation with a level shift below the
        ! lowest eigenvalue
        call bisection_ah(red_space_hess, grad_norm, red_space_basis, &
                          red_space_hess_eigvals, red_space_hess_right_eigvecs, &
                          red_space_hess_left_eigvecs, trust_radius, solution, &
                          red_space_solution, mu, settings, error)
        if (error /= 0) then
            write(stderr, *) "test_bisection_ah failed: Produced error when alpha "// &
                "has to be increased."
            test_bisection_ah = .false.
        end if
        if (abs(norm2(solution) - trust_radius) > tol) then
            write(stderr, *) "test_bisection_ah failed: Solution does not respect "// &
                "trust radius when alpha has to be increased."
            test_bisection_ah = .false.
        end if
        if (norm2(matmul(red_space_basis, matmul(red_space_hess, red_space_solution) - &
                         mu * red_space_solution) + grad) > tol) then
            write(stderr, *) "test_bisection_ah failed: Step does not satisfy "// &
                "level-shifted Newton equation when alpha has to be increased."
            test_bisection_ah = .false.
        end if
        if (mu >= minval(red_space_hess_eigvals)) then
            write(stderr, *) "test_bisection_ah failed: Level shift not below "// &
                "lowest reduced space Hessian eigenvalue when alpha has to be "// &
                "increased."
            test_bisection_ah = .false.
        end if

        ! let the gradient vanish, so that the lowest augmented Hessian eigenvector,
        ! the eigenvector of the negative lowest reduced space Hessian eigenvalue,
        ! which lies along the gradient direction, so that it is not the hard case, has
        ! no level-shift component, and check that an error is returned and an error
        ! message is printed
        grad_norm = 0.0_rp
        call setup_error_logging(settings, context)
        call bisection_ah(red_space_hess, grad_norm, red_space_basis, &
                          red_space_hess_eigvals, red_space_hess_right_eigvecs, &
                          red_space_hess_left_eigvecs, trust_radius, solution, &
                          red_space_solution, mu, settings, error)
        if (error /= 1) then
            write(stderr, *) "test_bisection_ah failed: No error returned when "// &
                "lowest augmented Hessian eigenvector has no level-shift component."
            test_bisection_ah = .false.
        end if
        if (len_trim(context%log_message) == 0) then
            write(stderr, *) "test_bisection_ah failed: No error message printed "// &
                "when lowest augmented Hessian eigenvector has no level-shift "// &
                "component."
            test_bisection_ah = .false.
        end if

        ! set up the reduced space near a minimum
        call setup_reduced_space(near_minimum, grad, grad_norm, red_space_basis, &
                                 red_space_hess)
        call ref_symm_mat_diag(red_space_hess, red_space_hess_eigvals, &
                               red_space_hess_right_eigvecs)
        red_space_hess_left_eigvecs = red_space_hess_right_eigvecs

        ! perform bisection and determine whether routine correctly falls back to the
        ! Newton step since the minimum is closer than the target trust radius and no
        ! level shift is necessary
        call bisection_ah(red_space_hess, grad_norm, red_space_basis, &
                          red_space_hess_eigvals, red_space_hess_right_eigvecs, &
                          red_space_hess_left_eigvecs, trust_radius, solution, &
                          red_space_solution, mu, settings, error)
        if (error /= 0) then
            write(stderr, *) "test_bisection_ah failed: Produced error instead of "// &
                "falling back to Newton step."
            test_bisection_ah = .false.
        end if
        if (abs(mu) > tol) then
            write(stderr, *) "test_bisection_ah failed: Level shift not zero for "// &
                "Newton step fallback."
            test_bisection_ah = .false.
        end if
        newton_red_space_solution = -grad_norm * matmul( &
            red_space_hess_right_eigvecs, &
            red_space_hess_left_eigvecs(1, :) / red_space_hess_eigvals)
        newton_solution = matmul(red_space_basis, newton_red_space_solution)
        if (any(abs(red_space_solution - newton_red_space_solution) > tol)) then
            write(stderr, *) "test_bisection_ah failed: Newton step fallback "// &
                "reduced space solution not correct."
            test_bisection_ah = .false.
        end if
        if (any(abs(solution - newton_solution) > tol)) then
            write(stderr, *) "test_bisection_ah failed: Newton step fallback full "// &
                "space solution not correct."
            test_bisection_ah = .false.
        end if

        ! set up hard case
        call setup_hard_case([5.0_rp, -2.0_rp, 3.0_rp], red_space_basis, grad_norm, &
                             red_space_hess)
        trust_radius = 0.6_rp
        call ref_symm_mat_diag(red_space_hess, red_space_hess_eigvals, &
                               red_space_hess_right_eigvecs)
        red_space_hess_left_eigvecs = red_space_hess_right_eigvecs

        ! test hard case: lowest reduced space Hessian eigenvalue is negative and its
        ! eigenvector has no component along the gradient direction, so no level shift
        ! can reproduce the trust region solution and it must be constructed directly
        ! from the eigendecomposition
        call bisection_ah(red_space_hess, grad_norm, red_space_basis, &
                          red_space_hess_eigvals, red_space_hess_right_eigvecs, &
                          red_space_hess_left_eigvecs, trust_radius, solution, &
                          red_space_solution, mu, settings, error)
        if (error /= 0) then
            write(stderr, *) "test_bisection_ah failed: Produced error for hard case."
            test_bisection_ah = .false.
        end if
        if (abs(norm2(solution) - trust_radius) > tol) then
            write(stderr, *) "test_bisection_ah failed: Hard case solution does "// &
                "not respect trust radius."
            test_bisection_ah = .false.
        end if
        if (abs(mu - minval(red_space_hess_eigvals)) > tol) then
            write(stderr, *) "test_bisection_ah failed: Hard case level shift not "// &
                "equal to lowest reduced space Hessian eigenvalue."
            test_bisection_ah = .false.
        end if
        grad_coupled_component = grad_norm / (minval(red_space_hess_eigvals) - &
                                              maxval(red_space_hess_eigvals))
        if (any(abs(red_space_solution - [ &
            grad_coupled_component, sqrt(trust_radius**2 - grad_coupled_component**2), &
            0.0_rp]) > tol)) then
            write(stderr, *) "test_bisection_ah failed: Hard case reduced space "// &
                "solution not correct."
            test_bisection_ah = .false.
        end if
        if (any(abs(solution - matmul(red_space_basis, red_space_solution)) > tol)) then
            write(stderr, *) "test_bisection_ah failed: Hard case full space "// &
                "solution not correct."
            test_bisection_ah = .false.
        end if

        ! test the hard case with a trust radius below the step at the crossover point,
        ! which then lies outside the trust region, so that the bisection has to start
        ! there and find the level-shifted Newton step on the trust region boundary,
        ! which has no component along the lowest eigenvector and a level shift below
        ! the lowest eigenvalue, here the step along the gradient direction with the
        ! level shift -5
        trust_radius = 0.1_rp
        grad = 0.0_rp
        grad(1) = grad_norm
        call bisection_ah(red_space_hess, grad_norm, red_space_basis, &
                          red_space_hess_eigvals, red_space_hess_right_eigvecs, &
                          red_space_hess_left_eigvecs, trust_radius, solution, &
                          red_space_solution, mu, settings, error)
        if (error /= 0) then
            write(stderr, *) "test_bisection_ah failed: Produced error for hard "// &
                "case with crossover step outside trust region."
            test_bisection_ah = .false.
        end if
        if (abs(norm2(solution) - trust_radius) > tol) then
            write(stderr, *) "test_bisection_ah failed: Hard case solution with "// &
                "crossover step outside trust region does not respect trust radius."
            test_bisection_ah = .false.
        end if
        if (norm2(matmul(red_space_basis, matmul(red_space_hess, red_space_solution) - &
                         mu * red_space_solution) + grad) > tol) then
            write(stderr, *) "test_bisection_ah failed: Hard case solution with "// &
                "crossover step outside trust region does not satisfy "// &
                "level-shifted Newton equation."
            test_bisection_ah = .false.
        end if
        if (abs(mu + 5.0_rp) > tol) then
            write(stderr, *) "test_bisection_ah failed: Hard case level shift with "// &
                "crossover step outside trust region not that of the step along "// &
                "gradient direction."
            test_bisection_ah = .false.
        end if
        if (any(abs(red_space_solution - [-trust_radius, 0.0_rp, 0.0_rp]) > tol)) then
            write(stderr, *) "test_bisection_ah failed: Hard case solution with "// &
                "crossover step outside trust region not along gradient direction."
            test_bisection_ah = .false.
        end if

        ! set up hard case with degenerate lowest eigenvalues
        call setup_hard_case([5.0_rp, -2.0_rp, -2.0_rp], red_space_basis, grad_norm, &
                             red_space_hess)
        trust_radius = 0.5_rp
        call ref_symm_mat_diag(red_space_hess, red_space_hess_eigvals, &
                               red_space_hess_right_eigvecs)
        red_space_hess_left_eigvecs = red_space_hess_right_eigvecs

        ! test hard case with a degenerate lowest eigenvalue: the eigenvector spanning
        ! the degenerate subspace used to fill the trust radius is not uniquely
        ! defined, so only invariant properties of the solution are checked rather than
        ! exact reduced space solution components
        call bisection_ah(red_space_hess, grad_norm, red_space_basis, &
                          red_space_hess_eigvals, red_space_hess_right_eigvecs, &
                          red_space_hess_left_eigvecs, trust_radius, solution, &
                          red_space_solution, mu, settings, error)
        if (error /= 0) then
            write(stderr, *) "test_bisection_ah failed: Produced error for "// &
                "degenerate hard case."
            test_bisection_ah = .false.
        end if
        if (any(ieee_is_nan(solution)) .or. any(ieee_is_nan(red_space_solution))) then
            write(stderr, *) "test_bisection_ah failed: Degenerate hard case "// &
                "solution contains NaN."
            test_bisection_ah = .false.
        end if
        if (abs(norm2(solution) - trust_radius) > tol) then
            write(stderr, *) "test_bisection_ah failed: Degenerate hard case "// &
                "solution does not respect trust radius."
            test_bisection_ah = .false.
        end if
        if (abs(mu - minval(red_space_hess_eigvals)) > tol) then
            write(stderr, *) "test_bisection_ah failed: Degenerate hard case level "// &
                "shift not equal to lowest reduced space Hessian eigenvalue."
            test_bisection_ah = .false.
        end if
        grad_coupled_component = grad_norm / (minval(red_space_hess_eigvals) - &
                                              maxval(red_space_hess_eigvals))
        if (abs(red_space_solution(1) - grad_coupled_component) > tol) then
            write(stderr, *) "test_bisection_ah failed: Degenerate hard case "// &
                "component along non-degenerate eigenvector not correct."
            test_bisection_ah = .false.
        end if
        if (abs(norm2(red_space_solution(2:3)) - &
                sqrt(trust_radius**2 - grad_coupled_component**2)) > tol) then
            write(stderr, *) "test_bisection_ah failed: Degenerate hard case fill "// &
                "magnitude within degenerate eigenspace not correct."
            test_bisection_ah = .false.
        end if
        if (any(abs(solution - matmul(red_space_basis, red_space_solution)) > tol)) then
            write(stderr, *) "test_bisection_ah failed: Degenerate hard case full "// &
                "space solution not correct."
            test_bisection_ah = .false.
        end if

    end function test_bisection_ah

    logical(c_bool) function test_bisection_mu() bind(C)
        !
        ! this function tests the level-shift bisection subroutine
        !
        use opentrustregion, only: solver_settings_type, bisection_mu, &
                                   bisection_imaginary_eigval_warning_msg
        use, intrinsic :: ieee_arithmetic, only: ieee_is_nan

        character(len=*), parameter :: case_names(4) = &
            [character(len=29) :: "symmetric, indefinite", "asymmetric, indefinite", &
             "symmetric, positive-definite", "asymmetric, positive-definite"]
        character(len=*), parameter :: hard_case_names(3) = &
            [character(len=68) :: "symmetric Hessian", &
             "asymmetric Hessian with right eigenvector along gradient direction", &
             "asymmetric Hessian with right eigenvector against gradient direction"]
        character(len=*), parameter :: newton_case_names(2) = &
            [character(len=38) :: "at the trust radius", &
             "exceeding the trust radius by one ulp"]
        type(solver_settings_type) :: settings
        real(rp) :: red_space_basis(n_param, n_trial), grad(n_param), grad_norm, &
                    red_space_hess(n_trial, n_trial), &
                    red_space_hess_left_eigvecs(n_trial, n_trial), &
                    red_space_hess_right_eigvecs(n_trial, n_trial), solution(n_param), &
                    red_space_solution(n_trial), trust_radius, mu, &
                    grad_coupled_component, expected_red_space_solution(n_trial), &
                    eigvec_overlap, rotation(n_trial, n_trial)
        complex(rp) :: red_space_hess_eigvals(n_trial)
        integer(ip) :: error, i_case
        logical :: positive_definite
        type(test_context_type), target :: context

        ! assume tests pass
        test_bisection_mu = .true.

        ! setup settings object
        call setup_settings(settings, context)

        ! perform bisection in the reduced space at a point with strong negative
        ! curvature and near a minimum, for the reduced space Hessian and for this
        ! Hessian made asymmetric, check whether error has occured and determine
        ! whether the solution is a level-shifted Newton step which reaches the trust
        ! radius in the indefinite case and the Newton step within it in the positive
        ! definite case, and whether it is returned in the full space
        trust_radius = 1.0_rp
        do i_case = 1, size(case_names)
            positive_definite = i_case > 2
            call setup_reduced_space( &
                merge(near_minimum, negative_curvature_point, positive_definite), &
                grad, grad_norm, red_space_basis, red_space_hess)
            if (mod(i_case, 2) == 0) &
                red_space_hess(1, 2) = red_space_hess(1, 2) + 0.5_rp
            call diagonalize(mod(i_case, 2) == 1)
            call bisection_mu(red_space_hess, grad_norm, red_space_basis, &
                              red_space_hess_eigvals, red_space_hess_right_eigvecs, &
                              red_space_hess_left_eigvecs, trust_radius, solution, &
                              red_space_solution, mu, settings, error)
            if (error /= 0) then
                write(stderr, *) "test_bisection_mu failed: Produced error for "// &
                    trim(case_names(i_case))//" Hessian."
                test_bisection_mu = .false.
            end if
            if ((positive_definite .and. norm2(solution) > trust_radius) .or. &
                (.not. positive_definite .and. &
                 abs(norm2(solution) - trust_radius) > tol)) then
                write(stderr, *) "test_bisection_mu failed: Solution does not "// &
                    "respect trust radius for "//trim(case_names(i_case))//" Hessian."
                test_bisection_mu = .false.
            end if
            if (positive_definite .and. abs(mu) > tol) then
                write(stderr, *) "test_bisection_mu failed: Level shift does not "// &
                    "vanish for "//trim(case_names(i_case))//" Hessian."
                test_bisection_mu = .false.
            end if
            if (norm2( &
                matmul(red_space_basis, matmul(red_space_hess, red_space_solution) - &
                       mu * red_space_solution) + grad) > tol) then
                write(stderr, *) "test_bisection_mu failed: Step does not satisfy "// &
                    "level-shifted Newton equation for "//trim(case_names(i_case))// &
                    " Hessian."
                test_bisection_mu = .false.
            end if
            if (any(abs(solution - matmul(red_space_basis, red_space_solution)) > &
                    tol)) then
                write(stderr, *) "test_bisection_mu failed: Full space solution "// &
                    "not correct for "//trim(case_names(i_case))//" Hessian."
                test_bisection_mu = .false.
            end if
        end do

        ! set up hard case
        call setup_hard_case([5.0_rp, -2.0_rp, 3.0_rp], red_space_basis, grad_norm, &
                             red_space_hess)
        trust_radius = 0.6_rp

        ! test hard case: lowest reduced space Hessian eigenvalue is negative and its
        ! eigenvector has no component along the gradient direction, so no level shift
        ! can reproduce the trust region solution and the remaining trust radius has to
        ! be filled with this eigenvector, for the symmetric Hessian and for this
        ! Hessian made asymmetric by coupling the gradient direction to the right
        ! eigenvector of the lowest eigenvalue, (1, 1, 0) / sqrt(2), while its left
        ! eigenvector remains orthogonal to the gradient for the hard case to arise, so
        ! that the hard case is detected from the left eigenvector while the remaining
        ! trust radius is filled with the right eigenvector, which has to point against
        ! the gradient whichever sign it is returned with
        do i_case = 1, size(hard_case_names)
            if (i_case > 1) red_space_hess(1, 2) = -7.0_rp
            call diagonalize(i_case == 1)
            if (i_case > 1) red_space_hess_right_eigvecs( &
                :, minloc(red_space_hess_eigvals%re, dim=1)) = &
                merge(1.0_rp, -1.0_rp, i_case == 2) * [1.0_rp, 1.0_rp, 0.0_rp] / &
                sqrt(2.0_rp)
            call bisection_mu(red_space_hess, grad_norm, red_space_basis, &
                              red_space_hess_eigvals, red_space_hess_right_eigvecs, &
                              red_space_hess_left_eigvecs, trust_radius, solution, &
                              red_space_solution, mu, settings, error)
            if (error /= 0) then
                write(stderr, *) "test_bisection_mu failed: Produced error for "// &
                    "hard case with "//trim(hard_case_names(i_case))//"."
                test_bisection_mu = .false.
            end if
            if (abs(norm2(solution) - trust_radius) > tol) then
                write(stderr, *) "test_bisection_mu failed: Hard case solution "// &
                    "does not respect trust radius for "// &
                    trim(hard_case_names(i_case))//"."
                test_bisection_mu = .false.
            end if
            if (abs(mu - minval(red_space_hess_eigvals%re)) > tol) then
                write(stderr, *) "test_bisection_mu failed: Hard case level shift "// &
                    "not equal to lowest reduced space Hessian eigenvalue for "// &
                    trim(hard_case_names(i_case))//"."
                test_bisection_mu = .false.
            end if
            ! the step at the lowest eigenvalue, which for the asymmetric Hessian is
            ! the minimum norm solution (-1, 1, 0) / 14 orthogonal to the right
            ! eigenvector, completed to the trust radius along the lowest eigenvector
            if (i_case == 1) then
                grad_coupled_component = grad_norm / &
                                         (minval(red_space_hess_eigvals%re) - &
                                          maxval(red_space_hess_eigvals%re))
                expected_red_space_solution = [ &
                    grad_coupled_component, &
                    sqrt(trust_radius**2 - grad_coupled_component**2), 0.0_rp]
            else
                expected_red_space_solution = &
                    grad_norm * [-1.0_rp, 1.0_rp, 0.0_rp] / 14.0_rp - &
                    sqrt(trust_radius**2 - grad_norm**2 / 98.0_rp) * &
                    [1.0_rp, 1.0_rp, 0.0_rp] / sqrt(2.0_rp)
            end if
            if (any(abs(red_space_solution - expected_red_space_solution) > tol)) then
                write(stderr, *) "test_bisection_mu failed: Hard case reduced "// &
                    "space solution not correct for "//trim(hard_case_names(i_case))// &
                    "."
                test_bisection_mu = .false.
            end if
            if (any(abs(solution - matmul(red_space_basis, red_space_solution)) > &
                    tol)) then
                write(stderr, *) "test_bisection_mu failed: Hard case full space "// &
                    "solution not correct for "//trim(hard_case_names(i_case))//"."
                test_bisection_mu = .false.
            end if
        end do

        ! let the trust radius be smaller than the step at the lowest eigenvalue of the
        ! symmetric hard case, so that the bisection has to start from this singular
        ! level shift and find the level-shifted Newton step on the trust region
        ! boundary, which has no component along the lowest eigenvector and a level
        ! shift below the lowest eigenvalue, here the step along the gradient direction
        ! with the level shift -5
        red_space_hess(1, 2) = 0.0_rp
        trust_radius = 0.1_rp
        grad = 0.0_rp
        grad(1) = grad_norm
        call diagonalize(.true.)
        call bisection_mu(red_space_hess, grad_norm, red_space_basis, &
                          red_space_hess_eigvals, red_space_hess_right_eigvecs, &
                          red_space_hess_left_eigvecs, trust_radius, solution, &
                          red_space_solution, mu, settings, error)
        if (error /= 0) then
            write(stderr, *) "test_bisection_mu failed: Produced error for hard "// &
                "case with step at lowest eigenvalue outside trust region."
            test_bisection_mu = .false.
        end if
        if (abs(norm2(solution) - trust_radius) > tol) then
            write(stderr, *) "test_bisection_mu failed: Hard case solution with "// &
                "step at lowest eigenvalue outside trust region does not respect "// &
                "trust radius."
            test_bisection_mu = .false.
        end if
        if (norm2(matmul(red_space_basis, matmul(red_space_hess, red_space_solution) - &
                         mu * red_space_solution) + grad) > tol) then
            write(stderr, *) "test_bisection_mu failed: Hard case solution with "// &
                "step at lowest eigenvalue outside trust region does not satisfy "// &
                "level-shifted Newton equation."
            test_bisection_mu = .false.
        end if
        if (abs(mu + 5.0_rp) > tol) then
            write(stderr, *) "test_bisection_mu failed: Hard case level shift with "// &
                "step at lowest eigenvalue outside trust region not that of the "// &
                "step along gradient direction."
            test_bisection_mu = .false.
        end if
        if (any(abs(red_space_solution - [-trust_radius, 0.0_rp, 0.0_rp]) > tol)) then
            write(stderr, *) "test_bisection_mu failed: Hard case solution with "// &
                "step at lowest eigenvalue outside trust region not along gradient "// &
                "direction."
            test_bisection_mu = .false.
        end if

        ! rotate the eigenvector of the lowest eigenvalue of the symmetric hard case
        ! slightly towards the gradient direction, so that the solution norm has a pole
        ! at the lowest eigenvalue but the step at a level shift just below it, which
        ! already has a component along this eigenvector, lies within the trust radius,
        ! and check that the rest of the trust radius is filled along this eigenvector
        ! against the gradient with the level shift equal to the lowest eigenvalue
        eigvec_overlap = 1e-9_rp
        rotation = reshape([sqrt(1.0_rp - eigvec_overlap**2), -eigvec_overlap, 0.0_rp, &
                            eigvec_overlap, sqrt(1.0_rp - eigvec_overlap**2), 0.0_rp, &
                            0.0_rp, 0.0_rp, 1.0_rp], [n_trial, n_trial])
        red_space_hess = matmul(matmul(rotation, red_space_hess), transpose(rotation))
        trust_radius = 0.6_rp
        call diagonalize(.true.)
        call bisection_mu(red_space_hess, grad_norm, red_space_basis, &
                          red_space_hess_eigvals, red_space_hess_right_eigvecs, &
                          red_space_hess_left_eigvecs, trust_radius, solution, &
                          red_space_solution, mu, settings, error)
        if (error /= 0) then
            write(stderr, *) "test_bisection_mu failed: Produced error for lowest "// &
                "eigenvector nearly orthogonal to gradient direction."
            test_bisection_mu = .false.
        end if
        if (abs(norm2(solution) - trust_radius) > tol) then
            write(stderr, *) "test_bisection_mu failed: Solution does not respect "// &
                "trust radius for lowest eigenvector nearly orthogonal to gradient "// &
                "direction."
            test_bisection_mu = .false.
        end if
        if (abs(mu + 2.0_rp) > tol) then
            write(stderr, *) "test_bisection_mu failed: Level shift not equal to "// &
                "lowest eigenvalue for lowest eigenvector nearly orthogonal to "// &
                "gradient direction."
            test_bisection_mu = .false.
        end if
        if (dot_product(red_space_solution, rotation(:, 2)) >= 0.0_rp) then
            write(stderr, *) "test_bisection_mu failed: Step not against gradient "// &
                "along lowest eigenvector nearly orthogonal to gradient direction."
            test_bisection_mu = .false.
        end if

        ! set up hard case with degenerate lowest eigenvalues
        call setup_hard_case([5.0_rp, -2.0_rp, -2.0_rp], red_space_basis, grad_norm, &
                             red_space_hess)
        trust_radius = 0.5_rp
        call diagonalize(.true.)

        ! test hard case with a degenerate lowest eigenvalue: the eigenvector spanning
        ! the degenerate subspace used to fill the trust radius is not uniquely
        ! defined, so only invariant properties of the solution are checked rather than
        ! exact reduced space solution components
        call bisection_mu(red_space_hess, grad_norm, red_space_basis, &
                          red_space_hess_eigvals, red_space_hess_right_eigvecs, &
                          red_space_hess_left_eigvecs, trust_radius, solution, &
                          red_space_solution, mu, settings, error)
        if (error /= 0) then
            write(stderr, *) "test_bisection_mu failed: Produced error for "// &
                "degenerate hard case."
            test_bisection_mu = .false.
        end if
        if (any(ieee_is_nan(solution)) .or. any(ieee_is_nan(red_space_solution))) then
            write(stderr, *) "test_bisection_mu failed: Degenerate hard case "// &
                "solution contains NaN."
            test_bisection_mu = .false.
        end if
        if (abs(norm2(solution) - trust_radius) > tol) then
            write(stderr, *) "test_bisection_mu failed: Degenerate hard case "// &
                "solution does not respect trust radius."
            test_bisection_mu = .false.
        end if
        if (abs(mu - minval(red_space_hess_eigvals%re)) > tol) then
            write(stderr, *) "test_bisection_mu failed: Degenerate hard case level "// &
                "shift not equal to lowest reduced space Hessian eigenvalue."
            test_bisection_mu = .false.
        end if
        grad_coupled_component = grad_norm / (minval(red_space_hess_eigvals%re) - &
                                              maxval(red_space_hess_eigvals%re))
        if (abs(red_space_solution(1) - grad_coupled_component) > tol) then
            write(stderr, *) "test_bisection_mu failed: Degenerate hard case "// &
                "component along non-degenerate eigenvector not correct."
            test_bisection_mu = .false.
        end if
        if (abs(norm2(red_space_solution(2:3)) - &
                sqrt(trust_radius**2 - grad_coupled_component**2)) > tol) then
            write(stderr, *) "test_bisection_mu failed: Degenerate hard case fill "// &
                "magnitude within degenerate eigenspace not correct."
            test_bisection_mu = .false.
        end if
        if (any(abs(solution - matmul(red_space_basis, red_space_solution)) > tol)) then
            write(stderr, *) "test_bisection_mu failed: Degenerate hard case full "// &
                "space solution not correct."
            test_bisection_mu = .false.
        end if

        ! let the Newton step of a positive-definite Hessian with a small lowest
        ! eigenvalue reach the trust radius exactly or exceed it by the smallest
        ! representable amount and check that it is taken without a level shift,
        ! directly in the first case, since a bisection towards a vanishing level shift
        ! could not converge, and in the second case after a bisection which resolves
        ! the vanishing level shift absolutely, since the bracket cannot shrink
        ! relative to it within the bisection iterations
        red_space_hess = 0.0_rp
        red_space_hess(1, 1) = 1e-3_rp
        red_space_hess(2, 2) = 3.0_rp
        red_space_hess(3, 3) = 4.0_rp
        call diagonalize(.true.)
        do i_case = 1, 2
            trust_radius = grad_norm / red_space_hess(1, 1)
            if (i_case == 2) trust_radius = nearest(trust_radius, -1.0_rp)
            call bisection_mu(red_space_hess, grad_norm, red_space_basis, &
                              red_space_hess_eigvals, red_space_hess_right_eigvecs, &
                              red_space_hess_left_eigvecs, trust_radius, solution, &
                              red_space_solution, mu, settings, error)
            if (error /= 0) then
                write(stderr, *) "test_bisection_mu failed: Produced error for "// &
                    "Newton step "//trim(newton_case_names(i_case))//"."
                test_bisection_mu = .false.
            end if
            if (abs(mu) > tol) then
                write(stderr, *) "test_bisection_mu failed: Level shift does not "// &
                    "vanish for Newton step "//trim(newton_case_names(i_case))//"."
                test_bisection_mu = .false.
            end if
            if (any(abs(red_space_solution - [-trust_radius, 0.0_rp, 0.0_rp]) > &
                    tol * trust_radius)) then
                write(stderr, *) "test_bisection_mu failed: Newton step "// &
                    trim(newton_case_names(i_case))//" not taken."
                test_bisection_mu = .false.
            end if
        end do

        ! let the lowest eigenvalue of a positive semidefinite Hessian vanish with an
        ! eigenvector that has a component along the gradient direction, so that the
        ! level-shifted Hessian is singular at a vanishing level shift, at which the
        ! solution norm has a pole, and check that the solution is the level-shifted
        ! Newton step on the trust region boundary with a negative level shift
        red_space_hess = 0.0_rp
        red_space_hess(1:2, 1:2) = 1.0_rp
        red_space_hess(3, 3) = 3.0_rp
        trust_radius = 0.5_rp
        call diagonalize(.true.)
        call bisection_mu(red_space_hess, grad_norm, red_space_basis, &
                          red_space_hess_eigvals, red_space_hess_right_eigvecs, &
                          red_space_hess_left_eigvecs, trust_radius, solution, &
                          red_space_solution, mu, settings, error)
        if (error /= 0) then
            write(stderr, *) "test_bisection_mu failed: Produced error for "// &
                "vanishing lowest eigenvalue."
            test_bisection_mu = .false.
        end if
        if (abs(norm2(solution) - trust_radius) > tol) then
            write(stderr, *) "test_bisection_mu failed: Solution does not respect "// &
                "trust radius for vanishing lowest eigenvalue."
            test_bisection_mu = .false.
        end if
        if (norm2(matmul(red_space_basis, matmul(red_space_hess, red_space_solution) - &
                         mu * red_space_solution) + grad) > tol) then
            write(stderr, *) "test_bisection_mu failed: Step does not satisfy "// &
                "level-shifted Newton equation for vanishing lowest eigenvalue."
            test_bisection_mu = .false.
        end if
        if (mu >= 0.0_rp) then
            write(stderr, *) "test_bisection_mu failed: Level shift not negative "// &
                "for vanishing lowest eigenvalue."
            test_bisection_mu = .false.
        end if

        ! let the lowest eigenvalues of an asymmetric Hessian be a complex pair, which
        ! leads to no pole in the solution norm, so that the trust radius cannot be
        ! reached, and check that the level shift is the real part of these eigenvalues
        ! and that a warning is printed
        red_space_hess = 0.0_rp
        red_space_hess(1, 1) = 5.0_rp
        red_space_hess(2:3, 2:3) = reshape([-2.0_rp, -1.0_rp, &
                                            1.0_rp, -2.0_rp], [2, 2])
        trust_radius = 0.6_rp
        call diagonalize(.false.)
        context%log_message = ""
        call bisection_mu(red_space_hess, grad_norm, red_space_basis, &
                          red_space_hess_eigvals, red_space_hess_right_eigvecs, &
                          red_space_hess_left_eigvecs, trust_radius, solution, &
                          red_space_solution, mu, settings, error)
        if (error /= 0) then
            write(stderr, *) "test_bisection_mu failed: Produced error for complex "// &
                "lowest eigenvalues."
            test_bisection_mu = .false.
        end if
        if (abs(mu + 2.0_rp) > tol) then
            write(stderr, *) "test_bisection_mu failed: Level shift not the real "// &
                "part of complex lowest eigenvalues."
            test_bisection_mu = .false.
        end if
        if (norm2(solution) >= trust_radius) then
            write(stderr, *) "test_bisection_mu failed: Trust radius filled for "// &
                "complex lowest eigenvalues."
            test_bisection_mu = .false.
        end if
        if (index(context%log_message, " "//bisection_imaginary_eigval_warning_msg) == &
            0) then
            write(stderr, *) "test_bisection_mu failed: No warning printed for "// &
                "complex lowest eigenvalues."
            test_bisection_mu = .false.
        end if

    contains

        subroutine diagonalize(symmetric)
            !
            ! this subroutine diagonalizes the reduced space Hessian as a symmetric or
            ! as a general matrix
            !
            logical, intent(in) :: symmetric

            real(rp) :: eigvals_re(n_trial), eigvals_im(n_trial)

            if (symmetric) then
                call ref_symm_mat_diag(red_space_hess, eigvals_re, &
                                       red_space_hess_right_eigvecs)
                eigvals_im = 0.0_rp
                red_space_hess_left_eigvecs = red_space_hess_right_eigvecs
            else
                call ref_general_mat_diag(red_space_hess, eigvals_re, eigvals_im, &
                                          red_space_hess_right_eigvecs, &
                                          red_space_hess_left_eigvecs)
            end if
            red_space_hess_eigvals = cmplx(eigvals_re, eigvals_im, kind=rp)

        end subroutine diagonalize

    end function test_bisection_mu

    logical(c_bool) function test_bracket() bind(C)
        !
        ! this function tests the bracketing function
        !
        use opentrustregion, only: solver_settings_type, obj_func_type, bracket, &
                                   error_obj_func

        type(solver_settings_type) :: settings
        procedure(obj_func_type), pointer :: obj_func_funptr
        real(rp) :: n
        integer(ip) :: error, i, i_call, n_calls, n_obj_func
        character(len=20) :: call_number
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

        ! define procedure pointer
        obj_func_funptr => bracket_obj_func

        ! start the counter of objective function evaluations, which the fault cases
        ! below check
        n_obj_func = 0

        ! bracket the minimum of one-dimensional functions which drive the bracketing
        ! through its parabolic extrapolation steps, its golden ratio steps and its
        ! error conditions and determine whether the expected multiplier or an error
        ! is returned
        do i = 1, n_cases
            context%func_case = i
            call setup_error_logging(settings, context)
            n = bracket(obj_func_funptr, direction, 0.0_rp, 1.0_rp, settings, &
                        n_obj_func, error)
            if (expect_error(i)) then
                if (error /= 1) then
                    write(stderr, *) "test_bracket failed: No error returned for "// &
                        "the function "//trim(case_names(i))//"."
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
        n = bracket(obj_func_funptr, direction, 1.0_rp, 0.0_rp, settings, n_obj_func, &
                    error)
        if (error /= 0) then
            write(stderr, *) "test_bracket failed: Produced error for swapped bounds."
            test_bracket = .false.
        end if
        if (abs(n - expected_n(1)) > tol * abs(expected_n(1))) then
            write(stderr, *) "test_bracket failed: Returned multiplier not correct "// &
                "for swapped bounds."
            test_bracket = .false.
        end if

        ! let every call of the objective function fail in turn for every function
        ! that is bracketed without error, the calls being counted in an unfaulted
        ! run, which has to report them, and check that the error is reported with its
        ! origin and that the bracketing stops at the failing call, which is still
        ! counted
        do i = 1, n_cases
            if (expect_error(i)) cycle
            context%func_case = i
            context%fail_call = 0
            context%n_calls = 0
            n_obj_func = 0
            n = bracket(obj_func_funptr, direction, 0.0_rp, 1.0_rp, settings, &
                        n_obj_func, error)
            n_calls = context%n_calls
            if (n_obj_func /= n_calls) then
                write(stderr, *) "test_bracket failed: Reported number of "// &
                    "objective function evaluations wrong for the function "// &
                    trim(case_names(i))//"."
                test_bracket = .false.
            end if
            do i_call = 1, n_calls
                write(call_number, '(I0)') i_call
                context%fail_call = i_call
                context%n_calls = 0
                n_obj_func = 0
                n = bracket(obj_func_funptr, direction, 0.0_rp, 1.0_rp, settings, &
                            n_obj_func, error)
                if (error /= error_obj_func + 1) then
                    write(stderr, *) "test_bracket failed: Did not report the "// &
                        "error of the objective function failing at call "// &
                        trim(call_number)//" with its origin for the function "// &
                        trim(case_names(i))//"."
                    test_bracket = .false.
                end if
                if (context%n_calls /= i_call) then
                    write(stderr, *) "test_bracket failed: Did not stop at the "// &
                        "objective function failing at call "//trim(call_number)// &
                        " for the function "//trim(case_names(i))//"."
                    test_bracket = .false.
                end if
                if (n_obj_func /= i_call) then
                    write(stderr, *) "test_bracket failed: Reported number of "// &
                        "objective function evaluations wrong for the objective "// &
                        "function failing at call "//trim(call_number)// &
                        " for the function "//trim(case_names(i))//"."
                    test_bracket = .false.
                end if
            end do
        end do
        context%fail_call = 0

    end function test_bracket

    logical(c_bool) function test_extend_matrix() bind(C)
        !
        ! this function tests the subroutine for extending a matrix
        !
        use opentrustregion, only: extend_matrix

        real(rp), allocatable :: matrix(:, :)
        real(rp) :: initial_matrix(2, 2), expected(3, 3), vector1(3), vector2(3)

        ! assume tests pass
        test_extend_matrix = .true.

        ! generate a matrix and the vectors to be added, which share their last element
        ! as those of the reduced-space Hessian do, so that the corner does not depend
        ! on which of them fills it
        call random_number(initial_matrix)
        call random_number(vector1)
        call random_number(vector2)
        vector2(3) = vector1(3)
        matrix = initial_matrix

        ! the first vector becomes the new column and the second the new row
        expected(:2, :2) = initial_matrix
        expected(:, 3) = vector1
        expected(3, :) = vector2

        ! call routine and determine if dimensions and values of resulting matrix match
        call extend_matrix(matrix, vector1, vector2)
        if (size(matrix, 1) /= 3 .or. size(matrix, 2) /= 3) then
            write(stderr, *) "test_extend_matrix failed: Incorrect matrix "// &
                "dimensions after extending."
            test_extend_matrix = .false.
            return
        end if
        if (norm2(matrix - expected) > tol) then
            write(stderr, *) "test_extend_matrix failed: Incorrect matrix values "// &
                "after extending."
            test_extend_matrix = .false.
        end if

        ! deallocate matrix
        deallocate(matrix)

    end function test_extend_matrix

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

        ! call routine and determine if dimensions of resulting matrix match and
        ! whether the initial matrix is extended by the column
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

    logical(c_bool) function test_mat_diag() bind(C)
        !
        ! this function tests the subroutine which determines the eigenvalues and left
        ! and right eigenvectors of a symmetric or a general matrix
        !
        use opentrustregion, only: solver_settings_type, mat_diag

        type(solver_settings_type) :: settings
        real(rp) :: matrix(n_param, n_param), right_eigvecs(n_param, n_param), &
                    left_eigvecs(n_param, n_param), ref_eigvals(n_param), &
                    ref_eigvecs(n_param, n_param)
        complex(rp) :: eigvals(n_param)
        integer(ip) :: error, i
        type(test_context_type), target :: context

        ! assume tests pass
        test_mat_diag = .true.

        ! setup settings object
        call setup_settings(settings, context)

        ! diagonalize a random symmetric matrix as symmetric and determine whether the
        ! eigenvalues are real and the left eigenvectors are the right ones, which are
        ! eigenvectors
        matrix = generate_random_symm_matrix(n_param)
        call ref_symm_mat_diag(matrix, ref_eigvals, ref_eigvecs)
        call mat_diag(matrix, .true., eigvals, right_eigvecs, left_eigvecs, settings, &
                      error)
        if (error /= 0) then
            write(stderr, *) "test_mat_diag failed: Produced error for symmetric "// &
                "matrix."
            test_mat_diag = .false.
        end if
        if (any(abs(eigvals%re - ref_eigvals) > tol) .or. &
            any(abs(eigvals%im) > 0.0_rp)) then
            write(stderr, *) "test_mat_diag failed: Incorrect eigenvalues for "// &
                "symmetric matrix."
            test_mat_diag = .false.
        end if
        if (norm2(matmul(matrix, right_eigvecs) - &
                  right_eigvecs * spread(eigvals%re, 1, n_param)) > tol .or. &
            norm2(left_eigvecs - right_eigvecs) > 0.0_rp) then
            write(stderr, *) "test_mat_diag failed: Incorrect eigenvectors for "// &
                "symmetric matrix."
            test_mat_diag = .false.
        end if

        ! diagonalize a random upper triangular matrix, whose eigenvalues are its
        ! diagonal elements and which a symmetric diagonalization would take as the
        ! symmetric matrix of its upper triangle, as general and determine whether
        ! these eigenvalues and the right eigenvectors are returned
        matrix = generate_random_upper_triangular_matrix(n_param)
        call mat_diag(matrix, .false., eigvals, right_eigvecs, left_eigvecs, settings, &
                      error)
        if (error /= 0) then
            write(stderr, *) "test_mat_diag failed: Produced error for general matrix."
            test_mat_diag = .false.
        end if
        do i = 1, n_param
            if (minval(abs(eigvals%re - matrix(i, i))) > tol) then
                write(stderr, *) "test_mat_diag failed: Incorrect eigenvalues for "// &
                    "general matrix."
                test_mat_diag = .false.
                exit
            end if
        end do
        if (norm2(matmul(matrix, right_eigvecs) - &
                  right_eigvecs * spread(eigvals%re, 1, n_param)) > tol) then
            write(stderr, *) "test_mat_diag failed: Incorrect right eigenvectors "// &
                "for general matrix."
            test_mat_diag = .false.
        end if

    end function test_mat_diag

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

    logical(c_bool) function test_general_mat_diag() bind(C)
        !
        ! this function tests the function for determining the eigenvalues and left and
        ! right eigenvectors for a square matrix
        !
        use opentrustregion, only: solver_settings_type, general_mat_diag

        type(solver_settings_type) :: settings
        real(rp) :: matrix(n_param, n_param), right_eigvecs(n_param, n_param), &
                    left_eigvecs(n_param, n_param), eigval_diag(n_param, n_param), &
                    rotation(2, 2), rotation_right_eigvecs(2, 2), &
                    rotation_left_eigvecs(2, 2)
        complex(rp) :: eigvals(n_param), rotation_eigvals(2), right_eigvec(2), &
                       left_eigvec(2), defective_eigvals(4)
        real(rp) :: defective_matrix(4, 4), defective_right_eigvecs(4, 4), &
                    defective_left_eigvecs(4, 4)
        integer(ip) :: error, i, i_case
        type(test_context_type), target :: context
        character(len=*), parameter :: defective_case_names(2) = &
            [character(len=29) :: "defective real eigenvalue", &
             "defective complex eigenvalues"]

        ! assume tests pass
        test_general_mat_diag = .true.

        ! setup settings object
        call setup_settings(settings, context)

        ! diagonalize a random upper triangular matrix, whose eigenvalues are its real
        ! diagonal elements, and determine whether the left and right eigenvectors are
        ! eigenvectors and biorthonormal
        matrix = generate_random_upper_triangular_matrix(n_param)
        call general_mat_diag(matrix, eigvals, right_eigvecs, left_eigvecs, settings, &
                              error)
        if (error /= 0) then
            write(stderr, *) "test_general_mat_diag failed: Produced error."
            test_general_mat_diag = .false.
        end if
        eigval_diag = 0.0_rp
        do i = 1, n_param
            eigval_diag(i, i) = eigvals(i)%re
            if (minval(abs(eigvals%re - matrix(i, i))) > tol) then
                write(stderr, *) "test_general_mat_diag failed: Diagonal element "// &
                    "of triangular matrix not among the eigenvalues."
                test_general_mat_diag = .false.
            end if
        end do
        if (norm2(matmul(transpose(left_eigvecs), right_eigvecs) - &
                  identity_matrix(n_param)) > tol) then
            write(stderr, *) "test_general_mat_diag failed: Left and right "// &
                "eigenvectors are not biorthonormal."
            test_general_mat_diag = .false.
        end if
        if (norm2(matmul(matrix, right_eigvecs) - &
                  matmul(right_eigvecs, eigval_diag)) > tol) then
            write(stderr, *) "test_general_mat_diag failed: Incorrect right "// &
                "eigenvectors or eigenvalues for matrix."
            test_general_mat_diag = .false.
        end if
        if (norm2(matmul(transpose(left_eigvecs), matrix) - &
                  matmul(eigval_diag, transpose(left_eigvecs))) > tol) then
            write(stderr, *) "test_general_mat_diag failed: Incorrect left "// &
                "eigenvectors or eigenvalues for matrix."
            test_general_mat_diag = .false.
        end if

        ! diagonalize a scaled rotation, whose eigenvalues are a complex conjugate pair
        ! with the cosine as real and the sine as imaginary part, sheared by a
        ! similarity transformation, which keeps the eigenvalues but lets the left and
        ! right eigenvectors differ, and determine whether both parts are returned
        rotation = matmul(matmul(reshape([1.0_rp, 0.0_rp, &
                                          0.5_rp, 1.0_rp], [2, 2]), &
                                 reshape([0.6_rp, 0.8_rp, &
                                          -0.8_rp, 0.6_rp], [2, 2])), &
                          reshape([1.0_rp, 0.0_rp, &
                                   -0.5_rp, 1.0_rp], [2, 2]))
        call general_mat_diag(rotation, rotation_eigvals, rotation_right_eigvecs, &
                              rotation_left_eigvecs, settings, error)
        if (error /= 0) then
            write(stderr, *) "test_general_mat_diag failed: Produced error for "// &
                "complex eigenvalues."
            test_general_mat_diag = .false.
        end if
        if (any(abs(rotation_eigvals%re - 0.6_rp) > tol) .or. &
            any(abs(abs(rotation_eigvals%im) - 0.8_rp) > tol) .or. &
            abs(sum(rotation_eigvals%im)) > tol) then
            write(stderr, *) "test_general_mat_diag failed: Incorrect complex "// &
                "eigenvalues."
            test_general_mat_diag = .false.
        end if

        ! determine whether the complex left and right eigenvectors of the first
        ! eigenvalue, whose real and imaginary parts are stored in consecutive columns,
        ! are eigenvectors and biorthonormal, comparing such that a NaN fails
        right_eigvec = cmplx(rotation_right_eigvecs(:, 1), &
                             rotation_right_eigvecs(:, 2), kind=rp)
        left_eigvec = cmplx(rotation_left_eigvecs(:, 1), rotation_left_eigvecs(:, 2), &
                            kind=rp)
        if (.not. abs(dot_product(left_eigvec, right_eigvec) - 1.0_rp) <= tol) then
            write(stderr, *) "test_general_mat_diag failed: Left and right "// &
                "eigenvectors of complex eigenvalue are not biorthonormal."
            test_general_mat_diag = .false.
        end if
        if (.not. all(abs(matmul(rotation, right_eigvec) - &
                          rotation_eigvals(1) * right_eigvec) <= tol)) then
            write(stderr, *) "test_general_mat_diag failed: Incorrect right "// &
                "eigenvector for complex eigenvalue."
            test_general_mat_diag = .false.
        end if
        if (.not. all(abs(matmul(conjg(left_eigvec), rotation) - &
                          rotation_eigvals(1) * conjg(left_eigvec)) <= tol)) then
            write(stderr, *) "test_general_mat_diag failed: Incorrect left "// &
                "eigenvector for complex eigenvalue."
            test_general_mat_diag = .false.
        end if

        ! diagonalize matrices with a Jordan block of a real eigenvalue and of the
        ! scaled rotation, whose left and right eigenvectors of the defective
        ! eigenvalues are orthogonal and cannot be scaled to biorthonormality, and
        ! check that an error is returned and an error message is printed
        do i_case = 1, 2
            defective_matrix = 0.0_rp
            if (i_case == 1) then
                defective_matrix(1, 2) = 1.0_rp
                defective_matrix(3, 3) = 1.0_rp
                defective_matrix(4, 4) = 2.0_rp
            else
                defective_matrix(1:2, 1:2) = reshape([0.6_rp, 0.8_rp, &
                                                      -0.8_rp, 0.6_rp], [2, 2])
                defective_matrix(3:4, 3:4) = defective_matrix(1:2, 1:2)
                defective_matrix(1:2, 3:4) = identity_matrix(2_ip)
            end if
            call setup_error_logging(settings, context)
            call general_mat_diag(defective_matrix, defective_eigvals, &
                                  defective_right_eigvecs, defective_left_eigvecs, &
                                  settings, error)
            if (error /= 1) then
                write(stderr, *) "test_general_mat_diag failed: No error returned "// &
                    "for "//trim(defective_case_names(i_case))//"."
                test_general_mat_diag = .false.
            end if
            if (len_trim(context%log_message) == 0) then
                write(stderr, *) "test_general_mat_diag failed: No error message "// &
                    "printed for "//trim(defective_case_names(i_case))//"."
                test_general_mat_diag = .false.
            end if
        end do

    end function test_general_mat_diag

    logical(c_bool) function test_mat_min_eig() bind(C)
        !
        ! this function tests the subroutine for determining the minimum eigenvalue and
        ! corresponding eigenvector of a matrix, for both the symmetric and the
        ! non-symmetric case
        !
        use opentrustregion, only: solver_settings_type, mat_min_eig

        type(solver_settings_type) :: settings
        real(rp) :: matrix(n_param, n_param), eigval, eigvec(n_param), &
                    ref_eigvals(n_param), ref_eigvecs(n_param, n_param), ref_min_eigval
        integer(ip) :: error, i, i_case
        type(test_context_type), target :: context

        ! assume tests pass
        test_mat_min_eig = .true.

        ! setup settings object
        call setup_settings(settings, context)

        ! determine whether the lowest eigenvalue and an eigenvector belonging to it
        ! are returned for a random symmetric matrix and for a random upper triangular
        ! matrix, whose lowest eigenvalue is its lowest diagonal element
        do i_case = 1, 2
            if (i_case == 1) then
                matrix = generate_random_symm_matrix(n_param)
                call ref_symm_mat_diag(matrix, ref_eigvals, ref_eigvecs)
                ref_min_eigval = ref_eigvals(1)
            else
                matrix = generate_random_upper_triangular_matrix(n_param)
                ref_min_eigval = minval([(matrix(i, i), i=1, n_param)])
            end if
            call mat_min_eig(matrix, i_case == 1, eigval, eigvec, settings, error)
            if (error /= 0) then
                write(stderr, *) "test_mat_min_eig failed: Produced error for "// &
                    trim(merge("symmetric    ", "non-symmetric", i_case == 1))// &
                    " matrix."
                test_mat_min_eig = .false.
            end if
            if (abs(eigval - ref_min_eigval) > tol) then
                write(stderr, *) "test_mat_min_eig failed: Incorrect minimum "// &
                    "eigenvalue for "//trim(merge("symmetric    ", "non-symmetric", &
                                                  i_case == 1))//" matrix."
                test_mat_min_eig = .false.
            end if
            if (norm2(matmul(matrix, eigvec) - eigval * eigvec) > tol) then
                write(stderr, *) "test_mat_min_eig failed: Incorrect eigenvector "// &
                    "corresponding to minimum eigenvalue for "// &
                    trim(merge("symmetric    ", "non-symmetric", i_case == 1))// &
                    " matrix."
                test_mat_min_eig = .false.
            end if
            if (abs(norm2(eigvec) - 1.0_rp) > tol) then
                write(stderr, *) "test_mat_min_eig failed: Eigenvector "// &
                    "corresponding to minimum eigenvalue not normalized for "// &
                    trim(merge("symmetric    ", "non-symmetric", i_case == 1))// &
                    " matrix."
                test_mat_min_eig = .false.
            end if
        end do

    end function test_mat_min_eig

    logical(c_bool) function test_mat_lowest_eigpairs() bind(C)
        !
        ! this function tests the subroutine which returns the lowest eigenpairs of a
        ! matrix
        !
        use opentrustregion, only: solver_settings_type, mat_lowest_eigpairs, &
                                   complex_eigval_warning_msg

        integer(ip), parameter :: n_target = 2

        type(solver_settings_type) :: settings
        real(rp) :: matrix(n_param, n_param), eigvals(n_target), &
                    eigvecs(n_param, n_target), ref_eigvals(n_param), &
                    ref_eigvecs(n_param, n_param), diag_copy(n_param)
        integer(ip) :: error, i, i_case
        type(test_context_type), target :: context

        ! assume tests pass
        test_mat_lowest_eigpairs = .true.

        ! setup settings object
        call setup_settings(settings, context)

        ! determine whether the lowest eigenvalues are returned in ascending order with
        ! their eigenvectors for a random symmetric matrix and, ordered by their real
        ! parts, for a random upper triangular matrix, whose eigenvalues are its
        ! diagonal elements
        do i_case = 1, 2
            if (i_case == 1) then
                matrix = generate_random_symm_matrix(n_param)
                call ref_symm_mat_diag(matrix, ref_eigvals, ref_eigvecs)
            else
                matrix = generate_random_upper_triangular_matrix(n_param)
                diag_copy = [(matrix(i, i), i=1, n_param)]
                do i = 1, n_target
                    ref_eigvals(i) = minval(diag_copy)
                    diag_copy(minloc(diag_copy, dim=1)) = huge(1.0_rp)
                end do
            end if
            call mat_lowest_eigpairs(matrix, i_case == 1, eigvals, eigvecs, settings, &
                                     error)
            if (error /= 0) then
                write(stderr, *) "test_mat_lowest_eigpairs failed: Produced error "// &
                    "for "//trim(merge("symmetric    ", "non-symmetric", &
                                       i_case == 1))//" matrix."
                test_mat_lowest_eigpairs = .false.
            end if
            if (any(abs(eigvals - ref_eigvals(:n_target)) > tol)) then
                write(stderr, *) "test_mat_lowest_eigpairs failed: Incorrect "// &
                    "lowest eigenvalues for "// &
                    trim(merge("symmetric    ", "non-symmetric", i_case == 1))// &
                    " matrix."
                test_mat_lowest_eigpairs = .false.
            end if
            if (norm2(matmul(matrix, eigvecs) - &
                      eigvecs * spread(eigvals, 1, n_param)) > tol) then
                write(stderr, *) "test_mat_lowest_eigpairs failed: Returned "// &
                    "eigenvectors are not eigenvectors for "// &
                    trim(merge("symmetric    ", "non-symmetric", i_case == 1))// &
                    " matrix."
                test_mat_lowest_eigpairs = .false.
            end if
        end do
        if (index(context%log_message, " "//complex_eigval_warning_msg) /= 0) then
            write(stderr, *) "test_mat_lowest_eigpairs failed: Warning for a "// &
                "complex conjugate pair printed for real eigenvalues."
            test_mat_lowest_eigpairs = .false.
        end if

        ! determine whether the complex conjugate pair of lowest eigenvalues of a
        ! rotation in the plane of the first two unit vectors, below the positive
        ! diagonal elements of the other unit vectors, is returned by its real part and
        ! normalized real vectors in this plane, the pair's invariant subspace, with a
        ! warning
        matrix = 0.0_rp
        matrix(1, 2) = -1.0_rp
        matrix(2, 1) = 1.0_rp
        do i = 3, n_param
            matrix(i, i) = real(i, kind=rp)
        end do
        context%log_message = ""
        call mat_lowest_eigpairs(matrix, .false., eigvals, eigvecs, settings, error)
        if (error /= 0) then
            write(stderr, *) "test_mat_lowest_eigpairs failed: Produced error for "// &
                "complex conjugate pair."
            test_mat_lowest_eigpairs = .false.
        end if
        if (any(abs(eigvals) > tol)) then
            write(stderr, *) "test_mat_lowest_eigpairs failed: Real part of "// &
                "complex conjugate pair not returned."
            test_mat_lowest_eigpairs = .false.
        end if
        if (any(abs(norm2(eigvecs, dim=1) - 1.0_rp) > tol)) then
            write(stderr, *) "test_mat_lowest_eigpairs failed: Vectors of complex "// &
                "conjugate pair not normalized."
            test_mat_lowest_eigpairs = .false.
        end if
        if (any(abs(eigvecs(3:, :)) > tol)) then
            write(stderr, *) "test_mat_lowest_eigpairs failed: Vectors of complex "// &
                "conjugate pair not in its invariant subspace."
            test_mat_lowest_eigpairs = .false.
        end if
        if (index(context%log_message, " "//complex_eigval_warning_msg) == 0) then
            write(stderr, *) "test_mat_lowest_eigpairs failed: Warning not printed "// &
                "for complex conjugate pair."
            test_mat_lowest_eigpairs = .false.
        end if

    end function test_mat_lowest_eigpairs

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

    logical(c_bool) function test_lowest_h_diag_unit_vectors() bind(C)
        !
        ! this function tests the subroutine which fills columns with unit vectors
        ! along the lowest Hessian diagonal elements, both when any sign is accepted
        ! and when only negative curvature is wanted
        !
        use opentrustregion, only: lowest_h_diag_unit_vectors

        integer(ip), parameter :: n_param = 4, n_vectors = 3
        real(rp), parameter :: h_diag(n_param) = [1.0_rp, -2.0_rp, 3.0_rp, -4.0_rp]
        integer(ip), parameter :: ref_idx(n_vectors) = [4_ip, 2_ip, 1_ip]

        real(rp) :: unit_vectors(n_param, n_vectors), expected(n_param)
        integer(ip) :: i

        ! assume tests pass
        test_lowest_h_diag_unit_vectors = .true.

        ! every column is filled when the sign is not restricted, in order of
        ! increasing diagonal element
        call lowest_h_diag_unit_vectors(h_diag, .false., unit_vectors)
        do i = 1, n_vectors
            expected = 0.0_rp
            expected(ref_idx(i)) = 1.0_rp
            if (norm2(unit_vectors(:, i) - expected) > tol) then
                write(stderr, *) "test_lowest_h_diag_unit_vectors failed: "// &
                    "Incorrect unit vector when any sign is accepted."
                test_lowest_h_diag_unit_vectors = .false.
            end if
        end do

        ! only the two negative elements produce a vector when negative curvature is
        ! required, leaving the trailing column vanishing
        call lowest_h_diag_unit_vectors(h_diag, .true., unit_vectors)
        do i = 1, count(h_diag < 0.0_rp)
            expected = 0.0_rp
            expected(ref_idx(i)) = 1.0_rp
            if (norm2(unit_vectors(:, i) - expected) > tol) then
                write(stderr, *) "test_lowest_h_diag_unit_vectors failed: "// &
                    "Incorrect unit vector when negative curvature is required."
                test_lowest_h_diag_unit_vectors = .false.
            end if
        end do
        if (norm2(unit_vectors(:, count(h_diag < 0.0_rp) + 1:)) > tol) then
            write(stderr, *) "test_lowest_h_diag_unit_vectors failed: Column "// &
                "without a negative diagonal element does not vanish."
            test_lowest_h_diag_unit_vectors = .false.
        end if

    end function test_lowest_h_diag_unit_vectors

    logical(c_bool) function test_generate_trial_vectors() bind(C)
        !
        ! this function tests the function which generates trial vectors for the
        ! Davidson procedure
        !
        use opentrustregion, only: solver_settings_type, generate_trial_vectors, &
                                   lin_dep_trial_vector_warning_msg, error_project, &
                                   error_get_extra_trial_vectors

        character(len=*), parameter :: case_names(5) = [ &
            character(len=48) :: "Hessian with only positive diagonal elements", &
            "Hessian with a negative diagonal element", "more supplied extra trial "// &
            "vectors than fit", "extra and random trial vectors beyond the space", &
            "gradient along negative Hessian diagonal element"]
        integer(ip), parameter :: &
            n_extra_ref(5) = [0_ip, 1_ip, n_param - 1, n_param - 2, 0_ip], &
            n_trial_ref(5) = [3_ip, 4_ip, n_param, n_param, 3_ip]

        type(solver_settings_type) :: settings
        real(rp), allocatable :: red_space_basis(:, :)
        real(rp) :: grad(n_param), h_diag(n_param), grad_norm, &
                    unit_vectors(n_param, n_param), projected_extra(n_param)
        integer(ip) :: error, i_case, n_leading
        type(test_context_type), target :: context
        type(hartmann6d_fault_context_type), target :: fault_context

        ! assume tests pass
        test_generate_trial_vectors = .true.

        ! setup settings object
        call setup_settings(settings, context)
        settings%n_random_trial_vectors = 2

        ! generate gradient and positive Hessian diagonal elements
        call random_number(grad)
        grad_norm = norm2(grad)
        call random_number(h_diag)
        h_diag = h_diag + 1.0_rp

        ! the extra trial vectors of every case are leading unit vectors
        unit_vectors = identity_matrix(n_param)

        do i_case = 1, size(case_names)
            select case (i_case)
            case (2)
                ! make the first Hessian diagonal element negative, so that the unit
                ! vector along it follows the gradient direction
                h_diag(1) = -h_diag(1)
            case (3)
                ! supply the scaled leading unit vectors, which replace the unit vector
                ! along the negative element, and request one more than fits alongside
                ! the gradient direction, so that they complete the space without
                ! random trial vectors and none of them has to be dropped
                settings%n_extra_trial_vectors = n_param
                settings%n_random_trial_vectors = 0
                settings%get_extra_trial_vectors => mock_get_extra_trial_vectors
            case (4)
                ! request fewer extra trial vectors than fit alongside the gradient
                ! direction, but more random ones than fit alongside these, so that the
                ! random trial vectors give way rather than the extra ones
                settings%n_extra_trial_vectors = n_param - 2
                settings%n_random_trial_vectors = 2
            case (5)
                ! let the gradient point along the negative element, so that the unit
                ! vector along it is linearly dependent on the gradient direction
                settings%n_extra_trial_vectors = 1
                settings%get_extra_trial_vectors => null()
                grad = 0.0_rp
                call random_number(grad(1))
                grad(1) = grad(1) + 1.0_rp
                grad_norm = norm2(grad)
            end select

            ! generate trial vectors and determine whether function returns an
            ! orthonormal space of the normalized gradient, followed by the vectors
            ! which together with it span the expected extra trial vectors and by as
            ! many of the requested random trial vectors as fit alongside them
            n_leading = 1 + n_extra_ref(i_case)
            context%log_message = ""
            red_space_basis = &
                generate_trial_vectors(grad, grad_norm, h_diag, settings, error)
            if (error /= 0) then
                write(stderr, *) "test_generate_trial_vectors failed: Produced "// &
                    "error for "//trim(case_names(i_case))//"."
                test_generate_trial_vectors = .false.
            end if
            if (size(red_space_basis, 2) /= n_trial_ref(i_case)) then
                write(stderr, *) "test_generate_trial_vectors failed: Incorrect "// &
                    "number of vectors for "//trim(case_names(i_case))//"."
                test_generate_trial_vectors = .false.
                cycle
            end if
            if (norm2(matmul(transpose(red_space_basis), red_space_basis) - &
                      identity_matrix(size(red_space_basis, 2, kind=ip))) > tol) then
                write(stderr, *) "test_generate_trial_vectors failed: Reduced "// &
                    "space basis not orthonormal for "//trim(case_names(i_case))//"."
                test_generate_trial_vectors = .false.
            end if
            if (any(abs(red_space_basis(:, 1) - grad / grad_norm) > tol)) then
                write(stderr, *) "test_generate_trial_vectors failed: First vector "// &
                    "is not the normalized gradient for "//trim(case_names(i_case))//"."
                test_generate_trial_vectors = .false.
            end if
            associate (leading => red_space_basis(:, :n_leading), &
                       extra => unit_vectors(:, :n_extra_ref(i_case)))
                if (norm2(extra - matmul(leading, &
                                         matmul(transpose(leading), extra))) > tol) then
                    write(stderr, *) "test_generate_trial_vectors failed: Leading "// &
                        "vectors do not span the extra trial vectors for "// &
                        trim(case_names(i_case))//"."
                    test_generate_trial_vectors = .false.
                end if
            end associate
            if (i_case == 3 .and. index( &
                context%log_message, " "//lin_dep_trial_vector_warning_msg) /= 0) then
                write(stderr, *) "test_generate_trial_vectors failed: Supplied "// &
                    "extra trial vector dropped although no more were requested "// &
                    "than fit alongside the gradient direction."
                test_generate_trial_vectors = .false.
            end if
        end do

        ! supply an extra trial vector along the first unit vector with a projection
        ! which removes the component along the symmetric combination of the first two
        ! unit vectors and a gradient without this component, as the projection
        ! requires of it, and check that the extra trial vector enters the space
        ! projected, along the difference of the first two unit vectors, and that no
        ! trial vector has this component
        call setup_settings(settings, context)
        settings%project => mock_project_out_pair
        settings%get_extra_trial_vectors => mock_get_extra_trial_vectors
        settings%n_extra_trial_vectors = 1
        settings%n_random_trial_vectors = 2
        call random_number(grad)
        grad(1:2) = grad(1:2) - 0.5_rp * sum(grad(1:2))
        grad_norm = norm2(grad)
        projected_extra = 0.0_rp
        projected_extra(1:2) = [1.0_rp, -1.0_rp] / sqrt(2.0_rp)
        red_space_basis = &
            generate_trial_vectors(grad, grad_norm, h_diag, settings, error)
        if (error /= 0) then
            write(stderr, *) "test_generate_trial_vectors failed: Produced error "// &
                "with projection."
            test_generate_trial_vectors = .false.
        else if (norm2(projected_extra - matmul(red_space_basis(:, :2), matmul( &
            transpose(red_space_basis(:, :2)), projected_extra))) > tol .or. &
                 any(abs(sum(red_space_basis(1:2, :), dim=1)) > tol)) then
            write(stderr, *) "test_generate_trial_vectors failed: Trial vectors "// &
                "not projected."
            test_generate_trial_vectors = .false.
        end if

        ! check that the errors of a failing extra trial vector function and of a
        ! failing projection are reported with their origins and that the function
        ! result is allocated nevertheless, since the caller assigns it before it
        ! inspects the error
        call setup_settings(settings, fault_context)
        settings%get_extra_trial_vectors => faulty_get_extra_trial_vectors
        settings%project => faulty_project
        do i_case = 1, 2
            fault_context%n_callback_calls = 0
            fault_context%fault_callback = &
                merge(fault_get_extra_trial_vectors, fault_project, i_case == 1)
            fault_context%fault_call = 1
            red_space_basis = &
                generate_trial_vectors(grad, grad_norm, h_diag, settings, error)
            if (error /= merge(error_get_extra_trial_vectors, error_project, &
                               i_case == 1) + 1) then
                write(stderr, *) "test_generate_trial_vectors failed: Error of "// &
                    "failing "// &
                    trim(fault_callback_names(fault_context%fault_callback))// &
                    " not reported with its origin."
                test_generate_trial_vectors = .false.
            end if
            if (.not. allocated(red_space_basis)) then
                write(stderr, *) "test_generate_trial_vectors failed: Reduced "// &
                    "space basis not allocated for failing "// &
                    trim(fault_callback_names(fault_context%fault_callback))//"."
                test_generate_trial_vectors = .false.
            end if
        end do

    end function test_generate_trial_vectors

    logical(c_bool) function test_generate_random_trial_vectors() bind(C)
        !
        ! this function tests the function which generates random trial vectors for the
        ! Davidson procedure
        !
        use opentrustregion, only: solver_settings_type, &
                                   generate_random_trial_vectors, error_project

        integer(ip), parameter :: n_fill = 2

        type(solver_settings_type) :: settings
        real(rp), allocatable :: red_space_basis(:, :)
        real(rp) :: first_vector(n_param)
        integer(ip) :: error, i, j
        type(test_context_type), target :: context
        type(hartmann6d_fault_context_type), target :: fault_context

        ! assume tests pass
        test_generate_random_trial_vectors = .true.

        ! setup settings object
        call setup_settings(settings, context)

        ! allocate reduced space basis and set first normalized basis vector, which is
        ! the leading column the routine has to leave alone and orthogonalize the
        ! randomly filled trailing columns against
        allocate(red_space_basis(n_param, n_fill + 1))
        call random_number(red_space_basis(:, 1))
        red_space_basis(:, 1) = red_space_basis(:, 1) / norm2(red_space_basis(:, 1))
        first_vector = red_space_basis(:, 1)

        ! generate trial vectors and determine whether function returns orthonormal
        ! trial vectors and leaves the first basis vector unchanged
        call generate_random_trial_vectors(red_space_basis, n_fill, settings, error)
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
        ! starting from a leading basis vector without this component, and determine
        ! whether the generated vectors do not have this component
        settings%project => mock_project_out_pair
        red_space_basis(1:2, 1) = red_space_basis(1:2, 1) - &
                                  0.5_rp * sum(red_space_basis(1:2, 1))
        red_space_basis(:, 1) = red_space_basis(:, 1) / norm2(red_space_basis(:, 1))
        call generate_random_trial_vectors(red_space_basis, n_fill, settings, error)
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

        ! check that the error of a projection function which fails at its first call
        ! is reported with its origin
        call setup_settings(settings, fault_context)
        settings%project => faulty_project
        fault_context%fault_callback = fault_project
        fault_context%fault_call = 1
        call generate_random_trial_vectors(red_space_basis, n_fill, settings, error)
        if (error /= error_project + 1) then
            write(stderr, *) "test_generate_random_trial_vectors failed: Did not "// &
                "report the error of a failing projection function with its origin."
            test_generate_random_trial_vectors = .false.
        end if
        deallocate(red_space_basis)

        ! fill a trial vector behind the first unit vector with a projection onto the
        ! first unit vector, so that every attempt produces a linearly dependent
        ! vector, and determine whether an error is returned and an error message is
        ! printed once the maximum number of attempts is reached
        call setup_settings(settings, context)
        allocate(red_space_basis(n_param, 2))
        red_space_basis = 0.0_rp
        red_space_basis(1, 1) = 1.0_rp
        settings%project => mock_project_first_component
        call setup_error_logging(settings, context)
        call generate_random_trial_vectors(red_space_basis, 1_ip, settings, error)
        if (error /= 1) then
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

    logical(c_bool) function test_block_davidson() bind(C)
        !
        ! this function tests the Davidson procedure for both a single tracked
        ! eigenpair and a larger block
        !
        use opentrustregion, only: hess_x_type, stability_settings_type, &
                                   block_davidson, error_hess_x, error_conv_check

        integer(ip), parameter :: n_block = 3
        character(len=*), parameter :: jacobi_davidson_case_names(5) = &
            [character(len=26) :: "for a single eigenpair", "for a block of two", &
             "with a residual tolerance", "by the Davidson method", &
             "before the requested start"]
        real(rp), parameter :: res_tol = 1e-8_rp, asymm_perturb = 0.1_rp
        real(rp) :: well_separated(n_param, n_param), ref_eigvals(n_param), &
                    eigvals(n_block), eigvecs(n_param, n_block), h_diag(n_param), &
                    ref_eigvecs(n_param, n_param), ref_eigvals_re(n_param), &
                    ref_eigvals_im(n_param), left_eigvecs(n_param, n_param)
        real(rp), allocatable :: trial_space(:, :)
        integer(ip) :: error, i, j, min_idx, i_case
        logical :: sorted(n_param)
        logical :: converged
        procedure(hess_x_type), pointer :: hess_x_funptr
        type(stability_settings_type) :: settings
        type(stability_state_type), target :: state
        type(quadratic_context_type), target :: context
        type(hartmann6d_fault_context_type), target :: fault_context

        ! assume tests pass
        test_block_davidson = .true.

        ! construct symmetric matrix with well separated eigenvalues
        do i = 1, n_param
            do j = 1, n_param
                well_separated(i, j) = 1.0_rp / real(i + j, kind=rp)
            end do
            well_separated(i, i) = well_separated(i, i) + real(i, kind=rp)
        end do
        context%hess = well_separated
        h_diag = [(context%hess(i, i), i=1, n_param)]

        ! obtain reference eigenvalues
        call ref_symm_mat_diag(context%hess, ref_eigvals, ref_eigvecs)

        ! setup settings object
        call setup_settings(settings, context)

        ! set Hessian linear transformation function pointer
        hess_x_funptr => quadratic_hess_x

        ! check block of eigenpairs with a residual threshold, which takes precedence
        ! over a convergence check, since the state only holds the lowest eigenpair
        settings%conv_check => mock_stability_conv_check
        context%stability_state => state
        context%seen_residual_rms = -1.0_rp
        trial_space = unit_vectors(n_block)
        call block_davidson(hess_x_funptr, h_diag, trial_space, .true., n_block, &
                            200_ip, eigvals, eigvecs, converged, settings, state, &
                            error, res_tol=res_tol)
        if (error /= 0) then
            write(stderr, *) "test_block_davidson failed: Produced error for block "// &
                "of eigenpairs."
            test_block_davidson = .false.
        end if
        if (context%seen_residual_rms >= 0.0_rp) then
            write(stderr, *) "test_block_davidson failed: Convergence check called "// &
                "despite a residual threshold."
            test_block_davidson = .false.
        end if
        settings%conv_check => null()
        context%stability_state => null()
        if (.not. converged) then
            write(stderr, *) "test_block_davidson failed: Block of eigenpairs did "// &
                "not converge."
            test_block_davidson = .false.
        end if
        do i = 1, n_block
            if (abs(eigvals(i) - ref_eigvals(i)) > tol) then
                write(stderr, *) "test_block_davidson failed: Block of eigenpairs "// &
                    "does not return the lowest eigenvalues."
                test_block_davidson = .false.
            end if
            if (norm2(matmul(context%hess, eigvecs(:, i)) - &
                      eigvals(i) * eigvecs(:, i)) > res_tol) then
                write(stderr, *) "test_block_davidson failed: Block of eigenpairs "// &
                    "does not return converged eigenvectors."
                test_block_davidson = .false.
            end if
            if (abs(norm2(eigvecs(:, i)) - 1.0_rp) > tol) then
                write(stderr, *) "test_block_davidson failed: Block of eigenpairs "// &
                    "does not return normalized eigenvectors."
                test_block_davidson = .false.
            end if
        end do

        ! check single eigenpair, with a preconditioner which has to be level-shifted
        ! to the stability threshold
        settings%conv_tol = res_tol
        settings%stability_thresh = -2e-2_rp
        settings%precond => mu_recording_precond
        context%recorded_mu = huge(1.0_rp)
        trial_space = unit_vectors(1_ip)
        call block_davidson(hess_x_funptr, h_diag, trial_space, .true., 1_ip, 200_ip, &
                            eigvals(:1), eigvecs(:, :1), converged, settings, state, &
                            error)
        if (abs(context%recorded_mu - settings%stability_thresh) > tol) then
            write(stderr, *) "test_block_davidson failed: Preconditioner not "// &
                "level-shifted to the stability threshold."
            test_block_davidson = .false.
        end if
        if (error /= 0) then
            write(stderr, *) "test_block_davidson failed: Produced error for "// &
                "single eigenpair."
            test_block_davidson = .false.
        end if
        if (.not. converged) then
            write(stderr, *) "test_block_davidson failed: Single eigenpair did not "// &
                "converge."
            test_block_davidson = .false.
        end if
        if (abs(eigvals(1) - ref_eigvals(1)) > tol) then
            write(stderr, *) "test_block_davidson failed: Single eigenpair does "// &
                "not return the lowest eigenvalue."
            test_block_davidson = .false.
        end if

        ! check that the state holds the iteration, at most one per parameter since
        ! every iteration adds a vector, the lowest eigenvalue and the RMS of its
        ! residual
        if (state%iteration < 2 .or. state%iteration > n_param) then
            write(stderr, *) "test_block_davidson failed: Iteration of single "// &
                "eigenpair not in the state."
            test_block_davidson = .false.
        end if
        if (abs(state%min_eigval - eigvals(1)) > tol .or. abs( &
            state%residual_rms - &
            norm2(matmul(context%hess, eigvecs(:, 1)) - eigvals(1) * eigvecs(:, 1)) / &
            sqrt(real(n_param, kind=rp))) > tol) then
            write(stderr, *) "test_block_davidson failed: Lowest eigenvalue or the "// &
                "RMS of its residual of single eigenpair not in the state."
            test_block_davidson = .false.
        end if

        ! check that a convergence check takes the place of the convergence threshold,
        ! which every residual satisfies here, so that only the convergence check
        ! demands further iterations until the eigenpair is converged, and that it
        ! reads the lowest eigenpair of the current iteration from the state
        call setup_settings(settings, context)
        settings%conv_tol = huge(1.0_rp)
        settings%conv_check => mock_stability_conv_check
        context%stability_state => state
        context%seen_residual_rms = -1.0_rp
        trial_space = unit_vectors(1_ip)
        call block_davidson(hess_x_funptr, h_diag, trial_space, .true., 1_ip, 200_ip, &
                            eigvals(:1), eigvecs(:, :1), converged, settings, state, &
                            error)
        if (error /= 0) then
            write(stderr, *) "test_block_davidson failed: Produced error with "// &
                "convergence check."
            test_block_davidson = .false.
        end if
        if (.not. converged .or. state%iteration < 2 .or. &
            norm2(matmul(context%hess, eigvecs(:, 1)) - eigvals(1) * eigvecs(:, 1)) / &
            sqrt(real(n_param, kind=rp)) >= tol * abs(eigvals(1))) then
            write(stderr, *) "test_block_davidson failed: Convergence check not "// &
                "used instead of the convergence threshold."
            test_block_davidson = .false.
        end if
        if (abs(context%seen_residual_rms - state%residual_rms) > tol * &
            state%residual_rms .or. state%residual_rms <= 0.0_rp) then
            write(stderr, *) "test_block_davidson failed: Convergence check did "// &
                "not read the lowest eigenpair of the current iteration from the state."
            test_block_davidson = .false.
        end if
        context%stability_state => null()

        ! check that a trial space which grows to span the full space with a
        ! convergence threshold that no residual satisfies converges, since its
        ! eigenpair is exact
        call setup_settings(settings, context)
        settings%conv_tol = 0.0_rp
        trial_space = unit_vectors(1_ip)
        call block_davidson(hess_x_funptr, h_diag, trial_space, .true., 1_ip, 200_ip, &
                            eigvals(:1), eigvecs(:, :1), converged, settings, state, &
                            error)
        if (error /= 0) then
            write(stderr, *) "test_block_davidson failed: Produced error when "// &
                "trial space spans the full space."
            test_block_davidson = .false.
        end if
        if (.not. converged .or. size(trial_space, 2) /= n_param) then
            write(stderr, *) "test_block_davidson failed: Did not converge once "// &
                "the trial space spans the full space."
            test_block_davidson = .false.
        end if
        if (abs(eigvals(1) - ref_eigvals(1)) > tol) then
            write(stderr, *) "test_block_davidson failed: Lowest eigenvalue not "// &
                "returned when trial space spans the full space."
            test_block_davidson = .false.
        end if

        ! check that a preconditioner which only returns the vector already spanning
        ! the trial space leaves nothing to expand it by, which ends the iterations
        ! unconverged without an error and without adding the linearly dependent vector
        call setup_settings(settings, context)
        settings%conv_tol = 1e-8_rp
        settings%precond => mock_precond_first_unit_vector
        trial_space = unit_vectors(1_ip)
        call block_davidson(hess_x_funptr, h_diag, trial_space, .true., 1_ip, 200_ip, &
                            eigvals(:1), eigvecs(:, :1), converged, settings, state, &
                            error)
        if (error /= 0) then
            write(stderr, *) "test_block_davidson failed: Produced error for "// &
                "linearly dependent trial vector."
            test_block_davidson = .false.
        end if
        if (converged) then
            write(stderr, *) "test_block_davidson failed: Converged although the "// &
                "trial space cannot be expanded."
            test_block_davidson = .false.
        end if
        if (size(trial_space, 2) /= 1) then
            write(stderr, *) "test_block_davidson failed: Expanded trial space by "// &
                "linearly dependent trial vector."
            test_block_davidson = .false.
        end if

        ! check that the Jacobi-Davidson method, switched on from the first iteration,
        ! solves the correction equations for a single tracked eigenpair, which takes
        ! Hessian linear transformations beyond the one of every trial vector that the
        ! Davidson method takes once they are solved to a relative tolerance below one
        ! from the second iteration, and that it is used neither for a block of two
        ! eigenpairs, whose trial space only spans the full space after this
        ! iteration, nor with a residual tolerance, nor by the Davidson method, nor
        ! before the iteration at which it is requested to start
        do i_case = 1, size(jacobi_davidson_case_names)
            call setup_settings(settings, context)
            state%n_hess_x = 0
            settings%conv_tol = res_tol
            settings%diag_solver = "jacobi-davidson"
            settings%jacobi_davidson_start = 0
            if (i_case == 4) settings%diag_solver = "davidson"
            if (i_case == 5) settings%jacobi_davidson_start = 200
            if (i_case == 2) then
                trial_space = unit_vectors(2_ip)
                call block_davidson(hess_x_funptr, h_diag, trial_space, .true., 2_ip, &
                                    200_ip, eigvals(:2), eigvecs(:, :2), converged, &
                                    settings, state, error)
            else if (i_case == 3) then
                trial_space = unit_vectors(1_ip)
                call block_davidson(hess_x_funptr, h_diag, trial_space, .true., 1_ip, &
                                    200_ip, eigvals(:1), eigvecs(:, :1), converged, &
                                    settings, state, error, res_tol=res_tol)
            else
                trial_space = unit_vectors(1_ip)
                call block_davidson(hess_x_funptr, h_diag, trial_space, .true., 1_ip, &
                                    200_ip, eigvals(:1), eigvecs(:, :1), converged, &
                                    settings, state, error)
            end if
            if (error /= 0 .or. .not. converged) then
                write(stderr, *) "test_block_davidson failed: Produced error or "// &
                    "did not converge with Jacobi-Davidson settings "// &
                    trim(jacobi_davidson_case_names(i_case))//"."
                test_block_davidson = .false.
            end if
            if ((state%n_hess_x > size(trial_space, 2)) .neqv. i_case == 1) then
                write(stderr, *) "test_block_davidson failed: Jacobi-Davidson "// &
                    "correction equations "// &
                    trim(merge("not solved", "solved    ", i_case == 1))//" "// &
                    trim(jacobi_davidson_case_names(i_case))//"."
                test_block_davidson = .false.
            end if
        end do

        ! construct non-symmetric transformation
        do i = 1, n_param
            do j = i + 1, n_param
                context%hess(i, j) = context%hess(i, j) + &
                                     asymm_perturb / real(i + j, kind=rp)
                context%hess(j, i) = context%hess(j, i) - &
                                     asymm_perturb / real(i + j, kind=rp)
            end do
        end do

        ! obtain reference eigenvalues, which are sorted by their real parts since they
        ! are returned in no particular order
        call ref_general_mat_diag(context%hess, ref_eigvals_re, ref_eigvals_im, &
                                  ref_eigvecs, left_eigvecs)
        sorted = .false.
        do i = 1, n_param
            min_idx = minloc(ref_eigvals_re, dim=1, mask=.not. sorted)
            ref_eigvals(i) = ref_eigvals_re(min_idx)
            sorted(min_idx) = .true.
        end do

        ! check block of eigenpairs
        call setup_settings(settings, context)
        trial_space = unit_vectors(n_block)
        call block_davidson(hess_x_funptr, h_diag, trial_space, .false., n_block, &
                            200_ip, eigvals, eigvecs, converged, settings, state, &
                            error, res_tol=res_tol)
        if (error /= 0) then
            write(stderr, *) "test_block_davidson failed: Produced error for "// &
                "non-symmetric block of eigenpairs."
            test_block_davidson = .false.
        end if
        if (.not. converged) then
            write(stderr, *) "test_block_davidson failed: Non-symmetric block of "// &
                "eigenpairs did not converge."
            test_block_davidson = .false.
        end if
        do i = 1, n_block
            if (abs(eigvals(i) - ref_eigvals(i)) > tol) then
                write(stderr, *) "test_block_davidson failed: Non-symmetric block "// &
                    "of eigenpairs does not return the eigenvalues with the lowest "// &
                    "real parts."
                test_block_davidson = .false.
            end if
            if (norm2(matmul(context%hess, eigvecs(:, i)) - &
                      eigvals(i) * eigvecs(:, i)) > res_tol) then
                write(stderr, *) "test_block_davidson failed: Non-symmetric block "// &
                    "of eigenpairs does not return converged eigenvectors."
                test_block_davidson = .false.
            end if
        end do

        ! test behavior when stopping at instability by building the well-separated
        ! symmetric matrix and forcing its first diagonal entry to be clearly unstable,
        ! letting a single iteration with a tight residual tolerance distinguish the
        ! two settings unambiguously
        context%hess = well_separated
        context%hess(1, 1) = -1.0_rp
        h_diag = [(context%hess(i, i), i=1, n_param)]
        call setup_settings(settings, context)
        settings%conv_tol = 1e-8_rp

        ! when not stopping at instability, a single iteration must not report
        ! convergence, since the residual is far from the tight tolerance and the
        ! eigenvalue is not consulted on its own
        trial_space = unit_vectors(1_ip)
        call block_davidson(hess_x_funptr, h_diag, trial_space, .true., 1_ip, 1_ip, &
                            eigvals(:1), eigvecs(:, :1), converged, settings, state, &
                            error)
        if (error /= 0) then
            write(stderr, *) "test_block_davidson failed: Produced error when not "// &
                "stopping at instability."
            test_block_davidson = .false.
        end if
        if (converged) then
            write(stderr, *) "test_block_davidson failed: Converged after a single "// &
                "iteration with an unconverged residual although not stopping at "// &
                "instability."
            test_block_davidson = .false.
        end if

        ! when stopping at instability, the same single iteration must instead report
        ! convergence immediately, since the (exactly known) eigenvalue is already
        ! below the instability threshold
        settings%stop_on_instability = .true.
        trial_space = unit_vectors(1_ip)
        call block_davidson(hess_x_funptr, h_diag, trial_space, .true., 1_ip, 1_ip, &
                            eigvals(:1), eigvecs(:, :1), converged, settings, state, &
                            error)
        if (error /= 0) then
            write(stderr, *) "test_block_davidson failed: Produced error when "// &
                "stopping at instability."
            test_block_davidson = .false.
        end if
        if (.not. converged) then
            write(stderr, *) "test_block_davidson failed: Did not converge on an "// &
                "established instability when stopping at instability."
            test_block_davidson = .false.
        end if
        if (abs(eigvals(1) + 1.0_rp) > tol) then
            write(stderr, *) "test_block_davidson failed: Incorrect Ritz value "// &
                "when stopping at instability."
            test_block_davidson = .false.
        end if

        ! a stability threshold below the eigenvalue does not establish an instability,
        ! so that the same single iteration must not report convergence
        settings%stability_thresh = -2.0_rp
        trial_space = unit_vectors(1_ip)
        call block_davidson(hess_x_funptr, h_diag, trial_space, .true., 1_ip, 1_ip, &
                            eigvals(:1), eigvecs(:, :1), converged, settings, state, &
                            error)
        if (error /= 0) then
            write(stderr, *) "test_block_davidson failed: Produced error when "// &
                "stopping at instability below the eigenvalue."
            test_block_davidson = .false.
        end if
        if (converged) then
            write(stderr, *) "test_block_davidson failed: Converged when stopping "// &
                "at instability with a stability threshold below the eigenvalue."
            test_block_davidson = .false.
        end if

        ! check that the errors of a failing Hessian linear transformation and a
        ! failing convergence check are reported with their origins
        fault_context%hess = hartmann6d_hessian(saddle_point)
        hess_x_funptr => faulty_hess_x
        do i_case = 1, 2
            call setup_settings(settings, fault_context)
            settings%conv_check => faulty_stability_conv_check
            fault_context%stability_state => state
            fault_context%n_callback_calls = 0
            fault_context%fault_callback = &
                merge(fault_hess_x, fault_stability_conv_check, i_case == 1)
            fault_context%fault_call = 1
            trial_space = unit_vectors(1_ip)
            call block_davidson(hess_x_funptr, &
                                [(fault_context%hess(i, i), i=1, n_param)], &
                                trial_space, .true., 1_ip, 200_ip, eigvals(:1), &
                                eigvecs(:, :1), converged, settings, state, error)
            if (error /= merge(error_hess_x, error_conv_check, i_case == 1) + 1) then
                write(stderr, *) "test_block_davidson failed: Error of failing "// &
                    trim(merge("Hessian linear transformation", &
                               "convergence check            ", i_case == 1))// &
                    " not reported with its origin."
                test_block_davidson = .false.
            end if
        end do

    contains

        function unit_vectors(n_vectors) result(vectors)
            !
            ! this function returns the first unit vectors of the given number
            !
            integer(ip), intent(in) :: n_vectors
            real(rp), allocatable :: vectors(:, :)

            integer(ip) :: i_vector

            allocate(vectors(n_param, n_vectors))
            vectors = 0.0_rp
            do i_vector = 1, n_vectors
                vectors(i_vector, i_vector) = 1.0_rp
            end do

        end function unit_vectors

    end function test_block_davidson

    logical(c_bool) function test_get_stability_trial_space() bind(C)
        !
        ! this function tests the subroutine which generates the initial trial space
        ! for the stability check
        !
        use opentrustregion, only: stability_settings_type, get_stability_trial_space, &
                                   error_get_extra_trial_vectors, error_project

        character(len=*), parameter :: case_names(4) = [ &
            character(len=30) :: "Hessian diagonal leading block", "supplied extra "// &
            "trial vectors", "no extra trial vectors", "vanishing extra trial vectors"]
        integer(ip), parameter :: n_trial_ref(4) = [3_ip, n_param, 2_ip, 1_ip]

        type(stability_settings_type) :: settings
        real(rp) :: h_diag(n_param), h_diag_copy(n_param), projected_extra(n_param)
        real(rp), allocatable :: red_space_basis(:, :)
        integer(ip) :: error, i, i_case, ref_idx(2), leading_idx(2)
        type(hartmann6d_context_type), target :: context
        type(hartmann6d_fault_context_type), target :: fault_context

        ! assume tests pass
        test_get_stability_trial_space = .true.

        ! setup settings object
        call setup_settings(settings, context)
        settings%n_extra_trial_vectors = 2
        settings%n_random_trial_vectors = 1

        ! take the Hessian diagonal elements at a minimum, which are positive, and
        ! independently determine the two lowest
        context%hess = hartmann6d_hessian(minimum1)
        h_diag = [(context%hess(i, i), i=1, n_param)]
        h_diag_copy = h_diag
        do i = 1, 2
            ref_idx(i) = minloc(h_diag_copy, dim=1)
            h_diag_copy(ref_idx(i)) = huge(1.0_rp)
        end do

        do i_case = 1, size(case_names)
            select case (i_case)
            case (2)
                ! supply extra trial vectors along the scaled first two unit vectors
                ! and request as many random trial vectors as there are parameters, so
                ! that the trial space is capped at the full space
                settings%n_extra_trial_vectors = 2
                settings%n_random_trial_vectors = n_param
                settings%get_extra_trial_vectors => mock_get_extra_trial_vectors
            case (3)
                ! request no extra trial vectors, so that the trial space consists of
                ! the random trial vectors alone
                settings%n_extra_trial_vectors = 0
                settings%n_random_trial_vectors = 2
                settings%get_extra_trial_vectors => null()
            case (4)
                ! supply extra trial vectors which all vanish, so that the trial space
                ! consists of the random trial vector alone
                settings%n_extra_trial_vectors = 2
                settings%n_random_trial_vectors = 1
                settings%get_extra_trial_vectors => mock_vanishing_extra_trial_vectors
            end select

            ! generate trial space and check its size and orthonormality
            context%log_message = ""
            call get_stability_trial_space(h_diag, red_space_basis, settings, error)
            if (error /= 0) then
                write(stderr, *) "test_get_stability_trial_space failed: Produced "// &
                    "error for "//trim(case_names(i_case))//"."
                test_get_stability_trial_space = .false.
                cycle
            end if
            if (size(red_space_basis, 2) /= n_trial_ref(i_case)) then
                write(stderr, *) "test_get_stability_trial_space failed: Incorrect "// &
                    "number of trial vectors for "//trim(case_names(i_case))//"."
                test_get_stability_trial_space = .false.
                cycle
            end if
            if (norm2(matmul(transpose(red_space_basis), red_space_basis) - &
                      identity_matrix(size(red_space_basis, 2, kind=ip))) > tol) then
                write(stderr, *) "test_get_stability_trial_space failed: Trial "// &
                    "space not orthonormal for "//trim(case_names(i_case))//"."
                test_get_stability_trial_space = .false.
            end if

            ! the leading vectors are the unit vectors along the lowest Hessian
            ! diagonal elements or along the supplied extra trial vectors, up to their
            ! sign
            if (any(i_case == [1, 2])) then
                leading_idx = merge([1_ip, 2_ip], ref_idx, i_case == 2)
                if (any(abs(abs([(red_space_basis(leading_idx(i), i), i=1, 2)]) - &
                            1.0_rp) > tol)) then
                    write(stderr, *) "test_get_stability_trial_space failed: "// &
                        "Leading vectors not along the expected unit vectors for "// &
                        trim(case_names(i_case))//"."
                    test_get_stability_trial_space = .false.
                end if
            end if

        end do

        ! with no random trial vectors requested either, dropping every extra trial
        ! vector leaves nothing to start the Davidson iterations from, which returns an
        ! error and prints an error message
        settings%get_extra_trial_vectors => mock_vanishing_extra_trial_vectors
        settings%n_extra_trial_vectors = 2
        settings%n_random_trial_vectors = 0
        call setup_error_logging(settings, context)
        call get_stability_trial_space(h_diag, red_space_basis, settings, error)
        if (error /= 1) then
            write(stderr, *) "test_get_stability_trial_space failed: No error "// &
                "returned when every trial vector is dropped."
            test_get_stability_trial_space = .false.
        end if
        if (len_trim(context%log_message) == 0) then
            write(stderr, *) "test_get_stability_trial_space failed: No error "// &
                "message printed when every trial vector is dropped."
            test_get_stability_trial_space = .false.
        end if

        ! supply an extra trial vector along the first unit vector with a projection
        ! which removes the component along the symmetric combination of the first two
        ! unit vectors and check that it leads the trial space projected, along the
        ! difference of the first two unit vectors, and that no trial vector has this
        ! component
        call setup_settings(settings, context)
        settings%project => mock_project_out_pair
        settings%get_extra_trial_vectors => mock_get_extra_trial_vectors
        settings%n_extra_trial_vectors = 1
        settings%n_random_trial_vectors = 2
        projected_extra = 0.0_rp
        projected_extra(1:2) = [1.0_rp, -1.0_rp] / sqrt(2.0_rp)
        call get_stability_trial_space(h_diag, red_space_basis, settings, error)
        if (error /= 0) then
            write(stderr, *) "test_get_stability_trial_space failed: Produced "// &
                "error with projection."
            test_get_stability_trial_space = .false.
        else if (abs(abs(dot_product(red_space_basis(:, 1), projected_extra)) - &
                     1.0_rp) > tol .or. &
                 any(abs(sum(red_space_basis(1:2, :), dim=1)) > tol)) then
            write(stderr, *) "test_get_stability_trial_space failed: Trial vectors "// &
                "not projected."
            test_get_stability_trial_space = .false.
        end if

        ! check that the errors of a failing extra trial vector function and of a
        ! failing projection are reported with their origins
        call setup_settings(settings, fault_context)
        settings%get_extra_trial_vectors => faulty_get_extra_trial_vectors
        settings%project => faulty_project
        do i_case = 1, 2
            fault_context%n_callback_calls = 0
            fault_context%fault_callback = &
                merge(fault_get_extra_trial_vectors, fault_project, i_case == 1)
            fault_context%fault_call = 1
            call get_stability_trial_space(h_diag, red_space_basis, settings, error)
            if (error /= merge(error_get_extra_trial_vectors, error_project, &
                               i_case == 1) + 1) then
                write(stderr, *) "test_get_stability_trial_space failed: Error of "// &
                    "failing "// &
                    trim(fault_callback_names(fault_context%fault_callback))// &
                    " not reported with its origin."
                test_get_stability_trial_space = .false.
            end if
        end do

    end function test_get_stability_trial_space

    logical(c_bool) function test_orthogonalize_trial_vectors() bind(C)
        !
        ! this function tests the subroutine which orthogonalizes trial vectors
        !
        use opentrustregion, only: stability_settings_type, &
                                   orthogonalize_trial_vectors, &
                                   lin_dep_trial_vector_warning_msg

        type(stability_settings_type) :: settings
        real(rp), allocatable :: red_space_basis(:, :)
        real(rp) :: vectors(n_param, 2), weights(2)
        integer(ip) :: error
        type(test_context_type), target :: context

        ! assume tests pass
        test_orthogonalize_trial_vectors = .true.

        ! setup settings object
        call setup_settings(settings, context)

        ! allocate reduced space basis with 4 vectors, where the first and third are
        ! random, the second vanishes and the fourth is linearly dependent on the
        ! random ones, so that the third vector has to take the place of the second
        allocate(red_space_basis(n_param, 4))
        call random_number(vectors)
        call random_number(weights)
        red_space_basis(:, 1) = vectors(:, 1)
        red_space_basis(:, 2) = 0.0_rp
        red_space_basis(:, 3) = vectors(:, 2)
        red_space_basis(:, 4) = matmul(vectors, weights)

        ! orthogonalize trial vectors and check that an orthonormal basis of the random
        ! vectors remains and that a warning is printed for the removal of the
        ! linearly dependent vector
        call orthogonalize_trial_vectors(red_space_basis, settings, error)
        if (error /= 0) then
            write(stderr, *) "test_orthogonalize_trial_vectors failed: Produced error."
            test_orthogonalize_trial_vectors = .false.
        end if
        if (size(red_space_basis, 2) /= 2) then
            write(stderr, *) "test_orthogonalize_trial_vectors failed: Linearly "// &
                "dependent or vanishing vector was not removed."
            test_orthogonalize_trial_vectors = .false.
        end if
        if (norm2(matmul(transpose(red_space_basis), red_space_basis) - &
                  identity_matrix(size(red_space_basis, 2, kind=ip))) > tol) then
            write(stderr, *) "test_orthogonalize_trial_vectors failed: Vectors are "// &
                "not orthonormal."
            test_orthogonalize_trial_vectors = .false.
        end if
        if (norm2(vectors - matmul(red_space_basis, matmul(transpose(red_space_basis), &
                                                           vectors))) > tol) then
            write(stderr, *) "test_orthogonalize_trial_vectors failed: Vectors do "// &
                "not span the random vectors."
            test_orthogonalize_trial_vectors = .false.
        end if
        if (index(context%log_message, " "//lin_dep_trial_vector_warning_msg) == 0) then
            write(stderr, *) "test_orthogonalize_trial_vectors failed: No warning "// &
                "printed when removing linearly dependent vector."
            test_orthogonalize_trial_vectors = .false.
        end if
        deallocate(red_space_basis)

        ! a block whose vectors all vanish leaves no vector space behind, which is a
        ! legitimate outcome for a caller that has no direction to contribute
        allocate(red_space_basis(n_param, 2))
        red_space_basis = 0.0_rp

        ! orthogonalize trial vectors and check that the vanishing vectors are removed
        call orthogonalize_trial_vectors(red_space_basis, settings, error)
        if (error /= 0) then
            write(stderr, *) "test_orthogonalize_trial_vectors failed: Produced "// &
                "error for a block of vanishing vectors."
            test_orthogonalize_trial_vectors = .false.
        end if
        if (size(red_space_basis, 2) /= 0) then
            write(stderr, *) "test_orthogonalize_trial_vectors failed: Block of "// &
                "vanishing vectors did not produce an empty vector space."
            test_orthogonalize_trial_vectors = .false.
        end if
        deallocate(red_space_basis)

    end function test_orthogonalize_trial_vectors

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

        ! generate vector to be orthogonalized, scaled to a norm between the threshold
        ! of the routine for a numerically zero vector and the one for linear
        ! dependence, so that the check for linear dependence has to be relative to
        ! the size of the vector, and orthonormal space spanned by eigenvectors of a
        ! symmetric matrix
        call random_number(vector)
        vector = 1e-15_rp * vector / norm2(vector)
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
        if (error /= 1) then
            write(stderr, *) "test_gram_schmidt failed: No error returned during "// &
                "orthogonalization for zero vector."
            test_gram_schmidt = .false.
        end if
        if (len_trim(context%log_message) == 0) then
            write(stderr, *) "test_gram_schmidt failed: No error message printed "// &
                "for zero vector."
            test_gram_schmidt = .false.
        end if

        ! orthonormalize a vector which deviates only slightly from the first vector of
        ! an orthonormal space, so that a single orthogonalization leaves overlaps
        ! above the threshold of the routine and it has to orthogonalize again, and
        ! check that the result is orthonormal to the space within this threshold
        call ref_symm_mat_diag(generate_random_symm_matrix(n_param), eigvals, eigvecs)
        space = eigvecs(:, :2)
        call random_number(vector)
        vector = space(:, 1) + 1e-8_rp * vector
        call gram_schmidt(vector, space, settings, error)
        if (error /= 0) then
            write(stderr, *) "test_gram_schmidt failed: Produced error for vector "// &
                "close to the space."
            test_gram_schmidt = .false.
        end if
        if (maxval(abs(matmul(vector, space))) > 1e-14_rp) then
            write(stderr, *) "test_gram_schmidt failed: Vector close to the space "// &
                "not orthogonalized to it."
            test_gram_schmidt = .false.
        end if
        if (abs(norm2(vector) - 1.0_rp) > tol) then
            write(stderr, *) "test_gram_schmidt failed: Vector close to the space "// &
                "not normalized."
            test_gram_schmidt = .false.
        end if

        ! define linearly dependent vector
        vector = space(:, 1)

        ! perform Gram-Schmidt orthogonalization and determine if function correctly
        ! returns the linear dependence without printing an error message, since the
        ! callers handle it
        call setup_error_logging(settings, context)
        call gram_schmidt(vector, space, settings, error)
        if (error /= error_gram_schmidt_lin_dep) then
            write(stderr, *) "test_gram_schmidt failed: Linear dependence not "// &
                "returned for linearly dependent vector."
            test_gram_schmidt = .false.
        end if
        if (len_trim(context%log_message) /= 0) then
            write(stderr, *) "test_gram_schmidt failed: Error message printed for "// &
                "linearly dependent vector."
            test_gram_schmidt = .false.
        end if

        ! define vector in space that is already complete
        call random_number(vector_small)
        space_small = identity_matrix(2_ip)

        ! perform Gram-Schmidt orthogonalization and determine if function correctly
        ! returns an error and prints an error message
        call setup_error_logging(settings, context)
        call gram_schmidt(vector_small, space_small, settings, error)
        if (error /= 1) then
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

    logical(c_bool) function test_level_shifted_diag_precond() bind(C)
        !
        ! this function tests the subroutine that constructs the level-shifted diagonal
        ! preconditioner
        !
        use opentrustregion, only: solver_settings_type, level_shifted_diag_precond, &
                                   error_precond, error_project

        real(rp) :: vector(3), mu, h_diag(3), precond_vector(3)
        type(solver_settings_type) :: settings
        integer(ip) :: error, i_case
        type(test_context_type), target :: context
        type(hartmann6d_fault_context_type), target :: fault_context

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

        ! check that the errors of a failing custom preconditioner and of a failing
        ! projection of the default preconditioner are reported with their origins
        call setup_settings(settings, fault_context)
        do i_case = 1, 2
            if (i_case == 1) then
                settings%precond => faulty_precond
            else
                settings%precond => null()
                settings%project => faulty_project
            end if
            fault_context%n_callback_calls = 0
            fault_context%fault_callback = &
                merge(fault_precond, fault_project, i_case == 1)
            fault_context%fault_call = 1
            call level_shifted_diag_precond(vector, mu, h_diag, precond_vector, &
                                            settings, error)
            if (error /= merge(error_precond, error_project, i_case == 1) + 1) then
                write(stderr, *) "test_level_shifted_diag_precond failed: Error of "// &
                    "failing "// &
                    trim(fault_callback_names(fault_context%fault_callback))// &
                    " not reported with its origin."
                test_level_shifted_diag_precond = .false.
            end if
        end do

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
        integer(ip) :: error, i_case
        type(test_context_type), target :: context
        type(hartmann6d_fault_context_type), target :: fault_context

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
            write(stderr, *) "test_abs_diag_precond failed: Produced error for "// &
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
            write(stderr, *) "test_abs_diag_precond failed: Produced error for "// &
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
            write(stderr, *) "test_abs_diag_precond failed: Produced error for "// &
                "custom preconditioner."
            test_abs_diag_precond = .false.
        end if
        if (any(abs(precond_vector - [0.0_rp, 0.0_rp, 0.0_rp]) > tol)) then
            write(stderr, *) "test_abs_diag_precond failed: Returned "// &
                "preconditioned vector not correct for custom preconditioner."
            test_abs_diag_precond = .false.
        end if

        ! check that the errors of a failing custom preconditioner and of a failing
        ! projection of the default preconditioner are reported with their origins
        call setup_settings(settings, fault_context)
        do i_case = 1, 2
            if (i_case == 1) then
                settings%precond => faulty_precond
            else
                settings%precond => null()
                settings%project => faulty_project
            end if
            fault_context%n_callback_calls = 0
            fault_context%fault_callback = &
                merge(fault_precond, fault_project, i_case == 1)
            fault_context%fault_call = 1
            call abs_diag_precond(vector, h_diag, precond_vector, settings, error)
            if (error /= merge(error_precond, error_project, i_case == 1) + 1) then
                write(stderr, *) "test_abs_diag_precond failed: Error of failing "// &
                    trim(fault_callback_names(fault_context%fault_callback))// &
                    " not reported with its origin."
                test_abs_diag_precond = .false.
            end if
        end do

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
        integer(ip) :: error, n_hess_x
        type(hartmann6d_context_type), target :: context
        type(hartmann6d_fault_context_type), target :: fault_context

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
        ! Hessian linear transformation of the projected vector are returned and
        ! whether the reported number of Hessian linear transformations agrees with the
        ! calls
        n_hess_x = 0
        context%n_hess_x_calls = 0
        call jacobi_davidson_correction(hess_x_funptr, vector, solution, eigval, &
                                        corr_vector, hess_vector, settings, n_hess_x, &
                                        error)
        if (error /= 0) then
            write(stderr, *) "test_jacobi_davidson_correction failed: Produced "// &
                "error for valid input."
            test_jacobi_davidson_correction = .false.
        end if
        proj_vector = vector - dot_product(vector, solution) * solution
        if (any(abs(hess_vector - matmul(context%hess, proj_vector)) > tol)) then
            write(stderr, *) "test_jacobi_davidson_correction failed: Returned "// &
                "Hessian linear transformation wrong for valid input."
            test_jacobi_davidson_correction = .false.
        end if
        expected_corr_vector = matmul(context%hess, proj_vector) - eigval * proj_vector
        expected_corr_vector = expected_corr_vector - &
                               dot_product(expected_corr_vector, solution) * solution
        if (any(abs(corr_vector - expected_corr_vector) > tol)) then
            write(stderr, *) "test_jacobi_davidson_correction failed: Returned "// &
                "correction vector wrong for valid input."
            test_jacobi_davidson_correction = .false.
        end if
        if (.not. check_call_counts(context, "jacobi_davidson_correction", &
                                    "for valid input", n_hess_x=n_hess_x)) &
            test_jacobi_davidson_correction = .false.

        ! calculate Jacobi-Davidson correction with a Hessian linear transformation
        ! which fails and check that its error is reported with its origin and that
        ! the failing call is still counted
        fault_context%hess = context%hess
        settings%context => fault_context
        hess_x_funptr => faulty_hess_x
        fault_context%fault_callback = fault_hess_x
        fault_context%fault_call = 1
        n_hess_x = 0
        call jacobi_davidson_correction(hess_x_funptr, vector, solution, eigval, &
                                        corr_vector, hess_vector, settings, n_hess_x, &
                                        error)
        if (error /= error_hess_x + 1) then
            write(stderr, *) "test_jacobi_davidson_correction failed: Did not "// &
                "report the error of a failing Hessian linear transformation with "// &
                "its origin."
            test_jacobi_davidson_correction = .false.
        end if
        if (.not. check_call_counts(fault_context, "jacobi_davidson_correction", &
                                    "for failing Hessian linear transformation", &
                                    n_hess_x=n_hess_x)) &
            test_jacobi_davidson_correction = .false.

    end function test_jacobi_davidson_correction

    logical(c_bool) function test_minres() bind(C)
        !
        ! this function tests the minimum residual method subroutine
        !
        use opentrustregion, only: solver_settings_type, hess_x_type, minres, &
                                   verbosity_warning, minres_max_iter_warning_msg

        type(solver_settings_type) :: settings
        procedure(hess_x_type), pointer :: hess_x_funptr
        real(rp), dimension(n_param) :: rhs, solution, vector, hess_vector, &
                                        corr_vector, guess
        real(rp) :: mu
        real(rp), parameter :: minres_tol = 1e-14_rp, large_r_tol = 1e-8_rp
        integer(ip), parameter :: n_large = 50
        real(rp), dimension(n_large) :: large_rhs, large_solution, large_vector, &
                                        large_hess_vector, large_residual
        integer(ip) :: error, i, n_hess_x, n_iter
        type(hartmann6d_context_type), target :: context
        type(quadratic_context_type), target :: quadratic_context

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
        ! solved, whether Hessian linear transformation is correct and whether the
        ! Hessian linear transformations of the iterations are counted, if rhs and
        ! solution are orthogonal (as in Jacobi-Davidson), the final vector will be
        ! orthogonal to the solution vector and consequently the Hessian linear
        ! transformation of the projected vector is equivalent to the Hessian linear
        ! transformation of the vector itself
        n_hess_x = 0
        context%n_hess_x_calls = 0
        call minres(-rhs, hess_x_funptr, solution, mu, minres_tol, vector, &
                    hess_vector, settings, n_hess_x, error)
        if (error /= 0) then
            write(stderr, *) "test_minres failed: Produced error for correction "// &
                "equation."
            test_minres = .false.
        end if
        if (context%n_hess_x_calls < 1 .or. n_hess_x /= context%n_hess_x_calls) then
            write(stderr, *) "test_minres failed: Hessian linear transformations "// &
                "of the iterations not counted for correction equation."
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
                    hess_vector, settings, n_hess_x, error)
        if (error /= 0) then
            write(stderr, *) "test_minres failed: Produced error for vanishing rhs."
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

        ! run minimum residual method from a random initial guess orthogonal to the
        ! solution, as the correction equation requires, and check if the equation is
        ! solved
        rhs = matmul(context%hess, solution) - mu * solution
        call random_number(guess)
        guess = guess - dot_product(guess, solution) * solution
        call minres(-rhs, hess_x_funptr, solution, mu, minres_tol, vector, &
                    hess_vector, settings, n_hess_x, error, guess=guess)
        if (error /= 0) then
            write(stderr, *) "test_minres failed: Produced error for initial guess."
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
        ! returned with only its own Hessian linear transformation, which is counted
        call random_number(guess)
        guess = guess - dot_product(guess, solution) * solution
        corr_vector = matmul(context%hess, guess) - mu * guess
        rhs = -(corr_vector - dot_product(corr_vector, solution) * solution)
        n_hess_x = 0
        context%n_hess_x_calls = 0
        call minres(-rhs, hess_x_funptr, solution, mu, minres_tol, vector, &
                    hess_vector, settings, n_hess_x, error, guess=guess)
        if (error /= 0) then
            write(stderr, *) "test_minres failed: Produced error for initial guess "// &
                "which solves the equation."
            test_minres = .false.
        end if
        if (context%n_hess_x_calls /= 1) then
            write(stderr, *) "test_minres failed: Initial guess which solves the "// &
                "equation not accepted after a single Hessian linear transformation."
            test_minres = .false.
        end if
        if (n_hess_x /= 1) then
            write(stderr, *) "test_minres failed: Hessian linear transformation of "// &
                "initial guess which solves the equation not counted."
            test_minres = .false.
        end if
        if (any(abs(vector - guess) > tol)) then
            write(stderr, *) "test_minres failed: Initial guess which solves the "// &
                "equation not returned."
            test_minres = .false.
        end if

        ! allow only a single iteration and check that the approximate solution is
        ! returned without an error and a warning is printed when the iteration limit
        ! is reached
        call setup_error_logging(settings, context)
        settings%verbose = verbosity_warning
        call minres(-rhs, hess_x_funptr, solution, mu, minres_tol, vector, &
                    hess_vector, settings, n_hess_x, error, max_iter=1_ip)
        if (error /= 0) then
            write(stderr, *) "test_minres failed: Produced error when iteration "// &
                "limit is reached."
            test_minres = .false.
        end if
        if (index(context%log_message, " "//minres_max_iter_warning_msg) == 0) then
            write(stderr, *) "test_minres failed: No warning printed when "// &
                "iteration limit is reached."
            test_minres = .false.
        end if
        corr_vector = matmul(context%hess, rhs) - mu * rhs
        corr_vector = corr_vector - dot_product(corr_vector, solution) * solution
        if (any(abs(vector + dot_product(rhs, corr_vector) / &
                    dot_product(corr_vector, corr_vector) * rhs) > tol)) then
            write(stderr, *) "test_minres failed: Returned solution is not the "// &
                "iterate of the single iteration when iteration limit is reached."
            test_minres = .false.
        end if
        if (any(abs(hess_vector - matmul(context%hess, vector)) > tol)) then
            write(stderr, *) "test_minres failed: Returned Hessian linear "// &
                "transformation wrong when iteration limit is reached."
            test_minres = .false.
        end if

        ! run minimum residual method from an initial guess for a vanishing right hand
        ! side and check that the solution and its Hessian linear transformation vanish
        rhs = 0.0_rp
        call minres(-rhs, hess_x_funptr, solution, mu, minres_tol, vector, &
                    hess_vector, settings, n_hess_x, error, guess=guess)
        if (error /= 0) then
            write(stderr, *) "test_minres failed: Produced error for initial guess "// &
                "and vanishing rhs."
            test_minres = .false.
        end if
        if (sum(abs(vector)) > tol) then
            write(stderr, *) "test_minres failed: Returned solution is not zero "// &
                "for initial guess and vanishing rhs."
            test_minres = .false.
        end if
        if (sum(abs(hess_vector)) > tol) then
            write(stderr, *) "test_minres failed: Returned Hessian linear "// &
                "transformation of solution is not zero for initial guess and "// &
                "vanishing rhs."
            test_minres = .false.
        end if

        ! solve the correction equation for a diagonal Hessian with a vanishing
        ! tolerance, which only the tests for tolerances below machine precision can
        ! satisfy, and check that it is solved without reaching the iteration limit,
        ! which only prints a warning
        context%hess = 0.0_rp
        do i = 1, n_param
            context%hess(i, i) = real(i, kind=rp)
        end do
        solution = 0.0_rp
        solution(1) = 1.0_rp
        mu = context%hess(1, 1)
        call random_number(rhs)
        rhs(1) = 0.0_rp
        call setup_error_logging(settings, context)
        settings%verbose = verbosity_warning
        n_hess_x = 0
        call minres(-rhs, hess_x_funptr, solution, mu, 0.0_rp, vector, hess_vector, &
                    settings, n_hess_x, error)
        if (error /= 0) then
            write(stderr, *) "test_minres failed: Produced error for vanishing "// &
                "tolerance."
            test_minres = .false.
        end if
        if (index(context%log_message, minres_max_iter_warning_msg) /= 0) then
            write(stderr, *) "test_minres failed: Iteration limit reached for "// &
                "vanishing tolerance."
            test_minres = .false.
        end if
        corr_vector = vector - dot_product(vector, solution) * solution
        corr_vector = matmul(context%hess, corr_vector) - mu * corr_vector
        corr_vector = corr_vector - dot_product(corr_vector, solution) * solution
        if (sum(abs(corr_vector + rhs)) > tol) then
            write(stderr, *) "test_minres failed: Returned solution does not solve "// &
                "Jacobi-Davidson correction equation for vanishing tolerance."
            test_minres = .false.
        end if

        ! solve the same correction equation with the iteration limit set to the
        ! number of iterations it took, one Hessian linear transformation each, and
        ! check that converging in the last allowed iteration does not print the
        ! warning for the iteration limit
        context%log_message = ""
        n_iter = n_hess_x
        call minres(-rhs, hess_x_funptr, solution, mu, 0.0_rp, vector, hess_vector, &
                    settings, n_hess_x, error, max_iter=n_iter)
        if (error /= 0) then
            write(stderr, *) "test_minres failed: Produced error for convergence "// &
                "in last allowed iteration."
            test_minres = .false.
        end if
        if (index(context%log_message, minres_max_iter_warning_msg) /= 0) then
            write(stderr, *) "test_minres failed: Iteration limit reached for "// &
                "convergence in last allowed iteration."
            test_minres = .false.
        end if

        ! solve the correction equation for the first unit vector and a diagonal
        ! Hessian with vanishing first element, which is therefore its Rayleigh
        ! quotient, and vanishing second element, so that the projected shifted Hessian
        ! is singular along the second coordinate, and with its other elements spread
        ! between 1 and 2, so that the iterations converge in all other coordinates
        ! long before the Krylov space is exhausted, for a right hand side with a large
        ! component along the second coordinate, which no step can remove, and check
        ! that the iterations stop at the least-squares solution, which reduces the
        ! residual in all other coordinates, before the step along the second
        ! coordinate grows by orders of magnitude, as it does when the iterations
        ! continue
        allocate(quadratic_context%hess(n_large, n_large))
        quadratic_context%hess = 0.0_rp
        do i = 3, n_large
            quadratic_context%hess(i, i) = &
                1.0_rp + real(i - 3, kind=rp) / real(n_large - 3, kind=rp)
        end do
        large_solution = 0.0_rp
        large_solution(1) = 1.0_rp
        large_rhs = 1.0_rp
        large_rhs(1) = 0.0_rp
        large_rhs(2) = 1e4_rp
        call setup_settings(settings, quadratic_context)
        hess_x_funptr => quadratic_hess_x
        call minres(-large_rhs, hess_x_funptr, large_solution, 0.0_rp, large_r_tol, &
                    large_vector, large_hess_vector, settings, n_hess_x, error)
        if (error /= 0) then
            write(stderr, *) "test_minres failed: Produced error for singular "// &
                "correction equation."
            test_minres = .false.
        end if
        large_residual = matmul(quadratic_context%hess, large_vector) + large_rhs
        if (norm2(large_residual(3:)) > 1e-1_rp * norm2(large_rhs(3:))) then
            write(stderr, *) "test_minres failed: Returned solution does not "// &
                "reduce the residual outside the singular direction for singular "// &
                "correction equation."
            test_minres = .false.
        end if
        if (abs(large_vector(2)) > 1e2_rp * large_rhs(2)) then
            write(stderr, *) "test_minres failed: Iterations did not stop at the "// &
                "least-squares solution of singular correction equation."
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
        integer(ip) :: error, i, n_calls, n_hess_x

        ! late Jacobi-Davidson iteration, whose tolerance 3**-30 for the correction
        ! equations lies far below the test tolerance, so that they are solved exactly
        integer(ip), parameter :: jacobi_davidson_iter = 31
        type(hartmann6d_context_type), target :: context
        type(hartmann6d_fault_context_type), target :: fault_context

        ! assume tests pass
        test_add_trial_vector = .true.

        ! start the counter of Hessian linear transformations
        n_hess_x = 0

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
        n_hess_x = 0
        call add_trial_vector(residual, level_shift, h_diag, 0_ip, solution, 0.0_rp, &
                              hess_x_funptr, red_space_basis, h_basis, settings, &
                              n_hess_x, error)
        if (error /= 0) then
            write(stderr, *) "test_add_trial_vector failed: Produced error for "// &
                "Davidson."
            test_add_trial_vector = .false.
        end if
        if (size(red_space_basis, 2) /= 3) then
            write(stderr, *) "test_add_trial_vector failed: Reduced space basis "// &
                "not extended by one vector for Davidson."
            test_add_trial_vector = .false.
        else if (size(h_basis, 2) /= 3) then
            write(stderr, *) "test_add_trial_vector failed: Hessian linear "// &
                "transformations of reduced space basis not extended by one vector "// &
                "for Davidson."
            test_add_trial_vector = .false.
        else
            expected_vec = residual / (h_diag - level_shift)
            expected_vec = expected_vec - &
                           matmul(initial_basis, matmul(expected_vec, initial_basis))
            expected_vec = expected_vec / norm2(expected_vec)
            if (any(abs(red_space_basis(:, 3) - expected_vec) > tol)) then
                write(stderr, *) "test_add_trial_vector failed: New trial vector "// &
                    "is not the orthonormalized preconditioned residual for Davidson."
                test_add_trial_vector = .false.
            end if
            if (any(abs(h_basis(:, 3) - matmul(context%hess, red_space_basis(:, 3))) > &
                    tol)) then
                write(stderr, *) "test_add_trial_vector failed: Hessian linear "// &
                    "transformation of new trial vector wrong for Davidson."
                test_add_trial_vector = .false.
            end if
            if (.not. check_call_counts(context, "add_trial_vector", "for Davidson", &
                                        n_hess_x=n_hess_x)) &
                test_add_trial_vector = .false.
        end if

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
        n_hess_x = 0
        call add_trial_vector(residual, level_shift, h_diag, jacobi_davidson_iter, &
                              solution, eigval, hess_x_funptr, red_space_basis, &
                              h_basis, settings, n_hess_x, error)
        if (error /= 0) then
            write(stderr, *) "test_add_trial_vector failed: Produced error for "// &
                "Jacobi-Davidson."
            test_add_trial_vector = .false.
        end if
        if (size(red_space_basis, 2) /= 3) then
            write(stderr, *) "test_add_trial_vector failed: Reduced space basis "// &
                "not extended by one vector for Jacobi-Davidson."
            test_add_trial_vector = .false.
        else if (size(h_basis, 2) /= 3) then
            write(stderr, *) "test_add_trial_vector failed: Hessian linear "// &
                "transformations of reduced space basis not extended by one vector "// &
                "for Jacobi-Davidson."
            test_add_trial_vector = .false.
        else
            if (any(abs(red_space_basis(:, 3) - expected_vec) > tol)) then
                write(stderr, *) "test_add_trial_vector failed: New trial vector "// &
                    "is not the orthonormalized solution of the correction "// &
                    "equations for Jacobi-Davidson."
                test_add_trial_vector = .false.
            end if
            if (any(abs(h_basis(:, 3) - matmul(context%hess, red_space_basis(:, 3))) > &
                    tol)) then
                write(stderr, *) "test_add_trial_vector failed: Hessian linear "// &
                    "transformation of new trial vector wrong for Jacobi-Davidson."
                test_add_trial_vector = .false.
            end if
            if (.not. check_call_counts(context, "add_trial_vector", &
                                        "for Jacobi-Davidson", n_hess_x=n_hess_x)) &
                test_add_trial_vector = .false.
        end if

        ! repeat in the first Jacobi-Davidson iteration, whose tolerance of 1 for the
        ! correction equations the minimum residual method always meets after its first
        ! iteration, whose solution lies along the right hand side, and determine
        ! whether the reduced space basis is extended by the orthonormalized residual
        red_space_basis = initial_basis
        h_basis = matmul(context%hess, initial_basis)
        call add_trial_vector(residual, level_shift, h_diag, 1_ip, solution, eigval, &
                              hess_x_funptr, red_space_basis, h_basis, settings, &
                              n_hess_x, error)
        if (error /= 0) then
            write(stderr, *) "test_add_trial_vector failed: Produced error for "// &
                "first Jacobi-Davidson iteration."
            test_add_trial_vector = .false.
        end if
        if (size(red_space_basis, 2) /= 3) then
            write(stderr, *) "test_add_trial_vector failed: Reduced space basis "// &
                "not extended by one vector for first Jacobi-Davidson iteration."
            test_add_trial_vector = .false.
        else
            expected_vec = residual - &
                           matmul(initial_basis, matmul(residual, initial_basis))
            expected_vec = expected_vec / norm2(expected_vec)
            if (abs(abs(dot_product(red_space_basis(:, 3), expected_vec)) - 1.0_rp) > &
                tol) then
                write(stderr, *) "test_add_trial_vector failed: New trial vector "// &
                    "is not the orthonormalized residual for first Jacobi-Davidson "// &
                    "iteration."
                test_add_trial_vector = .false.
            end if
        end if

        ! repeat with a slightly asymmetric Hessian linear transformation, whose new
        ! linear transformation then no longer respects Hessian symmetry with respect
        ! to the existing basis, and determine whether it is recalculated and counted
        red_space_basis = initial_basis
        h_basis = matmul(context%hess, initial_basis)
        hess_x_funptr => hess_x_fun_asymmetric
        context%n_hess_x_calls = 0
        n_hess_x = 0
        call add_trial_vector(residual, level_shift, h_diag, jacobi_davidson_iter, &
                              solution, eigval, hess_x_funptr, red_space_basis, &
                              h_basis, settings, n_hess_x, error)
        if (error /= 0) then
            write(stderr, *) "test_add_trial_vector failed: Produced error for "// &
                "asymmetric Hessian linear transformation."
            test_add_trial_vector = .false.
        end if
        if (size(red_space_basis, 2) /= 3) then
            write(stderr, *) "test_add_trial_vector failed: Reduced space basis "// &
                "not extended by one vector for asymmetric Hessian linear "// &
                "transformation."
            test_add_trial_vector = .false.
        else if (size(h_basis, 2) /= 3) then
            write(stderr, *) "test_add_trial_vector failed: Hessian linear "// &
                "transformations of reduced space basis not extended by one vector "// &
                "for asymmetric Hessian linear transformation."
            test_add_trial_vector = .false.
        else
            if (.not. check_call_counts( &
                context, "add_trial_vector", "for asymmetric Hessian linear "// &
                "transformation", n_hess_x=n_hess_x)) test_add_trial_vector = .false.
            call hess_x_funptr(red_space_basis(:, 3), expected_vec, error, &
                               settings%context)
            if (any(abs(h_basis(:, 3) - expected_vec) > tol)) then
                write(stderr, *) "test_add_trial_vector failed: Linear "// &
                    "transformation of new trial vector not recalculated for "// &
                    "asymmetric Hessian linear transformation."
                test_add_trial_vector = .false.
            end if
        end if

        ! let the recalculation, the last call of the asymmetric Hessian linear
        ! transformation, fail and check that its error is returned with its origin
        ! and the failing call is counted
        fault_context%hess = context%hess
        settings%context => fault_context
        red_space_basis = initial_basis
        h_basis = matmul(context%hess, initial_basis)
        fault_context%n_callback_calls = 0
        call add_trial_vector(residual, level_shift, h_diag, jacobi_davidson_iter, &
                              solution, eigval, hess_x_funptr, red_space_basis, &
                              h_basis, settings, n_hess_x, error)
        n_calls = fault_context%n_callback_calls(fault_hess_x)
        red_space_basis = initial_basis
        h_basis = matmul(context%hess, initial_basis)
        fault_context%fault_callback = fault_hess_x
        fault_context%fault_call = n_calls
        fault_context%n_callback_calls = 0
        fault_context%n_hess_x_calls = 0
        n_hess_x = 0
        call add_trial_vector(residual, level_shift, h_diag, jacobi_davidson_iter, &
                              solution, eigval, hess_x_funptr, red_space_basis, &
                              h_basis, settings, n_hess_x, error)
        if (error /= error_hess_x + 1) then
            write(stderr, *) "test_add_trial_vector failed: Did not report the "// &
                "error of a failing recalculated Hessian linear transformation "// &
                "with its origin."
            test_add_trial_vector = .false.
        end if
        if (.not. check_call_counts( &
            fault_context, "add_trial_vector", "for failing recalculated Hessian "// &
            "linear transformation", n_hess_x=n_hess_x)) test_add_trial_vector = .false.
        settings%context => context
        hess_x_funptr => hess_x_fun

        ! precondition a residual along the first basis vector with a constant Hessian
        ! diagonal so that the new vector is linearly dependent on the basis and
        ! determine whether this is returned without extending the basis
        red_space_basis = initial_basis
        h_basis = matmul(context%hess, initial_basis)
        residual = initial_basis(:, 1)
        h_diag = 1.0_rp
        call add_trial_vector(residual, 0.0_rp, h_diag, 0_ip, solution, 0.0_rp, &
                              hess_x_funptr, red_space_basis, h_basis, settings, &
                              n_hess_x, error)
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

        ! choose a Hessian which couples the first unit vector, the solution, only to
        ! the second unit vector, which is an eigenvector of the Hessian projected onto
        ! the complement of the solution, so that the solution of the correction
        ! equations lies along the second unit vector, and determine whether its linear
        ! dependence on the basis of the first two unit vectors is returned without
        ! extending the basis
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
        call add_trial_vector(residual, 0.0_rp, h_diag, jacobi_davidson_iter, &
                              solution, eigval, hess_x_funptr, red_space_basis, &
                              h_basis, settings, n_hess_x, error)
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

        ! add trial vector with a Hessian linear transformation which fails at its
        ! first call and check that its error is returned with its origin and the
        ! failing call is counted
        red_space_basis = initial_basis
        h_basis = matmul(context%hess, initial_basis)
        call random_number(residual)
        fault_context%hess = context%hess
        settings%context => fault_context
        hess_x_funptr => faulty_hess_x
        fault_context%fault_callback = fault_hess_x
        fault_context%fault_call = 1
        fault_context%n_callback_calls = 0
        fault_context%n_hess_x_calls = 0
        n_hess_x = 0
        call add_trial_vector(residual, 0.0_rp, h_diag, 0_ip, solution, 0.0_rp, &
                              hess_x_funptr, red_space_basis, h_basis, settings, &
                              n_hess_x, error)
        if (error /= error_hess_x + 1) then
            write(stderr, *) "test_add_trial_vector failed: Did not report the "// &
                "error of a failing Hessian linear transformation with its origin."
            test_add_trial_vector = .false.
        end if
        if (.not. check_call_counts(fault_context, "add_trial_vector", &
                                    "for failing Hessian linear transformation", &
                                    n_hess_x=n_hess_x)) test_add_trial_vector = .false.

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

        ! print row of results table with optional arguments, the micro iterations as
        ! Jacobi-Davidson and Davidson iterations
        call settings%print_results(1_ip, 2.0_rp, 3.0_rp, 4.0_rp, 5_ip, 2_ip, 7.0_rp, &
                                    8.0_rp)
        if (context%log_message /= "        1   |     2.00000000000000E+00   "// &
            "|   3.00E+00   |   4.00E+00  |   2 |   3  |    7.00E+00  |  8.00E+00") then
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

        ! check that a message longer than a line is split into several lines which
        ! together reproduce the message, each line is logged with a leading space,
        ! the logging function counts its calls in the host context
        context%log_message = ""
        context%n_calls = 0
        call print_message(settings, repeat("This is a long test message. ", 8), &
                           verbosity_warning)
        if (context%n_calls < 2) then
            write(stderr, *) "test_print_message failed: Long message is not split "// &
                "into several lines."
            test_print_message = .false.
        end if
        if (context%log_message /= &
            " "//trim(repeat("This is a long test message. ", 8))) then
            write(stderr, *) "test_print_message failed: Lines of long message do "// &
                "not reproduce the message."
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
            trust_radius_expand_ratio, trust_radius_shrink_factor, &
            trust_radius_expand_factor, function_unchanged_warning_msg
        use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_positive_inf, &
                                                 ieee_quiet_nan

        ! cases with the reason for accepting or rejecting the step, the size of the
        ! first rotation, the norm of the step, the ratio of actual to predicted
        ! function change, whether the micro iterations converged, and the expected
        ! acceptance and trust radius, for a starting trust radius of one
        character(len=*), parameter :: case_names(7) = &
            [character(len=47) :: "micro iterations have not converged", "ratio is "// &
             "negative", "individual rotations are too large", "ratio is too small", &
             "ratio is acceptable", "ratio is too large at the trust region boundary", &
             "ratio is too large within the trust region"]
        real(rp), parameter :: &
            case_rotations(7) = &
                [0.3_rp, 0.3_rp, 1.0_rp, 0.3_rp, 0.3_rp, 0.3_rp, 0.3_rp], &
            case_step_norms(7) = &
                [0.5_rp, 0.5_rp, 0.5_rp, 0.5_rp, 0.5_rp, 1.0_rp, 0.5_rp], &
            case_ratios(7) = [ &
                1.0_rp, -1.0_rp, 1.0_rp, 0.9_rp * trust_radius_shrink_ratio, 0.5_rp * &
                (trust_radius_shrink_ratio + trust_radius_expand_ratio), 1.1_rp * &
                trust_radius_expand_ratio, 1.1_rp * trust_radius_expand_ratio], &
            case_trust_radii(7) = [ &
                trust_radius_shrink_factor, trust_radius_shrink_factor, &
                trust_radius_shrink_factor, trust_radius_shrink_factor, 1.0_rp, &
                trust_radius_expand_factor, 1.0_rp]
        logical, parameter :: &
            case_converged(7) = &
                [.false., .true., .true., .true., .true., .true., .true.], &
            case_accepted(7) = &
                [.false., .false., .false., .true., .true., .true., .true.]

        ! cases with the objective function value, the actual and the predicted
        ! function change, whether the micro iterations converged and whether maximum
        ! precision is reached, which requires that the predicted change can not be
        ! resolved relative to the objective function value and that the actual change
        ! can not be resolved either or the step is rejected
        character(len=*), parameter :: precision_case_names(9) = &
            [character(len=68) :: "vanishing function changes", &
             "changes below the resolution of a large function value", &
             "changes above the resolution of a small function value", "vanishing "// &
             "changes of a vanishing function value", "accepted step with only the "// &
             "predicted change below the resolution", "rejected step with only the "// &
             "predicted change below the resolution", "unconverged step with only "// &
             "the predicted change below the resolution", "rejected step with a "// &
             "huge function increase", "only the actual change below the resolution"]
        real(rp), parameter :: &
            precision_case_funcs(9) = [-1.0_rp, -2000.0_rp, -1.0_rp, 0.0_rp, -1.0_rp, &
                                       -1.0_rp, -1.0_rp, -1.0_rp, -1.0_rp], &
            precision_case_actual(9) = [0.0_rp, 1e-13_rp, 1e-13_rp, 0.0_rp, 1e-3_rp, &
                                        1e-3_rp, -1e-3_rp, 1e10_rp, 1e-20_rp], &
            precision_case_predicted(9) = [0.0_rp, 1e-13_rp, 1e-13_rp, 0.0_rp, &
                                           1e-20_rp, -1e-20_rp, -1e-20_rp, -1e-6_rp, &
                                           1e-3_rp]
        logical, parameter :: &
            precision_case_converged(9) = [.true., .true., .true., .true., .true., &
                                           .true., .false., .true., .true.], &
            precision_case_reached(9) = [.true., .true., .false., .true., .false., &
                                         .true., .true., .false., .false.]

        ! cases of an actual or predicted function change which is not finite, as for
        ! an objective function or a Hessian linear transformation which overflows at
        ! the step
        character(len=*), parameter :: non_finite_case_names(3) = &
            [character(len=23) :: "overflowing function", "function not a number", &
             "overflowing prediction"]
        real(rp) :: non_finite_actual(3), non_finite_predicted(3)

        logical :: accept_step
        real(rp) :: solution(3)
        type(solver_state_type) :: state
        type(solver_settings_type) :: settings
        type(test_context_type), target :: context
        integer(ip) :: i

        ! assume tests pass
        test_accept_trust_region_step = .true.

        ! setup settings object
        call setup_settings(settings, context)

        solution = [0.3_rp, 0.3_rp, 0.3_rp]

        ! check for every case whether the step is accepted and the trust radius is
        ! reduced, kept or expanded, without maximum precision being reported as
        ! reached, the flag starts set so that a flag which is not cleared is detected
        do i = 1, size(case_names)
            state%trust_radius = 1.0_rp
            state%max_precision_reached = .true.
            solution(1) = case_rotations(i)
            accept_step = accept_trust_region_step(solution, case_step_norms(i), &
                                                   -1.0_rp, case_ratios(i), 1.0_rp, &
                                                   case_converged(i), settings, state)
            if (accept_step .neqv. case_accepted(i)) then
                write(stderr, *) "test_accept_trust_region_step failed: Step "// &
                    "acceptance incorrect when "//trim(case_names(i))//"."
                test_accept_trust_region_step = .false.
            end if
            if (abs(state%trust_radius - case_trust_radii(i)) > tol) then
                write(stderr, *) "test_accept_trust_region_step failed: Trust "// &
                    "radius not correctly updated when "//trim(case_names(i))//"."
                test_accept_trust_region_step = .false.
            end if
            if (state%max_precision_reached) then
                write(stderr, *) "test_accept_trust_region_step failed: Maximum "// &
                    "precision reported as reached when "//trim(case_names(i))//"."
                test_accept_trust_region_step = .false.
            end if
        end do

        ! check for every case of a function change which is not finite, with changes
        ! of the same sign which would otherwise accept the step, that the step is
        ! rejected and the trust radius reduced without maximum precision being
        ! reported
        non_finite_actual = [ieee_value(1.0_rp, ieee_positive_inf), &
                             ieee_value(1.0_rp, ieee_quiet_nan), 1.0_rp]
        non_finite_predicted = [1.0_rp, 1.0_rp, ieee_value(1.0_rp, ieee_positive_inf)]
        do i = 1, size(non_finite_case_names)
            state%trust_radius = 1.0_rp
            accept_step = accept_trust_region_step( &
                solution, 0.5_rp, -1.0_rp, non_finite_actual(i), &
                non_finite_predicted(i), .true., settings, state)
            if (accept_step) then
                write(stderr, *) "test_accept_trust_region_step failed: Step "// &
                    "accepted for "//trim(non_finite_case_names(i))//"."
                test_accept_trust_region_step = .false.
            end if
            if (abs(state%trust_radius - trust_radius_shrink_factor) > tol) then
                write(stderr, *) "test_accept_trust_region_step failed: Trust "// &
                    "radius not reduced for "//trim(non_finite_case_names(i))//"."
                test_accept_trust_region_step = .false.
            end if
            if (state%max_precision_reached) then
                write(stderr, *) "test_accept_trust_region_step failed: Maximum "// &
                    "precision reported as reached for "// &
                    trim(non_finite_case_names(i))//"."
                test_accept_trust_region_step = .false.
            end if
        end do

        ! check for every case whether maximum precision is reported as reached and,
        ! if so, whether the step is rejected with the trust radius left unchanged and
        ! a warning printed
        solution(1) = 0.3_rp
        do i = 1, size(precision_case_names)
            state%trust_radius = 1.0_rp
            context%log_message = ""
            accept_step = accept_trust_region_step( &
                solution, 0.5_rp, precision_case_funcs(i), precision_case_actual(i), &
                precision_case_predicted(i), precision_case_converged(i), settings, &
                state)
            if (state%max_precision_reached .neqv. precision_case_reached(i)) then
                write(stderr, *) "test_accept_trust_region_step failed: Maximum "// &
                    "precision detected incorrectly for "// &
                    trim(precision_case_names(i))//"."
                test_accept_trust_region_step = .false.
            end if
            if (.not. precision_case_reached(i)) cycle
            if (accept_step) then
                write(stderr, *) "test_accept_trust_region_step failed: Step "// &
                    "accepted at maximum precision for "// &
                    trim(precision_case_names(i))//"."
                test_accept_trust_region_step = .false.
            end if
            if (abs(state%trust_radius - 1.0_rp) > tol) then
                write(stderr, *) "test_accept_trust_region_step failed: Trust "// &
                    "radius changed at maximum precision for "// &
                    trim(precision_case_names(i))//"."
                test_accept_trust_region_step = .false.
            end if
            if (index(context%log_message, " "//function_unchanged_warning_msg) == 0) &
                then
                write(stderr, *) "test_accept_trust_region_step failed: No warning "// &
                    "printed at maximum precision for "// &
                    trim(precision_case_names(i))//"."
                test_accept_trust_region_step = .false.
            end if
        end do

    end function test_accept_trust_region_step

    logical(c_bool) function test_init_defaults() bind(C)
        !
        ! this function tests the subroutine that resolves settings defaults that
        ! depend on other settings
        !
        use opentrustregion, only: solver_settings_type, init_defaults, &
                                   verbosity_info, verbosity_debug

        ! the settings of the cases for the symmetry of the Hessian of the nested
        ! stability check settings
        character(len=*), parameter :: hess_symm_case_names(2) = &
            [character(len=44) :: "non-symmetric Hessian of the solver", &
             "non-symmetric Hessian of the nested settings"]
        type(solver_settings_type) :: settings
        type(test_context_type), target :: context, stability_context
        integer(ip) :: i_case

        ! assume tests pass
        test_init_defaults = .true.

        ! check that nested stability check settings without callback functions, a
        ! convergence check or a context of their own, with a verbosity below the
        ! solver's and an unset number of random trial vectors inherit the solver's
        ! preconditioning, projection and logging functions, context and verbosity,
        ! stop as soon as an instability is certain and use a single random trial
        ! vector
        call setup_settings(settings, context)
        settings%precond => identity_precond
        settings%project => identity_project
        settings%verbose = verbosity_debug
        settings%stability_settings%stop_on_instability = .false.
        settings%stability_settings%n_random_trial_vectors = -1
        call init_defaults(settings)
        if (.not. associated(settings%stability_settings%precond, &
                             identity_precond) .or. &
            .not. associated(settings%stability_settings%project, &
                             identity_project) .or. &
            .not. associated(settings%stability_settings%logger, logger)) then
            write(stderr, *) "test_init_defaults failed: Nested stability check "// &
                "settings did not inherit the solver's callback functions."
            test_init_defaults = .false.
        end if
        if (.not. associated(settings%stability_settings%context, context)) then
            write(stderr, *) "test_init_defaults failed: Nested stability check "// &
                "settings did not inherit the solver's context."
            test_init_defaults = .false.
        end if
        if (settings%stability_settings%verbose /= verbosity_debug) then
            write(stderr, *) "test_init_defaults failed: Nested stability check "// &
                "settings did not inherit the solver's higher verbosity."
            test_init_defaults = .false.
        end if
        if (.not. settings%stability_settings%stop_on_instability) then
            write(stderr, *) "test_init_defaults failed: Nested stability check "// &
                "settings without convergence check do not stop on instability."
            test_init_defaults = .false.
        end if
        if (settings%stability_settings%n_random_trial_vectors /= 1) then
            write(stderr, *) "test_init_defaults failed: Unset number of random "// &
                "trial vectors of nested stability check settings not set to one."
            test_init_defaults = .false.
        end if

        ! check that nested stability check settings with callback functions, a
        ! convergence check and a context of their own, a logging function while the
        ! solver has none, a verbosity above the solver's and no random trial vectors
        ! keep all of these and do not stop as soon as an instability is certain
        call setup_settings(settings, context)
        settings%precond => mock_precond
        settings%project => mock_project
        settings%logger => null()
        settings%verbose = verbosity_info
        settings%stability_settings%precond => identity_precond
        settings%stability_settings%project => identity_project
        settings%stability_settings%logger => logger
        settings%stability_settings%conv_check => mock_stability_conv_check
        settings%stability_settings%context => stability_context
        settings%stability_settings%verbose = verbosity_debug
        settings%stability_settings%stop_on_instability = .false.
        settings%stability_settings%n_random_trial_vectors = 0
        call init_defaults(settings)
        if (.not. associated(settings%stability_settings%precond, &
                             identity_precond) .or. &
            .not. associated(settings%stability_settings%project, &
                             identity_project) .or. &
            .not. associated(settings%stability_settings%logger, logger)) then
            write(stderr, *) "test_init_defaults failed: Callback functions of "// &
                "nested stability check settings replaced."
            test_init_defaults = .false.
        end if
        if (.not. associated(settings%stability_settings%context, stability_context)) &
            then
            write(stderr, *) "test_init_defaults failed: Context of nested "// &
                "stability check settings replaced."
            test_init_defaults = .false.
        end if
        if (settings%stability_settings%verbose /= verbosity_debug) then
            write(stderr, *) "test_init_defaults failed: Higher verbosity of "// &
                "nested stability check settings replaced."
            test_init_defaults = .false.
        end if
        if (settings%stability_settings%stop_on_instability) then
            write(stderr, *) "test_init_defaults failed: Nested stability check "// &
                "settings with convergence check stop on instability."
            test_init_defaults = .false.
        end if
        if (settings%stability_settings%n_random_trial_vectors /= 0) then
            write(stderr, *) "test_init_defaults failed: Vanishing number of "// &
                "random trial vectors of nested stability check settings replaced."
            test_init_defaults = .false.
        end if

        ! check that nested stability check settings without a convergence check of
        ! their own treat the Hessian as non-symmetric if the solver or they do so and
        ! then do not stop as soon as an instability is certain, which only the Ritz
        ! values of a symmetric Hessian establish
        do i_case = 1, size(hess_symm_case_names)
            call setup_settings(settings, context)
            settings%stability_settings%stop_on_instability = .false.
            select case (i_case)
            case (1)
                settings%hess_symm = .false.
            case (2)
                settings%stability_settings%hess_symm = .false.
            end select
            call init_defaults(settings)
            if (settings%stability_settings%hess_symm) then
                write(stderr, *) "test_init_defaults failed: Nested stability "// &
                    "check settings treat the Hessian as symmetric for "// &
                    trim(hess_symm_case_names(i_case))//"."
                test_init_defaults = .false.
            end if
            if (settings%stability_settings%stop_on_instability) then
                write(stderr, *) "test_init_defaults failed: Nested stability "// &
                    "check settings stop on instability for "// &
                    trim(hess_symm_case_names(i_case))//"."
                test_init_defaults = .false.
            end if
        end do

    end function test_init_defaults

    logical(c_bool) function test_solver_sanity_check() bind(C)
        !
        ! this function tests the subroutine which performs a sanity check for the
        ! solver
        !
        use opentrustregion, only: &
            solver_settings_type, solver_sanity_check, project_warning_msg, &
            random_trial_vector_warning_msg, subsystem_solver_options

        ! the subsystem solvers in mixed case, settings which the sanity check rejects
        ! and the settings that carry a projection
        character(len=*), parameter :: mixed_case_subsystem_solvers(5) = &
            [character(len=18) :: "Davidson_LS", "Davidson_AH", "Jacobi-Davidson_LS", &
             "Jacobi-Davidson_AH", "TCG"]
        character(len=*), parameter :: error_case_names(7) = &
            [character(len=55) :: "vanishing number of parameters", &
             "negative number of parameters", "vanishing number of micro iterations", &
             "negative number of random trial vectors", &
             "extra trial vector function without extra trial vectors", &
             "unknown subsystem solver", "rejected nested stability check settings"]
        character(len=*), parameter :: projection_case_names(2) = &
            [character(len=42) :: "solver settings", &
             "solver settings and inheriting nested ones"]
        type(solver_settings_type) :: settings
        integer(ip) :: n_param, error, i_solver, i_case, i_warning
        logical :: davidson
        type(test_context_type), target :: context

        ! assume tests pass
        test_solver_sanity_check = .true.

        ! check that every subsystem solver passes when given in mixed case, which is
        ! converted to lowercase, and that the number of random trial vectors is
        ! reduced to half the number of parameters with a warning for the
        ! Davidson-based subsystem solvers, which use them, and left unchanged without
        ! any message for the others
        do i_solver = 1, size(subsystem_solver_options)
            davidson = index(subsystem_solver_options(i_solver), "davidson") > 0
            call setup_settings(settings, context)
            settings%subsystem_solver = mixed_case_subsystem_solvers(i_solver)
            settings%n_random_trial_vectors = 3
            call solver_sanity_check(settings, 4_ip, error)
            if (error /= 0) then
                write(stderr, *) "test_solver_sanity_check failed: Error returned "// &
                    "for "//trim(subsystem_solver_options(i_solver))// &
                    " subsystem solver."
                test_solver_sanity_check = .false.
            end if
            if (settings%subsystem_solver /= subsystem_solver_options(i_solver)) then
                write(stderr, *) "test_solver_sanity_check failed: Not converted "// &
                    "to lowercase for "//trim(subsystem_solver_options(i_solver))// &
                    " subsystem solver."
                test_solver_sanity_check = .false.
            end if
            if (settings%n_random_trial_vectors /= merge(2_ip, 3_ip, davidson)) then
                write(stderr, *) "test_solver_sanity_check failed: Number of "// &
                    "random trial vectors wrong for "// &
                    trim(subsystem_solver_options(i_solver))//" subsystem solver."
                test_solver_sanity_check = .false.
            end if
            if ((davidson .and. &
                 index(context%log_message, " "//random_trial_vector_warning_msg) == &
                 0) .or. (.not. davidson .and. len_trim(context%log_message) > 0)) then
                write(stderr, *) "test_solver_sanity_check failed: Warning for the "// &
                    "number of random trial vectors wrong for "// &
                    trim(subsystem_solver_options(i_solver))//" subsystem solver."
                test_solver_sanity_check = .false.
            end if
        end do

        ! check that every rejected setting returns an error and prints an error
        ! message
        do i_case = 1, size(error_case_names)
            call setup_settings(settings, context)
            call setup_error_logging(settings, context)
            n_param = 4
            select case (i_case)
            case (1)
                n_param = 0
            case (2)
                n_param = -1
            case (3)
                settings%n_micro = 0
            case (4)
                settings%n_random_trial_vectors = -1
            case (5)
                settings%get_extra_trial_vectors => mock_get_extra_trial_vectors
                settings%n_extra_trial_vectors = 0
            case (6)
                settings%subsystem_solver = "unknown"
            case (7)
                ! the nested settings print the error message of their own check
                call setup_error_logging(settings%stability_settings, context)
                settings%stability_settings%diag_solver = "unknown"
            end select
            call solver_sanity_check(settings, n_param, error)
            if (error /= 1) then
                write(stderr, *) "test_solver_sanity_check failed: No error "// &
                    "returned for "//trim(error_case_names(i_case))//"."
                test_solver_sanity_check = .false.
            end if
            if (len_trim(context%log_message) == 0) then
                write(stderr, *) "test_solver_sanity_check failed: No error "// &
                    "message printed for "//trim(error_case_names(i_case))//"."
                test_solver_sanity_check = .false.
            end if
        end do

        ! check that a non-symmetric Hessian is accepted by the level-shifted Davidson
        ! subsystem solver and rejected with an error message by every other
        do i_solver = 1, size(subsystem_solver_options)
            call setup_settings(settings, context)
            call setup_error_logging(settings, context)
            settings%subsystem_solver = subsystem_solver_options(i_solver)
            settings%hess_symm = .false.
            call solver_sanity_check(settings, 4_ip, error)
            if (subsystem_solver_options(i_solver) == "davidson_ls") then
                if (error /= 0) then
                    write(stderr, *) "test_solver_sanity_check failed: "// &
                        "Non-symmetric Hessian rejected for davidson_ls subsystem "// &
                        "solver."
                    test_solver_sanity_check = .false.
                end if
            else
                if (error /= 1) then
                    write(stderr, *) "test_solver_sanity_check failed: "// &
                        "Non-symmetric Hessian not rejected for "// &
                        trim(subsystem_solver_options(i_solver))//" subsystem solver."
                    test_solver_sanity_check = .false.
                end if
                if (len_trim(context%log_message) == 0) then
                    write(stderr, *) "test_solver_sanity_check failed: No error "// &
                        "message printed for non-symmetric Hessian with "// &
                        trim(subsystem_solver_options(i_solver))//" subsystem solver."
                    test_solver_sanity_check = .false.
                end if
            end if
        end do

        ! check that a warning is printed exactly once when a projection is passed,
        ! also when the nested stability check settings carry it together with the
        ! logging function, as they do once they inherit both, so that their check
        ! prints the warning as well
        do i_case = 1, size(projection_case_names)
            call setup_settings(settings, context)
            if (i_case == 2) then
                call setup_settings(settings%stability_settings, context)
                settings%stability_settings%project => mock_project_out_pair
            end if
            settings%project => mock_project_out_pair
            call solver_sanity_check(settings, 4_ip, error)
            if (error /= 0) then
                write(stderr, *) "test_solver_sanity_check failed: Error returned "// &
                    "when a projection is passed with "// &
                    trim(projection_case_names(i_case))//"."
                test_solver_sanity_check = .false.
            end if
            i_warning = index(context%log_message, " "//project_warning_msg)
            if (i_warning == 0) then
                write(stderr, *) "test_solver_sanity_check failed: Warning not "// &
                    "printed when a projection is passed with "// &
                    trim(projection_case_names(i_case))//"."
                test_solver_sanity_check = .false.
            end if
            if (i_warning > 0 .and. index(context%log_message(i_warning + 1:), &
                                          " "//project_warning_msg) /= 0) then
                write(stderr, *) "test_solver_sanity_check failed: Warning printed "// &
                    "twice when a projection is passed with "// &
                    trim(projection_case_names(i_case))//"."
                test_solver_sanity_check = .false.
            end if
        end do

    end function test_solver_sanity_check

    logical(c_bool) function test_stability_sanity_check() bind(C)
        !
        ! this function tests the subroutine which performs a sanity check for the
        ! stability check
        !
        use opentrustregion, only: &
            stability_settings_type, stability_sanity_check, project_warning_msg, &
            random_trial_vector_warning_msg, diag_solver_options, &
            standalone_n_random_trial_vectors

        ! the diagonalization solvers in mixed case and settings which the sanity check
        ! rejects
        character(len=*), parameter :: mixed_case_diag_solvers(2) = &
            [character(len=15) :: "Davidson", "Jacobi-Davidson"]
        character(len=*), parameter :: error_case_names(9) = &
            [character(len=55) :: "vanishing number of parameters", &
             "extra trial vector function without extra trial vectors", &
             "vanishing number of trial vectors", &
             "no random and a negative number of extra trial vectors", "no extra "// &
             "and no random trial vectors that fit", "unknown diagonalization solver", &
             "non-symmetric Hessian with Jacobi-Davidson solver", "stopping at an "// &
             "instability with non-symmetric Hessian", "positive stability threshold"]
        character(len=*), parameter :: accepted_case_names(4) = &
            [character(len=49) :: "non-symmetric Hessian with Davidson solver", &
             "stopping at an instability with symmetric Hessian", &
             "no random but extra trial vectors", "vanishing stability threshold"]
        type(stability_settings_type) :: settings
        integer(ip) :: error, i_solver, i_case, n_param
        type(test_context_type), target :: context

        ! assume tests pass
        test_stability_sanity_check = .true.

        ! check that every diagonalization solver passes when given in mixed case,
        ! which is converted to lowercase, and that an unset number of random trial
        ! vectors resolves to the default of a stability check outside the solver
        ! without any message
        do i_solver = 1, size(diag_solver_options)
            call setup_settings(settings, context)
            settings%diag_solver = mixed_case_diag_solvers(i_solver)
            call stability_sanity_check(settings, &
                                        4_ip * standalone_n_random_trial_vectors, error)
            if (error /= 0) then
                write(stderr, *) "test_stability_sanity_check failed: Error "// &
                    "returned for "//trim(diag_solver_options(i_solver))// &
                    " diagonalization solver."
                test_stability_sanity_check = .false.
            end if
            if (settings%diag_solver /= diag_solver_options(i_solver)) then
                write(stderr, *) "test_stability_sanity_check failed: Not "// &
                    "converted to lowercase for "// &
                    trim(diag_solver_options(i_solver))//" diagonalization solver."
                test_stability_sanity_check = .false.
            end if
            if (settings%n_random_trial_vectors /= standalone_n_random_trial_vectors) &
                then
                write(stderr, *) "test_stability_sanity_check failed: Unset number "// &
                    "of random trial vectors not resolved for "// &
                    trim(diag_solver_options(i_solver))//" diagonalization solver."
                test_stability_sanity_check = .false.
            end if
            if (len_trim(context%log_message) > 0) then
                write(stderr, *) "test_stability_sanity_check failed: Message "// &
                    "printed for "//trim(diag_solver_options(i_solver))// &
                    " diagonalization solver."
                test_stability_sanity_check = .false.
            end if
        end do

        ! check that the number of random trial vectors is reduced to half the number
        ! of parameters with a warning
        call setup_settings(settings, context)
        settings%n_random_trial_vectors = 3
        call stability_sanity_check(settings, 4_ip, error)
        if (error /= 0) then
            write(stderr, *) "test_stability_sanity_check failed: Error returned "// &
                "when number of random trial vectors is reduced."
            test_stability_sanity_check = .false.
        end if
        if (settings%n_random_trial_vectors /= 2) then
            write(stderr, *) "test_stability_sanity_check failed: Number of random "// &
                "trial vectors not reduced."
            test_stability_sanity_check = .false.
        end if
        if (index(context%log_message, " "//random_trial_vector_warning_msg) == 0) then
            write(stderr, *) "test_stability_sanity_check failed: Warning not "// &
                "printed when number of random trial vectors is reduced."
            test_stability_sanity_check = .false.
        end if

        ! check that every rejected setting returns an error and prints an error
        ! message
        do i_case = 1, size(error_case_names)
            call setup_settings(settings, context)
            call setup_error_logging(settings, context)
            n_param = 4
            select case (i_case)
            case (1)
                n_param = 0
            case (2)
                settings%get_extra_trial_vectors => mock_get_extra_trial_vectors
                settings%n_extra_trial_vectors = 0
            case (3)
                settings%n_extra_trial_vectors = 0
                settings%n_random_trial_vectors = 0
            case (4)
                settings%n_extra_trial_vectors = -1
                settings%n_random_trial_vectors = 0
            case (5)
                ! a single parameter leaves no room for the requested random trial
                ! vectors once their number is reduced to half the number of
                ! parameters
                n_param = 1
                settings%n_extra_trial_vectors = 0
            case (6)
                settings%diag_solver = "unknown"
            case (7)
                settings%diag_solver = "jacobi-davidson"
                settings%hess_symm = .false.
            case (8)
                settings%hess_symm = .false.
                settings%stop_on_instability = .true.
            case (9)
                settings%stability_thresh = 1e-3_rp
            end select
            call stability_sanity_check(settings, n_param, error)
            if (error /= 1) then
                write(stderr, *) "test_stability_sanity_check failed: No error "// &
                    "returned for "//trim(error_case_names(i_case))//"."
                test_stability_sanity_check = .false.
            end if
            if (len_trim(context%log_message) == 0) then
                write(stderr, *) "test_stability_sanity_check failed: No error "// &
                    "message printed for "//trim(error_case_names(i_case))//"."
                test_stability_sanity_check = .false.
            end if
        end do

        ! check that settings close to rejected ones are accepted
        do i_case = 1, size(accepted_case_names)
            call setup_settings(settings, context)
            select case (i_case)
            case (1)
                settings%hess_symm = .false.
            case (2)
                settings%stop_on_instability = .true.
            case (3)
                settings%n_random_trial_vectors = 0
            case (4)
                settings%stability_thresh = 0.0_rp
            end select
            call stability_sanity_check(settings, 4_ip, error)
            if (error /= 0) then
                write(stderr, *) "test_stability_sanity_check failed: Error "// &
                    "returned for "//trim(accepted_case_names(i_case))//"."
                test_stability_sanity_check = .false.
            end if
            if (i_case == 3 .and. settings%n_random_trial_vectors /= 0) then
                write(stderr, *) "test_stability_sanity_check failed: Vanishing "// &
                    "number of random trial vectors replaced."
                test_stability_sanity_check = .false.
            end if
        end do

        ! check that a warning is printed when a projection is passed
        call setup_settings(settings, context)
        settings%project => mock_project_out_pair
        call stability_sanity_check(settings, 4_ip, error)
        if (error /= 0) then
            write(stderr, *) "test_stability_sanity_check failed: Error returned "// &
                "when a projection is passed."
            test_stability_sanity_check = .false.
        end if
        if (index(context%log_message, " "//project_warning_msg) == 0) then
            write(stderr, *) "test_stability_sanity_check failed: Warning not "// &
                "printed when a projection is passed."
            test_stability_sanity_check = .false.
        end if

    end function test_stability_sanity_check

    logical(c_bool) function test_level_shifted_davidson() bind(C)
        !
        ! this function tests the level-shifted Davidson subroutine
        !
        use opentrustregion, only: solver_settings_type, subsystem_solver_options, &
                                   obj_func_type, hess_x_type, level_shifted_davidson, &
                                   error_hess_x, error_obj_func, &
                                   trust_radius_expand_ratio, trust_radius_expand_factor

        integer(ip), parameter :: n_stagnation = 20
        character(len=*), parameter :: &
            growth_case_names(2) = [character(len=8) :: "diagonal", "rotated"], &
            stagnation_solvers(2) = &
                [character(len=18) :: "davidson_ls", "jacobi-davidson_ls"]
        type(solver_settings_type) :: settings
        type(solver_state_type), target :: state
        real(rp) :: func, grad_norm, mu, solution_norm, noise, residual_tol
        real(rp), dimension(n_param) :: grad, h_diag, solution, newton_solution, eigvals
        real(rp) :: rotation(n_param, n_param), h_diag_stagnation(n_stagnation), &
                    solution_stagnation(n_stagnation), double_well_grad(2), &
                    double_well_h_diag(2), double_well_solution(2)
        integer(ip) :: i, i_solver, i_case, n_micro_ref, error
        procedure(obj_func_type), pointer :: obj_func_funptr
        procedure(hess_x_type), pointer :: hess_x_funptr
        logical :: local_region
        type(hartmann6d_context_type), target :: context
        type(hartmann6d_fault_context_type), target :: fault_context
        type(quadratic_context_type), target :: quadratic_context
        type(double_well_state_context_type), target :: state_context

        ! assume tests pass
        test_level_shifted_davidson = .true.

        ! solve the subproblem near a minimum and near a saddle point with every
        ! Davidson-family subsystem solver
        do i_solver = 1, size(subsystem_solver_options)
            if (index(subsystem_solver_options(i_solver), "davidson") == 0) cycle
            call check_region("near minimum")
            call check_region("near saddle point")
        end do

        ! run level-shifted Davidson near the minimum with a trust radius shorter than
        ! the Newton step, so that the step reaches the trust region boundary, where
        ! the quadratic model predicts the function change well, and check that the
        ! trust radius is expanded
        call setup_settings(settings, context)
        settings%grad_noise = .false.
        call setup_hartmann6d_point(context, near_minimum, func, grad, h_diag)
        grad_norm = norm2(grad)
        obj_func_funptr => obj_func
        hess_x_funptr => hess_x_fun
        state%trust_radius = 1e-3_rp
        call level_shifted_davidson(func, grad, grad_norm, h_diag, n_param, &
                                    obj_func_funptr, hess_x_funptr, settings, .false., &
                                    state, solution, solution_norm, mu, error)
        if (error /= 0) then
            write(stderr, *) "test_level_shifted_davidson failed: Produced error "// &
                "near minimum with trust radius shorter than Newton step."
            test_level_shifted_davidson = .false.
        end if
        if (abs(norm2(solution) - 1e-3_rp) > tol) then
            write(stderr, *) "test_level_shifted_davidson failed: Solution does "// &
                "not lie at trust region boundary near minimum with trust radius "// &
                "shorter than Newton step."
            test_level_shifted_davidson = .false.
        end if
        if (ref_hartmann6d_ratio(context, solution) < trust_radius_expand_ratio) then
            write(stderr, *) "test_level_shifted_davidson failed: Step at trust "// &
                "region boundary near minimum not predicted well enough to expand "// &
                "the trust radius."
            test_level_shifted_davidson = .false.
        end if
        if (abs(state%trust_radius - trust_radius_expand_factor * 1e-3_rp) > tol) then
            write(stderr, *) "test_level_shifted_davidson failed: Trust radius not "// &
                "expanded for a well predicted step at the trust region boundary "// &
                "near minimum."
            test_level_shifted_davidson = .false.
        end if

        ! run level-shifted Jacobi-Davidson near a saddle point with a slightly
        ! asymmetric Hessian linear transformation and check whether the reported
        ! number of Hessian linear transformations agrees with the calls, including the
        ! ones recalculated when Hessian symmetry is violated
        call setup_settings(settings, context)
        settings%grad_noise = .false.
        settings%subsystem_solver = "jacobi-davidson_ls"
        settings%jacobi_davidson_start = 1
        call setup_hartmann6d_point(context, near_saddle_point, func, grad, h_diag)
        grad_norm = norm2(grad)
        obj_func_funptr => obj_func
        hess_x_funptr => hess_x_fun_asymmetric
        state%trust_radius = 0.4_rp
        state%n_hess_x = 0
        context%n_hess_x_calls = 0
        call level_shifted_davidson(func, grad, grad_norm, h_diag, n_param, &
                                    obj_func_funptr, hess_x_funptr, settings, .false., &
                                    state, solution, solution_norm, mu, error)
        if (error /= 0) then
            write(stderr, *) "test_level_shifted_davidson failed: Produced error "// &
                "with asymmetric Hessian linear transformation."
            test_level_shifted_davidson = .false.
        end if
        if (.not. check_call_counts(context, "level_shifted_davidson", &
                                    "with asymmetric Hessian linear transformation", &
                                    n_hess_x=state%n_hess_x)) &
            test_level_shifted_davidson = .false.

        ! run level-shifted Davidson, which accepts a Hessian declared not symmetric,
        ! on a quadratic model whose Hessian is lower triangular, so that it is not
        ! symmetric but has its diagonal elements as real eigenvalues, the lowest of
        ! them negative with the first unit vector as left eigenvector, along which the
        ! gradient has a component small enough that the level shift on the trust
        ! region boundary lies close to this eigenvalue, and check that the solution is
        ! a level-shifted Newton step of this Hessian to the reduction factor of the
        ! global region on the boundary of the trust region it started from with a
        ! level shift below its lowest eigenvalue, which the reduced space only reaches
        ! with the reduced Hessian from both sides of the Hessian and diagonalized as a
        ! non-symmetric matrix
        call setup_settings(settings, quadratic_context)
        settings%grad_noise = .false.
        settings%hess_symm = .false.
        settings%subsystem_solver = "davidson_ls"
        quadratic_context%hess = generate_random_symm_matrix(n_param)
        do i = 1, n_param
            quadratic_context%hess(:i - 1, i) = 0.0_rp
            quadratic_context%hess(i, i) = real(i, kind=rp)
        end do
        quadratic_context%hess(1, 1) = -1.0_rp
        call random_number(grad)
        grad = 0.1_rp * (grad + 0.5_rp)
        grad(1) = 0.01_rp
        quadratic_context%grad = grad
        h_diag = [(quadratic_context%hess(i, i), i=1, n_param)]
        obj_func_funptr => quadratic_obj_func
        hess_x_funptr => quadratic_hess_x
        state%trust_radius = 0.5_rp
        call level_shifted_davidson(0.0_rp, grad, norm2(grad), h_diag, n_param, &
                                    obj_func_funptr, hess_x_funptr, settings, .false., &
                                    state, solution, solution_norm, mu, error)
        if (error /= 0) then
            write(stderr, *) "test_level_shifted_davidson failed: Produced error "// &
                "for non-symmetric Hessian."
            test_level_shifted_davidson = .false.
        end if
        if (abs(solution_norm - 0.5_rp) > tol) then
            write(stderr, *) "test_level_shifted_davidson failed: Solution not on "// &
                "the trust region boundary for non-symmetric Hessian."
            test_level_shifted_davidson = .false.
        end if
        if (mu >= minval(h_diag)) then
            write(stderr, *) "test_level_shifted_davidson failed: Level shift not "// &
                "below the lowest eigenvalue for non-symmetric Hessian."
            test_level_shifted_davidson = .false.
        end if
        if (norm2(grad + matmul(quadratic_context%hess, solution) - mu * solution) > &
            settings%global_red_factor * norm2(grad)) then
            write(stderr, *) "test_level_shifted_davidson failed: Solution does "// &
                "not describe level-shifted Newton step of non-symmetric Hessian."
            test_level_shifted_davidson = .false.
        end if

        ! run level-shifted Davidson with a Hessian linear transformation and with an
        ! objective function which fail at their first call and check that their errors
        ! are reported with their origins and that a failing Hessian linear
        ! transformation or objective function is counted
        call setup_settings(settings, fault_context)
        settings%grad_noise = .false.
        call setup_hartmann6d_point(fault_context, near_saddle_point, func, grad, &
                                    h_diag)
        obj_func_funptr => faulty_obj_func
        hess_x_funptr => faulty_hess_x
        do i_case = 1, 2
            fault_context%fault_callback = &
                merge(fault_hess_x, fault_obj_func, i_case == 1)
            fault_context%fault_call = 1
            fault_context%n_callback_calls = 0
            fault_context%n_hess_x_calls = 0
            fault_context%n_obj_func_calls = 0
            state%trust_radius = 0.4_rp
            state%n_hess_x = 0
            state%n_obj_func = 0
            call level_shifted_davidson(func, grad, grad_norm, h_diag, n_param, &
                                        obj_func_funptr, hess_x_funptr, settings, &
                                        .false., state, solution, solution_norm, mu, &
                                        error)
            if (error /= merge(error_hess_x, error_obj_func, i_case == 1) + 1) then
                write(stderr, *) "test_level_shifted_davidson failed: Error of "// &
                    "failing "// &
                    trim(fault_callback_names(fault_context%fault_callback))// &
                    " not reported with its origin."
                test_level_shifted_davidson = .false.
            end if
            if (.not. check_call_counts( &
                fault_context, "level_shifted_davidson", "for failing "// &
                trim(fault_callback_names(fault_context%fault_callback)), &
                n_hess_x=state%n_hess_x, n_obj_func=state%n_obj_func)) &
                test_level_shifted_davidson = .false.
        end do

        ! run level-shifted Davidson on a quadratic model with a gradient too small for
        ! its step to change the objective function, which does not change at all, so
        ! that the first step reaches maximum precision, and check that the subproblem
        ! stops without error and reports this
        call setup_settings(settings, quadratic_context)
        settings%grad_noise = .false.
        quadratic_context%hess = identity_matrix(n_param)
        quadratic_context%grad = spread(1e-9_rp, 1, n_param)
        h_diag = 1.0_rp
        obj_func_funptr => constant_obj_func
        hess_x_funptr => quadratic_hess_x
        state%trust_radius = 0.4_rp
        state%max_precision_reached = .false.
        call level_shifted_davidson(1.0_rp, quadratic_context%grad, &
                                    norm2(quadratic_context%grad), h_diag, n_param, &
                                    obj_func_funptr, hess_x_funptr, settings, .false., &
                                    state, solution, solution_norm, mu, error)
        if (error /= 0) then
            write(stderr, *) "test_level_shifted_davidson failed: Produced error "// &
                "when the objective function does not change."
            test_level_shifted_davidson = .false.
        end if
        if (.not. state%max_precision_reached) then
            write(stderr, *) "test_level_shifted_davidson failed: Maximum "// &
                "precision not reached when the objective function does not change."
            test_level_shifted_davidson = .false.
        end if

        ! run level-shifted Davidson with an objective function which exceeds the
        ! current one, as when the objective function and the orbital update disagree,
        ! so that every step is rejected and the trust radius shrinks until the
        ! predicted change can not be resolved and maximum precision is reached, and
        ! check that it stops without error and reports this and that the micro
        ! iterations of all rejected steps are counted, every step evaluates the
        ! objective function after at least one micro iteration
        call setup_settings(settings, context)
        settings%grad_noise = .false.
        call setup_hartmann6d_point(context, near_saddle_point, func, grad, h_diag)
        obj_func_funptr => constant_obj_func
        hess_x_funptr => hess_x_fun
        state%trust_radius = 0.4_rp
        state%max_precision_reached = .false.
        state%n_obj_func = 0
        call level_shifted_davidson(0.5_rp, grad, norm2(grad), h_diag, n_param, &
                                    obj_func_funptr, hess_x_funptr, settings, .false., &
                                    state, solution, solution_norm, mu, error)
        if (error /= 0) then
            write(stderr, *) "test_level_shifted_davidson failed: Produced error "// &
                "when every step is rejected."
            test_level_shifted_davidson = .false.
        end if
        if (.not. state%max_precision_reached) then
            write(stderr, *) "test_level_shifted_davidson failed: Maximum "// &
                "precision not reached when every step is rejected."
            test_level_shifted_davidson = .false.
        end if
        if (state%micro_iteration < state%n_obj_func) then
            write(stderr, *) "test_level_shifted_davidson failed: Micro iterations "// &
                "of rejected steps not counted when every step is rejected."
            test_level_shifted_davidson = .false.
        end if

        ! let the residual stagnate with a quadratic model whose gradient only couples
        ! to the second unit vector and a preconditioner which returns the following
        ! unit vectors, which can never reduce the residual, and check that the
        ! Davidson method stops without sufficient residual reduction, so that its step
        ! is rejected, while the Jacobi-Davidson method switches to the correction
        ! equations after ten micro iterations without sufficient residual reduction,
        ! long before the reduced space reaches the dimension of the full parameter
        ! space or the Jacobi-Davidson method would otherwise be started
        do i_case = 1, 2
            call setup_settings(settings, quadratic_context)
            settings%grad_noise = .false.
            settings%n_random_trial_vectors = 0
            settings%precond => mock_precond_next_unit_vector
            settings%subsystem_solver = stagnation_solvers(i_case)
            quadratic_context%hess = identity_matrix(n_stagnation)
            quadratic_context%hess(1:2, 1:2) = reshape([2.0_rp, 1.0_rp, &
                                                        1.0_rp, 3.0_rp], [2, 2])
            quadratic_context%grad = spread(0.0_rp, 1, n_stagnation)
            quadratic_context%grad(1) = 1.0_rp
            quadratic_context%next_unit_vector = 2
            h_diag_stagnation = [(quadratic_context%hess(i, i), i=1, n_stagnation)]
            obj_func_funptr => quadratic_obj_func
            hess_x_funptr => quadratic_hess_x
            state%trust_radius = 1.0_rp
            state%n_obj_func = 0
            call level_shifted_davidson( &
                0.0_rp, quadratic_context%grad, 1.0_rp, h_diag_stagnation, &
                n_stagnation, obj_func_funptr, hess_x_funptr, settings, .false., &
                state, solution_stagnation, solution_norm, mu, error)
            if (error /= 0) then
                write(stderr, *) "test_level_shifted_davidson failed: Produced "// &
                    "error when residual stagnates with "// &
                    trim(stagnation_solvers(i_case))//" solver."
                test_level_shifted_davidson = .false.
            end if
            if (i_case == 1 .and. state%n_obj_func < 2) then
                write(stderr, *) "test_level_shifted_davidson failed: Davidson "// &
                    "method did not stop and have its step rejected when residual "// &
                    "stagnates."
                test_level_shifted_davidson = .false.
            end if
            if (i_case == 2 .and. &
                state%micro_iteration - state%jacobi_davidson_iteration /= 10) then
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
        ! floor, for a diagonal Hessian with only two distinct eigenvalues and without
        ! random trial vectors the preconditioned residuals stay in a small subspace so
        ! that a new trial vector becomes linearly dependent before the reduced space
        ! reaches the dimension of the full parameter space, while for the same Hessian
        ! in a random orthonormal basis the reduced space grows to this dimension, each
        ! trial vector costs one Hessian linear transformation
        do i_case = 1, 2
            call setup_settings(settings, quadratic_context)
            settings%grad_noise = .false.
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
            residual_tol = 1e2_rp * epsilon(1.0_rp) * 1e8_rp
            h_diag = [(quadratic_context%hess(i, i), i=1, n_param)]
            obj_func_funptr => quadratic_obj_func
            hess_x_funptr => quadratic_hess_x
            state%trust_radius = 1.0_rp
            state%n_hess_x = 0
            call level_shifted_davidson( &
                0.0_rp, quadratic_context%grad, norm2(quadratic_context%grad), h_diag, &
                n_param, obj_func_funptr, hess_x_funptr, settings, .false., state, &
                solution, solution_norm, mu, error)
            if (error /= 0) then
                write(stderr, *) "test_level_shifted_davidson failed: Produced "// &
                    "error when reduced space stops growing for the "// &
                    trim(growth_case_names(i_case))//" Hessian."
                test_level_shifted_davidson = .false.
            end if
            if (i_case == 1 .and. state%n_hess_x >= n_param) then
                write(stderr, *) "test_level_shifted_davidson failed: Reduced "// &
                    "space does not stop at a linearly dependent trial vector for "// &
                    "the diagonal Hessian."
                test_level_shifted_davidson = .false.
            end if
            if (i_case == 2 .and. state%n_hess_x /= n_param) then
                write(stderr, *) "test_level_shifted_davidson failed: Reduced "// &
                    "space does not grow to the dimension of the full parameter "// &
                    "space for the rotated Hessian."
                test_level_shifted_davidson = .false.
            end if
            if (abs(norm2(solution) - 1.0_rp) > residual_tol) then
                write(stderr, *) "test_level_shifted_davidson failed: Solution "// &
                    "does not lie at trust region boundary when reduced space "// &
                    "stops growing for the "//trim(growth_case_names(i_case))// &
                    " Hessian."
                test_level_shifted_davidson = .false.
            end if
            if (norm2(quadratic_context%grad + matmul( &
                quadratic_context%hess, solution) - mu * solution) > residual_tol) then
                write(stderr, *) "test_level_shifted_davidson failed: Solution "// &
                    "does not describe level-shifted Newton step when reduced "// &
                    "space stops growing for the "//trim(growth_case_names(i_case))// &
                    " Hessian."
                test_level_shifted_davidson = .false.
            end if
        end do

        ! run level-shifted Davidson with a perturbed gradient in the range of a
        ! projection on a diagonal Hessian that the projection preserves and the
        ! diagonal preconditioner inverts exactly, after a subproblem solved in either
        ! region, and check that the step deviates from the Newton step by roughly the
        ! size of the perturbation, which is the reduction factor of that region, stays
        ! in the range of the projection and needs no more micro iterations than
        ! without the perturbation
        do i_case = 1, 2
            local_region = i_case == 2
            call setup_settings(settings, quadratic_context)
            settings%grad_noise = .false.
            settings%project => mock_project_out_pair
            quadratic_context%hess = identity_matrix(n_param)
            do i = 1, n_param
                quadratic_context%hess(i, i) = real(max(i, 2), kind=rp)
            end do
            call random_number(grad)
            grad = grad + 0.5_rp
            grad(1:2) = grad(1:2) - 0.5_rp * sum(grad(1:2))
            quadratic_context%grad = grad
            grad_norm = norm2(grad)
            h_diag = [(quadratic_context%hess(i, i), i=1, n_param)]
            newton_solution = -grad / h_diag
            obj_func_funptr => quadratic_obj_func
            hess_x_funptr => quadratic_hess_x
            state%trust_radius = 10.0_rp
            call level_shifted_davidson(0.0_rp, grad, grad_norm, h_diag, n_param, &
                                        obj_func_funptr, hess_x_funptr, settings, &
                                        local_region, state, solution, solution_norm, &
                                        mu, error)
            n_micro_ref = state%micro_iteration
            settings%grad_noise = .true.
            noise = merge(settings%local_red_factor, settings%global_red_factor, &
                          local_region)
            state%trust_radius = 10.0_rp
            call level_shifted_davidson(0.0_rp, grad, grad_norm, h_diag, n_param, &
                                        obj_func_funptr, hess_x_funptr, settings, &
                                        local_region, state, solution, solution_norm, &
                                        mu, error)
            if (error /= 0) then
                write(stderr, *) "test_level_shifted_davidson failed: Produced "// &
                    "error for perturbed gradient in "// &
                    trim(merge("local ", "global", local_region))//" region."
                test_level_shifted_davidson = .false.
            end if
            if (norm2(solution - newton_solution) < &
                0.1_rp * noise * grad_norm / maxval(h_diag)) then
                write(stderr, *) "test_level_shifted_davidson failed: Solution "// &
                    "does not solve a perturbed system for perturbed gradient in "// &
                    trim(merge("local ", "global", local_region))//" region."
                test_level_shifted_davidson = .false.
            end if
            if (norm2(solution - newton_solution) > &
                10.0_rp * noise * grad_norm / minval(h_diag)) then
                write(stderr, *) "test_level_shifted_davidson failed: Solution "// &
                    "deviates from the Newton step by more than the perturbation "// &
                    "for perturbed gradient in "// &
                    trim(merge("local ", "global", local_region))//" region."
                test_level_shifted_davidson = .false.
            end if
            if (abs(solution(1) + solution(2)) > tol) then
                write(stderr, *) "test_level_shifted_davidson failed: Solution "// &
                    "leaves the range of the projection for perturbed gradient in "// &
                    trim(merge("local ", "global", local_region))//" region."
                test_level_shifted_davidson = .false.
            end if
            if (state%micro_iteration > n_micro_ref) then
                write(stderr, *) "test_level_shifted_davidson failed: Perturbed "// &
                    "gradient needs more micro iterations than the unperturbed one "// &
                    "in "//trim(merge("local ", "global", local_region))//" region."
                test_level_shifted_davidson = .false.
            end if
        end do

        ! solve the subproblem of a quadratic model with a random coupling, which the
        ! diagonal preconditioner does not invert, so that the residual decreases
        ! gradually, after a subproblem solved in the global region, with a local
        ! reduction factor far below the global one, and check that the perturbed
        ! gradient, whose perturbation is the global reduction factor, stops the micro
        ! iterations as soon as the residual falls below the perturbation, earlier
        ! than the unperturbed gradient, which is solved to the local reduction factor
        call setup_settings(settings, quadratic_context)
        settings%local_red_factor = 1e-12_rp
        quadratic_context%hess = 0.1_rp * generate_random_symm_matrix(n_stagnation)
        do i = 1, n_stagnation
            quadratic_context%hess(i, i) = quadratic_context%hess(i, i) + &
                                           real(i + 4, kind=rp)
        end do
        quadratic_context%grad = spread(0.0_rp, 1, n_stagnation)
        call random_number(quadratic_context%grad)
        quadratic_context%grad = quadratic_context%grad + 0.5_rp
        h_diag_stagnation = [(quadratic_context%hess(i, i), i=1, n_stagnation)]
        obj_func_funptr => quadratic_obj_func
        hess_x_funptr => quadratic_hess_x
        do i_case = 1, 2
            settings%grad_noise = i_case == 2
            state%trust_radius = 10.0_rp
            call level_shifted_davidson( &
                0.0_rp, quadratic_context%grad, norm2(quadratic_context%grad), &
                h_diag_stagnation, n_stagnation, obj_func_funptr, hess_x_funptr, &
                settings, .false., state, solution_stagnation, solution_norm, mu, error)
            if (error /= 0) then
                write(stderr, *) "test_level_shifted_davidson failed: Produced "// &
                    "error for "// &
                    trim(merge("unperturbed", "perturbed  ", i_case == 1))// &
                    " gradient on a coupled quadratic model."
                test_level_shifted_davidson = .false.
            end if
            if (i_case == 1) n_micro_ref = state%micro_iteration
        end do
        if (state%micro_iteration >= n_micro_ref) then
            write(stderr, *) "test_level_shifted_davidson failed: Micro iterations "// &
                "for perturbed gradient not stopped at the perturbation."
            test_level_shifted_davidson = .false.
        end if

        ! solve the subproblem of the double well function with a state holding the
        ! values of an earlier subproblem and a Hessian linear transformation which
        ! reads the state through the pointer in its context, and check that the
        ! initial linear transformations already see no micro iterations and the
        ! residual of the vanishing solution, which is the gradient
        call setup_settings(settings, state_context)
        settings%grad_noise = .false.
        call setup_double_well_state_context(state_context, state, func, &
                                             double_well_grad, double_well_h_diag)
        obj_func_funptr => double_well_obj_func
        hess_x_funptr => state_reading_hess_x
        call level_shifted_davidson(func, double_well_grad, norm2(double_well_grad), &
                                    double_well_h_diag, 2_ip, obj_func_funptr, &
                                    hess_x_funptr, settings, .false., state, &
                                    double_well_solution, solution_norm, mu, error)
        if (error /= 0) then
            write(stderr, *) "test_level_shifted_davidson failed: Produced error "// &
                "with a state of an earlier subproblem."
            test_level_shifted_davidson = .false.
        end if
        if (.not. check_reset_state_seen(state_context, "level_shifted_davidson", &
                                         "Hessian linear transformation")) &
            test_level_shifted_davidson = .false.

    contains

        subroutine check_region(region)
            !
            ! this subroutine solves the subproblem in a region with the current
            ! subsystem solver, the Jacobi-Davidson ones switching after the first
            ! micro iteration of every step, and checks the solution, the level shift,
            ! the trust region, the state and the counters
            !
            character(len=*), intent(in) :: region

            character(len=60) :: case_name
            real(rp) :: ratio
            logical :: jacobi_davidson

            case_name = region//" with "//trim(subsystem_solver_options(i_solver))// &
                        " solver"
            jacobi_davidson = index(subsystem_solver_options(i_solver), "jacobi") > 0

            ! set up the Hartmann 6D function at the point
            call setup_settings(settings, context)
            call setup_hartmann6d_point(context, merge( &
                near_minimum, near_saddle_point, region == "near minimum"), func, &
                                        grad, h_diag)
            grad_norm = norm2(grad)
            obj_func_funptr => obj_func
            hess_x_funptr => hess_x_fun
            settings%grad_noise = .false.
            settings%subsystem_solver = subsystem_solver_options(i_solver)
            settings%jacobi_davidson_start = 1
            state%trust_radius = 0.4_rp
            state%n_hess_x = 0
            state%n_obj_func = 0
            context%n_hess_x_calls = 0
            context%n_obj_func_calls = 0

            ! run level-shifted Davidson
            call level_shifted_davidson(func, grad, grad_norm, h_diag, n_param, &
                                        obj_func_funptr, hess_x_funptr, settings, &
                                        .false., state, solution, solution_norm, mu, &
                                        error)
            if (error /= 0) then
                write(stderr, *) "test_level_shifted_davidson failed: Produced "// &
                    "error "//trim(case_name)//"."
                test_level_shifted_davidson = .false.
            end if

            ! the micro iterations of every step after its first one are
            ! Jacobi-Davidson iterations for the Jacobi-Davidson solvers
            if (state%jacobi_davidson_iteration /= merge( &
                state%micro_iteration - state%n_obj_func, 0_ip, jacobi_davidson)) then
                write(stderr, *) "test_level_shifted_davidson failed: "// &
                    "Jacobi-Davidson iterations wrong "//trim(case_name)//"."
                test_level_shifted_davidson = .false.
            end if

            select case (region)
            case ("near minimum")
                if (abs(mu) > tol) then
                    write(stderr, *) "test_level_shifted_davidson failed: Level "// &
                        "shift is not zero "//trim(case_name)//"."
                    test_level_shifted_davidson = .false.
                end if
                if (norm2(grad + matmul(context%hess, solution)) > &
                    settings%local_red_factor * grad_norm) then
                    write(stderr, *) "test_level_shifted_davidson failed: Solution "// &
                        "does not describe Newton step "//trim(case_name)//"."
                    test_level_shifted_davidson = .false.
                end if
            case ("near saddle point")
                if (mu >= 0.0_rp) then
                    write(stderr, *) "test_level_shifted_davidson failed: Level "// &
                        "shift is not negative "//trim(case_name)//"."
                    test_level_shifted_davidson = .false.
                end if
                if (norm2(grad + matmul(context%hess, solution) - mu * solution) > &
                    settings%global_red_factor * grad_norm) then
                    write(stderr, *) "test_level_shifted_davidson failed: Solution "// &
                        "does not describe level-shifted Newton step "// &
                        trim(case_name)//"."
                    test_level_shifted_davidson = .false.
                end if
            end select

            ! the relative residual is that of the returned level-shifted Newton step
            if (abs(state%relative_residual - norm2(grad + matmul( &
                context%hess, solution) - mu * solution) / grad_norm) > tol) then
                write(stderr, *) "test_level_shifted_davidson failed: Relative "// &
                    "residual wrong "//trim(case_name)//"."
                test_level_shifted_davidson = .false.
            end if

            ! the step stays within the trust region near the minimum and lies at its
            ! boundary near the saddle point
            ratio = ref_hartmann6d_ratio(context, solution)
            if (region == "near minimum" .and. &
                norm2(solution) > ref_step_trust_radius(state%trust_radius, ratio, &
                                                        norm2(solution)) + tol) then
                write(stderr, *) "test_level_shifted_davidson failed: Solution "// &
                    "does not stay within trust region "//trim(case_name)//"."
                test_level_shifted_davidson = .false.
            end if
            if (region == "near saddle point" .and. abs( &
                norm2(solution) - ref_step_trust_radius(state%trust_radius, ratio, &
                                                        norm2(solution))) > tol) then
                write(stderr, *) "test_level_shifted_davidson failed: Solution "// &
                    "does not lie at trust region boundary "//trim(case_name)//"."
                test_level_shifted_davidson = .false.
            end if

            ! the counters agree with the calls
            if (.not. check_call_counts(context, "level_shifted_davidson", &
                                        trim(case_name), n_hess_x=state%n_hess_x, &
                                        n_obj_func=state%n_obj_func)) &
                test_level_shifted_davidson = .false.

        end subroutine check_region

    end function test_level_shifted_davidson

    logical(c_bool) function test_truncated_conjugate_gradient() bind(C)
        !
        ! this function tests the truncated conjugate gradient subroutine
        !
        use opentrustregion, only: obj_func_type, hess_x_type, solver_settings_type, &
                                   truncated_conjugate_gradient, error_hess_x, &
                                   error_obj_func, function_unchanged_warning_msg, &
                                   precond_floor

        real(rp) :: func, ratio, solution_norm, metric_norm
        real(rp), dimension(n_param) :: grad, h_diag, solution
        integer(ip) :: i, i_case, error
        procedure(obj_func_type), pointer :: obj_func_funptr
        procedure(hess_x_type), pointer :: hess_x_funptr
        type(solver_settings_type) :: settings
        type(solver_state_type), target :: state
        type(hartmann6d_context_type), target :: context
        type(hartmann6d_fault_context_type), target :: fault_context
        type(quadratic_context_type), target :: quadratic_context

        ! assume tests pass
        test_truncated_conjugate_gradient = .true.

        ! setup settings object
        call setup_settings(settings, context)

        ! initialize variables
        state%trust_radius = 0.4_rp
        obj_func_funptr => obj_func
        hess_x_funptr => hess_x_fun

        ! start in quadratic region near minimum
        context%vars = near_minimum
        func = hartmann6d_func(context%vars)
        call hartmann6d_gradient(context%vars, grad)
        context%hess = hartmann6d_hessian(context%vars)
        h_diag = [(context%hess(i, i), i=1, size(h_diag))]

        ! run truncated conjugate gradient with micro iterations left in the state by
        ! a previous subproblem, check whether the solution stays within the trust
        ! region and reduces the function value, whether its norm in the metric of the
        ! absolute diagonal preconditioner is returned, whether the micro iterations,
        ! one per Hessian linear transformation, and the converged residual relative to
        ! the initial one are in the state and whether the reported numbers of Hessian
        ! linear transformations and objective function evaluations agree with the
        ! calls
        state%micro_iteration = 100
        state%jacobi_davidson_iteration = 100
        state%n_hess_x = 0
        state%n_obj_func = 0
        context%n_hess_x_calls = 0
        context%n_obj_func_calls = 0
        call truncated_conjugate_gradient(func, grad, h_diag, n_param, &
                                          obj_func_funptr, hess_x_funptr, settings, &
                                          state, solution, solution_norm, error)
        if (error /= 0) then
            write(stderr, *) "test_truncated_conjugate_gradient failed: Produced "// &
                "error near minimum."
            test_truncated_conjugate_gradient = .false.
        end if
        if (.not. check_call_counts(context, "truncated_conjugate_gradient", &
                                    "near minimum", n_hess_x=state%n_hess_x, &
                                    n_obj_func=state%n_obj_func)) &
            test_truncated_conjugate_gradient = .false.
        if (state%micro_iteration < 1 .or. &
            state%micro_iteration /= context%n_hess_x_calls .or. &
            state%jacobi_davidson_iteration /= 0 .or. &
            state%relative_residual <= 0.0_rp .or. state%relative_residual > 1e-3_rp) &
            then
            write(stderr, *) "test_truncated_conjugate_gradient failed: Micro "// &
                "iterations or converged relative residual not in the state near "// &
                "minimum."
            test_truncated_conjugate_gradient = .false.
        end if
        if (hartmann6d_func(context%vars + solution) >= func) then
            write(stderr, *) "test_truncated_conjugate_gradient failed: Solution "// &
                "does not reduce function value near minimum."
            test_truncated_conjugate_gradient = .false.
        end if
        ratio = (hartmann6d_func(context%vars + solution) - func) / &
                dot_product(solution, grad + 0.5_rp * matmul(context%hess, solution))
        metric_norm = sqrt(dot_product(solution, &
                                       solution / max(abs(h_diag), precond_floor)))
        if (metric_norm > &
            ref_step_trust_radius(state%trust_radius, ratio, metric_norm) + tol) then
            write(stderr, *) "test_truncated_conjugate_gradient failed: Solution "// &
                "does not stay within trust region near minimum."
            test_truncated_conjugate_gradient = .false.
        end if
        if (abs(solution_norm - metric_norm) > tol) then
            write(stderr, *) "test_truncated_conjugate_gradient failed: Norm of "// &
                "the solution in the metric of the absolute diagonal "// &
                "preconditioner not returned near minimum."
            test_truncated_conjugate_gradient = .false.
        end if

        ! run truncated conjugate gradient again with an identity preconditioner, which
        ! makes the metric of the trust region Euclidean, and check that the norm of
        ! the solution is measured with it
        settings%precond => identity_precond
        state%trust_radius = 0.4_rp
        call truncated_conjugate_gradient(func, grad, h_diag, n_param, &
                                          obj_func_funptr, hess_x_funptr, settings, &
                                          state, solution, solution_norm, error)
        if (error /= 0) then
            write(stderr, *) "test_truncated_conjugate_gradient failed: Produced "// &
                "error near minimum with preconditioner."
            test_truncated_conjugate_gradient = .false.
        end if
        if (abs(solution_norm - norm2(solution)) > tol) then
            write(stderr, *) "test_truncated_conjugate_gradient failed: Norm of "// &
                "the solution not measured with the preconditioner near minimum."
            test_truncated_conjugate_gradient = .false.
        end if
        settings%precond => null()

        ! run truncated conjugate gradient with a trust radius twice as long as the
        ! first solution in the metric of the trust region but shorter than its
        ! Euclidean norm, so that the step stays inside the trust region, and check
        ! that the trust radius is not expanded for this step although the quadratic
        ! model is accurate
        state%trust_radius = 2.0_rp * metric_norm
        call truncated_conjugate_gradient(func, grad, h_diag, n_param, &
                                          obj_func_funptr, hess_x_funptr, settings, &
                                          state, solution, solution_norm, error)
        if (error /= 0) then
            write(stderr, *) "test_truncated_conjugate_gradient failed: Produced "// &
                "error for a step inside the trust region."
            test_truncated_conjugate_gradient = .false.
        end if
        if (abs(state%trust_radius - 2.0_rp * metric_norm) > tol) then
            write(stderr, *) "test_truncated_conjugate_gradient failed: Trust "// &
                "radius expanded for a step inside the trust region."
            test_truncated_conjugate_gradient = .false.
        end if

        ! run truncated conjugate gradient with a trust radius half as long as the
        ! first solution, so that the step reaches the trust region boundary although
        ! the curvature is positive, and check whether the solution lies at the trust
        ! region boundary and reduces the function value
        state%trust_radius = 0.5_rp * metric_norm
        call truncated_conjugate_gradient(func, grad, h_diag, n_param, &
                                          obj_func_funptr, hess_x_funptr, settings, &
                                          state, solution, solution_norm, error)
        if (error /= 0) then
            write(stderr, *) "test_truncated_conjugate_gradient failed: Produced "// &
                "error near minimum with trust radius shorter than step."
            test_truncated_conjugate_gradient = .false.
        end if
        if (hartmann6d_func(context%vars + solution) >= func) then
            write(stderr, *) "test_truncated_conjugate_gradient failed: Solution "// &
                "does not reduce function value near minimum with trust radius "// &
                "shorter than step."
            test_truncated_conjugate_gradient = .false.
        end if
        ratio = (hartmann6d_func(context%vars + solution) - func) / &
                dot_product(solution, grad + 0.5_rp * matmul(context%hess, solution))
        metric_norm = sqrt(dot_product(solution, &
                                       solution / max(abs(h_diag), precond_floor)))
        if (abs(metric_norm - ref_step_trust_radius(state%trust_radius, ratio, &
                                                    metric_norm)) > tol) then
            write(stderr, *) "test_truncated_conjugate_gradient failed: Solution "// &
                "does not lie at trust region boundary near minimum with trust "// &
                "radius shorter than step."
            test_truncated_conjugate_gradient = .false.
        end if

        ! start near saddle point
        context%vars = near_saddle_point
        func = hartmann6d_func(context%vars)
        call hartmann6d_gradient(context%vars, grad)
        context%hess = hartmann6d_hessian(context%vars)
        h_diag = [(context%hess(i, i), i=1, size(h_diag))]
        state%trust_radius = 0.4_rp

        ! run truncated conjugate gradient, check whether the solution lies at the
        ! trust region boundary and reduces the function value
        call truncated_conjugate_gradient(func, grad, h_diag, n_param, &
                                          obj_func_funptr, hess_x_funptr, settings, &
                                          state, solution, solution_norm, error)
        if (error /= 0) then
            write(stderr, *) "test_truncated_conjugate_gradient failed: Produced "// &
                "error near saddle point."
            test_truncated_conjugate_gradient = .false.
        end if
        if (hartmann6d_func(context%vars + solution) >= func) then
            write(stderr, *) "test_truncated_conjugate_gradient failed: Solution "// &
                "does not reduce function value near saddle point."
            test_truncated_conjugate_gradient = .false.
        end if
        ratio = (hartmann6d_func(context%vars + solution) - func) / &
                dot_product(solution, grad + 0.5_rp * matmul(context%hess, solution))
        metric_norm = sqrt(dot_product(solution, &
                                       solution / max(abs(h_diag), precond_floor)))
        if (abs(metric_norm - ref_step_trust_radius(state%trust_radius, ratio, &
                                                    metric_norm)) > tol) then
            write(stderr, *) "test_truncated_conjugate_gradient failed: Solution "// &
                "does not lie at trust region boundary near saddle point."
            test_truncated_conjugate_gradient = .false.
        end if

        ! run truncated conjugate gradient with a Hessian linear transformation and
        ! with an objective function which fail at their first call and check that
        ! their errors are reported with their origins and that a failing Hessian
        ! linear transformation or objective function is counted
        settings%context => fault_context
        fault_context%vars = context%vars
        fault_context%hess = context%hess
        obj_func_funptr => faulty_obj_func
        hess_x_funptr => faulty_hess_x
        do i_case = 1, 2
            fault_context%fault_callback = &
                merge(fault_hess_x, fault_obj_func, i_case == 1)
            fault_context%fault_call = 1
            fault_context%n_callback_calls = 0
            fault_context%n_hess_x_calls = 0
            fault_context%n_obj_func_calls = 0
            state%trust_radius = 0.4_rp
            state%n_hess_x = 0
            state%n_obj_func = 0
            call truncated_conjugate_gradient( &
                func, grad, h_diag, n_param, obj_func_funptr, hess_x_funptr, settings, &
                state, solution, solution_norm, error)
            if (error /= merge(error_hess_x, error_obj_func, i_case == 1) + 1) then
                write(stderr, *) "test_truncated_conjugate_gradient failed: Error "// &
                    "of failing "// &
                    trim(fault_callback_names(fault_context%fault_callback))// &
                    " not reported with its origin."
                test_truncated_conjugate_gradient = .false.
            end if
            if (.not. check_call_counts( &
                fault_context, "truncated_conjugate_gradient", "for failing "// &
                trim(fault_callback_names(fault_context%fault_callback)), &
                n_hess_x=state%n_hess_x, n_obj_func=state%n_obj_func)) &
                test_truncated_conjugate_gradient = .false.
        end do
        fault_context%fault_callback = 0

        ! run truncated conjugate gradient on a quadratic model with a diagonal
        ! Hessian, which the absolute diagonal preconditioner inverts exactly, so that
        ! the first micro iteration solves the system it is given, and a gradient along
        ! the first parameter, without and with the gradient perturbation, and check
        ! that the solution only leaves the first parameter if the gradient is
        ! perturbed, and then by the reduction factor of the local region, which the
        ! residual of the unperturbed Newton equation shows
        call setup_settings(settings, quadratic_context)
        quadratic_context%hess = identity_matrix(n_param)
        do i = 1, n_param
            quadratic_context%hess(i, i) = real(i, kind=rp)
        end do
        quadratic_context%grad = spread(0.0_rp, 1, n_param)
        quadratic_context%grad(1) = 0.1_rp
        h_diag = [(quadratic_context%hess(i, i), i=1, n_param)]
        obj_func_funptr => quadratic_obj_func
        hess_x_funptr => quadratic_hess_x
        do i_case = 1, 2
            settings%grad_noise = i_case == 2
            state%trust_radius = 10.0_rp
            call truncated_conjugate_gradient( &
                0.0_rp, quadratic_context%grad, h_diag, n_param, obj_func_funptr, &
                hess_x_funptr, settings, state, solution, solution_norm, error)
            if (error /= 0) then
                write(stderr, *) "test_truncated_conjugate_gradient failed: "// &
                    "Produced error "// &
                    trim(merge("without", "with   ", i_case == 1))// &
                    " gradient perturbation."
                test_truncated_conjugate_gradient = .false.
            end if
            if ((maxval(abs(solution(2:))) > 0.0_rp) .neqv. settings%grad_noise) then
                write(stderr, *) "test_truncated_conjugate_gradient failed: "// &
                    "Solution "//trim(merge("left  ", "stayed", i_case == 1))// &
                    " along the gradient "// &
                    trim(merge("without", "with   ", i_case == 1))// &
                    " gradient perturbation."
                test_truncated_conjugate_gradient = .false.
            end if
            if (settings%grad_noise .and. abs(norm2( &
                matmul(quadratic_context%hess, solution) + quadratic_context%grad) - &
                settings%local_red_factor * norm2(quadratic_context%grad)) > tol) then
                write(stderr, *) "test_truncated_conjugate_gradient failed: "// &
                    "Gradient not perturbed by the reduction factor of the local "// &
                    "region."
                test_truncated_conjugate_gradient = .false.
            end if
        end do

        ! run truncated conjugate gradient on the quadratic model with an identity
        ! preconditioner, so that the micro iterations converge gradually, and with a
        ! gradient perturbation of a tenth of the gradient norm, which a projection
        ! that removes every random vector leaves out, so that only the floor the
        ! perturbation sets on the convergence of the micro iterations remains, and
        ! check that they stop once the relative residual falls below the size of the
        ! perturbation instead of converging further
        settings%precond => identity_precond
        settings%project => vanishing_project
        settings%local_red_factor = 0.1_rp
        quadratic_context%grad = 1.0_rp
        state%trust_radius = 10.0_rp
        call truncated_conjugate_gradient( &
            0.0_rp, quadratic_context%grad, h_diag, n_param, obj_func_funptr, &
            hess_x_funptr, settings, state, solution, solution_norm, error)
        if (error /= 0) then
            write(stderr, *) "test_truncated_conjugate_gradient failed: Produced "// &
                "error with a large gradient perturbation."
            test_truncated_conjugate_gradient = .false.
        end if
        if (state%relative_residual > settings%local_red_factor) then
            write(stderr, *) "test_truncated_conjugate_gradient failed: Micro "// &
                "iterations stopped above the relative size of the gradient "// &
                "perturbation."
            test_truncated_conjugate_gradient = .false.
        end if
        if (state%relative_residual <= 1e-3_rp) then
            write(stderr, *) "test_truncated_conjugate_gradient failed: Micro "// &
                "iterations converged below the relative size of the gradient "// &
                "perturbation."
            test_truncated_conjugate_gradient = .false.
        end if

        ! run truncated conjugate gradient on a quadratic model with a gradient too
        ! small for its step to change the objective function, which does not change
        ! at all, so that the first step reaches maximum precision, and check that the
        ! subproblem stops without error, reports this with a warning and returns no
        ! norm for the step that is not taken
        call setup_settings(settings, quadratic_context)
        quadratic_context%hess = identity_matrix(n_param)
        quadratic_context%grad = spread(1e-9_rp, 1, n_param)
        h_diag = 1.0_rp
        obj_func_funptr => constant_obj_func
        hess_x_funptr => quadratic_hess_x
        state%trust_radius = 0.4_rp
        call truncated_conjugate_gradient( &
            1.0_rp, quadratic_context%grad, h_diag, n_param, obj_func_funptr, &
            hess_x_funptr, settings, state, solution, solution_norm, error)
        if (error /= 0) then
            write(stderr, *) "test_truncated_conjugate_gradient failed: Produced "// &
                "error when the objective function does not change."
            test_truncated_conjugate_gradient = .false.
        end if
        if (.not. state%max_precision_reached) then
            write(stderr, *) "test_truncated_conjugate_gradient failed: Maximum "// &
                "precision not reached when the objective function does not change."
            test_truncated_conjugate_gradient = .false.
        end if
        if (abs(solution_norm) > 0.0_rp) then
            write(stderr, *) "test_truncated_conjugate_gradient failed: Norm of "// &
                "the solution not vanishing for a step that is not taken when "// &
                "function value does not change."
            test_truncated_conjugate_gradient = .false.
        end if
        if (index(quadratic_context%log_message, &
                  " "//function_unchanged_warning_msg) == 0) then
            write(stderr, *) "test_truncated_conjugate_gradient failed: Warning "// &
                "not printed when function value does not change."
            test_truncated_conjugate_gradient = .false.
        end if

        ! run truncated conjugate gradient on the quadratic model with an objective
        ! function which exceeds the current one, as when the objective function and
        ! the orbital update disagree, so that every step is rejected and the trust
        ! radius shrinks until the predicted change can not be resolved, and check
        ! that the subproblem stops without error at maximum precision and returns no
        ! norm for the step that is not taken
        quadratic_context%grad = spread(0.1_rp, 1, n_param)
        state%trust_radius = 0.4_rp
        call truncated_conjugate_gradient( &
            0.5_rp, quadratic_context%grad, h_diag, n_param, obj_func_funptr, &
            hess_x_funptr, settings, state, solution, solution_norm, error)
        if (error /= 0) then
            write(stderr, *) "test_truncated_conjugate_gradient failed: Produced "// &
                "error when every step is rejected."
            test_truncated_conjugate_gradient = .false.
        end if
        if (.not. state%max_precision_reached) then
            write(stderr, *) "test_truncated_conjugate_gradient failed: Maximum "// &
                "precision not reached when every step is rejected."
            test_truncated_conjugate_gradient = .false.
        end if
        if (abs(solution_norm) > 0.0_rp) then
            write(stderr, *) "test_truncated_conjugate_gradient failed: Norm of "// &
                "the solution not vanishing for a step that is not taken when "// &
                "every step is rejected."
            test_truncated_conjugate_gradient = .false.
        end if
        hess_x_funptr => hess_x_fun

        ! start near the minimum with a trust radius far beyond the step, so that the
        ! step ends inside the trust region, and let the objective function reject this
        ! first step, so that the trust radius shrinks without a new evaluation while
        ! the step still lies within it before the step is cut back to the shrunk trust
        ! region boundary, evaluated again and accepted, and check that the objective
        ! function is evaluated and counted twice and that the returned step lies at
        ! the trust region boundary
        call setup_settings(settings, fault_context)
        fault_context%vars = near_minimum
        fault_context%hess = hartmann6d_hessian(near_minimum)
        func = hartmann6d_func(near_minimum)
        call hartmann6d_gradient(near_minimum, grad)
        h_diag = [(fault_context%hess(i, i), i=1, size(h_diag))]
        obj_func_funptr => first_step_rejecting_obj_func
        state%trust_radius = 10.0_rp
        state%n_obj_func = 0
        fault_context%n_callback_calls = 0
        call truncated_conjugate_gradient(func, grad, h_diag, n_param, &
                                          obj_func_funptr, hess_x_funptr, settings, &
                                          state, solution, solution_norm, error)
        if (error /= 0) then
            write(stderr, *) "test_truncated_conjugate_gradient failed: Produced "// &
                "error when step inside trust region is rejected."
            test_truncated_conjugate_gradient = .false.
        end if
        if (fault_context%n_callback_calls(fault_obj_func) /= 2) then
            write(stderr, *) "test_truncated_conjugate_gradient failed: Objective "// &
                "function not evaluated once more for cut back step when step "// &
                "inside trust region is rejected."
            test_truncated_conjugate_gradient = .false.
        end if
        if (state%n_obj_func /= 2) then
            write(stderr, *) "test_truncated_conjugate_gradient failed: Objective "// &
                "function evaluations not counted when step inside trust region is "// &
                "rejected."
            test_truncated_conjugate_gradient = .false.
        end if
        ratio = (hartmann6d_func(near_minimum + solution) - func) / dot_product( &
            solution, grad + 0.5_rp * matmul(fault_context%hess, solution))
        metric_norm = sqrt(dot_product(solution, &
                                       solution / max(abs(h_diag), precond_floor)))
        if (abs(solution_norm - metric_norm) > tol) then
            write(stderr, *) "test_truncated_conjugate_gradient failed: Norm of "// &
                "the cut back step not returned when step inside trust region is "// &
                "rejected."
            test_truncated_conjugate_gradient = .false.
        end if
        if (abs(metric_norm - ref_step_trust_radius(state%trust_radius, ratio, &
                                                    metric_norm)) > tol) then
            write(stderr, *) "test_truncated_conjugate_gradient failed: Cut back "// &
                "step does not lie at trust region boundary when step inside trust "// &
                "region is rejected."
            test_truncated_conjugate_gradient = .false.
        end if

        ! let the evaluation of the cut back step fail and check that its error is
        ! reported with its origin
        fault_context%fault_callback = fault_obj_func
        fault_context%fault_call = 2
        fault_context%n_callback_calls = 0
        state%trust_radius = 10.0_rp
        call truncated_conjugate_gradient(func, grad, h_diag, n_param, &
                                          obj_func_funptr, hess_x_funptr, settings, &
                                          state, solution, solution_norm, error)
        if (error /= error_obj_func + 1) then
            write(stderr, *) "test_truncated_conjugate_gradient failed: Did not "// &
                "report the error of the objective function failing for the cut "// &
                "back step with its origin."
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

    logical(c_bool) function test_string_in() bind(C)
        !
        ! this function tests the function that checks whether a string is contained in
        ! another string
        !
        use opentrustregion, only: string_in

        character(len=*), parameter :: input = "OpenTrustRegion123!"

        ! assume tests pass
        test_string_in = .true.

        ! test if substring is in string
        if (.not. string_in("Trust", input)) then
            write(stderr, *) "test_string_in failed: Substring not correctly "// &
                "detected in string."
            test_string_in = .false.
        end if
        if (string_in("trust", input)) then
            write(stderr, *) "test_string_in failed: Substring detected in string "// &
                "even though it is not present."
            test_string_in = .false.
        end if

    end function test_string_in

    logical(c_bool) function test_perturb_vector() bind(C)
        !
        ! this function tests the subroutine which perturbs a vector by a random vector
        ! of a given size relative to it
        !
        use opentrustregion, only: solver_settings_type, perturb_vector, error_project

        real(rp), parameter :: noise = 1e-3_rp
        type(solver_settings_type) :: settings
        real(rp) :: vector(n_param), perturbed_vector(n_param)
        integer(ip) :: error, i_case
        type(test_context_type), target :: context
        type(hartmann6d_fault_context_type), target :: fault_context

        ! assume tests pass
        test_perturb_vector = .true.

        ! setup settings object
        call setup_settings(settings, context)

        ! perturb a random vector in the range of a projection without and with this
        ! projection and check that the perturbation has the given size relative to
        ! the vector and, with the projection, stays in its range
        do i_case = 1, 2
            if (i_case == 2) settings%project => mock_project_out_pair
            call random_number(vector)
            vector(1:2) = vector(1:2) - 0.5_rp * sum(vector(1:2))
            perturbed_vector = vector
            call perturb_vector(perturbed_vector, noise, settings, error)
            if (error /= 0) then
                write(stderr, *) "test_perturb_vector failed: Produced error "// &
                    trim(merge("without", "with   ", i_case == 1))//" projection."
                test_perturb_vector = .false.
            end if
            if (abs(norm2(perturbed_vector - vector) - noise * norm2(vector)) > tol) &
                then
                write(stderr, *) "test_perturb_vector failed: Perturbation does "// &
                    "not have the given relative size "// &
                    trim(merge("without", "with   ", i_case == 1))//" projection."
                test_perturb_vector = .false.
            end if
        end do
        if (abs(perturbed_vector(1) + perturbed_vector(2)) > tol) then
            write(stderr, *) "test_perturb_vector failed: Perturbation leaves the "// &
                "range of the projection."
            test_perturb_vector = .false.
        end if

        ! perturb a vector with a projection which removes every random vector and
        ! check that the vector is left unperturbed
        settings%project => vanishing_project
        perturbed_vector = vector
        call perturb_vector(perturbed_vector, noise, settings, error)
        if (error /= 0) then
            write(stderr, *) "test_perturb_vector failed: Produced error when the "// &
                "projection removes every random vector."
            test_perturb_vector = .false.
        end if
        if (maxval(abs(perturbed_vector - vector)) > 0.0_rp) then
            write(stderr, *) "test_perturb_vector failed: Vector perturbed when "// &
                "the projection removes every random vector."
            test_perturb_vector = .false.
        end if

        ! check that the error of a failing projection is reported with its origin
        call setup_settings(settings, fault_context)
        settings%project => faulty_project
        fault_context%fault_callback = fault_project
        fault_context%fault_call = 1
        call perturb_vector(perturbed_vector, noise, settings, error)
        if (error /= error_project + 1) then
            write(stderr, *) "test_perturb_vector failed: Error of failing "// &
                "projection not reported with its origin."
            test_perturb_vector = .false.
        end if

    end function test_perturb_vector

    logical(c_bool) function test_noise_level() bind(C)
        !
        ! this function tests the function which returns the relative size of the
        ! gradient perturbation
        !
        use opentrustregion, only: solver_settings_type, noise_level

        type(solver_settings_type) :: settings

        ! assume tests pass
        test_noise_level = .true.

        ! check that the perturbation is the reduction factor of the region with the
        ! gradient perturbation, which differ between the regions for the default
        ! settings, and vanishes without it
        settings%grad_noise = .true.
        if (abs(noise_level(settings, .true.) - settings%local_red_factor) > &
            0.0_rp .or. &
            abs(noise_level(settings, .false.) - settings%global_red_factor) > 0.0_rp) &
            then
            write(stderr, *) "test_noise_level failed: Perturbation not the "// &
                "reduction factor of the region."
            test_noise_level = .false.
        end if
        settings%grad_noise = .false.
        if (abs(noise_level(settings, .true.)) > 0.0_rp .or. &
            abs(noise_level(settings, .false.)) > 0.0_rp) then
            write(stderr, *) "test_noise_level failed: Perturbation does not "// &
                "vanish without gradient perturbation."
            test_noise_level = .false.
        end if

    end function test_noise_level

    logical(c_bool) function test_micro_conv_tol() bind(C)
        !
        ! this function tests the function which returns the residual norm at which the
        ! micro iterations stop
        !
        use opentrustregion, only: solver_settings_type, micro_conv_tol, &
                                   level_shift_local_thres, residual_norm_floor

        ! cases with the level shift, the gradient perturbation and the initial
        ! residual norm
        character(len=*), parameter :: case_names(5) = &
            [character(len=48) :: "vanishing level shift", "level shift at the "// &
             "threshold of the local region", "large negative level shift", &
             "gradient perturbation above the reduction factor", &
             "vanishing initial residual norm"]
        real(rp), parameter :: &
            case_level_shifts(5) = &
                [0.0_rp, level_shift_local_thres, -1.0_rp, 0.0_rp, 0.0_rp], &
            case_noises(5) = [0.0_rp, 0.0_rp, 0.0_rp, 0.5_rp, 0.0_rp], &
            case_residual_norms(5) = [2.0_rp, 2.0_rp, 2.0_rp, 2.0_rp, 1e-20_rp]
        type(solver_settings_type) :: settings
        real(rp) :: ref_conv_tols(5)
        integer(ip) :: i

        ! assume tests pass
        test_micro_conv_tol = .true.

        ! the stopping tolerance is the reduction factor of the local region inside the
        ! threshold of the level shift and that of the global region outside, which
        ! differ for the default settings, the gradient perturbation where it is
        ! larger, relative to the initial residual norm and never below the residual
        ! norm floor
        ref_conv_tols = [ &
            2.0_rp * settings%local_red_factor, 2.0_rp * settings%global_red_factor, &
            2.0_rp * settings%global_red_factor, 1.0_rp, residual_norm_floor]
        do i = 1, size(case_names)
            if (abs(micro_conv_tol(case_level_shifts(i), case_noises(i), &
                                   case_residual_norms(i), settings) - &
                    ref_conv_tols(i)) > 0.0_rp) then
                write(stderr, *) "test_micro_conv_tol failed: Stopping tolerance "// &
                    "wrong for "//trim(case_names(i))//"."
                test_micro_conv_tol = .false.
            end if
        end do

    end function test_micro_conv_tol

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
