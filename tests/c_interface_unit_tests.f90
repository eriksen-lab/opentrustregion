! Copyright (C) 2025- Jonas Greiner
!
! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.

module c_interface_unit_tests

    use opentrustregion, only: rp, ip, stderr
    use c_interface, only: c_rp, c_ip, update_orbs_c_type, hess_x_c_type, &
                           obj_func_c_type, precond_c_type, project_c_type, &
                           conv_check_c_type, logger_c_type
    use test_reference, only: tol, tol_c, n_param, n_param_c
    use, intrinsic :: iso_c_binding, only: c_bool, c_ptr, c_loc, c_funptr, c_funloc, &
                                           c_char, c_associated, c_null_ptr, &
                                           c_null_char, c_null_funptr, c_f_pointer

    implicit none

    ! create function pointers to ensure that routines comply with interface
    procedure(update_orbs_c_type), pointer :: mock_update_orbs_ptr => mock_update_orbs
    procedure(update_orbs_c_type), pointer :: mock_update_orbs_no_hess_x_ptr => &
        mock_update_orbs_no_hess_x
    procedure(hess_x_c_type), pointer :: mock_hess_x_ptr => mock_hess_x
    procedure(obj_func_c_type), pointer :: mock_obj_func_ptr => mock_obj_func
    procedure(precond_c_type), pointer :: mock_precond_ptr => mock_precond
    procedure(project_c_type), pointer :: mock_project_ptr => mock_project
    procedure(conv_check_c_type), pointer :: mock_conv_check_ptr => mock_conv_check
    procedure(logger_c_type), pointer :: mock_logger_ptr => mock_logger
    procedure(precond_c_type), pointer :: mock_stability_precond_ptr => &
        mock_stability_precond
    procedure(project_c_type), pointer :: mock_stability_project_ptr => &
        mock_stability_project
    procedure(logger_c_type), pointer :: mock_stability_logger_ptr => &
        mock_stability_logger

contains

    function mock_update_orbs(kappa, func, grad, h_diag, hess_x_c_funptr, context_c) &
        result(error) bind(C)
        !
        ! this function is a test function for the orbital update C function
        !
        use test_reference, only: check_host_context_c, host_context_error

        real(c_rp), intent(in) :: kappa(*)
        real(c_rp), intent(out) :: func, grad(*), h_diag(*)
        type(c_funptr), intent(inout) :: hess_x_c_funptr
        type(c_ptr), intent(in), value :: context_c
        integer(c_ip) :: error

        ! check host context
        call check_host_context_c(context_c)

        func = sum(kappa(:n_param))

        grad(:n_param) = 2 * kappa(:n_param)

        h_diag(:n_param) = 3 * kappa(:n_param)

        hess_x_c_funptr = c_funloc(mock_hess_x)

        error = host_context_error(context_c)

    end function mock_update_orbs

    function mock_update_orbs_no_hess_x(kappa, func, grad, h_diag, hess_x_c_funptr, &
                                        context_c) result(error) bind(C)
        !
        ! this function is a test function for an orbital update C function which
        ! reports success but does not provide a Hessian linear transformation
        !
        real(c_rp), intent(in) :: kappa(*)
        real(c_rp), intent(out) :: func, grad(*), h_diag(*)
        type(c_funptr), intent(inout) :: hess_x_c_funptr
        type(c_ptr), intent(in), value :: context_c
        integer(c_ip) :: error

        error = mock_update_orbs(kappa, func, grad, h_diag, hess_x_c_funptr, context_c)
        hess_x_c_funptr = c_null_funptr

    end function mock_update_orbs_no_hess_x

    function mock_hess_x(x, hess_x, context_c) result(error) bind(C)
        !
        ! this function is a test function for the Hessian linear transformation C
        ! function
        !
        use test_reference, only: check_host_context_c, host_context_error

        real(c_rp), intent(in) :: x(*)
        real(c_rp), intent(out) :: hess_x(*)
        type(c_ptr), intent(in), value :: context_c
        integer(c_ip) :: error

        ! check host context
        call check_host_context_c(context_c)

        hess_x(:n_param) = 4 * x(:n_param)

        error = host_context_error(context_c)

    end function mock_hess_x

    function mock_obj_func(kappa, func, context_c) result(error) bind(C)
        !
        ! this function is a test function for the C objective function
        !
        use test_reference, only: check_host_context_c, host_context_error

        real(c_rp), intent(in) :: kappa(*)
        real(c_rp), intent(out) :: func
        type(c_ptr), intent(in), value :: context_c
        integer(c_ip) :: error

        ! check host context
        call check_host_context_c(context_c)

        func = sum(kappa(:n_param))

        error = host_context_error(context_c)

    end function mock_obj_func

    function mock_precond(residual, mu, precond_residual, context_c) result(error) &
        bind(C)
        !
        ! this function is a test function for the C preconditioner function
        !
        use test_reference, only: check_host_context_c, host_context_error

        real(c_rp), intent(in) :: residual(*), mu
        real(c_rp), intent(out) :: precond_residual(*)
        type(c_ptr), intent(in), value :: context_c
        integer(c_ip) :: error

        ! check host context
        call check_host_context_c(context_c)

        precond_residual(:n_param) = mu * residual(:n_param)

        error = host_context_error(context_c)

    end function mock_precond

    function mock_project(vector, context_c) result(error) bind(C)
        !
        ! this function is a test function for the C projection function
        !
        use test_reference, only: check_host_context_c, host_context_error

        real(c_rp), intent(inout), target :: vector(*)
        type(c_ptr), intent(in), value :: context_c
        integer(c_ip) :: error

        ! check host context
        call check_host_context_c(context_c)

        vector(:n_param) = 2 * vector(:n_param)

        error = host_context_error(context_c)

    end function mock_project

    function mock_conv_check(converged, context_c) result(error) bind(C)
        !
        ! this function is a test function for the convergence check function
        !
        use test_reference, only: check_host_context_c, host_context_error

        logical(c_bool), intent(out) :: converged
        type(c_ptr), intent(in), value :: context_c
        integer(c_ip) :: error

        ! check host context
        call check_host_context_c(context_c)

        converged = .true.

        error = host_context_error(context_c)

    end function mock_conv_check

    subroutine mock_logger(message_c, context_c) bind(C)
        !
        ! this subroutine is a test subroutine for the C logging function
        !
        use test_reference, only: host_context_type, check_host_context_c

        character(kind=c_char), intent(in) :: message_c(*)
        type(c_ptr), intent(in), value :: context_c

        character(len=4) :: message
        type(host_context_type), pointer :: context

        ! check host context
        call check_host_context_c(context_c)

        ! record call with null-terminated test message in host context
        message = transfer(message_c(1:4), message)
        if (message == "test" .and. message_c(5) == c_null_char .and. &
            c_associated(context_c)) then
            call c_f_pointer(context_c, context)
            context%logger_called = .true.
        end if

    end subroutine mock_logger

    function mock_stability_precond(residual, mu, precond_residual, context_c) &
        result(error) bind(C)
        !
        ! this function is a test function for a C preconditioner function the nested
        ! stability check settings provide of their own, which behaves like the
        ! solver's but is a different function
        !
        real(c_rp), intent(in) :: residual(*), mu
        real(c_rp), intent(out) :: precond_residual(*)
        type(c_ptr), intent(in), value :: context_c
        integer(c_ip) :: error

        error = mock_precond(residual, mu, precond_residual, context_c)

    end function mock_stability_precond

    function mock_stability_project(vector, context_c) result(error) bind(C)
        !
        ! this function is a test function for a C projection function the nested
        ! stability check settings provide of their own, which behaves like the
        ! solver's but is a different function
        !
        real(c_rp), intent(inout), target :: vector(*)
        type(c_ptr), intent(in), value :: context_c
        integer(c_ip) :: error

        error = mock_project(vector, context_c)

    end function mock_stability_project

    subroutine mock_stability_logger(message_c, context_c) bind(C)
        !
        ! this subroutine is a test subroutine for a C logging function the nested
        ! stability check settings provide of their own, which behaves like the
        ! solver's but is a different function
        !
        character(kind=c_char), intent(in) :: message_c(*)
        type(c_ptr), intent(in), value :: context_c

        call mock_logger(message_c, context_c)

    end subroutine mock_stability_logger

    logical(c_bool) function test_solver_c_wrapper() bind(C)
        !
        ! this function tests the C wrapper for the solver
        !
        use c_interface, only: solver_settings_type_c, solver, solver_c_wrapper
        use opentrustregion, only: standard_solver => solver, solver_settings_type, &
                                   default_stability_settings
        use opentrustregion_mock, only: &
            mock_solver, test_passed, mock_error, mock_n_update_orbs, mock_n_hess_x, &
            mock_stability_n_hess_x, mock_solver_n_nested_calls, &
            received_solver_settings, received_stability_callbacks
        use test_reference, only: get_reference_solver_values, stability_host_context, &
                                  host_context, arm_host_context_c, &
                                  host_context_reached, unset_callbacks, &
                                  ref_solver_settings, operator(/=)

        type(c_funptr) :: update_orbs_c_funptr, obj_func_c_funptr
        type(solver_settings_type_c) :: settings
        type(solver_settings_type) :: expected_settings
        procedure(precond_c_type), pointer :: expected_precond
        procedure(project_c_type), pointer :: expected_project
        procedure(logger_c_type), pointer :: expected_logger
        type(c_ptr) :: expected_context
        integer(ip) :: expected_nested_calls
        integer(c_ip) :: error
        integer(ip) :: icase
        character(len=34), parameter :: case_names(4) = &
            [character(len=34) :: "without nested context", "with nested context", &
             "with nested callback functions", "with uninitialized nested settings"]

        ! assume tests pass
        test_solver_c_wrapper = .true.

        ! inject mock function
        solver => mock_solver

        ! get C function pointers to Fortran functions
        update_orbs_c_funptr = c_funloc(mock_update_orbs)
        obj_func_c_funptr = c_funloc(mock_obj_func)

        ! run with nested stability check settings that provide neither callback
        ! functions nor a context of their own, that provide a context, that provide
        ! a context and callback functions and that provide both but were not
        ! initialized, the callback bundle of the nested settings has to hold the
        ! nested settings' own callback functions and context where these were provided
        ! and initialized and the solver's otherwise
        do icase = 1, size(case_names)
            ! associate optional settings with the reference values and callback
            ! functions
            call get_reference_solver_values(settings)
            call unset_callbacks(settings%stability_settings)
            settings%stability_settings%context = c_null_ptr
            settings%precond = c_funloc(mock_precond)
            settings%project = c_funloc(mock_project)
            settings%conv_check = c_funloc(mock_conv_check)
            settings%logger = c_funloc(mock_logger)

            ! set host contexts and the nested settings' own callback functions
            call arm_host_context_c(settings%context)
            if (icase >= 2) &
                settings%stability_settings%context = c_loc(stability_host_context)
            if (icase >= 3) then
                settings%stability_settings%precond = c_funloc(mock_stability_precond)
                settings%stability_settings%project = c_funloc(mock_stability_project)
                settings%stability_settings%logger = c_funloc(mock_stability_logger)
            end if
            if (icase == 4) settings%stability_settings%initialized = .false.

            ! set expected settings, nested callback functions and context
            expected_settings = ref_solver_settings
            if (icase == 4) &
                expected_settings%stability_settings = default_stability_settings
            if (icase == 3) then
                expected_precond => mock_stability_precond
                expected_project => mock_stability_project
                expected_logger => mock_stability_logger
            else
                expected_precond => mock_precond
                expected_project => mock_project
                expected_logger => mock_logger
            end if
            if (icase == 2 .or. icase == 3) then
                expected_context = c_loc(stability_host_context)
                expected_nested_calls = mock_solver_n_nested_calls
            else
                expected_context = c_loc(host_context)
                expected_nested_calls = 0
            end if

            ! clear the result of the mock so that a missing call is detected
            test_passed = .false.

            ! call solver
            error = solver_c_wrapper(update_orbs_c_funptr, obj_func_c_funptr, &
                                     n_param_c, settings)

            ! check if logging subroutine was correctly called
            if (.not. host_context%logger_called) then
                test_solver_c_wrapper = .false.
                write(stderr, *) "test_solver_c_wrapper failed: Called logging "// &
                    "subroutine wrong "//trim(case_names(icase))//"."
            end if

            ! check if optional settings are correctly passed
            if (received_solver_settings /= expected_settings) then
                test_solver_c_wrapper = .false.
                write(stderr, *) "test_solver_c_wrapper failed: Passed optional "// &
                    "settings associated with wrong values "// &
                    trim(case_names(icase))//"."
            end if

            ! check the callback bundle of the nested settings
            if (.not. associated(received_stability_callbacks%precond, &
                                 expected_precond)) then
                test_solver_c_wrapper = .false.
                write(stderr, *) "test_solver_c_wrapper failed: Preconditioner of "// &
                    "internal stability check wrong "//trim(case_names(icase))//"."
            end if
            if (.not. associated(received_stability_callbacks%project, &
                                 expected_project)) then
                test_solver_c_wrapper = .false.
                write(stderr, *) "test_solver_c_wrapper failed: Projection of "// &
                    "internal stability check wrong "//trim(case_names(icase))//"."
            end if
            if (.not. associated(received_stability_callbacks%logger, &
                                 expected_logger)) then
                test_solver_c_wrapper = .false.
                write(stderr, *) "test_solver_c_wrapper failed: Logging function "// &
                    "of internal stability check wrong "//trim(case_names(icase))//"."
            end if
            if (.not. c_associated(received_stability_callbacks%host_context, &
                                   expected_context)) then
                test_solver_c_wrapper = .false.
                write(stderr, *) "test_solver_c_wrapper failed: Context of "// &
                    "internal stability check wrong "//trim(case_names(icase))//"."
            end if

            ! check if output variables are as expected
            if (error /= int(mock_error, kind=c_ip)) then
                test_solver_c_wrapper = .false.
                write(stderr, *) "test_solver_c_wrapper failed: Returned error "// &
                    "code wrong "//trim(case_names(icase))//"."
            end if

            ! check if output fields are written back with the values set by the solver
            if (settings%max_precision_reached) then
                test_solver_c_wrapper = .false.
                write(stderr, *) "test_solver_c_wrapper failed: Returned maximum "// &
                    "precision reached flag wrong "//trim(case_names(icase))//"."
            end if
            if (settings%n_update_orbs /= int(mock_n_update_orbs, kind=c_ip)) then
                test_solver_c_wrapper = .false.
                write(stderr, *) "test_solver_c_wrapper failed: Returned number of "// &
                    "orbital updates wrong "//trim(case_names(icase))//"."
            end if
            if (settings%n_hess_x /= int(mock_n_hess_x, kind=c_ip)) then
                test_solver_c_wrapper = .false.
                write(stderr, *) "test_solver_c_wrapper failed: Returned number of "// &
                    "Hessian linear transformations wrong "//trim(case_names(icase))// &
                    "."
            end if
            if (settings%stability_settings%n_hess_x /= &
                int(mock_stability_n_hess_x, kind=c_ip)) then
                test_solver_c_wrapper = .false.
                write(stderr, *) "test_solver_c_wrapper failed: Returned number of "// &
                    "Hessian linear transformations of internal stability check "// &
                    "wrong "//trim(case_names(icase))//"."
            end if

            ! check that the callback functions of the internal stability check received
            ! the context of the nested settings exactly when it was provided and
            ! initialized
            if (stability_host_context%n_calls /= expected_nested_calls) then
                test_solver_c_wrapper = .false.
                write(stderr, *) "test_solver_c_wrapper failed: Callback functions "// &
                    "of internal stability check received context of nested "// &
                    "settings wrongly "//trim(case_names(icase))//"."
            end if

            ! check that the host context reached the callback functions unchanged
            test_solver_c_wrapper = test_solver_c_wrapper .and. &
                                    host_context_reached("solver_c_wrapper")

            ! check if test has passed
            test_solver_c_wrapper = test_solver_c_wrapper .and. test_passed
        end do

        ! restore the procedure pointer so later tests do not inherit the mock
        solver => standard_solver

    end function test_solver_c_wrapper

    logical(c_bool) function test_stability_check_c_wrapper() bind(C)
        !
        ! this function tests the C wrapper for the stability check
        !
        use c_interface, only: stability_settings_type_c, stability_check, &
                               stability_check_c_wrapper
        use opentrustregion, only: standard_stability_check => stability_check
        use opentrustregion_mock, only: mock_stability_check, test_passed, mock_error, &
                                        mock_stability_check_n_hess_x
        use test_reference, only: get_reference_stability_values, host_context, &
                                  arm_host_context_c, host_context_reached

        type(c_funptr) :: hess_x_c_funptr
        real(c_rp) :: h_diag(n_param)
        real(c_rp), target :: kappa(n_param)
        type(stability_settings_type_c) :: settings
        logical(c_bool) :: stable
        type(c_ptr) :: kappa_c_ptr
        integer(c_ip) :: error
        integer(ip) :: icase
        character(len=26), parameter :: case_names(2) = &
            [character(len=26) :: "without returned direction", &
             "with returned direction"]

        ! assume tests pass
        test_stability_check_c_wrapper = .true.

        ! inject mock function
        stability_check => mock_stability_check

        ! get C function pointers to Fortran functions
        hess_x_c_funptr = c_funloc(mock_hess_x)

        ! initialize Hessian diagonal
        h_diag = 3.0_c_rp

        ! run once without and once with a returned direction
        do icase = 1, size(case_names)
            ! associate optional settings with the reference values and callback
            ! functions
            call get_reference_stability_values(settings)
            settings%precond = c_funloc(mock_precond)
            settings%project = c_funloc(mock_project)
            settings%logger = c_funloc(mock_logger)

            ! set host context
            call arm_host_context_c(settings%context)

            ! associate returned direction pointer
            kappa = 0.0_c_rp
            if (icase == 1) then
                kappa_c_ptr = c_null_ptr
            else
                kappa_c_ptr = c_loc(kappa)
            end if

            ! clear the result of the mock so that a missing call is detected
            test_passed = .false.

            ! call stability check
            stable = .false.
            error = stability_check_c_wrapper(h_diag, hess_x_c_funptr, n_param_c, &
                                              stable, settings, kappa_c_ptr)

            ! check if logging subroutine was correctly called
            if (.not. host_context%logger_called) then
                test_stability_check_c_wrapper = .false.
                write(stderr, *) "test_stability_check_c_wrapper failed: Called "// &
                    "logging subroutine wrong "//trim(case_names(icase))//"."
            end if

            ! check if output variables are as expected
            if (.not. stable) then
                test_stability_check_c_wrapper = .false.
                write(stderr, *) "test_stability_check_c_wrapper failed: Returned "// &
                    "stability boolean wrong "//trim(case_names(icase))//"."
            end if
            if (error /= int(mock_error, kind=c_ip)) then
                test_stability_check_c_wrapper = .false.
                write(stderr, *) "test_stability_check_c_wrapper failed: Returned "// &
                    "error code wrong "//trim(case_names(icase))//"."
            end if
            if (icase == 2 .and. any(abs(kappa - 1.0_c_rp) > tol_c)) then
                test_stability_check_c_wrapper = .false.
                write(stderr, *) "test_stability_check_c_wrapper failed: Returned "// &
                    "direction wrong."
            end if

            ! check if output field is written back with the value set by the
            ! stability check
            if (settings%n_hess_x /= int(mock_stability_check_n_hess_x, kind=c_ip)) then
                test_stability_check_c_wrapper = .false.
                write(stderr, *) "test_stability_check_c_wrapper failed: Returned "// &
                    "number of Hessian linear transformations wrong "// &
                    trim(case_names(icase))//"."
            end if

            ! check that the host context reached the callback functions unchanged
            test_stability_check_c_wrapper = &
                test_stability_check_c_wrapper .and. &
                host_context_reached("stability_check_c_wrapper")

            ! check if test has passed
            test_stability_check_c_wrapper = test_stability_check_c_wrapper .and. &
                                             test_passed
        end do

        ! restore the procedure pointer so later tests do not inherit the mock
        stability_check => standard_stability_check

    end function test_stability_check_c_wrapper

    logical(c_bool) function test_store_optional_c_callbacks() bind(C)
        !
        ! this function tests the subroutine that stores the optional C callback
        ! functions and the host context of C settings in a callback bundle
        !
        use c_interface, only: c_callbacks_type, store_optional_c_callbacks
        use test_reference, only: host_context

        type(c_callbacks_type) :: callbacks

        ! assume test passes
        test_store_optional_c_callbacks = .true.

        ! store the callback functions and host context of initialized settings
        call store_optional_c_callbacks( &
            callbacks, .true._c_bool, c_funloc(mock_precond), c_funloc(mock_project), &
            c_funloc(mock_conv_check), c_funloc(mock_logger), c_loc(host_context))
        if (.not. (associated(callbacks%precond, mock_precond) .and. &
                   associated(callbacks%project, mock_project) .and. &
                   associated(callbacks%conv_check, mock_conv_check) .and. &
                   associated(callbacks%logger, mock_logger))) then
            test_store_optional_c_callbacks = .false.
            write(stderr, *) "test_store_optional_c_callbacks failed: Callback "// &
                "functions of initialized settings not stored."
        end if
        if (.not. c_associated(callbacks%host_context, c_loc(host_context))) then
            test_store_optional_c_callbacks = .false.
            write(stderr, *) "test_store_optional_c_callbacks failed: Host context "// &
                "of initialized settings not stored."
        end if

        ! store initialized settings that provide nothing, the callback functions and
        ! host context already in the bundle have to be kept
        call store_optional_c_callbacks(callbacks, .true._c_bool, c_null_funptr, &
                                        c_null_funptr, c_null_funptr, c_null_funptr, &
                                        c_null_ptr)
        if (.not. (associated(callbacks%precond, mock_precond) .and. &
                   associated(callbacks%project, mock_project) .and. &
                   associated(callbacks%conv_check, mock_conv_check) .and. &
                   associated(callbacks%logger, mock_logger))) then
            test_store_optional_c_callbacks = .false.
            write(stderr, *) "test_store_optional_c_callbacks failed: Callback "// &
                "functions in the bundle replaced by ones that were not provided."
        end if
        if (.not. c_associated(callbacks%host_context, c_loc(host_context))) then
            test_store_optional_c_callbacks = .false.
            write(stderr, *) "test_store_optional_c_callbacks failed: Host context "// &
                "in the bundle replaced by one that was not provided."
        end if

        ! store the callback functions and host context of settings that were not
        ! initialized in an empty bundle, nothing may be stored
        callbacks = c_callbacks_type()
        call store_optional_c_callbacks( &
            callbacks, .false._c_bool, c_funloc(mock_precond), c_funloc(mock_project), &
            c_funloc(mock_conv_check), c_funloc(mock_logger), c_loc(host_context))
        if (associated(callbacks%precond) .or. associated(callbacks%project) .or. &
            associated(callbacks%conv_check) .or. associated(callbacks%logger)) then
            test_store_optional_c_callbacks = .false.
            write(stderr, *) "test_store_optional_c_callbacks failed: Callback "// &
                "functions of settings that were not initialized stored."
        end if
        if (c_associated(callbacks%host_context)) then
            test_store_optional_c_callbacks = .false.
            write(stderr, *) "test_store_optional_c_callbacks failed: Host context "// &
                "of settings that were not initialized stored."
        end if

    end function test_store_optional_c_callbacks

    logical(c_bool) function test_update_orbs_f_wrapper() bind(C)
        !
        ! this function tests the Fortran wrapper for the orbital update
        !
        use opentrustregion, only: update_orbs_type, hess_x_type
        use c_interface, only: c_callbacks_type, update_orbs_f_wrapper
        use test_reference, only: check_update_orbs_funptr, host_context, &
                                  arm_host_context_c, host_context_reached

        procedure(update_orbs_type), pointer :: update_orbs_funptr
        type(c_callbacks_type), target :: callbacks, stability_callbacks
        class(*), pointer :: context, foreign_context
        real(rp), target :: foreign_context_target
        procedure(hess_x_type), pointer :: hess_x_funptr
        real(rp) :: kappa(n_param), func, grad(n_param), h_diag(n_param)
        integer(ip) :: error

        ! inject mock C function and hand the bundle to the wrapper as its context
        callbacks%update_orbs => mock_update_orbs
        call arm_host_context_c(callbacks%host_context)
        context => callbacks

        ! attach the bundle of an internal stability check, which has to receive the
        ! returned Hessian linear transformation as well
        callbacks%stability => stability_callbacks

        ! get pointer to subroutine
        update_orbs_funptr => update_orbs_f_wrapper

        ! test orbital update wrapper
        test_update_orbs_f_wrapper = check_update_orbs_funptr( &
            update_orbs_funptr, "update_orbs_f_wrapper", "", context)

        ! check that the wrapper handed the host context to the C function
        test_update_orbs_f_wrapper = test_update_orbs_f_wrapper .and. logical( &
            host_context_reached("update_orbs_f_wrapper"), kind=c_bool)

        ! check that the returned Hessian linear transformation is handed to the bundle
        ! of the internal stability check
        if (.not. associated(stability_callbacks%hess_x, mock_hess_x)) then
            test_update_orbs_f_wrapper = .false.
            write(stderr, *) "test_update_orbs_f_wrapper failed: Hessian linear "// &
                "transformation not handed to internal stability check."
        end if

        ! an orbital update that succeeds without providing a Hessian linear
        ! transformation is an error
        callbacks%update_orbs => mock_update_orbs_no_hess_x
        kappa = 1.0_rp
        call update_orbs_f_wrapper(kappa, func, grad, h_diag, hess_x_funptr, error, &
                                   context)
        if (error /= 1) then
            test_update_orbs_f_wrapper = .false.
            write(stderr, *) "test_update_orbs_f_wrapper failed: Did not report a "// &
                "missing Hessian linear transformation."
        end if

        ! check that an error of the C function is passed on rather than replaced by
        ! the one of the missing Hessian linear transformation
        host_context%mock_error = 2
        call update_orbs_f_wrapper(kappa, func, grad, h_diag, hess_x_funptr, error, &
                                   context)
        if (error /= 2) then
            test_update_orbs_f_wrapper = .false.
            write(stderr, *) "test_update_orbs_f_wrapper failed: Did not pass on "// &
                "the error of the C function."
        end if
        host_context%mock_error = 0

        ! a context this module did not create is reported rather than dereferenced
        foreign_context_target = 0.0_rp
        foreign_context => foreign_context_target
        kappa = 1.0_rp
        call update_orbs_f_wrapper(kappa, func, grad, h_diag, hess_x_funptr, error, &
                                   foreign_context)
        if (error /= 1) then
            test_update_orbs_f_wrapper = .false.
            write(stderr, *) "test_update_orbs_f_wrapper failed: Did not report an "// &
                "invalid context."
        end if

    end function test_update_orbs_f_wrapper

    logical(c_bool) function test_hess_x_f_wrapper() bind(C)
        !
        ! this function tests the Fortran wrapper for the Hessian linear transformation
        !
        use opentrustregion, only: hess_x_type
        use c_interface, only: c_callbacks_type, hess_x_f_wrapper
        use test_reference, only: check_hess_x_funptr, host_context, &
                                  arm_host_context_c, host_context_reached

        procedure(hess_x_type), pointer :: hess_x_funptr
        type(c_callbacks_type), target :: callbacks
        class(*), pointer :: context, foreign_context
        real(rp), target :: foreign_context_target
        real(rp) :: x(n_param), hess_x(n_param)
        integer(ip) :: error

        ! inject mock C function and hand the bundle to the wrapper as its context
        callbacks%hess_x => mock_hess_x
        call arm_host_context_c(callbacks%host_context)
        context => callbacks

        ! get pointer to subroutine
        hess_x_funptr => hess_x_f_wrapper

        ! test Hessian linear transformation wrapper
        test_hess_x_f_wrapper = &
            check_hess_x_funptr(hess_x_funptr, "hess_x_f_wrapper", "", context)

        ! check that the wrapper handed the host context to the C function
        test_hess_x_f_wrapper = test_hess_x_f_wrapper .and. logical( &
            host_context_reached("hess_x_f_wrapper"), kind=c_bool)

        ! check that an error of the C function is passed on
        host_context%mock_error = 2
        x = 1.0_rp
        call hess_x_f_wrapper(x, hess_x, error, context)
        if (error /= 2) then
            test_hess_x_f_wrapper = .false.
            write(stderr, *) "test_hess_x_f_wrapper failed: Did not pass on the "// &
                "error of the C function."
        end if
        host_context%mock_error = 0

        ! a context this module did not create is reported rather than dereferenced
        x = 1.0_rp
        foreign_context_target = 0.0_rp
        foreign_context => foreign_context_target
        call hess_x_f_wrapper(x, hess_x, error, foreign_context)
        if (error /= 1) then
            test_hess_x_f_wrapper = .false.
            write(stderr, *) "test_hess_x_f_wrapper failed: Did not report an "// &
                "invalid context."
        end if

    end function test_hess_x_f_wrapper

    logical(c_bool) function test_obj_func_f_wrapper() bind(C)
        !
        ! this function tests the Fortran wrapper for the objective function
        !
        use opentrustregion, only: obj_func_type
        use c_interface, only: c_callbacks_type, obj_func_f_wrapper
        use test_reference, only: check_obj_func_funptr, host_context, &
                                  arm_host_context_c, host_context_reached

        procedure(obj_func_type), pointer :: obj_func_funptr
        type(c_callbacks_type), target :: callbacks
        class(*), pointer :: context, foreign_context
        real(rp), target :: foreign_context_target
        real(rp) :: kappa(n_param), func
        integer(ip) :: error

        ! inject mock C function and hand the bundle to the wrapper as its context
        callbacks%obj_func => mock_obj_func
        call arm_host_context_c(callbacks%host_context)
        context => callbacks

        ! get pointer to subroutine
        obj_func_funptr => obj_func_f_wrapper

        ! test objective function wrapper
        test_obj_func_f_wrapper = &
            check_obj_func_funptr(obj_func_funptr, "obj_func_f_wrapper", "", context)

        ! check that the wrapper handed the host context to the C function
        test_obj_func_f_wrapper = test_obj_func_f_wrapper .and. logical( &
            host_context_reached("obj_func_f_wrapper"), kind=c_bool)

        ! check that an error of the C function is passed on
        host_context%mock_error = 2
        kappa = 1.0_rp
        func = obj_func_f_wrapper(kappa, error, context)
        if (error /= 2) then
            test_obj_func_f_wrapper = .false.
            write(stderr, *) "test_obj_func_f_wrapper failed: Did not pass on the "// &
                "error of the C function."
        end if
        host_context%mock_error = 0

        ! a context this module did not create is reported rather than dereferenced
        foreign_context_target = 0.0_rp
        foreign_context => foreign_context_target
        kappa = 1.0_rp
        func = obj_func_f_wrapper(kappa, error, foreign_context)
        if (error /= 1) then
            test_obj_func_f_wrapper = .false.
            write(stderr, *) "test_obj_func_f_wrapper failed: Did not report an "// &
                "invalid context."
        end if

    end function test_obj_func_f_wrapper

    logical(c_bool) function test_precond_f_wrapper() bind(C)
        !
        ! this function tests the Fortran wrapper for the preconditioner function
        !
        use opentrustregion, only: precond_type
        use c_interface, only: c_callbacks_type, precond_f_wrapper
        use test_reference, only: check_precond_funptr, host_context, &
                                  arm_host_context_c, host_context_reached

        procedure(precond_type), pointer :: precond_funptr
        type(c_callbacks_type), target :: callbacks
        class(*), pointer :: context, foreign_context
        real(rp), target :: foreign_context_target
        real(rp) :: residual(n_param), precond_residual(n_param)
        integer(ip) :: error

        ! inject mock C function and hand the bundle to the wrapper as its context
        callbacks%precond => mock_precond
        call arm_host_context_c(callbacks%host_context)
        context => callbacks

        ! get pointer to subroutine
        precond_funptr => precond_f_wrapper

        ! test preconditioner wrapper
        test_precond_f_wrapper = &
            check_precond_funptr(precond_funptr, "precond_f_wrapper", "", context)

        ! check that the wrapper handed the host context to the C function
        test_precond_f_wrapper = test_precond_f_wrapper .and. logical( &
            host_context_reached("precond_f_wrapper"), kind=c_bool)

        ! check that an error of the C function is passed on
        host_context%mock_error = 2
        residual = 1.0_rp
        call precond_f_wrapper(residual, 1.0_rp, precond_residual, error, context)
        if (error /= 2) then
            test_precond_f_wrapper = .false.
            write(stderr, *) "test_precond_f_wrapper failed: Did not pass on the "// &
                "error of the C function."
        end if
        host_context%mock_error = 0

        ! a context this module did not create is reported rather than dereferenced
        foreign_context_target = 0.0_rp
        foreign_context => foreign_context_target
        residual = 1.0_rp
        call precond_f_wrapper(residual, 1.0_rp, precond_residual, error, &
                               foreign_context)
        if (error /= 1) then
            test_precond_f_wrapper = .false.
            write(stderr, *) "test_precond_f_wrapper failed: Did not report an "// &
                "invalid context."
        end if

    end function test_precond_f_wrapper

    logical(c_bool) function test_project_f_wrapper() bind(C)
        !
        ! this function tests the Fortran wrapper for the projection function
        !
        use opentrustregion, only: project_type
        use c_interface, only: c_callbacks_type, project_f_wrapper
        use test_reference, only: check_project_funptr, host_context, &
                                  arm_host_context_c, host_context_reached

        procedure(project_type), pointer :: project_funptr
        type(c_callbacks_type), target :: callbacks
        class(*), pointer :: context, foreign_context
        real(rp), target :: foreign_context_target
        real(rp) :: vector(n_param)
        integer(ip) :: error

        ! inject mock C function and hand the bundle to the wrapper as its context
        callbacks%project => mock_project
        call arm_host_context_c(callbacks%host_context)
        context => callbacks

        ! get pointer to subroutine
        project_funptr => project_f_wrapper

        ! test projection wrapper
        test_project_f_wrapper = &
            check_project_funptr(project_funptr, "project_f_wrapper", "", context)

        ! check that the wrapper handed the host context to the C function
        test_project_f_wrapper = test_project_f_wrapper .and. logical( &
            host_context_reached("project_f_wrapper"), kind=c_bool)

        ! check that an error of the C function is passed on
        host_context%mock_error = 2
        vector = 1.0_rp
        call project_f_wrapper(vector, error, context)
        if (error /= 2) then
            test_project_f_wrapper = .false.
            write(stderr, *) "test_project_f_wrapper failed: Did not pass on the "// &
                "error of the C function."
        end if
        host_context%mock_error = 0

        ! a context this module did not create is reported rather than dereferenced
        foreign_context_target = 0.0_rp
        foreign_context => foreign_context_target
        vector = 1.0_rp
        call project_f_wrapper(vector, error, foreign_context)
        if (error /= 1) then
            test_project_f_wrapper = .false.
            write(stderr, *) "test_project_f_wrapper failed: Did not report an "// &
                "invalid context."
        end if

    end function test_project_f_wrapper

    logical(c_bool) function test_conv_check_f_wrapper() bind(C)
        !
        ! this function tests the Fortran wrapper for the convergence check function
        !
        use opentrustregion, only: conv_check_type
        use c_interface, only: c_callbacks_type, conv_check_f_wrapper
        use test_reference, only: check_conv_check_funptr, host_context, &
                                  arm_host_context_c, host_context_reached

        procedure(conv_check_type), pointer :: conv_check_funptr
        type(c_callbacks_type), target :: callbacks
        class(*), pointer :: context, foreign_context
        real(rp), target :: foreign_context_target
        logical :: converged
        integer(ip) :: error

        ! inject mock C function and hand the bundle to the wrapper as its context
        callbacks%conv_check => mock_conv_check
        call arm_host_context_c(callbacks%host_context)
        context => callbacks

        ! get pointer to subroutine
        conv_check_funptr => conv_check_f_wrapper

        ! test convergence check wrapper
        test_conv_check_f_wrapper = check_conv_check_funptr( &
            conv_check_funptr, "conv_check_f_wrapper", "", context)

        ! check that the wrapper handed the host context to the C function
        test_conv_check_f_wrapper = test_conv_check_f_wrapper .and. logical( &
            host_context_reached("conv_check_f_wrapper"), kind=c_bool)

        ! check that an error of the C function is passed on
        host_context%mock_error = 2
        converged = conv_check_f_wrapper(error, context)
        if (error /= 2) then
            test_conv_check_f_wrapper = .false.
            write(stderr, *) "test_conv_check_f_wrapper failed: Did not pass on "// &
                "the error of the C function."
        end if
        host_context%mock_error = 0

        ! a context this module did not create is reported rather than dereferenced
        foreign_context_target = 0.0_rp
        foreign_context => foreign_context_target
        converged = conv_check_f_wrapper(error, foreign_context)
        if (error /= 1) then
            test_conv_check_f_wrapper = .false.
            write(stderr, *) "test_conv_check_f_wrapper failed: Did not report an "// &
                "invalid context."
        end if

    end function test_conv_check_f_wrapper

    logical(c_bool) function test_logger_f_wrapper() bind(C)
        !
        ! this function tests the Fortran wrapper for the logging function
        !
        use c_interface, only: c_callbacks_type, logger_f_wrapper
        use test_reference, only: host_context, arm_host_context_c, host_context_reached

        type(c_callbacks_type), target :: callbacks
        class(*), pointer :: context, foreign_context
        real(rp), target :: foreign_context_target

        ! assume tests pass
        test_logger_f_wrapper = .true.

        ! inject mock C function and hand the bundle to the wrapper as its context
        callbacks%logger => mock_logger
        call arm_host_context_c(callbacks%host_context)
        context => callbacks

        ! call subroutine
        call logger_f_wrapper("test", context)

        ! check if logging test boolean is as expected
        if (.not. host_context%logger_called) then
            test_logger_f_wrapper = .false.
            write(stderr, *) "test_logger_f_wrapper failed: Returned logging "// &
                "subroutine wrong."
        end if

        ! check that the wrapper handed the host context to the C function
        test_logger_f_wrapper = test_logger_f_wrapper .and. logical( &
            host_context_reached("logger_f_wrapper"), kind=c_bool)

        ! a context this module did not create is dropped rather than dereferenced,
        ! the logger has no error channel so it must simply not be called
        foreign_context_target = 0.0_rp
        foreign_context => foreign_context_target
        host_context%logger_called = .false.
        call logger_f_wrapper("test", foreign_context)
        if (host_context%logger_called) then
            test_logger_f_wrapper = .false.
            write(stderr, *) "test_logger_f_wrapper failed: Called logging "// &
                "subroutine with an invalid context."
        end if

    end function test_logger_f_wrapper

    logical(c_bool) function test_init_solver_settings_c() bind(C)
        !
        ! this function tests that the solver settings initialization routine correctly
        ! initializes all settings to their default values
        !
        use c_interface, only: solver_settings_type_c, init_solver_settings_c
        use opentrustregion, only: default_solver_settings
        use test_reference, only: operator(/=), host_context, callbacks_unset

        type(solver_settings_type_c) :: settings

        ! assume test passes
        test_init_solver_settings_c = .true.

        ! set callback functions and host contexts which the initialization has to
        ! discard
        settings%precond = c_funloc(mock_precond)
        settings%project = c_funloc(mock_project)
        settings%conv_check = c_funloc(mock_conv_check)
        settings%logger = c_funloc(mock_logger)
        settings%context = c_loc(host_context)
        settings%stability_settings%precond = c_funloc(mock_precond)
        settings%stability_settings%project = c_funloc(mock_project)
        settings%stability_settings%logger = c_funloc(mock_logger)
        settings%stability_settings%context = c_loc(host_context)

        ! initialize settings
        call init_solver_settings_c(settings)

        ! check function pointers and host contexts
        if (.not. callbacks_unset(settings)) then
            write(stderr, *) "test_init_solver_settings_c failed: Function "// &
                "pointers not discarded."
            test_init_solver_settings_c = .false.
        end if
        if (c_associated(settings%context) .or. &
            c_associated(settings%stability_settings%context)) then
            write(stderr, *) "test_init_solver_settings_c failed: Host contexts "// &
                "not discarded."
            test_init_solver_settings_c = .false.
        end if

        ! check settings
        if (settings /= default_solver_settings) then
            write(stderr, *) "test_init_solver_settings_c failed: Settings not "// &
                "initialized correctly."
            test_init_solver_settings_c = .false.
        end if

    end function test_init_solver_settings_c

    logical(c_bool) function test_init_stability_settings_c() bind(C)
        !
        ! this function tests that the stability check settings initialization routine
        ! correctly initializes all settings to their default values
        !
        use c_interface, only: stability_settings_type_c, init_stability_settings_c
        use opentrustregion, only: default_stability_settings
        use test_reference, only: operator(/=), host_context, callbacks_unset

        type(stability_settings_type_c) :: settings

        ! assume test passes
        test_init_stability_settings_c = .true.

        ! set callback functions and host context which the initialization has to
        ! discard
        settings%precond = c_funloc(mock_precond)
        settings%project = c_funloc(mock_project)
        settings%logger = c_funloc(mock_logger)
        settings%context = c_loc(host_context)

        ! initialize settings
        call init_stability_settings_c(settings)

        ! check function pointers and host context
        if (.not. callbacks_unset(settings)) then
            write(stderr, *) "test_init_stability_settings_c failed: Function "// &
                "pointers not discarded."
            test_init_stability_settings_c = .false.
        end if
        if (c_associated(settings%context)) then
            write(stderr, *) "test_init_stability_settings_c failed: Host context "// &
                "not discarded."
            test_init_stability_settings_c = .false.
        end if

        ! check settings
        if (settings /= default_stability_settings) then
            write(stderr, *) "test_init_stability_settings_c failed: Settings not "// &
                "initialized correctly."
            test_init_stability_settings_c = .false.
        end if

    end function test_init_stability_settings_c

    logical(c_bool) function test_assign_solver_f_c() bind(C)
        !
        ! this function tests that the function that converts solver settings from C to
        ! Fortran correctly perform this conversion
        !
        use c_interface, only: solver_settings_type_c, assignment(=)
        use opentrustregion, only: solver_settings_type, default_solver_settings
        use test_reference, only: ref_solver_settings, get_reference_solver_values, &
                                  operator(/=), callbacks_unset, unset_callbacks, &
                                  callbacks_wrapped

        type(solver_settings_type_c) :: settings_c
        type(solver_settings_type) :: settings, expected_settings
        integer(ip) :: i
        character(len=21), parameter :: logical_names(3) = &
            [character(len=21) :: "stability", "line_search", "max_precision_reached"]

        ! assume test passes
        test_assign_solver_f_c = .true.

        ! initialize the C settings with the reference values and callback functions
        call get_reference_solver_values(settings_c)

        ! convert to Fortran settings
        settings = settings_c

        ! check that every callback function, including those of the nested stability
        ! check settings, is converted to its wrapper
        if (.not. callbacks_wrapped(settings)) then
            write(stderr, *) "test_assign_solver_f_c failed: Callback functions "// &
                "not converted to their wrappers."
            test_assign_solver_f_c = .false.
        end if

        ! check against reference values
        if (settings /= ref_solver_settings) then
            write(stderr, *) "test_assign_solver_f_c failed: Settings not "// &
                "converted correctly."
            test_assign_solver_f_c = .false.
        end if

        ! convert again with only one of the logicals other than initialized set at a
        ! time, since the logicals cannot all differ from each other and from their
        ! default values, so that a logical that is swapped with another or not
        ! converted is detected
        do i = 1, size(logical_names)
            call get_reference_solver_values(settings_c, &
                                             trim(logical_names(i))//c_null_char)
            settings_c%initialized = .true.
            settings_c%stability_settings%initialized = .true.
            expected_settings = ref_solver_settings
            expected_settings%stability = logical_names(i) == "stability"
            expected_settings%line_search = logical_names(i) == "line_search"
            expected_settings%max_precision_reached = &
                logical_names(i) == "max_precision_reached"
            settings = settings_c
            if (settings /= expected_settings) then
                write(stderr, *) "test_assign_solver_f_c failed: Settings with "// &
                    "only "//trim(logical_names(i))//" set not converted correctly."
                test_assign_solver_f_c = .false.
            end if
        end do

        ! convert initialized C settings without callback functions and check that no
        ! callback functions are associated
        call unset_callbacks(settings_c)
        settings = settings_c
        if (.not. callbacks_unset(settings)) then
            write(stderr, *) "test_assign_solver_f_c failed: Function pointers "// &
                "associated for callback functions that were not provided."
            test_assign_solver_f_c = .false.
        end if

        ! convert C settings with callback functions that were not initialized, the
        ! custom values and callback functions they still carry have to be replaced by
        ! the default settings
        call get_reference_solver_values(settings_c)
        settings_c%initialized = .false.
        settings = settings_c
        if (.not. callbacks_unset(settings)) then
            write(stderr, *) "test_assign_solver_f_c failed: Function pointers "// &
                "converted for settings that were not initialized."
            test_assign_solver_f_c = .false.
        end if
        if (settings /= default_solver_settings) then
            write(stderr, *) "test_assign_solver_f_c failed: Settings that were "// &
                "not initialized not converted to default values."
            test_assign_solver_f_c = .false.
        end if

    end function test_assign_solver_f_c

    logical(c_bool) function test_assign_stability_f_c() bind(C)
        !
        ! this function tests that the function that converts stability check settings
        ! from C to Fortran correctly performs this conversion
        !
        use c_interface, only: stability_settings_type_c, assignment(=)
        use opentrustregion, only: stability_settings_type, default_stability_settings
        use test_reference, only: ref_stability_settings, &
                                  get_reference_stability_values, operator(/=), &
                                  callbacks_unset, unset_callbacks, callbacks_wrapped

        type(stability_settings_type_c) :: settings_c
        type(stability_settings_type) :: settings

        ! assume test passes
        test_assign_stability_f_c = .true.

        ! initialize the C settings with the reference values and callback functions
        call get_reference_stability_values(settings_c)

        ! convert to Fortran settings
        settings = settings_c

        ! check that every callback function is converted to its wrapper
        if (.not. callbacks_wrapped(settings)) then
            write(stderr, *) "test_assign_stability_f_c failed: Callback functions "// &
                "not converted to their wrappers."
            test_assign_stability_f_c = .false.
        end if

        ! check against reference values
        if (settings /= ref_stability_settings) then
            write(stderr, *) "test_assign_stability_f_c failed: Settings not "// &
                "converted correctly."
            test_assign_stability_f_c = .false.
        end if

        ! convert initialized C settings without callback functions and check that no
        ! callback functions are associated
        call unset_callbacks(settings_c)
        settings = settings_c
        if (.not. callbacks_unset(settings)) then
            write(stderr, *) "test_assign_stability_f_c failed: Function pointers "// &
                "associated for callback functions that were not provided."
            test_assign_stability_f_c = .false.
        end if

        ! convert C settings with callback functions that were not initialized, the
        ! custom values and callback functions they still carry have to be replaced by
        ! the default settings
        call get_reference_stability_values(settings_c)
        settings_c%initialized = .false.
        settings = settings_c
        if (.not. callbacks_unset(settings)) then
            write(stderr, *) "test_assign_stability_f_c failed: Function pointers "// &
                "converted for settings that were not initialized."
            test_assign_stability_f_c = .false.
        end if
        if (settings /= default_stability_settings) then
            write(stderr, *) "test_assign_stability_f_c failed: Settings that were "// &
                "not initialized not converted to default values."
            test_assign_stability_f_c = .false.
        end if

    end function test_assign_stability_f_c

    logical(c_bool) function test_assign_solver_c_f() bind(C)
        !
        ! this function tests that the function that converts solver settings from
        ! Fortran to C correctly performs this conversion
        !
        use opentrustregion, only: solver_settings_type
        use c_interface, only: solver_settings_type_c, assignment(=)
        use test_reference, only: ref_solver_settings, operator(/=), callbacks_unset

        type(solver_settings_type) :: settings
        type(solver_settings_type_c) :: settings_c
        integer(ip) :: i
        character(len=21), parameter :: logical_names(3) = &
            [character(len=21) :: "stability", "line_search", "max_precision_reached"]

        ! assume test passes
        test_assign_solver_c_f = .true.

        ! convert Fortran settings with the reference values to C settings
        settings = ref_solver_settings
        settings_c = settings

        ! check that no callback function pointers are associated
        if (.not. callbacks_unset(settings_c)) then
            write(stderr, *) "test_assign_solver_c_f failed: Callback function "// &
                "pointers associated."
            test_assign_solver_c_f = .false.
        end if

        ! check the converted C settings against the reference values
        if (settings_c /= ref_solver_settings) then
            write(stderr, *) "test_assign_solver_c_f failed: Settings not "// &
                "converted correctly."
            test_assign_solver_c_f = .false.
        end if

        ! convert again with only one of the logicals other than initialized set at a
        ! time, since the logicals cannot all differ from each other and from their
        ! default values, so that a logical that is swapped with another or not
        ! converted is detected
        do i = 1, size(logical_names)
            settings = ref_solver_settings
            settings%stability = logical_names(i) == "stability"
            settings%line_search = logical_names(i) == "line_search"
            settings%max_precision_reached = logical_names(i) == "max_precision_reached"
            settings_c = settings
            if (settings_c /= settings) then
                write(stderr, *) "test_assign_solver_c_f failed: Settings with "// &
                    "only "//trim(logical_names(i))//" set not converted correctly."
                test_assign_solver_c_f = .false.
            end if
        end do

    end function test_assign_solver_c_f

    logical(c_bool) function test_assign_stability_c_f() bind(C)
        !
        ! this function tests that the function that converts stability check settings
        ! from Fortran to C correctly performs this conversion
        !
        use opentrustregion, only: stability_settings_type
        use c_interface, only: stability_settings_type_c, assignment(=)
        use test_reference, only: ref_stability_settings, operator(/=), callbacks_unset

        type(stability_settings_type) :: settings
        type(stability_settings_type_c) :: settings_c

        ! assume test passes
        test_assign_stability_c_f = .true.

        ! convert Fortran settings with the reference values to C settings
        settings = ref_stability_settings
        settings_c = settings

        ! check that no callback function pointers are associated
        if (.not. callbacks_unset(settings_c)) then
            write(stderr, *) "test_assign_stability_c_f failed: Callback function "// &
                "pointers associated."
            test_assign_stability_c_f = .false.
        end if

        ! check the converted C settings against the reference values
        if (settings_c /= ref_stability_settings) then
            write(stderr, *) "test_assign_stability_c_f failed: Settings not "// &
                "converted correctly."
            test_assign_stability_c_f = .false.
        end if

    end function test_assign_stability_c_f

    logical(c_bool) function test_character_to_c() bind(C)
        !
        ! this function tests conversion of a Fortran character string to a C
        ! null-terminated character array
        !
        use c_interface, only: character_to_c
        use opentrustregion, only: kw_len

        character(len=*), parameter :: test_string = "test  "
        character(kind=c_char) :: char_c(kw_len + 1)
        integer :: n, i

        ! assume test passes
        test_character_to_c = .true.

        ! check that the array has the size of the keyword fields of the C settings,
        ! before it is assigned to an array of this size
        if (size(character_to_c(test_string)) /= kw_len + 1) then
            write(stderr, *) "test_character_to_c failed: Character array has "// &
                "wrong size."
            test_character_to_c = .false.
            return
        end if

        ! perform conversion
        char_c = character_to_c(test_string)

        ! check characters, trailing blanks are dropped
        n = len_trim(test_string)
        do i = 1, n
            if (char_c(i) /= test_string(i:i)) then
                write(stderr, *) "test_character_to_c failed: Character array "// &
                    "mismatch at character ", i, "."
                test_character_to_c = .false.
            end if
        end do

        ! check null terminator and that the remainder is filled with null characters
        if (any(char_c(n + 1:) /= c_null_char)) then
            write(stderr, *) "test_character_to_c failed: Character array not "// &
                "filled with null characters after the string."
            test_character_to_c = .false.
        end if

    end function test_character_to_c

    logical(c_bool) function test_character_from_c() bind(C)
        !
        ! this function tests conversion of a C null-terminated character array to a
        ! Fortran character string
        !
        use c_interface, only: character_from_c

        character(kind=c_char), parameter :: test_array(5) = &
            ["t", "e", "s", "t", c_null_char]
        character(len=:), allocatable :: char_f
        integer :: i

        ! assume test passes
        test_character_from_c = .true.

        ! perform conversion
        char_f = character_from_c(test_array)

        ! check equality
        if (len(char_f) /= size(test_array) - 1) then
            write(stderr, *) "test_character_from_c failed: Converted string has "// &
                "wrong length."
            test_character_from_c = .false.
            return
        end if

        ! check characters
        do i = 1, len(char_f)
            if (test_array(i) /= char_f(i:i)) then
                write(stderr, *) "test_character_from_c failed: String mismatch at "// &
                    "character ", i
                test_character_from_c = .false.
            end if
        end do

    end function test_character_from_c

end module c_interface_unit_tests
