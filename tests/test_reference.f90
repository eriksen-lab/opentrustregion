! Copyright (C) 2025- Jonas Greiner
!
! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.

module test_reference

    use opentrustregion, only: rp, ip, kw_len, stderr, solver_settings_type, &
                               stability_settings_type
    use c_interface, only: c_rp, c_ip
    use, intrinsic :: iso_c_binding, only: c_bool, c_char, c_null_char, c_funptr, &
                                           c_f_procpointer, c_associated, c_ptr, &
                                           c_null_ptr, c_loc, c_f_pointer

    implicit none

    ! tolerance for real comparisons
    real(rp), parameter :: tol = 1e-10_rp
    real(c_rp), parameter :: tol_c = real(tol, kind=c_rp)

    ! number of parameters
    integer(ip), parameter :: n_param = 3_ip
    integer(c_ip), bind(C, name="test_n_param") :: n_param_c = n_param

    ! reference settings, every field other than the logicals has a value distinct
    ! from every other field of the same kind and from its default value, so that a
    ! field read from or written to the wrong place or not converted at all is
    ! detected, the logicals, which cannot all differ from each other and from their
    ! default values, are told apart by tests that set them one at a time or flip
    ! them, the nested stability check settings of the reference solver settings are
    ! the reference stability check settings
    type(stability_settings_type), parameter :: ref_stability_settings = &
        stability_settings_type(precond=null(), project=null(), logger=null(), &
                                initialized=.true., conv_tol=5e-6_rp, &
                                n_random_trial_vectors=6, n_iter=60, &
                                jacobi_davidson_start=13, seed=34, verbose=2, &
                                n_hess_x=12, diag_solver="jacobi-davidson")
    type(solver_settings_type), parameter :: ref_solver_settings = &
        solver_settings_type(precond=null(), project=null(), conv_check=null(), &
                             logger=null(), stability=.true., line_search=.false., &
                             initialized=.true., max_precision_reached=.true., &
                             conv_tol=2e-3_rp, start_trust_radius=0.3_rp, &
                             global_red_factor=3e-2_rp, local_red_factor=4e-3_rp, &
                             n_random_trial_vectors=5, n_macro=300, n_micro=200, &
                             jacobi_davidson_start=10, seed=33, verbose=3, &
                             n_update_orbs=7, n_hess_x=11, subsystem_solver="tcg", &
                             stability_settings=ref_stability_settings)

    ! host context handed to the callback functions by the tests that exercise it,
    ! together with the bookkeeping those callback functions keep, so that a test can
    ! check the context reached all of them unchanged, and the error code the C mock
    ! callback functions report
    type :: host_context_type
        integer(ip) :: n_calls = 0
        logical :: logger_called = .false.
        integer(c_ip) :: mock_error = 0
    end type
    type(host_context_type), target :: host_context, stability_host_context

    logical :: host_context_armed = .false., host_context_wrong = .false., &
               host_context_missing = .false.

    ! checks and operations on all callback functions of settings, so that a new
    ! callback function only needs to be added here
    interface callbacks_unset
        module procedure solver_callbacks_unset
        module procedure stability_callbacks_unset
        module procedure solver_callbacks_unset_c
        module procedure stability_callbacks_unset_c
    end interface

    interface unset_callbacks
        module procedure unset_solver_callbacks_c
        module procedure unset_stability_callbacks_c
    end interface

    interface callbacks_wrapped
        module procedure solver_callbacks_wrapped
        module procedure stability_callbacks_wrapped
    end interface

    interface operator(==)
        module procedure equal_solver
        module procedure equal_stability
        module procedure equal_solver_c
        module procedure equal_stability_c
    end interface

    interface operator(/=)
        module procedure not_equal_solver
        module procedure not_equal_stability
        module procedure not_equal_solver_c
        module procedure not_equal_stability_c
    end interface

contains

    subroutine check_host_context(context)
        !
        ! this subroutine records the host context a callback function received so that
        ! a test can check it was handed the context it armed, unchanged, on every call
        !
        class(*), intent(in), pointer :: context

        ! a callback function reached without a context is only a failure while a test
        ! has armed one
        if (.not. associated(context)) then
            if (host_context_armed) host_context_missing = .true.
            return
        end if

        select type (context)
        class is (host_context_type)
            context%n_calls = context%n_calls + 1
        class default
            host_context_wrong = .true.
        end select

    end subroutine check_host_context

    subroutine check_host_context_c(context_c)
        !
        ! this subroutine records the host context a C callback function received so
        ! that a test can check it was handed the context it armed, unchanged, on every
        ! call
        !
        type(c_ptr), intent(in) :: context_c

        type(host_context_type), pointer :: context

        ! a callback function reached without a context is only a failure while a test
        ! has armed one
        if (.not. c_associated(context_c)) then
            if (host_context_armed) host_context_missing = .true.
            return
        end if

        if (c_associated(context_c, c_loc(host_context)) .or. &
            c_associated(context_c, c_loc(stability_host_context))) then
            call c_f_pointer(context_c, context)
            context%n_calls = context%n_calls + 1
        else
            host_context_wrong = .true.
        end if

    end subroutine check_host_context_c

    function host_context_error(context_c) result(error)
        !
        ! this function returns the error code a C callback function reports, which a
        ! test sets in the host context to check that the error is passed on, and no
        ! error for a missing or foreign context
        !
        type(c_ptr), intent(in) :: context_c
        integer(c_ip) :: error

        type(host_context_type), pointer :: context

        error = 0
        if (c_associated(context_c, c_loc(host_context)) .or. &
            c_associated(context_c, c_loc(stability_host_context))) then
            call c_f_pointer(context_c, context)
            error = context%mock_error
        end if

    end function host_context_error

    function ref_character_to_c(char_f) result(char_c)
        !
        ! this function converts a Fortran keyword to a C null-terminated character
        ! array of the size of the keyword fields of the C settings, whose remainder
        ! is filled with null characters, independently of the conversion routine
        ! under test
        !
        character(len=*), intent(in) :: char_f
        character(kind=c_char) :: char_c(kw_len + 1)

        integer(ip) :: i

        char_c = c_null_char
        do i = 1, len_trim(char_f)
            char_c(i) = char_f(i:i)
        end do

    end function ref_character_to_c

    function ref_character_from_c(char_c) result(char_f)
        !
        ! this function converts a C null-terminated character array to a Fortran
        ! character string, independently of the conversion routine under test
        !
        character(kind=c_char), intent(in) :: char_c(*)
        character(len=:), allocatable :: char_f

        integer(ip) :: i

        char_f = ""
        i = 1
        do while (char_c(i) /= c_null_char)
            char_f = char_f//char_c(i)
            i = i + 1
        end do

    end function ref_character_from_c

    subroutine reset_host_context()
        !
        ! this subroutine clears the bookkeeping of the callback functions and arms the
        ! host context for a test
        !
        host_context%n_calls = 0
        stability_host_context%n_calls = 0
        host_context%logger_called = .false.
        stability_host_context%logger_called = .false.
        host_context%mock_error = 0
        stability_host_context%mock_error = 0
        host_context_armed = .true.
        host_context_wrong = .false.
        host_context_missing = .false.

    end subroutine reset_host_context

    subroutine arm_host_context(settings, context)
        !
        ! this subroutine points a settings object at the host context the callback
        ! functions expect, by default the global one, and clears their bookkeeping
        !
        use opentrustregion, only: settings_type

        class(settings_type), intent(inout) :: settings
        class(host_context_type), intent(inout), target, optional :: context

        call reset_host_context()
        if (present(context)) then
            context%n_calls = 0
            settings%context => context
        else
            settings%context => host_context
        end if

    end subroutine arm_host_context

    subroutine arm_host_context_c(context_c)
        !
        ! this subroutine points a C callback bundle at the host context the callback
        ! functions expect and clears their bookkeeping
        !
        type(c_ptr), intent(out) :: context_c

        context_c = c_loc(host_context)
        call reset_host_context()

    end subroutine arm_host_context_c

    subroutine get_reference_solver_values(values_out, true_logical) bind(C)
        !
        ! this subroutine exports the reference solver settings as C settings whose
        ! fields are set one by one without the conversion routines under test, the
        ! callback function pointers and host contexts are set to distinct addresses,
        ! and if the name of a logical (prefixed by "stability_settings." for the
        ! nested settings) is given, which only the layout tests do, only this logical
        ! is set so that they can tell swapped logicals apart, otherwise the logicals
        ! take their reference values
        !
        use c_interface, only: solver_settings_type_c
        use, intrinsic :: iso_c_binding, only: c_intptr_t, c_null_funptr, c_null_ptr

        type(solver_settings_type_c), intent(out) :: values_out
        character(kind=c_char), intent(in), optional :: true_logical(*)

        character(len=:), allocatable :: name

        ! set callback function pointers and host context to distinct addresses
        values_out%precond = transfer(1_c_intptr_t, c_null_funptr)
        values_out%project = transfer(2_c_intptr_t, c_null_funptr)
        values_out%conv_check = transfer(3_c_intptr_t, c_null_funptr)
        values_out%logger = transfer(4_c_intptr_t, c_null_funptr)
        values_out%context = transfer(5_c_intptr_t, c_null_ptr)

        ! set logicals
        name = ""
        if (present(true_logical)) name = ref_character_from_c(true_logical)
        if (len(name) == 0) then
            values_out%stability = ref_solver_settings%stability
            values_out%line_search = ref_solver_settings%line_search
            values_out%initialized = ref_solver_settings%initialized
            values_out%max_precision_reached = ref_solver_settings%max_precision_reached
        else
            values_out%stability = name == "stability"
            values_out%line_search = name == "line_search"
            values_out%initialized = name == "initialized"
            values_out%max_precision_reached = name == "max_precision_reached"
        end if

        ! set reals
        values_out%conv_tol = ref_solver_settings%conv_tol
        values_out%start_trust_radius = ref_solver_settings%start_trust_radius
        values_out%global_red_factor = ref_solver_settings%global_red_factor
        values_out%local_red_factor = ref_solver_settings%local_red_factor

        ! set integers
        values_out%n_random_trial_vectors = ref_solver_settings%n_random_trial_vectors
        values_out%n_macro = ref_solver_settings%n_macro
        values_out%n_micro = ref_solver_settings%n_micro
        values_out%jacobi_davidson_start = ref_solver_settings%jacobi_davidson_start
        values_out%seed = ref_solver_settings%seed
        values_out%verbose = ref_solver_settings%verbose
        values_out%n_update_orbs = ref_solver_settings%n_update_orbs
        values_out%n_hess_x = ref_solver_settings%n_hess_x

        ! set keyword
        values_out%subsystem_solver = &
            ref_character_to_c(ref_solver_settings%subsystem_solver)

        ! set nested stability check settings
        call get_reference_stability_values(values_out%stability_settings)
        if (len(name) > 0) values_out%stability_settings%initialized = &
            name == "stability_settings.initialized"

    end subroutine get_reference_solver_values

    subroutine get_reference_stability_values(values_out)
        !
        ! this subroutine sets C stability check settings to the reference values field
        ! by field without the conversion routines under test, the callback function
        ! pointers and host context are set to distinct addresses which also differ
        ! from those of the solver settings
        !
        use c_interface, only: stability_settings_type_c
        use, intrinsic :: iso_c_binding, only: c_intptr_t, c_null_funptr, c_null_ptr

        type(stability_settings_type_c), intent(out) :: values_out

        ! set callback function pointers and host context to distinct addresses
        values_out%precond = transfer(6_c_intptr_t, c_null_funptr)
        values_out%project = transfer(7_c_intptr_t, c_null_funptr)
        values_out%logger = transfer(8_c_intptr_t, c_null_funptr)
        values_out%context = transfer(9_c_intptr_t, c_null_ptr)

        ! set logical, real and integers
        values_out%initialized = ref_stability_settings%initialized
        values_out%conv_tol = ref_stability_settings%conv_tol
        values_out%n_random_trial_vectors = &
            ref_stability_settings%n_random_trial_vectors
        values_out%n_iter = ref_stability_settings%n_iter
        values_out%jacobi_davidson_start = ref_stability_settings%jacobi_davidson_start
        values_out%seed = ref_stability_settings%seed
        values_out%verbose = ref_stability_settings%verbose
        values_out%n_hess_x = ref_stability_settings%n_hess_x

        ! set keyword
        values_out%diag_solver = ref_character_to_c(ref_stability_settings%diag_solver)

    end subroutine get_reference_stability_values

    subroutine reference_field(name_c, field_value, keyword_c) bind(C)
        !
        ! this subroutine returns the value of a field of the C solver settings
        ! exported by get_reference_solver_values by its name, prefixed by
        ! "stability_settings." for the nested settings, which also describe
        ! standalone stability check settings, a numeric, logical (as 1 or 0) or
        ! pointer (as its address) field in field_value and a keyword field as a
        ! null-terminated character array in keyword_c, so that the C and Python tests
        ! need no values of their own and compare what they read through their own
        ! declarations with what the Fortran declaration reads
        !
        use c_interface, only: solver_settings_type_c
        use, intrinsic :: iso_c_binding, only: c_intptr_t

        character(kind=c_char), intent(in) :: name_c(*)
        real(c_rp), intent(out) :: field_value
        character(kind=c_char), intent(out) :: keyword_c(kw_len + 1)

        type(solver_settings_type_c) :: values
        character(len=:), allocatable :: name

        ! export reference values and select field by name
        field_value = 0.0_c_rp
        keyword_c = c_null_char
        call get_reference_solver_values(values)
        name = ref_character_from_c(name_c)
        select case (name)
        case ("precond")
            field_value = address(values%precond)
        case ("project")
            field_value = address(values%project)
        case ("conv_check")
            field_value = address(values%conv_check)
        case ("logger")
            field_value = address(values%logger)
        case ("context")
            field_value = real(transfer(values%context, 0_c_intptr_t), kind=c_rp)
        case ("stability")
            field_value = merge(1.0_c_rp, 0.0_c_rp, values%stability)
        case ("line_search")
            field_value = merge(1.0_c_rp, 0.0_c_rp, values%line_search)
        case ("initialized")
            field_value = merge(1.0_c_rp, 0.0_c_rp, values%initialized)
        case ("max_precision_reached")
            field_value = merge(1.0_c_rp, 0.0_c_rp, values%max_precision_reached)
        case ("conv_tol")
            field_value = values%conv_tol
        case ("start_trust_radius")
            field_value = values%start_trust_radius
        case ("global_red_factor")
            field_value = values%global_red_factor
        case ("local_red_factor")
            field_value = values%local_red_factor
        case ("n_random_trial_vectors")
            field_value = values%n_random_trial_vectors
        case ("n_macro")
            field_value = values%n_macro
        case ("n_micro")
            field_value = values%n_micro
        case ("jacobi_davidson_start")
            field_value = values%jacobi_davidson_start
        case ("seed")
            field_value = values%seed
        case ("verbose")
            field_value = values%verbose
        case ("n_update_orbs")
            field_value = values%n_update_orbs
        case ("n_hess_x")
            field_value = values%n_hess_x
        case ("subsystem_solver")
            keyword_c = values%subsystem_solver
        case ("stability_settings.precond")
            field_value = address(values%stability_settings%precond)
        case ("stability_settings.project")
            field_value = address(values%stability_settings%project)
        case ("stability_settings.logger")
            field_value = address(values%stability_settings%logger)
        case ("stability_settings.context")
            field_value = real( &
                transfer(values%stability_settings%context, 0_c_intptr_t), kind=c_rp)
        case ("stability_settings.initialized")
            field_value = &
                merge(1.0_c_rp, 0.0_c_rp, values%stability_settings%initialized)
        case ("stability_settings.conv_tol")
            field_value = values%stability_settings%conv_tol
        case ("stability_settings.n_random_trial_vectors")
            field_value = values%stability_settings%n_random_trial_vectors
        case ("stability_settings.n_iter")
            field_value = values%stability_settings%n_iter
        case ("stability_settings.jacobi_davidson_start")
            field_value = values%stability_settings%jacobi_davidson_start
        case ("stability_settings.seed")
            field_value = values%stability_settings%seed
        case ("stability_settings.verbose")
            field_value = values%stability_settings%verbose
        case ("stability_settings.n_hess_x")
            field_value = values%stability_settings%n_hess_x
        case ("stability_settings.diag_solver")
            keyword_c = values%stability_settings%diag_solver
        case default
            write(stderr, *) "reference_field: unknown field "//name//"."
            field_value = huge(1.0_c_rp)
        end select

    contains

        real(c_rp) function address(funptr)
            !
            ! this function returns the address of a C function pointer
            !
            type(c_funptr), intent(in) :: funptr

            address = real(transfer(funptr, 0_c_intptr_t), kind=c_rp)

        end function address

    end subroutine reference_field

    subroutine unset_solver_callbacks_c(settings_c)
        !
        ! this subroutine unsets every callback function of C solver settings,
        ! including those of the nested stability check settings
        !
        use c_interface, only: solver_settings_type_c
        use, intrinsic :: iso_c_binding, only: c_null_funptr

        type(solver_settings_type_c), intent(inout) :: settings_c

        settings_c%precond = c_null_funptr
        settings_c%project = c_null_funptr
        settings_c%conv_check = c_null_funptr
        settings_c%logger = c_null_funptr
        call unset_stability_callbacks_c(settings_c%stability_settings)

    end subroutine unset_solver_callbacks_c

    subroutine unset_stability_callbacks_c(settings_c)
        !
        ! this subroutine unsets every callback function of C stability check settings
        !
        use c_interface, only: stability_settings_type_c
        use, intrinsic :: iso_c_binding, only: c_null_funptr

        type(stability_settings_type_c), intent(inout) :: settings_c

        settings_c%precond = c_null_funptr
        settings_c%project = c_null_funptr
        settings_c%logger = c_null_funptr

    end subroutine unset_stability_callbacks_c

    logical function host_context_reached(test_name, context)
        !
        ! this function checks that the callback functions all received the host
        ! context they were armed with, by default the global one, and disarms it again
        !
        character(len=*), intent(in) :: test_name
        class(host_context_type), intent(in), optional :: context

        integer(ip) :: n_calls

        ! assume test passes
        host_context_reached = .true.

        if (host_context_wrong) then
            host_context_reached = .false.
            write(stderr, *) "test_"//test_name//" failed: A callback function "// &
                "received a host context other than the one that was set."
        end if
        if (host_context_missing) then
            host_context_reached = .false.
            write(stderr, *) "test_"//test_name//" failed: A callback function was "// &
                "reached without the host context that was set."
        end if
        if (present(context)) then
            n_calls = context%n_calls
        else
            n_calls = host_context%n_calls
        end if
        if (n_calls == 0) then
            host_context_reached = .false.
            write(stderr, *) "test_"//test_name//" failed: No callback function "// &
                "received the host context that was set."
        end if

        host_context_armed = .false.

    end function host_context_reached

    logical function equal_solver(lhs, rhs)
        !
        ! this function overloads the comparison operator to compare solver settings to
        ! different solver settings
        !
        type(solver_settings_type), intent(in) :: lhs, rhs

        equal_solver = &
            (lhs%stability .eqv. rhs%stability) .and. &
            (lhs%line_search .eqv. rhs%line_search) .and. &
            (lhs%initialized .eqv. rhs%initialized) .and. &
            (lhs%max_precision_reached .eqv. rhs%max_precision_reached) .and. &
            abs(lhs%conv_tol - rhs%conv_tol) <= tol .and. &
            abs(lhs%start_trust_radius - rhs%start_trust_radius) <= tol .and. &
            abs(lhs%global_red_factor - rhs%global_red_factor) <= tol .and. &
            abs(lhs%local_red_factor - rhs%local_red_factor) <= tol .and. &
            lhs%n_random_trial_vectors == rhs%n_random_trial_vectors .and. &
            lhs%n_macro == rhs%n_macro .and. lhs%n_micro == rhs%n_micro .and. &
            lhs%jacobi_davidson_start == rhs%jacobi_davidson_start .and. &
            lhs%seed == rhs%seed .and. lhs%verbose == rhs%verbose .and. &
            lhs%n_update_orbs == rhs%n_update_orbs .and. &
            lhs%n_hess_x == rhs%n_hess_x .and. &
            lhs%subsystem_solver == rhs%subsystem_solver .and. &
            lhs%stability_settings == rhs%stability_settings

    end function equal_solver

    logical function not_equal_solver(lhs, rhs)
        !
        ! this function overloads the negated comparison operator to compare solver
        ! settings to different solver settings
        !
        type(solver_settings_type), intent(in) :: lhs, rhs

        not_equal_solver = .not. (lhs == rhs)

    end function not_equal_solver

    logical function equal_stability(lhs, rhs)
        !
        ! this function overloads the comparison operator to compare stability settings
        ! to different stability settings
        !
        type(stability_settings_type), intent(in) :: lhs, rhs

        equal_stability = &
            (lhs%initialized .eqv. rhs%initialized) .and. &
            abs(lhs%conv_tol - rhs%conv_tol) <= tol .and. &
            lhs%n_random_trial_vectors == rhs%n_random_trial_vectors .and. &
            lhs%n_iter == rhs%n_iter .and. &
            lhs%jacobi_davidson_start == rhs%jacobi_davidson_start .and. &
            lhs%seed == rhs%seed .and. lhs%verbose == rhs%verbose .and. &
            lhs%n_hess_x == rhs%n_hess_x .and. lhs%diag_solver == rhs%diag_solver

    end function equal_stability

    logical function not_equal_stability(lhs, rhs)
        !
        ! this function overloads the negated comparison operator to compare stability
        ! settings to different stability settings
        !
        type(stability_settings_type), intent(in) :: lhs, rhs

        not_equal_stability = .not. (lhs == rhs)

    end function not_equal_stability

    logical function equal_solver_c(lhs_c, rhs)
        !
        ! this function overloads the comparison operator to compare C solver settings
        ! to Fortran solver settings field by field, independently of the conversion
        ! routines under test
        !
        use c_interface, only: solver_settings_type_c

        type(solver_settings_type_c), intent(in) :: lhs_c
        type(solver_settings_type), intent(in) :: rhs

        equal_solver_c = &
            (lhs_c%stability .eqv. logical(rhs%stability, kind=c_bool)) .and. &
            (lhs_c%line_search .eqv. logical(rhs%line_search, kind=c_bool)) .and. &
            (lhs_c%initialized .eqv. logical(rhs%initialized, kind=c_bool)) .and. &
            (lhs_c%max_precision_reached .eqv. &
             logical(rhs%max_precision_reached, kind=c_bool)) .and. &
            abs(lhs_c%conv_tol - rhs%conv_tol) <= tol .and. &
            abs(lhs_c%start_trust_radius - rhs%start_trust_radius) <= tol .and. &
            abs(lhs_c%global_red_factor - rhs%global_red_factor) <= tol .and. &
            abs(lhs_c%local_red_factor - rhs%local_red_factor) <= tol .and. &
            lhs_c%n_random_trial_vectors == rhs%n_random_trial_vectors .and. &
            lhs_c%n_macro == rhs%n_macro .and. lhs_c%n_micro == rhs%n_micro .and. &
            lhs_c%jacobi_davidson_start == rhs%jacobi_davidson_start .and. &
            lhs_c%seed == rhs%seed .and. lhs_c%verbose == rhs%verbose .and. &
            lhs_c%n_update_orbs == rhs%n_update_orbs .and. &
            lhs_c%n_hess_x == rhs%n_hess_x .and. &
            ref_character_from_c(lhs_c%subsystem_solver) == &
            trim(rhs%subsystem_solver) .and. &
            lhs_c%stability_settings == rhs%stability_settings

    end function equal_solver_c

    logical function not_equal_solver_c(lhs_c, rhs)
        !
        ! this function overloads the negated comparison operator to compare C solver
        ! settings to Fortran solver settings
        !
        use c_interface, only: solver_settings_type_c

        type(solver_settings_type_c), intent(in) :: lhs_c
        type(solver_settings_type), intent(in) :: rhs

        not_equal_solver_c = .not. (lhs_c == rhs)

    end function not_equal_solver_c

    logical function equal_stability_c(lhs_c, rhs)
        !
        ! this function overloads the comparison operator to compare C stability
        ! settings to Fortran stability settings field by field, independently of the
        ! conversion routines under test
        !
        use c_interface, only: stability_settings_type_c

        type(stability_settings_type_c), intent(in) :: lhs_c
        type(stability_settings_type), intent(in) :: rhs

        equal_stability_c = &
            (lhs_c%initialized .eqv. logical(rhs%initialized, kind=c_bool)) .and. &
            abs(lhs_c%conv_tol - rhs%conv_tol) <= tol .and. &
            lhs_c%n_random_trial_vectors == rhs%n_random_trial_vectors .and. &
            lhs_c%n_iter == rhs%n_iter .and. &
            lhs_c%jacobi_davidson_start == rhs%jacobi_davidson_start .and. &
            lhs_c%seed == rhs%seed .and. lhs_c%verbose == rhs%verbose .and. &
            lhs_c%n_hess_x == rhs%n_hess_x .and. &
            ref_character_from_c(lhs_c%diag_solver) == trim(rhs%diag_solver)

    end function equal_stability_c

    logical function not_equal_stability_c(lhs_c, rhs)
        !
        ! this function overloads the negated comparison operator to compare C
        ! stability settings to Fortran stability settings
        !
        use c_interface, only: stability_settings_type_c

        type(stability_settings_type_c), intent(in) :: lhs_c
        type(stability_settings_type), intent(in) :: rhs

        not_equal_stability_c = .not. (lhs_c == rhs)

    end function not_equal_stability_c

    logical function solver_callbacks_unset(settings)
        !
        ! this function checks that no callback function of solver settings, including
        ! those of the nested stability check settings, is associated
        !
        type(solver_settings_type), intent(in) :: settings

        solver_callbacks_unset = .not. ( &
            associated(settings%precond) .or. associated(settings%project) .or. &
            associated(settings%conv_check) .or. associated(settings%logger)) .and. &
                                 stability_callbacks_unset(settings%stability_settings)

    end function solver_callbacks_unset

    logical function stability_callbacks_unset(settings)
        !
        ! this function checks that no callback function of stability check settings is
        ! associated
        !
        type(stability_settings_type), intent(in) :: settings

        stability_callbacks_unset = .not. (associated(settings%precond) .or. &
                                           associated(settings%project) .or. &
                                           associated(settings%logger))

    end function stability_callbacks_unset

    logical function solver_callbacks_unset_c(settings_c)
        !
        ! this function checks that no callback function of C solver settings, including
        ! those of the nested stability check settings, is associated
        !
        use c_interface, only: solver_settings_type_c

        type(solver_settings_type_c), intent(in) :: settings_c

        solver_callbacks_unset_c = &
            .not. (c_associated(settings_c%precond) .or. &
                   c_associated(settings_c%project) .or. &
                   c_associated(settings_c%conv_check) .or. &
                   c_associated(settings_c%logger)) .and. &
            stability_callbacks_unset_c(settings_c%stability_settings)

    end function solver_callbacks_unset_c

    logical function stability_callbacks_unset_c(settings_c)
        !
        ! this function checks that no callback function of C stability check settings
        ! is associated
        !
        use c_interface, only: stability_settings_type_c

        type(stability_settings_type_c), intent(in) :: settings_c

        stability_callbacks_unset_c = .not. (c_associated(settings_c%precond) .or. &
                                             c_associated(settings_c%project) .or. &
                                             c_associated(settings_c%logger))

    end function stability_callbacks_unset_c

    logical function solver_callbacks_wrapped(settings)
        !
        ! this function checks that every callback function of solver settings converted
        ! from C, including those of the nested stability check settings, is associated
        ! with its wrapper
        !
        use c_interface, only: precond_f_wrapper, project_f_wrapper, &
                               conv_check_f_wrapper, logger_f_wrapper

        type(solver_settings_type), intent(in) :: settings

        solver_callbacks_wrapped = &
            associated(settings%precond, precond_f_wrapper) .and. &
            associated(settings%project, project_f_wrapper) .and. &
            associated(settings%conv_check, conv_check_f_wrapper) .and. &
            associated(settings%logger, logger_f_wrapper) .and. &
            stability_callbacks_wrapped(settings%stability_settings)

    end function solver_callbacks_wrapped

    logical function stability_callbacks_wrapped(settings)
        !
        ! this function checks that every callback function of stability check
        ! settings converted from C is associated with its wrapper
        !
        use c_interface, only: precond_f_wrapper, project_f_wrapper, logger_f_wrapper

        type(stability_settings_type), intent(in) :: settings

        stability_callbacks_wrapped = &
            associated(settings%precond, precond_f_wrapper) .and. &
            associated(settings%project, project_f_wrapper) .and. &
            associated(settings%logger, logger_f_wrapper)

    end function stability_callbacks_wrapped

    logical(c_bool) function is_default_solver_settings(settings_c) bind(C)
        !
        ! this function checks whether C solver settings hold the default values
        ! without callback functions and host contexts, so that a C test can compare
        ! against them without values of its own
        !
        use opentrustregion, only: default_solver_settings
        use c_interface, only: solver_settings_type_c

        type(solver_settings_type_c), intent(in) :: settings_c

        is_default_solver_settings = &
            settings_c == default_solver_settings .and. &
            callbacks_unset(settings_c) .and. &
            .not. c_associated(settings_c%context) .and. &
            .not. c_associated(settings_c%stability_settings%context)

    end function is_default_solver_settings

    logical(c_bool) function is_default_stability_settings(settings_c) bind(C)
        !
        ! this function checks whether C stability check settings hold the default
        ! values without callback functions and host context, so that a C test can
        ! compare against them without values of its own
        !
        use opentrustregion, only: default_stability_settings
        use c_interface, only: stability_settings_type_c

        type(stability_settings_type_c), intent(in) :: settings_c

        is_default_stability_settings = settings_c == default_stability_settings .and. &
                                        callbacks_unset(settings_c) .and. &
                                        .not. c_associated(settings_c%context)

    end function is_default_stability_settings

    function check_update_orbs_funptr(update_orbs_funptr, test_name, message, context) &
        result(test_passed)
        !
        ! this function tests a provided orbital updating function pointer
        !
        use opentrustregion, only: update_orbs_type, hess_x_type

        procedure(update_orbs_type), intent(in), pointer :: update_orbs_funptr
        character(len=*), intent(in) :: test_name, message
        class(*), intent(in), pointer :: context
        logical :: test_passed

        real(rp), allocatable :: kappa(:), grad(:), h_diag(:)
        real(rp) :: func
        integer(ip) :: error
        procedure(hess_x_type), pointer :: hess_x_funptr

        ! assume tests pass
        test_passed = .true.

        ! check if function pointer is associated
        if (.not. associated(update_orbs_funptr)) then
            test_passed = .false.
            write(stderr, *) "test_"//test_name//" failed: Orbital updating "// &
                "function provided"//message//" not associated with value."
            return
        end if

        ! allocate arrays
        allocate(kappa(n_param), grad(n_param), h_diag(n_param))

        ! initialize orbital update
        kappa = 1.0_rp

        ! call orbital update
        call update_orbs_funptr(kappa, func, grad, h_diag, hess_x_funptr, error, &
                                context)

        ! check for error
        if (error /= 0) then
            test_passed = .false.
            write(stderr, *) "test_"//test_name//" failed: Error produced"//message//"."
        end if

        ! check objective function value
        if (abs(func - 3.0_rp) > tol) then
            test_passed = .false.
            write(stderr, *) "test_"//test_name// &
                " failed: Objective function value returned"//message//" wrong."
        end if

        ! check gradient
        if (any(abs(grad - 2.0_rp) > tol)) then
            test_passed = .false.
            write(stderr, *) "test_"//test_name//" failed: Gradient returned"// &
                message//" wrong."
        end if

        ! check Hessian diagonal
        if (any(abs(h_diag - 3.0_rp) > tol)) then
            test_passed = .false.
            write(stderr, *) "test_"//test_name// &
                " failed: Hessian diagonal returned"//message//" wrong."
        end if

        ! deallocate arrays
        deallocate(kappa, grad, h_diag)

        ! test returned Hessian linear transformation, the function pointer is only
        ! defined if the orbital update did not produce an error
        if (error == 0) then
            test_passed = test_passed .and. check_hess_x_funptr( &
                hess_x_funptr, test_name, &
                " by Hessian linear transformation function returned"//message, context)
        end if

    end function check_update_orbs_funptr

    function check_update_orbs_c_funptr(update_orbs_c_funptr, test_name, message, &
                                        context_c) result(test_passed)
        !
        ! this function tests a provided orbital updating C function pointer
        !
        use c_interface, only: update_orbs_c_type

        type(c_funptr), intent(in) :: update_orbs_c_funptr
        character(len=*), intent(in) :: test_name, message
        type(c_ptr), intent(in) :: context_c
        logical :: test_passed

        procedure(update_orbs_c_type), pointer :: update_orbs_funptr
        real(c_rp), allocatable :: kappa(:), grad(:), h_diag(:)
        real(c_rp) :: func
        integer(c_ip) :: error
        type(c_funptr) :: hess_x_c_funptr

        ! assume tests pass
        test_passed = .true.

        ! check if function pointer is associated
        if (.not. c_associated(update_orbs_c_funptr)) then
            test_passed = .false.
            write(stderr, *) "test_"//test_name//" failed: Orbital updating "// &
                "function provided"//message//" not associated with value."
            return
        end if

        ! convert to Fortran function pointer
        call c_f_procpointer(cptr=update_orbs_c_funptr, fptr=update_orbs_funptr)

        ! allocate arrays
        allocate(kappa(n_param), grad(n_param), h_diag(n_param))

        ! initialize orbital update
        kappa = 1.0_c_rp

        ! call orbital update
        error = &
            update_orbs_funptr(kappa, func, grad, h_diag, hess_x_c_funptr, context_c)

        ! check for error
        if (error /= 0) then
            test_passed = .false.
            write(stderr, *) "test_"//test_name//" failed: Error produced"//message//"."
        end if

        ! check objective function value
        if (abs(func - 3.0_c_rp) > tol_c) then
            test_passed = .false.
            write(stderr, *) "test_"//test_name// &
                " failed: Objective function value returned"//message//" wrong."
        end if

        ! check gradient
        if (any(abs(grad - 2.0_c_rp) > tol_c)) then
            test_passed = .false.
            write(stderr, *) "test_"//test_name//" failed: Gradient returned"// &
                message//" wrong."
        end if

        ! check Hessian diagonal
        if (any(abs(h_diag - 3.0_c_rp) > tol_c)) then
            test_passed = .false.
            write(stderr, *) "test_"//test_name// &
                " failed: Hessian diagonal returned"//message//" wrong."
        end if

        ! deallocate arrays
        deallocate(kappa, grad, h_diag)

        ! test returned Hessian linear transformation, the function pointer is only
        ! defined if the orbital update did not produce an error
        if (error == 0) then
            test_passed = test_passed .and. check_hess_x_c_funptr( &
                hess_x_c_funptr, test_name, " by Hessian linear transformation "// &
                "function returned"//message, context_c)
        end if

    end function check_update_orbs_c_funptr

    function check_hess_x_funptr(hess_x_funptr, test_name, message, context) &
        result(test_passed)
        !
        ! this function tests a provided Hessian linear transformation function pointer
        !
        use opentrustregion, only: hess_x_type

        procedure(hess_x_type), intent(in), pointer :: hess_x_funptr
        character(len=*), intent(in) :: test_name, message
        class(*), intent(in), pointer :: context
        logical :: test_passed

        real(rp), allocatable :: x(:), hess_x(:)
        integer(ip) :: error

        ! assume tests pass
        test_passed = .true.

        ! check if function pointer is associated
        if (.not. associated(hess_x_funptr)) then
            test_passed = .false.
            write(stderr, *) "test_"//test_name// &
                " failed: Hessian linear transformation function provided"//message// &
                " not associated with value."
            return
        end if

        ! allocate arrays
        allocate(x(n_param), hess_x(n_param))

        ! initialize trial vector
        x = 1.0_rp

        ! call Hessian linear transformation
        call hess_x_funptr(x, hess_x, error, context)

        ! check for error
        if (error /= 0) then
            test_passed = .false.
            write(stderr, *) "test_"//test_name//" failed: Error produced"//message//"."
        end if

        ! check Hessian linear transformation
        if (any(abs(hess_x - 4.0_rp) > tol)) then
            test_passed = .false.
            write(stderr, *) "test_"//test_name// &
                " failed: Hessian linear transformation returned"//message//" wrong."
        end if

        ! deallocate arrays
        deallocate(x, hess_x)

    end function check_hess_x_funptr

    function check_hess_x_c_funptr(hess_x_c_funptr, test_name, message, context_c) &
        result(test_passed)
        !
        ! this function tests a provided Hessian linear transformation C function
        ! pointer
        !
        use c_interface, only: hess_x_c_type

        type(c_funptr), intent(in) :: hess_x_c_funptr
        character(len=*), intent(in) :: test_name, message
        type(c_ptr), intent(in) :: context_c
        logical :: test_passed

        procedure(hess_x_c_type), pointer :: hess_x_funptr_c
        real(c_rp), allocatable :: x(:), hess_x(:)
        integer(c_ip) :: error

        ! assume tests pass
        test_passed = .true.

        ! check if function pointer is associated
        if (.not. c_associated(hess_x_c_funptr)) then
            test_passed = .false.
            write(stderr, *) "test_"//test_name// &
                " failed: Hessian linear transformation function provided"//message// &
                " not associated with value."
            return
        end if

        ! convert to Fortran function pointer
        call c_f_procpointer(cptr=hess_x_c_funptr, fptr=hess_x_funptr_c)

        ! allocate arrays
        allocate(x(n_param), hess_x(n_param))

        ! initialize trial vector
        x = 1.0_c_rp

        ! call Hessian linear transformation
        error = hess_x_funptr_c(x, hess_x, context_c)

        ! check for error
        if (error /= 0) then
            test_passed = .false.
            write(stderr, *) "test_"//test_name//" failed: Error produced"//message//"."
        end if

        ! check Hessian linear transformation
        if (any(abs(hess_x - 4.0_c_rp) > tol_c)) then
            test_passed = .false.
            write(stderr, *) "test_"//test_name// &
                " failed: Hessian linear transformation returned"//message//" wrong."
        end if

        ! deallocate arrays
        deallocate(x, hess_x)

    end function check_hess_x_c_funptr

    function check_obj_func_funptr(obj_func_funptr, test_name, message, context) &
        result(test_passed)
        !
        ! this function tests a provided objective function function pointer
        !
        use opentrustregion, only: obj_func_type

        procedure(obj_func_type), intent(in), pointer :: obj_func_funptr
        character(len=*), intent(in) :: test_name, message
        class(*), intent(in), pointer :: context
        logical :: test_passed

        real(rp), allocatable :: kappa(:)
        real(rp) :: func
        integer(ip) :: error

        ! assume tests pass
        test_passed = .true.

        ! check if function pointer is associated
        if (.not. associated(obj_func_funptr)) then
            test_passed = .false.
            write(stderr, *) "test_"//test_name//" failed: Objective function "// &
                "provided"//message//" not associated with value."
            return
        end if

        ! allocate arrays
        allocate(kappa(n_param))

        ! initialize orbital update
        kappa = 1.0_rp

        ! call objective function
        func = obj_func_funptr(kappa, error, context)

        ! check for error
        if (error /= 0) then
            write(stderr, *) "test_"//test_name//" failed: Error produced"//message//"."
            test_passed = .false.
        end if

        ! check objective function
        if (abs(func - 3.0_rp) > tol) then
            write(stderr, *) "test_"//test_name//" failed: Function value returned"// &
                message//" wrong."
            test_passed = .false.
        end if

        ! deallocate arrays
        deallocate(kappa)

    end function check_obj_func_funptr

    function check_obj_func_c_funptr(obj_func_c_funptr, test_name, message, context_c) &
        result(test_passed)
        !
        ! this function tests a provided objective function C function pointer
        !
        use c_interface, only: obj_func_c_type

        type(c_funptr), intent(in) :: obj_func_c_funptr
        character(len=*), intent(in) :: test_name, message
        type(c_ptr), intent(in) :: context_c
        logical :: test_passed

        procedure(obj_func_c_type), pointer :: obj_func_funptr
        real(c_rp), allocatable :: kappa(:)
        real(c_rp) :: func
        integer(c_ip) :: error

        ! assume tests pass
        test_passed = .true.

        ! check if function pointer is associated
        if (.not. c_associated(obj_func_c_funptr)) then
            test_passed = .false.
            write(stderr, *) "test_"//test_name//" failed: Objective function "// &
                "provided"//message//" not associated with value."
            return
        end if

        ! convert to Fortran function pointer
        call c_f_procpointer(cptr=obj_func_c_funptr, fptr=obj_func_funptr)

        ! allocate arrays
        allocate(kappa(n_param))

        ! initialize orbital update
        kappa = 1.0_c_rp

        ! call objective function
        error = obj_func_funptr(kappa, func, context_c)

        ! check for error
        if (error /= 0) then
            write(stderr, *) "test_"//test_name//" failed: Error produced"//message//"."
            test_passed = .false.
        end if

        ! check objective function
        if (abs(func - 3.0_c_rp) > tol_c) then
            write(stderr, *) "test_"//test_name//" failed: Function value returned"// &
                message//" wrong."
            test_passed = .false.
        end if

        ! deallocate arrays
        deallocate(kappa)

    end function check_obj_func_c_funptr

    function check_precond_funptr(precond_funptr, test_name, message, context) &
        result(test_passed)
        !
        ! this function tests a provided preconditioner function pointer
        !
        use opentrustregion, only: precond_type

        procedure(precond_type), intent(in), pointer :: precond_funptr
        character(len=*), intent(in) :: test_name, message
        class(*), intent(in), pointer :: context
        logical :: test_passed

        real(rp), allocatable :: residual(:), precond_residual(:)
        integer(ip) :: error

        ! assume tests pass
        test_passed = .true.

        ! check if function pointer is associated
        if (.not. associated(precond_funptr)) then
            test_passed = .false.
            write(stderr, *) "test_"//test_name//" failed: Preconditioner function "// &
                "provided"//message//" not associated with value."
            return
        end if

        ! allocate arrays
        allocate(residual(n_param), precond_residual(n_param))

        ! initialize residual
        residual = 1.0_rp

        ! call preconditioning subroutine
        call precond_funptr(residual, 5.0_rp, precond_residual, error, context)

        ! check for error
        if (error /= 0) then
            write(stderr, *) "test_"//test_name//" failed: Error produced"//message//"."
            test_passed = .false.
        end if

        ! check preconditioned residual
        if (any(abs(precond_residual - 5.0_rp) > tol)) then
            write(stderr, *) "test_"//test_name// &
                " failed: Preconditioned residual returned"//message//" wrong."
            test_passed = .false.
        end if

        ! deallocate arrays
        deallocate(residual, precond_residual)

    end function check_precond_funptr

    function check_precond_c_funptr(precond_c_funptr, test_name, message, context_c) &
        result(test_passed)
        !
        ! this function tests a provided preconditioner C function pointer
        !
        use c_interface, only: precond_c_type

        type(c_funptr), intent(in) :: precond_c_funptr
        character(len=*), intent(in) :: test_name, message
        type(c_ptr), intent(in) :: context_c
        logical :: test_passed

        procedure(precond_c_type), pointer :: precond_funptr
        real(c_rp), allocatable :: residual(:), precond_residual(:)
        integer(c_ip) :: error

        ! assume tests pass
        test_passed = .true.

        ! check if function pointer is associated
        if (.not. c_associated(precond_c_funptr)) then
            test_passed = .false.
            write(stderr, *) "test_"//test_name//" failed: Preconditioner function "// &
                "provided"//message//" not associated with value."
            return
        end if

        ! convert to Fortran function pointer
        call c_f_procpointer(cptr=precond_c_funptr, fptr=precond_funptr)

        ! allocate arrays
        allocate(residual(n_param), precond_residual(n_param))

        ! initialize residual
        residual = 1.0_c_rp

        ! call preconditioning function
        error = precond_funptr(residual, 5.0_c_rp, precond_residual, context_c)

        ! check for error
        if (error /= 0) then
            write(stderr, *) "test_"//test_name//" failed: Error produced"//message//"."
            test_passed = .false.
        end if

        ! check preconditioned residual
        if (any(abs(precond_residual - 5.0_c_rp) > tol_c)) then
            write(stderr, *) "test_"//test_name// &
                " failed: Preconditioned residual returned"//message//" wrong."
            test_passed = .false.
        end if

        ! deallocate arrays
        deallocate(residual, precond_residual)

    end function check_precond_c_funptr

    function check_project_funptr(project_funptr, test_name, message, context) &
        result(test_passed)
        !
        ! this function tests a provided projection function pointer
        !
        use opentrustregion, only: project_type

        procedure(project_type), intent(in), pointer :: project_funptr
        character(len=*), intent(in) :: test_name, message
        class(*), intent(in), pointer :: context
        logical :: test_passed

        real(rp), allocatable :: vector(:)
        integer(ip) :: error

        ! assume tests pass
        test_passed = .true.

        ! check if function pointer is associated
        if (.not. associated(project_funptr)) then
            test_passed = .false.
            write(stderr, *) "test_"//test_name//" failed: Project function "// &
                "provided"//message//" not associated with value."
            return
        end if

        ! allocate arrays
        allocate(vector(n_param))

        ! initialize vector
        vector = 1.0_rp

        ! call projection subroutine
        call project_funptr(vector, error, context)

        ! check for error
        if (error /= 0) then
            write(stderr, *) "test_"//test_name//" failed: Error produced"//message//"."
            test_passed = .false.
        end if

        ! check projected vector
        if (any(abs(vector - 2.0_rp) > tol)) then
            write(stderr, *) "test_"//test_name// &
                " failed: Projected vector returned"//message//" wrong."
            test_passed = .false.
        end if

        ! deallocate arrays
        deallocate(vector)

    end function check_project_funptr

    function check_project_c_funptr(project_c_funptr, test_name, message, context_c) &
        result(test_passed)
        !
        ! this function tests a provided projection C function pointer
        !
        use c_interface, only: project_c_type

        type(c_funptr), intent(in) :: project_c_funptr
        character(len=*), intent(in) :: test_name, message
        type(c_ptr), intent(in) :: context_c
        logical :: test_passed

        procedure(project_c_type), pointer :: project_funptr
        real(c_rp), allocatable :: vector(:)
        integer(c_ip) :: error

        ! assume tests pass
        test_passed = .true.

        ! check if function pointer is associated
        if (.not. c_associated(project_c_funptr)) then
            test_passed = .false.
            write(stderr, *) "test_"//test_name//" failed: Project function "// &
                "provided"//message//" not associated with value."
            return
        end if

        ! convert to Fortran function pointer
        call c_f_procpointer(cptr=project_c_funptr, fptr=project_funptr)

        ! allocate arrays
        allocate(vector(n_param))

        ! initialize vector
        vector = 1.0_c_rp

        ! call projection function
        error = project_funptr(vector, context_c)

        ! check for error
        if (error /= 0) then
            write(stderr, *) "test_"//test_name//" failed: Error produced"//message//"."
            test_passed = .false.
        end if

        ! check projected vector
        if (any(abs(vector - 2.0_c_rp) > tol_c)) then
            write(stderr, *) "test_"//test_name// &
                " failed: Projected vector returned"//message//" wrong."
            test_passed = .false.
        end if

        ! deallocate arrays
        deallocate(vector)

    end function check_project_c_funptr

    function check_conv_check_funptr(conv_check_funptr, test_name, message, context) &
        result(test_passed)
        !
        ! this function tests a provided convergence check function pointer
        !
        use opentrustregion, only: conv_check_type

        procedure(conv_check_type), intent(in), pointer :: conv_check_funptr
        character(len=*), intent(in) :: test_name, message
        class(*), intent(in), pointer :: context
        logical :: test_passed

        logical :: converged
        integer(ip) :: error

        ! assume tests pass
        test_passed = .true.

        ! check if function pointer is associated
        if (.not. associated(conv_check_funptr)) then
            test_passed = .false.
            write(stderr, *) "test_"//test_name//" failed: Convergence check "// &
                "function provided"//message//" not associated with value."
            return
        end if

        ! call convergence check function
        converged = conv_check_funptr(error, context)

        ! check for error
        if (error /= 0) then
            write(stderr, *) "test_"//test_name//" failed: Error produced"//message//"."
            test_passed = .false.
        end if

        ! check convergence logical
        if (.not. converged) then
            write(stderr, *) "test_"//test_name// &
                " failed: Convergence logical returned"//message//" wrong."
            test_passed = .false.
        end if

    end function check_conv_check_funptr

    function check_conv_check_c_funptr(conv_check_c_funptr, test_name, message, &
                                       context_c) result(test_passed)
        !
        ! this function tests a provided convergence check C function pointer
        !
        use c_interface, only: conv_check_c_type

        type(c_funptr), intent(in) :: conv_check_c_funptr
        character(len=*), intent(in) :: test_name, message
        type(c_ptr), intent(in) :: context_c
        logical :: test_passed

        procedure(conv_check_c_type), pointer :: conv_check_funptr
        logical(c_bool) :: converged
        integer(ip) :: error

        ! assume tests pass
        test_passed = .true.

        ! check if function pointer is associated
        if (.not. c_associated(conv_check_c_funptr)) then
            test_passed = .false.
            write(stderr, *) "test_"//test_name//" failed: Convergence check "// &
                "function provided"//message//" not associated with value."
            return
        end if

        ! convert to Fortran function pointer
        call c_f_procpointer(cptr=conv_check_c_funptr, fptr=conv_check_funptr)

        ! call convergence check function
        error = conv_check_funptr(converged, context_c)

        ! check for error
        if (error /= 0) then
            write(stderr, *) "test_"//test_name//" failed: Error produced"//message//"."
            test_passed = .false.
        end if

        ! check convergence logical
        if (.not. converged) then
            write(stderr, *) "test_"//test_name// &
                " failed: Convergence logical returned"//message//" wrong."
            test_passed = .false.
        end if

    end function check_conv_check_c_funptr

end module test_reference
