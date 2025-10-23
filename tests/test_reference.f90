! Copyright (C) 2025- Jonas Greiner
!
! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.

module test_reference

    use opentrustregion, only: rp, ip
    use c_interface, only: c_rp, c_ip
    use, intrinsic :: iso_c_binding, only: c_bool, c_null_funptr

    implicit none

    ! tolerance for real comparisons
    real(rp), parameter :: tol = 1e-10_rp
    real(c_rp), parameter :: tol_c = real(tol, kind=c_rp)

    ! derived types for solver settings
    type ref_settings_type
        logical :: stability, line_search, davidson, jacobi_davidson, &
                   prefer_jacobi_davidson
        real(rp) :: conv_tol, start_trust_radius, global_red_factor, local_red_factor
        integer(ip) :: n_random_trial_vectors, n_macro, n_micro, seed, verbose, n_iter
    end type

    type, bind(C) :: ref_settings_type_c
        logical(c_bool) :: stability, line_search, davidson, jacobi_davidson, &
                           prefer_jacobi_davidson
        real(c_rp) :: conv_tol, start_trust_radius, global_red_factor, local_red_factor
        integer(c_ip) :: n_random_trial_vectors, n_macro, n_micro, seed, verbose, n_iter
    end type

    ! general reference parameters
    logical, parameter :: stability = .false., line_search = .true., &
                          davidson = .false., jacobi_davidson = .false., &
                          prefer_jacobi_davidson = .true.
    real(rp), parameter :: conv_tol = 1e-3_rp, start_trust_radius = 0.2_rp, &
                           global_red_factor = 1e-2_rp, local_red_factor = 1e-3_rp
    integer(ip), parameter :: n_random_trial_vectors = 5, n_macro = 300, &
                              n_micro = 200, seed = 33, verbose = 3, n_iter = 50

    type(ref_settings_type) :: ref_settings = &
        ref_settings_type(stability = stability, line_search = line_search, &
                          davidson = davidson, jacobi_davidson = jacobi_davidson, &
                          prefer_jacobi_davidson = prefer_jacobi_davidson, &
                          conv_tol = conv_tol, &
                          start_trust_radius = start_trust_radius, &
                          global_red_factor = global_red_factor, &
                          local_red_factor = local_red_factor, &
                          n_random_trial_vectors = n_random_trial_vectors, &
                          n_macro = n_macro, n_micro = n_micro, seed = seed, &
                          verbose = verbose, n_iter = n_iter)

    type(ref_settings_type_c) :: ref_settings_c = &
        ref_settings_type_c(stability = logical(stability), &
                            line_search = logical(line_search), &
                            davidson = logical(davidson), &
                            jacobi_davidson = logical(jacobi_davidson), &
                            prefer_jacobi_davidson = logical(prefer_jacobi_davidson), &
                            conv_tol = real(conv_tol, kind=c_rp), &
                            start_trust_radius = real(start_trust_radius, kind=c_rp), &
                            global_red_factor = real(global_red_factor, kind=c_rp), &
                            local_red_factor = real(local_red_factor, kind=c_rp), &
                            n_random_trial_vectors = int(n_random_trial_vectors, &
                                                         kind=c_ip), &
                            n_macro = int(n_macro, kind=c_ip), &
                            n_micro = int(n_micro, kind=c_ip), &
                            seed = int(seed, kind=c_ip), &
                            verbose = int(verbose, kind=c_ip), &
                            n_iter = int(n_iter, kind=c_ip))

    interface assignment(=)
        module procedure assign_ref_to_solver_c
        module procedure assign_ref_to_stability_c
    end interface

    interface operator(==)
        module procedure equal_solver_to_ref
        module procedure equal_stability_to_ref
        module procedure equal_solver_to_ref_c
        module procedure equal_stability_to_ref_c
        module procedure equal_solver
        module procedure equal_stability
        module procedure equal_solver_c
        module procedure equal_stability_c
    end interface

    interface operator(/=)
        module procedure not_equal_solver_to_ref
        module procedure not_equal_stability_to_ref
        module procedure not_equal_solver_to_ref_c
        module procedure not_equal_stability_to_ref_c
        module procedure not_equal_solver
        module procedure not_equal_stability
        module procedure not_equal_solver_c
        module procedure not_equal_stability_c
    end interface

contains

    subroutine get_reference_values(ref_settings_out) bind(C)
        !
        ! this subroutine exports the reference values for tests
        !
        type(ref_settings_type_c), intent(out) :: ref_settings_out

        ref_settings_out = ref_settings_c

    end subroutine get_reference_values

    subroutine assign_ref_to_solver_c(lhs, rhs)
        !
        ! this subroutine overloads the assignment operator to set solver settings to 
        ! reference values
        !
        use c_interface, only: solver_settings_type_c

        type(solver_settings_type_c), intent(out) :: lhs
        type(ref_settings_type_c), intent(in) :: rhs

        ! unassociate function pointers
        lhs%precond = c_null_funptr
        lhs%conv_check = c_null_funptr
        lhs%logger = c_null_funptr

        ! set reference values
        lhs%stability = rhs%stability
        lhs%line_search = rhs%line_search
        lhs%davidson = rhs%davidson
        lhs%jacobi_davidson = rhs%jacobi_davidson
        lhs%prefer_jacobi_davidson = rhs%prefer_jacobi_davidson
        lhs%conv_tol = rhs%conv_tol
        lhs%start_trust_radius = rhs%start_trust_radius
        lhs%global_red_factor = rhs%global_red_factor
        lhs%local_red_factor  = rhs%local_red_factor
        lhs%n_random_trial_vectors = rhs%n_random_trial_vectors
        lhs%n_macro = rhs%n_macro
        lhs%n_micro = rhs%n_micro
        lhs%seed = rhs%seed
        lhs%verbose = rhs%verbose

        ! set initialization logical
        lhs%initialized = .true._c_bool

    end subroutine assign_ref_to_solver_c

    subroutine assign_ref_to_stability_c(lhs, rhs)
        !
        ! this subroutine overloads the assignment operator to set stability settings 
        ! to reference values
        !
        use c_interface, only: stability_settings_type_c

        type(stability_settings_type_c), intent(out) :: lhs
        type(ref_settings_type_c), intent(in) :: rhs

        ! unassociate function pointers
        lhs%precond = c_null_funptr
        lhs%logger = c_null_funptr

        ! set reference values
        lhs%jacobi_davidson = rhs%jacobi_davidson
        lhs%conv_tol = rhs%conv_tol
        lhs%n_random_trial_vectors = rhs%n_random_trial_vectors
        lhs%n_iter = rhs%n_iter
        lhs%verbose = rhs%verbose

        ! set initialization logical
        lhs%initialized = .true._c_bool

    end subroutine assign_ref_to_stability_c

    logical function equal_solver_to_ref(lhs, rhs)
        !
        ! this function overloads the comparison operator to compare solver settings 
        ! to reference values
        !
        use opentrustregion, only: solver_settings_type

        type(solver_settings_type), intent(in) :: lhs
        type(ref_settings_type), intent(in) :: rhs

        equal_solver_to_ref = (lhs%stability .eqv. rhs%stability) .and. &
            (lhs%line_search .eqv. rhs%line_search) .and. &
            (lhs%davidson .eqv. rhs%davidson) .and. &
            (lhs%jacobi_davidson .eqv. rhs%jacobi_davidson) .and. &
            (lhs%prefer_jacobi_davidson .eqv. rhs%prefer_jacobi_davidson) .and. &
            abs(lhs%conv_tol - rhs%conv_tol) <= tol .and. &
            abs(lhs%start_trust_radius - rhs%start_trust_radius) <= tol .and. &
            abs(lhs%global_red_factor - rhs%global_red_factor) <= tol .and. &
            abs(lhs%local_red_factor - rhs%local_red_factor) <= tol .and. &
            lhs%n_random_trial_vectors == rhs%n_random_trial_vectors .and. &
            lhs%n_macro == rhs%n_macro .and. lhs%n_micro == rhs%n_micro .and. &
            lhs%seed == rhs%seed .and. lhs%verbose == rhs%verbose

    end function equal_solver_to_ref

    logical function not_equal_solver_to_ref(lhs, rhs)
        !
        ! this function overloads the negated comparison operator to compare solver 
        ! settings to reference values
        !
        use opentrustregion, only: solver_settings_type

        type(solver_settings_type), intent(in) :: lhs
        type(ref_settings_type), intent(in) :: rhs

        not_equal_solver_to_ref = .not. (lhs == rhs)

    end function not_equal_solver_to_ref

    logical function equal_stability_to_ref(lhs, rhs)
        !
        ! this function overloads the comparison operator to compare stability settings 
        ! to reference values
        !
        use opentrustregion, only: stability_settings_type

        type(stability_settings_type), intent(in) :: lhs
        type(ref_settings_type), intent(in) :: rhs

        equal_stability_to_ref = (lhs%jacobi_davidson .eqv. rhs%jacobi_davidson) .and. &
            abs(lhs%conv_tol - rhs%conv_tol) <= tol .and. &
            lhs%n_random_trial_vectors == rhs%n_random_trial_vectors .and. &
            lhs%n_iter == rhs%n_iter .and. lhs%verbose == rhs%verbose

    end function equal_stability_to_ref

    logical function not_equal_stability_to_ref(lhs, rhs)
        !
        ! this function overloads the negated comparison operator to compare stability
        ! settings to reference values
        !
        use opentrustregion, only: stability_settings_type

        type(stability_settings_type), intent(in) :: lhs
        type(ref_settings_type), intent(in) :: rhs
        
        not_equal_stability_to_ref = .not. (lhs == rhs)

    end function not_equal_stability_to_ref

    logical function equal_solver_to_ref_c(lhs, rhs)
        !
        ! this function overloads the comparison operator to compare solver settings 
        ! to reference values
        !
        use c_interface, only: solver_settings_type_c

        type(solver_settings_type_c), intent(in) :: lhs
        type(ref_settings_type_c), intent(in) :: rhs
        
        equal_solver_to_ref_c = (lhs%stability .eqv. rhs%stability) .and. &
            (lhs%line_search .eqv. rhs%line_search) .and. &
            (lhs%davidson .eqv. rhs%davidson) .and. &
            (lhs%jacobi_davidson .eqv. rhs%jacobi_davidson) .and. &
            (lhs%prefer_jacobi_davidson .eqv. rhs%prefer_jacobi_davidson) .and. &
            abs(lhs%conv_tol - rhs%conv_tol) <= tol .and. &
            abs(lhs%start_trust_radius - rhs%start_trust_radius) <= tol .and. &
            abs(lhs%global_red_factor - rhs%global_red_factor) <= tol .and. &
            abs(lhs%local_red_factor - rhs%local_red_factor) <= tol .and. &
            lhs%n_random_trial_vectors == rhs%n_random_trial_vectors .and. &
            lhs%n_macro == rhs%n_macro .and. lhs%n_micro == rhs%n_micro .and. &
            lhs%seed == rhs%seed .and. lhs%verbose == rhs%verbose

    end function equal_solver_to_ref_c

    logical function not_equal_solver_to_ref_c(lhs, rhs)
        !
        ! this function overloads the negated comparison operator to compare solver 
        ! settings to reference values
        !
        use c_interface, only: solver_settings_type_c

        type(solver_settings_type_c), intent(in) :: lhs
        type(ref_settings_type_c), intent(in) :: rhs
        
        not_equal_solver_to_ref_c = .not. (lhs == rhs)

    end function not_equal_solver_to_ref_c

    logical function equal_stability_to_ref_c(lhs, rhs)
        !
        ! this function overloads the comparison operator to compare stability settings 
        ! to reference values
        !
        use c_interface, only: stability_settings_type_c

        type(stability_settings_type_c), intent(in) :: lhs
        type(ref_settings_type_c), intent(in) :: rhs
        
        equal_stability_to_ref_c = (lhs%jacobi_davidson .eqv. rhs%jacobi_davidson) &
            .and. abs(lhs%conv_tol - rhs%conv_tol) <= tol .and. &
            lhs%n_random_trial_vectors == rhs%n_random_trial_vectors .and. &
            lhs%n_iter == rhs%n_iter .and. lhs%verbose == rhs%verbose

    end function equal_stability_to_ref_c

    logical function not_equal_stability_to_ref_c(lhs, rhs)
        !
        ! this function overloads the negated comparison operator to compare stability
        ! settings to reference values
        !
        use c_interface, only: stability_settings_type_c

        type(stability_settings_type_c), intent(in) :: lhs
        type(ref_settings_type_c), intent(in) :: rhs
        
        not_equal_stability_to_ref_c = .not. (lhs == rhs)

    end function not_equal_stability_to_ref_c

    logical function equal_solver(lhs, rhs)
        !
        ! this function overloads the comparison operator to compare solver settings 
        ! to different solver settings
        !
        use opentrustregion, only: solver_settings_type

        type(solver_settings_type), intent(in) :: lhs, rhs
        
        equal_solver = (lhs%stability .eqv. rhs%stability) .and. &
            (lhs%line_search .eqv. rhs%line_search) .and. &
            (lhs%davidson .eqv. rhs%davidson) .and. &
            (lhs%jacobi_davidson .eqv. rhs%jacobi_davidson) .and. &
            (lhs%prefer_jacobi_davidson .eqv. rhs%prefer_jacobi_davidson) .and. &
            (lhs%initialized .eqv. rhs%initialized) .and. &
            abs(lhs%conv_tol - rhs%conv_tol) <= tol .and. &
            abs(lhs%start_trust_radius - rhs%start_trust_radius) <= tol .and. &
            abs(lhs%global_red_factor - rhs%global_red_factor) <= tol .and. &
            abs(lhs%local_red_factor - rhs%local_red_factor) <= tol .and. &
            lhs%n_random_trial_vectors == rhs%n_random_trial_vectors .and. &
            lhs%n_macro == rhs%n_macro .and. lhs%n_micro == rhs%n_micro .and. &
            lhs%seed == rhs%seed .and. lhs%verbose == rhs%verbose

    end function equal_solver

    logical function not_equal_solver(lhs, rhs)
        !
        ! this function overloads the negated comparison operator to compare solver 
        ! settings to different solver settings
        !
        use opentrustregion, only: solver_settings_type

        type(solver_settings_type), intent(in) :: lhs, rhs
        
        not_equal_solver = .not. (lhs == rhs)

    end function not_equal_solver

    logical function equal_stability(lhs, rhs)
        !
        ! this function overloads the comparison operator to compare stability settings 
        ! to different stability settings
        !
        use opentrustregion, only: stability_settings_type

        type(stability_settings_type), intent(in) :: lhs, rhs
        
        equal_stability = (lhs%jacobi_davidson .eqv. rhs%jacobi_davidson) .and. &
            abs(lhs%conv_tol - rhs%conv_tol) <= tol .and. &
            lhs%n_random_trial_vectors == rhs%n_random_trial_vectors .and. &
            lhs%n_iter == rhs%n_iter .and. lhs%verbose == rhs%verbose

    end function equal_stability

    logical function not_equal_stability(lhs, rhs)
        !
        ! this function overloads the negated comparison operator to compare stability
        ! settings to different stability settings
        !
        use opentrustregion, only: stability_settings_type

        type(stability_settings_type), intent(in) :: lhs, rhs
        
        not_equal_stability = .not. (lhs == rhs)

    end function not_equal_stability

    logical function equal_solver_c(lhs, rhs)
        !
        ! this function overloads the comparison operator to compare solver settings 
        ! to different solver settings
        !
        use c_interface, only: solver_settings_type_c

        type(solver_settings_type_c), intent(in) :: lhs, rhs
        
        equal_solver_c = (lhs%stability .eqv. rhs%stability) .and. &
            (lhs%line_search .eqv. rhs%line_search) .and. &
            (lhs%davidson .eqv. rhs%davidson) .and. &
            (lhs%jacobi_davidson .eqv. rhs%jacobi_davidson) .and. &
            (lhs%prefer_jacobi_davidson .eqv. rhs%prefer_jacobi_davidson) .and. &
            (lhs%initialized .eqv. rhs%initialized) .and. &
            abs(lhs%conv_tol - rhs%conv_tol) <= tol .and. &
            abs(lhs%start_trust_radius - rhs%start_trust_radius) <= tol .and. &
            abs(lhs%global_red_factor - rhs%global_red_factor) <= tol .and. &
            abs(lhs%local_red_factor - rhs%local_red_factor) <= tol .and. &
            lhs%n_random_trial_vectors == rhs%n_random_trial_vectors .and. &
            lhs%n_macro == rhs%n_macro .and. lhs%n_micro == rhs%n_micro .and. &
            lhs%seed == rhs%seed .and. lhs%verbose == rhs%verbose

    end function equal_solver_c

    logical function not_equal_solver_c(lhs, rhs)
        !
        ! this function overloads the negated comparison operator to compare solver 
        ! settings to different solver settings
        !
        use c_interface, only: solver_settings_type_c

        type(solver_settings_type_c), intent(in) :: lhs, rhs
        
        not_equal_solver_c = .not. (lhs == rhs)

    end function not_equal_solver_c

    logical function equal_stability_c(lhs, rhs)
        !
        ! this function overloads the comparison operator to compare stability settings 
        ! to different stability settings
        !
        use c_interface, only: stability_settings_type_c

        type(stability_settings_type_c), intent(in) :: lhs, rhs
        
        equal_stability_c = (lhs%jacobi_davidson .eqv. rhs%jacobi_davidson) .and. &
            abs(lhs%conv_tol - rhs%conv_tol) <= tol .and. &
            lhs%n_random_trial_vectors == rhs%n_random_trial_vectors .and. &
            lhs%n_iter == rhs%n_iter .and. lhs%verbose == rhs%verbose

    end function equal_stability_c

    logical function not_equal_stability_c(lhs, rhs)
        !
        ! this function overloads the negated comparison operator to compare stability
        ! settings to different stability settings
        !
        use c_interface, only: stability_settings_type_c

        type(stability_settings_type_c), intent(in) :: lhs, rhs
        
        not_equal_stability_c = .not. (lhs == rhs)

    end function not_equal_stability_c

end module test_reference
