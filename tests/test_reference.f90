! Copyright (C) 2025- Jonas Greiner
!
! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.

module test_reference

    use opentrustregion, only: rp, ip, kw_len
    use c_interface, only: c_rp, c_ip
    use, intrinsic :: iso_c_binding, only: c_bool, c_char

    implicit none

    ! tolerance for real comparisons
    real(rp), parameter :: tol = 1e-10_rp
    real(c_rp), parameter :: tol_c = real(tol, kind=c_rp)

    ! derived types for solver settings
    type ref_settings_type
        logical :: stability, line_search
        real(rp) :: conv_tol, start_trust_radius, global_red_factor, local_red_factor
        integer(ip) :: n_random_trial_vectors, n_macro, n_micro, &
                       jacobi_davidson_start, seed, verbose, n_iter
        character(kw_len, c_char) :: subsystem_solver, diag_solver
    end type

    type, bind(C) :: ref_settings_type_c
        logical(c_bool) :: stability, line_search
        real(c_rp) :: conv_tol, start_trust_radius, global_red_factor, local_red_factor
        integer(c_ip) :: n_random_trial_vectors, n_macro, n_micro, &
                         jacobi_davidson_start, seed, verbose, n_iter
        character(c_char) :: subsystem_solver(kw_len + 1), diag_solver(kw_len + 1)
    end type

    ! general reference parameters
    type(ref_settings_type) :: ref_settings = &
        ref_settings_type(stability = .true., line_search = .true., &
                          conv_tol = 1e-3_rp, start_trust_radius = 0.2_rp, &
                          global_red_factor = 1e-2_rp, local_red_factor = 1e-3_rp, &
                          n_random_trial_vectors = 5, n_macro = 300, n_micro = 200, &
                          jacobi_davidson_start = 10, seed = 33, verbose = 3, &
                          n_iter = 50, subsystem_solver = "tcg", &
                          diag_solver = "jacobi-davidson")

    interface assignment(=)
        module procedure assign_ref_to_solver
        module procedure assign_ref_to_stability
        module procedure assign_ref_to_solver_c
        module procedure assign_ref_to_stability_c
        module procedure assign_ref_to_ref_c
    end interface

    interface operator(==)
        module procedure equal_solver_to_ref
        module procedure equal_stability_to_ref
        module procedure equal_solver_c_to_ref
        module procedure equal_stability_c_to_ref
        module procedure equal_solver
        module procedure equal_stability
        module procedure equal_solver_c
        module procedure equal_stability_c
    end interface

    interface operator(/=)
        module procedure not_equal_solver_to_ref
        module procedure not_equal_stability_to_ref
        module procedure not_equal_solver_c_to_ref
        module procedure not_equal_stability_c_to_ref
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

        ref_settings_out = ref_settings

    end subroutine get_reference_values

    subroutine assign_ref_to_solver(lhs, rhs)
        !
        ! this subroutine overloads the assignment operator to set solver settings to 
        ! reference values
        !
        use opentrustregion, only: solver_settings_type

        type(solver_settings_type), intent(out) :: lhs
        type(ref_settings_type), intent(in) :: rhs

        ! unassociate function pointers
        lhs%precond => null()
        lhs%conv_check => null()
        lhs%logger => null()

        ! set reference values
        lhs%stability = rhs%stability
        lhs%line_search = rhs%line_search
        lhs%conv_tol = rhs%conv_tol
        lhs%start_trust_radius = rhs%start_trust_radius
        lhs%global_red_factor = rhs%global_red_factor
        lhs%local_red_factor  = rhs%local_red_factor
        lhs%n_random_trial_vectors = rhs%n_random_trial_vectors
        lhs%n_macro = rhs%n_macro
        lhs%n_micro = rhs%n_micro
        lhs%jacobi_davidson_start = rhs%jacobi_davidson_start
        lhs%seed = rhs%seed
        lhs%verbose = rhs%verbose
        lhs%subsystem_solver = rhs%subsystem_solver

        ! set initialization logical
        lhs%initialized = .true.

    end subroutine assign_ref_to_solver

    subroutine assign_ref_to_stability(lhs, rhs)
        !
        ! this subroutine overloads the assignment operator to set stability settings 
        ! to reference values
        !
        use opentrustregion, only: stability_settings_type

        type(stability_settings_type), intent(out) :: lhs
        type(ref_settings_type), intent(in) :: rhs

        ! unassociate function pointers
        lhs%precond => null()
        lhs%logger => null()

        ! set reference values
        lhs%conv_tol = rhs%conv_tol
        lhs%n_random_trial_vectors = rhs%n_random_trial_vectors
        lhs%n_iter = rhs%n_iter
        lhs%jacobi_davidson_start = rhs%jacobi_davidson_start
        lhs%seed = rhs%seed
        lhs%verbose = rhs%verbose
        lhs%diag_solver = rhs%diag_solver

        ! set initialization logical
        lhs%initialized = .true.

    end subroutine assign_ref_to_stability

    subroutine assign_ref_to_solver_c(lhs_c, rhs)
        !
        ! this subroutine overloads the assignment operator to set C solver settings to 
        ! reference values
        !
        use c_interface, only: solver_settings_type_c, assignment(=)
        use opentrustregion, only: solver_settings_type

        type(solver_settings_type_c), intent(out) :: lhs_c
        type(ref_settings_type), intent(in) :: rhs

        type(solver_settings_type) :: lhs

        lhs = rhs
        lhs_c = lhs

    end subroutine assign_ref_to_solver_c

    subroutine assign_ref_to_stability_c(lhs_c, rhs)
        !
        ! this subroutine overloads the assignment operator to set C stability settings 
        ! to reference values
        !
        use c_interface, only: stability_settings_type_c, assignment(=)
        use opentrustregion, only: stability_settings_type

        type(stability_settings_type_c), intent(out) :: lhs_c
        type(ref_settings_type), intent(in) :: rhs

        type(stability_settings_type) :: lhs

        lhs = rhs
        lhs_c = lhs

    end subroutine assign_ref_to_stability_c

    subroutine assign_ref_to_ref_c(lhs, rhs)
        !
        ! this subroutine overloads the assignment operator to convert reference values
        ! to their C counterpart
        !
        use c_interface, only: character_to_c

        type(ref_settings_type_c), intent(out) :: lhs
        type(ref_settings_type), intent(in) :: rhs

        lhs%stability = logical(rhs%stability, kind=c_bool)
        lhs%line_search = logical(rhs%line_search, kind=c_bool)
        lhs%conv_tol = real(rhs%conv_tol, kind=c_rp)
        lhs%start_trust_radius = real(rhs%start_trust_radius, kind=c_rp)
        lhs%global_red_factor = real(rhs%global_red_factor, kind=c_rp)
        lhs%local_red_factor  = real(rhs%local_red_factor, kind=c_rp)
        lhs%n_random_trial_vectors = int(rhs%n_random_trial_vectors, kind=c_ip)
        lhs%n_macro = int(rhs%n_macro, kind=c_ip)
        lhs%n_micro = int(rhs%n_micro, kind=c_ip)
        lhs%jacobi_davidson_start = int(rhs%jacobi_davidson_start, kind=c_ip)
        lhs%seed = int(rhs%seed, kind=c_ip)
        lhs%verbose = int(rhs%verbose, kind=c_ip)
        lhs%n_iter = int(rhs%n_iter, kind=c_ip)
        lhs%subsystem_solver = character_to_c(rhs%subsystem_solver)
        lhs%diag_solver = character_to_c(rhs%diag_solver)

    end subroutine assign_ref_to_ref_c

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
            abs(lhs%conv_tol - rhs%conv_tol) <= tol .and. &
            abs(lhs%start_trust_radius - rhs%start_trust_radius) <= tol .and. &
            abs(lhs%global_red_factor - rhs%global_red_factor) <= tol .and. &
            abs(lhs%local_red_factor - rhs%local_red_factor) <= tol .and. &
            lhs%n_random_trial_vectors == rhs%n_random_trial_vectors .and. &
            lhs%n_macro == rhs%n_macro .and. lhs%n_micro == rhs%n_micro .and. &
            lhs%jacobi_davidson_start == rhs%jacobi_davidson_start .and. &
            lhs%seed == rhs%seed .and. lhs%verbose == rhs%verbose .and. &
            lhs%subsystem_solver == rhs%subsystem_solver

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

        equal_stability_to_ref = abs(lhs%conv_tol - rhs%conv_tol) <= tol .and. &
            lhs%n_random_trial_vectors == rhs%n_random_trial_vectors .and. &
            lhs%n_iter == rhs%n_iter .and. &
            lhs%jacobi_davidson_start == rhs%jacobi_davidson_start .and. &
            lhs%seed == rhs%seed .and. lhs%verbose == rhs%verbose .and. &
            lhs%diag_solver == rhs%diag_solver

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

    logical function equal_solver_c_to_ref(lhs_c, rhs)
        !
        ! this function overloads the comparison operator to compare solver settings to 
        ! reference values
        !
        use c_interface, only: solver_settings_type_c, assignment(=)
        use opentrustregion, only: solver_settings_type

        type(solver_settings_type_c), intent(in) :: lhs_c
        type(ref_settings_type), intent(in) :: rhs

        type(solver_settings_type) :: lhs

        lhs = lhs_c
        equal_solver_c_to_ref = lhs == rhs

    end function equal_solver_c_to_ref

    logical function not_equal_solver_c_to_ref(lhs, rhs)
        !
        ! this function overloads the negated comparison operator to compare solver 
        ! settings to reference values
        !
        use c_interface, only: solver_settings_type_c

        type(solver_settings_type_c), intent(in) :: lhs
        type(ref_settings_type), intent(in) :: rhs
        
        not_equal_solver_c_to_ref = .not. (lhs == rhs)

    end function not_equal_solver_c_to_ref

    logical function equal_stability_c_to_ref(lhs_c, rhs)
        !
        ! this function overloads the comparison operator to compare stability settings 
        ! to reference values
        !
        use c_interface, only: stability_settings_type_c, assignment(=)
        use opentrustregion, only: stability_settings_type

        type(stability_settings_type_c), intent(in) :: lhs_c
        type(ref_settings_type), intent(in) :: rhs

        type(stability_settings_type) :: lhs
        
        lhs = lhs_c
        equal_stability_c_to_ref = lhs == rhs

    end function equal_stability_c_to_ref

    logical function not_equal_stability_c_to_ref(lhs, rhs)
        !
        ! this function overloads the negated comparison operator to compare stability
        ! settings to reference values
        !
        use c_interface, only: stability_settings_type_c

        type(stability_settings_type_c), intent(in) :: lhs
        type(ref_settings_type), intent(in) :: rhs
        
        not_equal_stability_c_to_ref = .not. (lhs == rhs)

    end function not_equal_stability_c_to_ref

    logical function equal_solver(lhs, rhs)
        !
        ! this function overloads the comparison operator to compare solver settings 
        ! to different solver settings
        !
        use opentrustregion, only: solver_settings_type

        type(solver_settings_type), intent(in) :: lhs, rhs
        
        equal_solver = (lhs%stability .eqv. rhs%stability) .and. &
            (lhs%line_search .eqv. rhs%line_search) .and. &
            (lhs%initialized .eqv. rhs%initialized) .and. &
            abs(lhs%conv_tol - rhs%conv_tol) <= tol .and. &
            abs(lhs%start_trust_radius - rhs%start_trust_radius) <= tol .and. &
            abs(lhs%global_red_factor - rhs%global_red_factor) <= tol .and. &
            abs(lhs%local_red_factor - rhs%local_red_factor) <= tol .and. &
            lhs%n_random_trial_vectors == rhs%n_random_trial_vectors .and. &
            lhs%n_macro == rhs%n_macro .and. lhs%n_micro == rhs%n_micro .and. &
            lhs%jacobi_davidson_start == rhs%jacobi_davidson_start .and. &
            lhs%seed == rhs%seed .and. lhs%verbose == rhs%verbose .and. &
            lhs%subsystem_solver == rhs%subsystem_solver

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
        
        equal_stability = abs(lhs%conv_tol - rhs%conv_tol) <= tol .and. &
            lhs%n_random_trial_vectors == rhs%n_random_trial_vectors .and. &
            lhs%n_iter == rhs%n_iter .and. &
            lhs%jacobi_davidson_start == rhs%jacobi_davidson_start .and. &
            lhs%seed == rhs%seed .and. lhs%verbose == rhs%verbose .and. &
            lhs%diag_solver == rhs%diag_solver

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

    logical function equal_solver_c(lhs_c, rhs)
        !
        ! this function overloads the comparison operator to compare solver settings 
        ! to different solver settings
        !
        use c_interface, only: solver_settings_type_c, assignment(=)
        use opentrustregion, only: solver_settings_type

        type(solver_settings_type_c), intent(in) :: lhs_c
        type(solver_settings_type), intent(in) :: rhs
        
        type(solver_settings_type) :: lhs

        lhs = lhs_c
        equal_solver_c = lhs == rhs

    end function equal_solver_c

    logical function not_equal_solver_c(lhs_c, rhs)
        !
        ! this function overloads the negated comparison operator to compare solver 
        ! settings to different solver settings
        !
        use c_interface, only: solver_settings_type_c
        use opentrustregion, only: solver_settings_type

        type(solver_settings_type_c), intent(in) :: lhs_c
        type(solver_settings_type), intent(in) :: rhs
        
        not_equal_solver_c = .not. (lhs_c == rhs)

    end function not_equal_solver_c

    logical function equal_stability_c(lhs_c, rhs)
        !
        ! this function overloads the comparison operator to compare stability settings 
        ! to different stability settings
        !
        use c_interface, only: stability_settings_type_c, assignment(=)
        use opentrustregion, only: stability_settings_type

        type(stability_settings_type_c), intent(in) :: lhs_c
        type(stability_settings_type), intent(in) :: rhs
        
        type(stability_settings_type) :: lhs

        lhs = lhs_c
        equal_stability_c = lhs == rhs

    end function equal_stability_c

    logical function not_equal_stability_c(lhs_c, rhs)
        !
        ! this function overloads the negated comparison operator to compare stability
        ! settings to different stability settings
        !
        use c_interface, only: stability_settings_type_c
        use opentrustregion, only: stability_settings_type

        type(stability_settings_type_c), intent(in) :: lhs_c
        type(stability_settings_type), intent(in) :: rhs
        
        not_equal_stability_c = .not. (lhs_c == rhs)

    end function not_equal_stability_c

end module test_reference
