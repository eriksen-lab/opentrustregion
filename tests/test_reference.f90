! Copyright (C) 2025- Jonas Greiner
!
! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.

module test_reference

    use opentrustregion, only: rp, ip, kw_len, stderr
    use c_interface, only: c_rp, c_ip
    use, intrinsic :: iso_c_binding, only: c_bool, c_char, c_funptr, c_f_procpointer

    implicit none

    ! tolerance for real comparisons
    real(rp), parameter :: tol = 1e-10_rp
    real(c_rp), parameter :: tol_c = real(tol, kind=c_rp)

    ! number of parameters
    integer(ip), parameter :: n_param = 3_ip
    integer(c_ip), parameter :: n_param_c = int(n_param, kind=c_ip)

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

    function test_update_orbs_funptr(update_orbs_funptr, test_name, message) &
        result(test_passed)
        !
        ! this function tests a provided orbital updating function pointer
        !
        use opentrustregion, only: update_orbs_type, hess_x_type

        procedure(update_orbs_type), intent(in), pointer :: update_orbs_funptr
        character(*), intent(in) :: test_name, message
        logical :: test_passed

        real(rp), allocatable :: kappa(:), grad(:), h_diag(:)
        real(rp) :: func
        integer(ip) :: error
        procedure(hess_x_type), pointer :: hess_x_funptr

        ! assume tests pass
        test_passed = .true.

        ! allocate arrays
        allocate(kappa(n_param), grad(n_param), h_diag(n_param))

        ! initialize orbital update
        kappa = 1.0_rp

        ! call orbital update
        call update_orbs_funptr(kappa, func, grad, h_diag, hess_x_funptr, error)

        ! check for error
        if (error /= 0) then
            test_passed = .false.
            write (stderr, *) "test_"//test_name//" failed: Error produced"//message// &
                "."
        end if

        ! check objective function value
        if (abs(func - 3.0_rp) > tol) then
            test_passed = .false.
            write (stderr, *) "test_"//test_name//" failed: Objective function "// &
                "value returned"//message//" wrong."
        end if

        ! check gradient
        if (any(abs(grad - 2.0_rp) > tol)) then
            test_passed = .false.
            write (stderr, *) "test_"//test_name//" failed: Gradient returned"// &
                message//" wrong."
        end if

        ! check Hessian diagonal
        if (any(abs(h_diag - 3.0_rp) > tol)) then
            test_passed = .false.
            write (stderr, *) "test_"//test_name//" failed: Hessian diagonal "// &
                "returned"//message//" wrong."
        end if

        ! deallocate arrays
        deallocate(kappa, grad, h_diag)

        ! test returned Hessian linear transformation
        test_passed = test_passed .and. &
            test_hess_x_funptr(hess_x_funptr, test_name, " by Hessian linear "// &
                               "transformation function returned"//message)

    end function test_update_orbs_funptr

    function test_update_orbs_c_funptr(update_orbs_c_funptr, test_name, message) &
        result(test_passed)
        !
        ! this function tests a provided orbital updating C function pointer
        !
        use c_interface, only: update_orbs_c_type

        type(c_funptr), intent(in) :: update_orbs_c_funptr
        character(*), intent(in) :: test_name, message
        logical :: test_passed

        procedure(update_orbs_c_type), pointer :: update_orbs_funptr
        real(c_rp), allocatable :: kappa(:), grad(:), h_diag(:)
        real(c_rp) :: func
        integer(c_ip) :: error
        type(c_funptr) :: hess_x_c_funptr

        ! assume tests pass
        test_passed = .true.

        ! convert to Fortran function pointer
        call c_f_procpointer(cptr=update_orbs_c_funptr, fptr=update_orbs_funptr)

        ! allocate arrays
        allocate(kappa(n_param), grad(n_param), h_diag(n_param))

        ! initialize orbital update
        kappa = 1.0_c_rp

        ! call orbital update
        error = update_orbs_funptr(kappa, func, grad, h_diag, hess_x_c_funptr)

        ! check for error
        if (error /= 0) then
            test_passed = .false.
            write (stderr, *) "test_"//test_name//" failed: Error produced"//message// &
                "."
        end if

        ! check objective function value
        if (abs(func - 3.0_c_rp) > tol_c) then
            test_passed = .false.
            write (stderr, *) "test_"//test_name//" failed: Objective function "// &
                "value returned"//message//" wrong."
        end if

        ! check gradient
        if (any(abs(grad - 2.0_c_rp) > tol_c)) then
            test_passed = .false.
            write (stderr, *) "test_"//test_name//" failed: Gradient returned"// &
                message//" wrong."
        end if

        ! check Hessian diagonal
        if (any(abs(h_diag - 3.0_c_rp) > tol_c)) then
            test_passed = .false.
            write (stderr, *) "test_"//test_name//" failed: Hessian diagonal "// &
                "returned"//message//" wrong."
        end if

        ! deallocate arrays
        deallocate(kappa, grad, h_diag)

        ! test returned Hessian linear transformation
        test_passed = test_passed .and. &
            test_hess_x_c_funptr(hess_x_c_funptr, test_name, " by Hessian linear "// &
                                 "transformation function returned"//message)

    end function test_update_orbs_c_funptr

    function test_hess_x_funptr(hess_x_funptr, test_name, message) result(test_passed)
        !
        ! this function tests a provided Hessian linear transformation function pointer
        !
        use opentrustregion, only: hess_x_type

        procedure(hess_x_type), intent(in), pointer :: hess_x_funptr
        character(*), intent(in) :: test_name, message
        logical :: test_passed

        real(rp), allocatable :: x(:), hess_x(:)
        integer(ip) :: error

        ! assume tests pass
        test_passed = .true.

        ! allocate arrays
        allocate(x(n_param), hess_x(n_param))

        ! initialize trial vector
        x = 1.0_rp

        ! call Hessian linear transformation
        call hess_x_funptr(x, hess_x, error)

        ! check for error
        if (error /= 0) then
            test_passed = .false.
            write (stderr, *) "test_"//test_name//" failed: Error produced"//message// &
                "."
        end if

        ! check Hessian linear transformation
        if (any(abs(hess_x - 4.0_rp) > tol)) then
            test_passed = .false.
            write (stderr, *) "test_"//test_name//" failed: Hessian linear "// &
                "transformation returned"//message//" wrong."
        end if

        ! deallocate arrays
        deallocate(x, hess_x)

    end function test_hess_x_funptr

    function test_hess_x_c_funptr(hess_x_c_funptr, test_name, message) &
        result(test_passed)
        !
        ! this function tests a provided Hessian linear transformation C function 
        ! pointer
        !
        use c_interface, only: hess_x_c_type

        type(c_funptr), intent(in) :: hess_x_c_funptr
        character(*), intent(in) :: test_name, message
        logical :: test_passed

        procedure(hess_x_c_type), pointer :: hess_x_funptr_c
        real(c_rp), allocatable :: x(:), hess_x(:)
        integer(c_ip) :: error

        ! assume tests pass
        test_passed = .true.

        ! convert to Fortran function pointer
        call c_f_procpointer(cptr=hess_x_c_funptr, fptr=hess_x_funptr_c)

        ! allocate arrays
        allocate(x(n_param), hess_x(n_param))

        ! initialize trial vector
        x = 1.0_c_rp

        ! call Hessian linear transformation
        error = hess_x_funptr_c(x, hess_x)

        ! check for error
        if (error /= 0) then
            test_passed = .false.
            write (stderr, *) "test_"//test_name//" failed: Error produced"//message// &
               "."
        end if

        ! check Hessian linear transformation
        if (any(abs(hess_x - 4.0_c_rp) > tol_c)) then
            test_passed = .false.
            write (stderr, *) "test_"//test_name//" failed: Hessian linear "// &
                "transformation returned"//message//" wrong."
        end if

        ! deallocate arrays
        deallocate(x, hess_x)

    end function test_hess_x_c_funptr

    function test_obj_func_funptr(obj_func_funptr, test_name, message) &
        result(test_passed)
        !
        ! this function tests a provided objective function function pointer
        !
        use opentrustregion, only: obj_func_type

        procedure(obj_func_type), intent(in), pointer :: obj_func_funptr
        character(*), intent(in) :: test_name, message
        logical :: test_passed

        real(rp), allocatable :: kappa(:)
        real(rp) :: func
        integer(ip) :: error

        ! assume tests pass
        test_passed = .true.

        ! allocate arrays
        allocate(kappa(n_param))

        ! initialize orbital update
        kappa = 1.0_rp

        ! call objective function
        func = obj_func_funptr(kappa, error)

        ! check for error
        if (error /= 0) then
            write (stderr, *) "test_"//test_name//" failed: Error produced"//message// &
                "."
            test_passed = .false.
        end if

        ! check objective function
        if (abs(func - 3.0_rp) > tol) then
            write (stderr, *) "test_"//test_name//" failed: Function value returned"// &
                message//" wrong."
            test_passed = .false.
        end if

        ! deallocate arrays
        deallocate(kappa)

    end function test_obj_func_funptr

    function test_obj_func_c_funptr(obj_func_c_funptr, test_name, message) &
        result(test_passed)
        !
        ! this function tests a provided objective function C function pointer
        !
        use c_interface, only: obj_func_c_type

        type(c_funptr), intent(in) :: obj_func_c_funptr
        character(*), intent(in) :: test_name, message
        logical :: test_passed

        procedure(obj_func_c_type), pointer :: obj_func_funptr
        real(c_rp), allocatable :: kappa(:)
        real(c_rp) :: func
        integer(c_ip) :: error

        ! assume tests pass
        test_passed = .true.

        ! convert to Fortran function pointer
        call c_f_procpointer(cptr=obj_func_c_funptr, fptr=obj_func_funptr)

        ! allocate arrays
        allocate(kappa(n_param))

        ! initialize orbital update
        kappa = 1.0_c_rp

        ! call objective function
        error = obj_func_funptr(kappa, func)

        ! check for error
        if (error /= 0) then
            write (stderr, *) "test_"//test_name//" failed: Error produced"//message// &
                "."
            test_passed = .false.
        end if

        ! check objective function
        if (abs(func - 3.0_c_rp) > tol_c) then
            write (stderr, *) "test_"//test_name//" failed: Function value returned"// &
                message//" wrong."
            test_passed = .false.
        end if

        ! deallocate arrays
        deallocate(kappa)

    end function test_obj_func_c_funptr

    function test_precond_funptr(precond_funptr, test_name, message) &
        result(test_passed)
        !
        ! this function tests a provided preconditioner function pointer
        !
        use opentrustregion, only: precond_type

        procedure(precond_type), intent(in), pointer :: precond_funptr
        character(*), intent(in) :: test_name, message
        logical :: test_passed
        
        real(rp), allocatable :: residual(:), precond_residual(:)
        integer(ip) :: error

        ! assume tests pass
        test_passed = .true.

        ! allocate arrays
        allocate(residual(n_param), precond_residual(n_param))

        ! initialize residual
        residual = 1.0_rp

        ! call preconditioning subroutine
        call precond_funptr(residual, 5.0_rp, precond_residual, error)

        ! check for error
        if (error /= 0) then
            write (stderr, *) "test_"//test_name//" failed: Error produced"//message// &
                "."
            test_passed = .false.
        end if

        ! check preconditioned residual
        if (any(abs(precond_residual - 5.0_rp) > tol)) then
            write (stderr, *) "test_"//test_name//" failed: Preconditioned "// &
                "residual returned"//message//" wrong."
            test_passed = .false.
        end if

        ! deallocate arrays
        deallocate(residual, precond_residual)

    end function test_precond_funptr

    function test_precond_c_funptr(precond_c_funptr, test_name, message) &
        result(test_passed)
        !
        ! this function tests a provided preconditioner C function pointer
        !
        use c_interface, only: precond_c_type

        type(c_funptr), intent(in) :: precond_c_funptr
        character(*), intent(in) :: test_name, message
        logical :: test_passed

        procedure(precond_c_type), pointer :: precond_funptr
        real(c_rp), allocatable :: residual(:), precond_residual(:)
        integer(c_ip) :: error

        ! assume tests pass
        test_passed = .true.

        ! convert to Fortran function pointer
        call c_f_procpointer(cptr=precond_c_funptr, fptr=precond_funptr)

        ! allocate arrays
        allocate(residual(n_param), precond_residual(n_param))

        ! initialize residual
        residual = 1.0_c_rp

        ! call preconditioning function
        error = precond_funptr(residual, 5.0_c_rp, precond_residual)

        ! check for error
        if (error /= 0) then
            write (stderr, *) "test_"//test_name//" failed: Error produced"//message// &
                "."
            test_passed = .false.
        end if

        ! check preconditioned residual
        if (any(abs(precond_residual - 5.0_c_rp) > tol_c)) then
            write (stderr, *) "test_"//test_name//" failed: Preconditioned "// &
                "residual returned"//message//" wrong."
            test_passed = .false.
        end if

        ! deallocate arrays
        deallocate(residual, precond_residual)

    end function test_precond_c_funptr

    function test_project_funptr(project_funptr, test_name, message) result(test_passed)
        !
        ! this function tests a provided projection function pointer
        !
        use opentrustregion, only: project_type

        procedure(project_type), intent(in), pointer :: project_funptr
        character(*), intent(in) :: test_name, message
        logical :: test_passed
        
        real(rp), allocatable :: vector(:)
        integer(ip) :: error

        ! assume tests pass
        test_passed = .true.

        ! allocate arrays
        allocate(vector(n_param))

        ! initialize vector
        vector = 1.0_rp

        ! call projection subroutine
        call project_funptr(vector, error)

        ! check for error
        if (error /= 0) then
            write (stderr, *) "test_"//test_name//" failed: Error produced"//message// &
                "."
            test_passed = .false.
        end if

        ! check projected vector
        if (any(abs(vector - 2.0_rp) > tol)) then
            write (stderr, *) "test_"//test_name//" failed: Projected vector "// &
                "returned"//message//" wrong."
            test_passed = .false.
        end if

        ! deallocate arrays
        deallocate(vector)

    end function test_project_funptr

    function test_project_c_funptr(project_c_funptr, test_name, message) &
        result(test_passed)
        !
        ! this function tests a provided projection C function pointer
        !
        use c_interface, only: project_c_type

        type(c_funptr), intent(in) :: project_c_funptr
        character(*), intent(in) :: test_name, message
        logical :: test_passed

        procedure(project_c_type), pointer :: project_funptr
        real(c_rp), allocatable :: vector(:)
        integer(c_ip) :: error

        ! assume tests pass
        test_passed = .true.

        ! convert to Fortran function pointer
        call c_f_procpointer(cptr=project_c_funptr, fptr=project_funptr)

        ! allocate arrays
        allocate(vector(n_param))

        ! initialize vector
        vector = 1.0_c_rp

        ! call projection function
        error = project_funptr(vector)

        ! check for error
        if (error /= 0) then
            write (stderr, *) "test_"//test_name//" failed: Error produced"//message// &
                "."
            test_passed = .false.
        end if

        ! check projected vector
        if (any(abs(vector - 2.0_c_rp) > tol_c)) then
            write (stderr, *) "test_"//test_name//" failed: Projected vector "// &
                "returned"//message//" wrong."
            test_passed = .false.
        end if

        ! deallocate arrays
        deallocate(vector)

    end function test_project_c_funptr

    function test_conv_check_funptr(conv_check_funptr, test_name, message) &
        result(test_passed)
        !
        ! this function tests a provided convergence check function pointer
        !
        use opentrustregion, only: conv_check_type

        procedure(conv_check_type), intent(in), pointer :: conv_check_funptr
        character(*), intent(in) :: test_name, message
        logical :: test_passed
        
        logical :: converged
        integer(ip) :: error

        ! assume tests pass
        test_passed = .true.

        ! call convergence check function
        converged = conv_check_funptr(error)

        ! check for error
        if (error /= 0) then
            write (stderr, *) "test_"//test_name//" failed: Error produced"//message// &
                "."
            test_passed = .false.
        end if

        ! check convergence logical
        if (.not. converged) then
            write (stderr, *) "test_"//test_name//" failed: Convergence logical "// &
                "returned"//message//" wrong."
            test_passed = .false.
        end if

    end function test_conv_check_funptr

    function test_conv_check_c_funptr(conv_check_c_funptr, test_name, message) &
        result(test_passed)
        !
        ! this function tests a provided convergence check C function pointer
        !
        use c_interface, only: conv_check_c_type

        type(c_funptr), intent(in) :: conv_check_c_funptr
        character(*), intent(in) :: test_name, message
        logical :: test_passed
        
        procedure(conv_check_c_type), pointer :: conv_check_funptr
        logical(c_bool) :: converged
        integer(ip) :: error

        ! assume tests pass
        test_passed = .true.

        ! convert to Fortran function pointer
        call c_f_procpointer(cptr=conv_check_c_funptr, fptr=conv_check_funptr)

        ! call convergence check function
        error = conv_check_funptr(converged)

        ! check for error
        if (error /= 0) then
            write (stderr, *) "test_"//test_name//" failed: Error produced"//message// &
                "."
            test_passed = .false.
        end if

        ! check convergence logical
        if (.not. converged) then
            write (stderr, *) "test_"//test_name//" failed: Convergence logical "// &
                "returned"//message//" wrong."
            test_passed = .false.
        end if

    end function test_conv_check_c_funptr

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
        lhs%project => null()
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
        lhs%project => null()
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
