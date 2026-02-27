! Copyright (C) 2025- Jonas Greiner
!
! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.

module c_interface

    use opentrustregion, only: rp, ip, kw_len, standard_solver => solver, &
                               standard_stability_check => stability_check, &
                               default_solver_settings, default_stability_settings, &
                               update_orbs_type, hess_x_type, obj_func_type, &
                               precond_type, project_type, conv_check_type, logger_type
    use, intrinsic :: iso_c_binding, only: c_double, c_int64_t, c_int32_t, c_bool, &
                                           c_ptr, c_funptr, c_f_pointer, &
                                           c_f_procpointer, c_associated, c_char, &
                                           c_null_char, c_null_funptr

    implicit none

    integer, parameter :: c_rp = c_double
#ifdef USE_ILP64
    integer, parameter :: c_ip = c_int64_t  ! 64-bit integers
    logical(c_bool), bind(C) :: ilp64 = .true._c_bool
#else
    integer, parameter :: c_ip = c_int32_t  ! 32-bit integers
    logical(c_bool), bind(C) :: ilp64 = .false._c_bool
#endif

    ! maximum keyword length
    integer(c_ip), bind(C) :: kw_len_c = kw_len

    ! define procedure pointer which will point to the Fortran procedures
    procedure(update_orbs_c_type), pointer :: update_orbs_before_wrapping => null()
    procedure(hess_x_c_type), pointer :: hess_x_before_wrapping => null()
    procedure(obj_func_c_type), pointer :: obj_func_before_wrapping => null()
    procedure(precond_c_type), pointer :: precond_before_wrapping => null()
    procedure(project_c_type), pointer :: project_before_wrapping => null()
    procedure(conv_check_c_type), pointer :: conv_check_before_wrapping => null()
    procedure(logger_c_type), pointer :: logger_before_wrapping => null()

    ! C-interoperable interfaces for the callback functions
    abstract interface
        function update_orbs_c_type(kappa_c, func_c, grad_c, h_diag_c, &
                                    hess_x_c_funptr) result(error) bind(C)
            import :: c_rp, c_funptr, c_ip

            real(c_rp), intent(in) :: kappa_c(*)
            real(c_rp), intent(out) :: func_c
            real(c_rp), intent(out) :: grad_c(*), h_diag_c(*)
            type(c_funptr), intent(out) :: hess_x_c_funptr
            integer(c_ip) :: error
        end function update_orbs_c_type
    end interface

    abstract interface
        function hess_x_c_type(x_c, hess_x_c) result(error) bind(C)
            import :: c_rp, c_ip

            real(c_rp), intent(in) :: x_c(*)
            real(c_rp), intent(out) :: hess_x_c(*)
            integer(c_ip) :: error
        end function hess_x_c_type
    end interface

    abstract interface
        function obj_func_c_type(kappa_c, func) result(error) bind(C)
            import :: c_rp, c_ip

            real(c_rp), intent(in) :: kappa_c(*)
            real(c_rp), intent(out) :: func
            integer(c_ip) :: error
        end function obj_func_c_type
    end interface

    abstract interface
        function precond_c_type(residual_c, mu_c, precond_residual_c) result(error) &
            bind(C)
            import :: c_rp, c_ip

            real(c_rp), intent(in) :: residual_c(*), mu_c
            real(c_rp), intent(out) :: precond_residual_c(*)
            integer(c_ip) :: error
        end function precond_c_type
    end interface

    abstract interface
        function project_c_type(vector_c) result(error) bind(C)
            import :: c_rp, c_ip

            real(c_rp), intent(inout), target :: vector_c(*)
            integer(c_ip) :: error
        end function project_c_type
    end interface

    abstract interface
        function conv_check_c_type(converged) result(error) bind(C)
            import :: c_bool, c_ip

            logical(c_bool), intent(out) :: converged
            integer(c_ip) :: error
        end function conv_check_c_type
    end interface

    abstract interface
        subroutine logger_c_type(message) bind(C)
            import :: c_char

            character(c_char), intent(in) :: message(*)
        end subroutine logger_c_type
    end interface

    ! derived type for solver settings
    type, bind(C) :: solver_settings_type_c
        type(c_funptr) :: precond, project, conv_check, logger
        logical(c_bool) :: stability, line_search, initialized
        real(c_rp) :: conv_tol, start_trust_radius, global_red_factor, local_red_factor
        integer(c_ip) :: n_random_trial_vectors, n_macro, n_micro, &
                         jacobi_davidson_start, seed, verbose
        character(c_char) :: subsystem_solver(kw_len + 1)
    end type

    type, bind(C) :: stability_settings_type_c
        type(c_funptr) :: precond, project, logger
        logical(c_bool) :: initialized
        real(c_rp) :: conv_tol
        integer(c_ip) :: n_random_trial_vectors, n_iter, jacobi_davidson_start, seed, &
                         verbose
        character(c_char) :: diag_solver(kw_len + 1)
    end type

    procedure(standard_solver), pointer :: solver => standard_solver
    procedure(standard_stability_check), pointer :: stability_check => &
        standard_stability_check

    ! create function pointers to ensure that routines comply with interface
    procedure(update_orbs_type), pointer :: update_orbs_f_wrapper_ptr => &
        update_orbs_f_wrapper
    procedure(hess_x_type), pointer :: hess_x_f_wrapper_ptr => hess_x_f_wrapper
    procedure(obj_func_type), pointer :: obj_func_f_wrapper_ptr => obj_func_f_wrapper
    procedure(precond_type), pointer :: precond_f_wrapper_ptr => precond_f_wrapper
    procedure(project_type), pointer :: project_f_wrapper_ptr => project_f_wrapper
    procedure(conv_check_type), pointer :: conv_check_f_wrapper_ptr => &
        conv_check_f_wrapper
    procedure(logger_type), pointer :: logger_f_wrapper_ptr => logger_f_wrapper

    ! interfaces for converting C settings to Fortran settings
    interface assignment(=)
        module procedure assign_solver_f_c
        module procedure assign_stability_f_c
        module procedure assign_solver_c_f
        module procedure assign_stability_c_f
    end interface

contains

    function solver_c_wrapper(update_orbs_c_funptr, obj_func_c_funptr, n_param_c, &
                              settings_c) result(error_c) bind(C, name="solver")
        !
        ! this function exposes a Fortran-implemented solver subroutine to C
        !
        use opentrustregion, only: solver_settings_type

        type(c_funptr), intent(in), value :: update_orbs_c_funptr, obj_func_c_funptr
        integer(c_ip), intent(in), value :: n_param_c
        type(solver_settings_type_c), intent(in), value :: settings_c
        integer(c_ip) :: error_c

        procedure(update_orbs_f_wrapper), pointer :: update_orbs
        procedure(obj_func_f_wrapper), pointer :: obj_func
        integer(ip) :: n_param, error
        type(solver_settings_type) :: settings

        ! associate the input C pointer to update_orbs subroutine to a Fortran
        ! procedure pointer
        call c_f_procpointer(cptr=update_orbs_c_funptr, &
                             fptr=update_orbs_before_wrapping)
        call c_f_procpointer(cptr=obj_func_c_funptr, fptr=obj_func_before_wrapping)

        ! associate procedure pointer to wrapper function
        update_orbs => update_orbs_f_wrapper
        obj_func => obj_func_f_wrapper

        ! convert dummy argument to Fortran kind
        n_param = int(n_param_c, kind=ip)

        ! convert settings
        settings = settings_c

        ! call solver
        call solver(update_orbs, obj_func, n_param, error, settings)

        ! convert return arguments to C kind
        error_c = int(error, kind=c_ip)

    end function solver_c_wrapper

    function stability_check_c_wrapper(h_diag_c, hess_x_c_funptr, n_param_c, &
                                       stable_c, settings_c, kappa_c_ptr) &                    
        result(error_c) bind(C, name="stability_check")
        !
        ! this function exposes a Fortran-implemented stability check subroutine to C
        !
        use opentrustregion, only: stability_settings_type

        real(c_rp), intent(in), target :: h_diag_c(*)
        integer(c_ip), intent(in), value :: n_param_c
        type(c_funptr), intent(in), value :: hess_x_c_funptr
        logical(c_bool), intent(out) :: stable_c
        type(stability_settings_type_c), intent(in), value :: settings_c
        type(c_ptr), intent(in), value :: kappa_c_ptr
        integer(c_ip) :: error_c

        real(rp), pointer :: h_diag_ptr(:), kappa_ptr(:)
        real(c_rp), pointer :: kappa_ptr_c(:)
        logical :: stable
        integer(ip) :: error
        procedure(hess_x_f_wrapper), pointer :: hess_x
        type(stability_settings_type) :: settings

        ! associate the input C pointer to update_orbs subroutine to a Fortran
        ! procedure pointer
        call c_f_procpointer(cptr=hess_x_c_funptr, fptr=hess_x_before_wrapping)

        ! associate procedure pointer to wrapper function
        hess_x => hess_x_f_wrapper

        ! convert C pointers
        if (rp == c_rp) then
            h_diag_ptr => h_diag_c(:n_param_c)
        else
            allocate(h_diag_ptr(n_param_c))
            h_diag_ptr = real(h_diag_c(:n_param_c), kind=rp)
        end if
        if (c_associated(kappa_c_ptr)) then
            if (rp == c_rp) then
                call c_f_pointer(cptr=kappa_c_ptr, fptr=kappa_ptr, shape=[n_param_c])
            else
                call c_f_pointer(cptr=kappa_c_ptr, fptr=kappa_ptr_c, shape=[n_param_c])
                allocate(kappa_ptr(n_param_c))
            end if
        end if

        ! convert settings
        settings = settings_c

        ! call stability check
        if (c_associated(kappa_c_ptr)) then
            call stability_check(h_diag_ptr, hess_x, stable, error, settings, &
                                 kappa=kappa_ptr)
            if (rp /= c_rp) then
                kappa_ptr_c = kappa_ptr
                deallocate(kappa_ptr)
            end if
        else
            call stability_check(h_diag_ptr, hess_x, stable, error, settings)
        end if
        if (rp /= c_rp) deallocate(h_diag_ptr)

        ! convert return arguments to C kind
        stable_c = logical(stable, kind=c_bool)
        error_c = int(error, kind=c_ip)

    end function stability_check_c_wrapper

    subroutine update_orbs_f_wrapper(kappa, func, grad, h_diag, hess_x, error)
        !
        ! this subroutine exposes a C-implemented orbital update function to Fortran
        !
        use opentrustregion, only: hess_x_type

        real(rp), intent(in), target :: kappa(:)
        real(rp), intent(out) :: func
        real(rp), intent(out), target :: grad(:), h_diag(:)
        procedure(hess_x_type), intent(out), pointer :: hess_x
        integer(ip), intent(out) :: error

        real(c_rp) :: func_c
        real(c_rp), pointer :: kappa_c(:), grad_c(:), h_diag_c(:)
        type(c_funptr) :: hess_x_c_funptr
        integer(c_ip) :: error_c

        ! convert arguments to C kind
        if (rp == c_rp) then
            kappa_c => kappa
            grad_c => grad
            h_diag_c => h_diag
        else
            allocate(kappa_c(size(kappa)))
            allocate(grad_c(size(kappa)))
            allocate(h_diag_c(size(kappa)))
            kappa_c = real(kappa, kind=c_rp)
        end if

        ! call update_orbs C function
        error_c = update_orbs_before_wrapping(kappa_c, func_c, grad_c, h_diag_c, &
                                              hess_x_c_funptr)

        ! convert arguments to Fortran kind
        func = real(func_c, kind=rp)
        error = int(error_c, kind=ip)
        if (rp /= c_rp) then
            grad = real(grad_c, kind=rp)
            h_diag = real(h_diag_c, kind=rp)
            deallocate(kappa_c)
            deallocate(grad_c)
            deallocate(h_diag_c)
        end if

        ! associate the input C pointer to hess_x function to a Fortran procedure
        ! pointer
        call c_f_procpointer(cptr=hess_x_c_funptr, fptr=hess_x_before_wrapping)

        ! associate procedure pointer to wrapper function
        hess_x => hess_x_f_wrapper

    end subroutine update_orbs_f_wrapper

    subroutine hess_x_f_wrapper(x, hess_x, error)
        !
        ! this subroutine exposes a C-implemented Hessian linear transformation to 
        ! Fortran
        !
        real(rp), intent(in), target :: x(:)
        real(rp), intent(out), target :: hess_x(:)
        integer(ip), intent(out) :: error

        real(c_rp), pointer :: x_c(:), hess_x_c(:)
        integer(c_ip) :: error_c

        ! convert arguments to C kind
        if (rp == c_rp) then
            x_c => x
            hess_x_c => hess_x
        else
            allocate(x_c(size(x)))
            allocate(hess_x_c(size(x)))
            x_c = real(x, kind=c_rp)
        end if

        ! call C function
        error_c = hess_x_before_wrapping(x_c, hess_x_c)

        ! convert arguments to Fortran kind
        error = int(error_c, kind=ip)
        if (rp /= c_rp) then
            hess_x = real(hess_x_c, kind=rp)
            deallocate(x_c)
            deallocate(hess_x_c)
        end if

    end subroutine hess_x_f_wrapper

    function obj_func_f_wrapper(kappa, error) result(obj_func)
        !
        ! this function exposes a C-implemented objective function to Fortran
        !
        real(rp), intent(in), target :: kappa(:)
        integer(ip), intent(out) :: error
        real(rp) :: obj_func

        real(c_rp) :: obj_func_c
        real(c_rp), pointer :: kappa_c(:)
        integer(c_ip) :: error_c

        ! convert arguments to C kind
        if (rp == c_rp) then
            kappa_c => kappa
        else
            allocate(kappa_c(size(kappa)))
            kappa_c = real(kappa, kind=c_rp)
        end if

        ! call obj_func C function
        error_c = obj_func_before_wrapping(kappa_c, obj_func_c)

        ! convert arguments to Fortran kind
        obj_func = real(obj_func_c, kind=rp)
        error = int(error_c, kind=ip)
        if (rp /= c_rp) then
            deallocate(kappa_c)
        end if

    end function obj_func_f_wrapper

    subroutine precond_f_wrapper(residual, mu, precond_residual, error)
        !
        ! this subroutine exposes a C-implemented preconditioner function to Fortran
        !
        real(rp), intent(in), target :: residual(:)
        real(rp), intent(in) :: mu
        real(rp), intent(out), target :: precond_residual(:)
        integer(ip), intent(out) :: error

        real(c_rp) :: mu_c
        real(c_rp), pointer :: residual_c(:), precond_residual_c(:)
        integer(c_ip) :: error_c

        ! convert arguments to C kind
        mu_c = real(mu, kind=c_rp)
        if (rp == c_rp) then
            residual_c => residual
            precond_residual_c => precond_residual
        else
            allocate(residual_c(size(residual)))
            allocate(precond_residual_c(size(residual)))
            residual_c = real(residual, kind=c_rp)
        end if

        ! call precond C function
        error_c = precond_before_wrapping(residual_c, mu_c, precond_residual_c)

        ! convert arguments to Fortran kind
        error = int(error_c, kind=ip)
        if (rp /= c_rp) then
            precond_residual = real(precond_residual_c, kind=rp)
            deallocate(residual_c)
            deallocate(precond_residual_c)
        end if

    end subroutine precond_f_wrapper

    subroutine project_f_wrapper(vector, error)
        !
        ! this subroutine exposes a C-implemented projection function to Fortran
        !
        real(rp), intent(inout), target :: vector(:)
        integer(ip), intent(out) :: error

        real(c_rp), pointer :: vector_c(:)
        integer(c_ip) :: error_c

        ! convert arguments to C kind
        if (rp == c_rp) then
            vector_c => vector
        else
            allocate(vector_c(size(vector)))
            vector_c = real(vector, kind=c_rp)
        end if

        ! call project C function
        error_c = project_before_wrapping(vector_c)

        ! convert arguments to Fortran kind
        error = int(error_c, kind=ip)
        if (rp /= c_rp) then
            vector = real(vector_c, kind=rp)
            deallocate(vector_c)
        end if

    end subroutine project_f_wrapper

    function conv_check_f_wrapper(error) result(converged)
        !
        ! this function exposes a C-implemented convergence check function to Fortran
        !
        integer(ip), intent(out) :: error
        logical :: converged

        integer(c_ip) :: error_c
        logical(c_bool) :: converged_c

        ! call conv_check C function
        error_c = conv_check_before_wrapping(converged_c)

        ! convert arguments to Fortran kind
        converged = logical(converged_c)
        error = int(error_c, kind=ip)

    end function conv_check_f_wrapper

    subroutine logger_f_wrapper(message)
        !
        ! this subroutine exposes a C-implemented logger function to Fortran
        !
        character(*), intent(in) :: message

        character(kind=c_char), allocatable :: message_c(:)
        integer(ip) :: message_len, i

        ! copy to C character
        message_len = len_trim(message)
        allocate(message_c(0:message_len))
        do i = 1, message_len
            message_c(i - 1) = message(i:i)
        end do
        
        ! append null terminator
        message_c(message_len) = c_null_char 

        ! call logging C function
        call logger_before_wrapping(message_c)
        deallocate(message_c)

    end subroutine logger_f_wrapper

    subroutine init_solver_settings_c(settings) bind(C, name="init_solver_settings")
        !
        ! this subroutine initializes the C solver settings
        !
        type(solver_settings_type_c), intent(inout) :: settings

        settings = default_solver_settings

    end subroutine init_solver_settings_c

    subroutine init_stability_settings_c(settings) &
        bind(C, name="init_stability_settings")
        !
        ! this subroutine initializes the C stability settings
        !
        type(stability_settings_type_c), intent(inout) :: settings

        settings = default_stability_settings

    end subroutine init_stability_settings_c

    subroutine assign_solver_f_c(settings, settings_c)
        !
        ! this subroutine converts solver settings from C to Fortran
        !
        use opentrustregion, only: solver_settings_type

        type(solver_settings_type), intent(out) :: settings
        type(solver_settings_type_c), intent(in) :: settings_c

        if (settings_c%initialized) then
            ! convert callback functions
            if (c_associated(settings_c%precond)) then
                call c_f_procpointer(cptr=settings_c%precond, &
                                     fptr=precond_before_wrapping)
                settings%precond => precond_f_wrapper
            else
                settings%precond => null()
            end if
            if (c_associated(settings_c%project)) then
                call c_f_procpointer(cptr=settings_c%project, &
                                     fptr=project_before_wrapping)
                settings%project => project_f_wrapper
            else
                settings%project => null()
            end if
            if (c_associated(settings_c%conv_check)) then
                call c_f_procpointer(cptr=settings_c%conv_check, &
                                     fptr=conv_check_before_wrapping)
                settings%conv_check => conv_check_f_wrapper
            else
                settings%conv_check => null()
            end if
            if (c_associated(settings_c%logger)) then
                call c_f_procpointer(cptr=settings_c%logger, &
                                     fptr=logger_before_wrapping)
                settings%logger => logger_f_wrapper
            else
                settings%logger => null()
            end if

            ! convert logicals
            settings%stability = logical(settings_c%stability)
            settings%line_search = logical(settings_c%line_search)

            ! convert reals
            settings%conv_tol = real(settings_c%conv_tol, kind=rp)
            settings%start_trust_radius = real(settings_c%start_trust_radius, kind=rp)
            settings%global_red_factor = real(settings_c%global_red_factor, kind=rp)
            settings%local_red_factor = real(settings_c%local_red_factor, kind=rp)

            ! convert integers
            settings%n_random_trial_vectors = int(settings_c%n_random_trial_vectors, &
                                                  kind=ip)
            settings%n_macro = int(settings_c%n_macro, kind=ip)
            settings%n_micro = int(settings_c%n_micro, kind=ip)
            settings%jacobi_davidson_start = int(settings_c%jacobi_davidson_start, &
                                                 kind=ip)
            settings%seed = int(settings_c%seed, kind=ip)
            settings%verbose = int(settings_c%verbose, kind=ip)

            ! convert characters
            settings%subsystem_solver = character_from_c(settings_c%subsystem_solver)

            ! set settings to initialized
            settings%initialized = .true.
        end if

    end subroutine assign_solver_f_c

    subroutine assign_stability_f_c(settings, settings_c)
        !
        ! this subroutine converts stability check settings from C to Fortran
        !
        use opentrustregion, only: stability_settings_type

        type(stability_settings_type), intent(out) :: settings
        type(stability_settings_type_c), intent(in) :: settings_c

        if (settings_c%initialized) then
            ! convert callback functions
            if (c_associated(settings_c%precond)) then
                call c_f_procpointer(cptr=settings_c%precond, &
                                     fptr=precond_before_wrapping)
                settings%precond => precond_f_wrapper
            else
                settings%precond => null()
            end if
            if (c_associated(settings_c%project)) then
                call c_f_procpointer(cptr=settings_c%project, &
                                     fptr=project_before_wrapping)
                settings%project => project_f_wrapper
            else
                settings%project => null()
            end if
            if (c_associated(settings_c%logger)) then
                call c_f_procpointer(cptr=settings_c%logger, &
                                     fptr=logger_before_wrapping)
                settings%logger => logger_f_wrapper
            else
                settings%logger => null()
            end if

            ! convert reals
            settings%conv_tol = real(settings_c%conv_tol, kind=rp)

            ! convert integers
            settings%n_random_trial_vectors = int(settings_c%n_random_trial_vectors, &
                                                  kind=ip)
            settings%n_iter = int(settings_c%n_iter, kind=ip)
            settings%jacobi_davidson_start = int(settings_c%jacobi_davidson_start, &
                                                 kind=ip)
            settings%seed = int(settings_c%seed, kind=ip)
            settings%verbose = int(settings_c%verbose, kind=ip)

            ! convert characters
            settings%diag_solver = character_from_c(settings_c%diag_solver)

            ! set settings to initialized
            settings%initialized = .true.
        end if

    end subroutine assign_stability_f_c

    subroutine assign_solver_c_f(settings_c, settings)
        !
        ! this subroutine converts solver settings from Fortran to C
        !
        use opentrustregion, only: solver_settings_type

        type(solver_settings_type_c), intent(out) :: settings_c
        type(solver_settings_type), intent(in) :: settings

        if (settings%initialized) then
            ! callback functions cannot be converted
            settings_c%precond = c_null_funptr
            settings_c%project = c_null_funptr
            settings_c%conv_check = c_null_funptr
            settings_c%logger = c_null_funptr

            ! convert logicals
            settings_c%stability = logical(settings%stability, kind=c_bool)
            settings_c%line_search = logical(settings%line_search, kind=c_bool)

            ! convert reals
            settings_c%conv_tol = real(settings%conv_tol, kind=c_rp)
            settings_c%start_trust_radius = real(settings%start_trust_radius, kind=c_rp)
            settings_c%global_red_factor = real(settings%global_red_factor, kind=c_rp)
            settings_c%local_red_factor = real(settings%local_red_factor, kind=c_rp)

            ! convert integers
            settings_c%n_random_trial_vectors = int(settings%n_random_trial_vectors, &
                                                    kind=c_ip)
            settings_c%n_macro = int(settings%n_macro, kind=c_ip)
            settings_c%n_micro = int(settings%n_micro, kind=c_ip)
            settings_c%jacobi_davidson_start = int(settings%jacobi_davidson_start, &
                                                   kind=c_ip)
            settings_c%seed = int(settings%seed, kind=c_ip)
            settings_c%verbose = int(settings%verbose, kind=c_ip)

            ! convert characters
            settings_c%subsystem_solver = character_to_c(settings%subsystem_solver)

            ! set settings to initialized
            settings_c%initialized = .true._c_bool
        end if

    end subroutine assign_solver_c_f

    subroutine assign_stability_c_f(settings_c, settings)
        !
        ! this subroutine converts stability settings from Fortran to C
        !
        use opentrustregion, only: stability_settings_type

        type(stability_settings_type_c), intent(out) :: settings_c
        type(stability_settings_type), intent(in) :: settings

        if (settings%initialized) then
            ! callback functions cannot be converted
            settings_c%precond = c_null_funptr
            settings_c%project = c_null_funptr
            settings_c%logger = c_null_funptr

            ! convert reals
            settings_c%conv_tol = real(settings%conv_tol, kind=c_rp)

            ! convert integers
            settings_c%n_random_trial_vectors = int(settings%n_random_trial_vectors, &
                                                    kind=c_ip)
            settings_c%n_iter = int(settings%n_iter, kind=c_ip)
            settings_c%jacobi_davidson_start = int(settings%jacobi_davidson_start, &
                                                   kind=c_ip)
            settings_c%seed = int(settings%seed, kind=c_ip)
            settings_c%verbose = int(settings%verbose, kind=c_ip)

            ! convert characters
            settings_c%diag_solver = character_to_c(settings%diag_solver)

            ! set settings to initialized
            settings_c%initialized = .true._c_bool
        end if

    end subroutine assign_stability_c_f

    function character_from_c(char_c) result(char_f)
        !
        ! this function converts a C null-terminated character array to a Fortran 
        ! character
        !
        character(c_char), intent(in) :: char_c(*)
        character(:), allocatable :: char_f

        integer(ip) :: n

        ! allocate Fortran character
        n = 1
        do while (char_c(n) /= c_null_char)
            n = n + 1
        end do
        allocate(character(n - 1) :: char_f)

        ! copy and convert each character
        char_f = transfer(char_c(1:n - 1), char_f)

    end function character_from_c

    function character_to_c(char_f) result(char_c)
        !
        ! this function converts a Fortran character string to a C null-terminated 
        ! character array
        !
        character(*), intent(in) :: char_f
        character(c_char), allocatable :: char_c(:)

        integer(ip) :: n

        ! allocate C null-terminated character array
        n = len_trim(char_f)
        allocate(char_c(n + 1))

        ! copy and convert each character
        char_c(1:n) = transfer(char_f(1:n), char_c(1:n))
        char_c(n + 1) = c_null_char

    end function character_to_c

end module c_interface
