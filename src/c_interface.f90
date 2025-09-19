! Copyright (C) 2025- Jonas Greiner
!
! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.

module c_interface

    use opentrustregion, only: rp, ip, standard_solver => solver, &
                               standard_stability_check => stability_check
    use, intrinsic :: iso_c_binding, only: c_long, c_double, c_bool, c_ptr, c_funptr, &
                                           c_f_pointer, c_f_procpointer, c_associated, &
                                           c_char, c_null_char

    implicit none

    ! define procedure pointer which will point to the Fortran procedures
    procedure(update_orbs_c_type), pointer :: update_orbs_before_wrapping => null()
    procedure(hess_x_c_type), pointer :: hess_x_before_wrapping => null()
    procedure(obj_func_c_type), pointer :: obj_func_before_wrapping => null()
    procedure(precond_c_type), pointer :: precond_before_wrapping => null()
    procedure(conv_check_c_type), pointer :: conv_check_before_wrapping => null()
    procedure(logger_c_type), pointer :: logger_before_wrapping => null()

    ! C-interoperable interfaces for the callback functions
    abstract interface
        function update_orbs_c_type(kappa_c, func_c, grad_c_ptr, h_diag_c_ptr, &
                                    hess_x_c_funptr) result(error) bind(C)
            use, intrinsic :: iso_c_binding, only: c_double, c_ptr, c_funptr, c_bool

            real(c_double), intent(in) :: kappa_c(*)
            real(c_double), intent(out) :: func_c
            type(c_ptr), intent(out) :: grad_c_ptr, h_diag_c_ptr
            type(c_funptr), intent(out) :: hess_x_c_funptr
            logical(c_bool) :: error
        end function update_orbs_c_type
    end interface

    abstract interface
        function hess_x_c_type(x_c, hess_x_c_ptr) result(error) bind(C)
            use, intrinsic :: iso_c_binding, only: c_double, c_ptr, c_bool

            real(c_double), intent(in) :: x_c(*)
            type(c_ptr), intent(out) :: hess_x_c_ptr
            logical(c_bool) :: error
        end function hess_x_c_type
    end interface

    abstract interface
        function obj_func_c_type(kappa_c, func) result(error) bind(C)
            use, intrinsic :: iso_c_binding, only: c_double, c_bool

            real(c_double), intent(in) :: kappa_c(*)
            real(c_double), intent(out) :: func
            logical(c_bool) :: error
        end function obj_func_c_type
    end interface

    abstract interface
        function precond_c_type(residual_c, mu_c, precond_residual_c_ptr) &
            result(error) bind(C)
            use, intrinsic :: iso_c_binding, only: c_double, c_ptr, c_bool

            real(c_double), intent(in) :: residual_c(*), mu_c
            type(c_ptr), intent(out) :: precond_residual_c_ptr
            logical(c_bool) :: error
        end function precond_c_type
    end interface

    abstract interface
        function conv_check_c_type(converged) result(error) bind(C)
            use, intrinsic :: iso_c_binding, only: c_bool

            logical(c_bool) :: error
            logical(c_bool), intent(out) :: converged
        end function conv_check_c_type
    end interface

    abstract interface
        subroutine logger_c_type(message) bind(C)
            use, intrinsic :: iso_c_binding, only: c_char

            character(kind=c_char), intent(in) :: message(*)
        end subroutine logger_c_type
    end interface

    ! interface for setting default values
    interface set_default_c_ptr
        module procedure set_default_c_ptr_integer, set_default_c_ptr_real, &
            set_default_c_ptr_logical
    end interface

    ! interfaces for solver and stability subroutines, different routines can be
    ! injected, for example for testing
    interface
        subroutine solver_type(update_orbs, obj_func, n_param, error, precond, &
                               conv_check, stability, line_search, davidson, &
                               jacobi_davidson, prefer_jacobi_davidson, conv_tol, &
                               n_random_trial_vectors, start_trust_radius, &
                               n_macro, n_micro, global_red_factor, local_red_factor, &
                               seed, verbose, logger)

            use opentrustregion, only: rp, ip, update_orbs_type, obj_func_type, &
                                       precond_type, conv_check_type, logger_type

            procedure(update_orbs_type), intent(in), pointer :: update_orbs
            procedure(obj_func_type), intent(in), pointer :: obj_func
            integer(ip), intent(in) :: n_param
            logical, intent(out) :: error
            procedure(precond_type), intent(in), pointer, optional :: precond
            procedure(conv_check_type), intent(in), pointer, optional :: conv_check
            logical, intent(in), optional :: stability, line_search, davidson, &
                                             jacobi_davidson, prefer_jacobi_davidson
            real(rp), intent(in), optional :: conv_tol, start_trust_radius, &
                                              global_red_factor, local_red_factor
            integer(ip), intent(in), optional :: n_random_trial_vectors, n_macro, &
                                                 n_micro, seed, verbose
            procedure(logger_type), intent(in), pointer, optional :: logger       

        end subroutine solver_type
    end interface

    interface
        subroutine stability_check_type(h_diag, hess_x, stable, kappa, error, precond, &
                                        jacobi_davidson, conv_tol, &
                                        n_random_trial_vectors, n_iter, verbose, logger)

            use opentrustregion, only: rp, ip, hess_x_type, precond_type, &
                                       conv_check_type, logger_type

            real(rp), intent(in) :: h_diag(:)
            procedure(hess_x_type), intent(in), pointer :: hess_x
            logical, intent(out) :: stable, error
            real(rp), intent(out) :: kappa(:)
            procedure(precond_type), intent(in), pointer, optional :: precond
            logical, intent(in), optional :: jacobi_davidson
            real(rp), intent(in), optional :: conv_tol
            integer(ip), intent(in), optional :: n_random_trial_vectors, n_iter, &
                                                 verbose
            procedure(logger_type), intent(in), pointer, optional :: logger   

        end subroutine stability_check_type
    end interface

    procedure(solver_type), pointer :: solver => standard_solver
    procedure(stability_check_type), pointer :: stability_check => &
        standard_stability_check

contains

    function solver_c_wrapper(update_orbs_c_funptr, obj_func_c_funptr, n_param_c, &
                              precond_c_funptr, conv_check_c_funptr, stability_c_ptr, &
                              line_search_c_ptr, davidson_c_ptr, &
                              jacobi_davidson_c_ptr, prefer_jacobi_davidson_c_ptr, &
                              conv_tol_c_ptr, n_random_trial_vectors_c_ptr, &
                              start_trust_radius_c_ptr, n_macro_c_ptr, n_micro_c_ptr, &
                              global_red_factor_c_ptr, local_red_factor_c_ptr, &
                              seed_c_ptr, verbose_c_ptr, logger_c_funptr) &
        result(error_c) bind(C, name="solver")
        !
        ! this subroutine wraps the solver subroutine to convert C variables to Fortran
        ! variables
        !
        use opentrustregion, only: solver_stability_default, &
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
                                   

        type(c_funptr), intent(in), value :: update_orbs_c_funptr, obj_func_c_funptr
        integer(c_long), intent(in), value :: n_param_c
        type(c_funptr), intent(in), value :: precond_c_funptr, conv_check_c_funptr, &
                                             logger_c_funptr
        type(c_ptr), intent(in), value :: stability_c_ptr, line_search_c_ptr, &
                                          davidson_c_ptr, jacobi_davidson_c_ptr, &
                                          prefer_jacobi_davidson_c_ptr, &
                                          conv_tol_c_ptr, &
                                          n_random_trial_vectors_c_ptr, &
                                          start_trust_radius_c_ptr, n_macro_c_ptr, &
                                          n_micro_c_ptr, global_red_factor_c_ptr, &
                                          local_red_factor_c_ptr, seed_c_ptr, &
                                          verbose_c_ptr
        logical(c_bool) :: error_c

        logical :: error, stability, line_search, davidson, jacobi_davidson, &
                   prefer_jacobi_davidson
        real(rp) :: conv_tol, start_trust_radius, global_red_factor, local_red_factor
        integer(ip) :: n_param, n_random_trial_vectors, n_macro, n_micro, seed, verbose
        procedure(update_orbs_c_wrapper), pointer :: update_orbs
        procedure(obj_func_c_wrapper), pointer :: obj_func
        procedure(precond_c_wrapper), pointer :: precond
        procedure(conv_check_c_wrapper), pointer :: conv_check
        procedure(logger_c_wrapper), pointer :: logger

        ! set optional arguments to default values
        stability = set_default_c_ptr(stability_c_ptr, solver_stability_default)
        line_search = set_default_c_ptr(line_search_c_ptr, solver_line_search_default)
        davidson = set_default_c_ptr(davidson_c_ptr, solver_davidson_default)
        jacobi_davidson = set_default_c_ptr(jacobi_davidson_c_ptr, &
                                            solver_jacobi_davidson_default)
        prefer_jacobi_davidson = set_default_c_ptr(prefer_jacobi_davidson_c_ptr, &
                                                  solver_prefer_jacobi_davidson_default)
        conv_tol = set_default_c_ptr(conv_tol_c_ptr, solver_conv_tol_default)
        n_random_trial_vectors = set_default_c_ptr(n_random_trial_vectors_c_ptr, &
                                                  solver_n_random_trial_vectors_default)
        start_trust_radius = set_default_c_ptr(start_trust_radius_c_ptr, &
                                               solver_start_trust_radius_default)
        n_macro = set_default_c_ptr(n_macro_c_ptr, solver_n_macro_default)
        n_micro = set_default_c_ptr(n_micro_c_ptr, solver_n_micro_default)
        global_red_factor = set_default_c_ptr(global_red_factor_c_ptr, &
                                              solver_global_red_factor_default)
        local_red_factor = set_default_c_ptr(local_red_factor_c_ptr, &
                                             solver_local_red_factor_default)
        seed = set_default_c_ptr(seed_c_ptr, solver_seed_default)
        verbose = set_default_c_ptr(verbose_c_ptr, solver_verbose_default)

        ! associate the input C pointer to update_orbs subroutine to a Fortran
        ! procedure pointer
        call c_f_procpointer(cptr=update_orbs_c_funptr, &
                             fptr=update_orbs_before_wrapping)
        call c_f_procpointer(cptr=obj_func_c_funptr, fptr=obj_func_before_wrapping)

        ! associate procedure pointer to wrapper function
        update_orbs => update_orbs_c_wrapper
        obj_func => obj_func_c_wrapper

        ! convert dummy argument to Fortran kind
        n_param = int(n_param_c, kind=ip)

        ! call solver with or without optional procedures
        if (c_associated(precond_c_funptr)) then
            call c_f_procpointer(cptr=precond_c_funptr, fptr=precond_before_wrapping)
            precond => precond_c_wrapper
        else
            precond => null()
        end if
        if (c_associated(conv_check_c_funptr)) then
            call c_f_procpointer(cptr=conv_check_c_funptr, &
                                 fptr=conv_check_before_wrapping)
            conv_check => conv_check_c_wrapper
        else
            conv_check => null()
        end if
        if (c_associated(logger_c_funptr)) then
            call c_f_procpointer(cptr=logger_c_funptr, fptr=logger_before_wrapping)
            logger => logger_c_wrapper
        else
            logger => null()
        end if
        call solver(update_orbs, obj_func, n_param, error, precond, conv_check, &
                    stability, line_search, davidson, jacobi_davidson, &
                    prefer_jacobi_davidson, conv_tol, n_random_trial_vectors, &
                    start_trust_radius, n_macro, n_micro, global_red_factor, &
                    local_red_factor, seed, verbose, logger)

        ! convert return arguments to C kind
        error_c = error

    end function solver_c_wrapper

    function stability_check_c_wrapper(h_diag_c, hess_x_c_funptr, n_param_c, stable_c, &
                                       kappa_c, precond_c_funptr, &
                                       jacobi_davidson_c_ptr, conv_tol_c_ptr, &
                                       n_random_trial_vectors_c_ptr, n_iter_c_ptr, &
                                       verbose_c_ptr, logger_c_funptr) result(error_c) &
        bind(C, name="stability_check")
        !
        ! this subroutine wraps the stability check subroutine to convert C variables
        ! to Fortran variables
        !
        use opentrustregion, only: stability_jacobi_davidson_default, &
                                   stability_conv_tol_default, &
                                   stability_n_random_trial_vectors_default, &
                                   stability_n_iter_default, &
                                   stability_verbose_default

        integer(c_long), intent(in), value :: n_param_c
        real(c_double), intent(in), dimension(n_param_c) :: h_diag_c
        type(c_funptr), intent(in), value :: hess_x_c_funptr
        logical(c_bool), intent(out) :: stable_c
        real(c_double), intent(out) :: kappa_c(n_param_c)
        type(c_funptr), intent(in), value :: precond_c_funptr, logger_c_funptr
        type(c_ptr), intent(in), value :: jacobi_davidson_c_ptr, conv_tol_c_ptr, &
                                          n_random_trial_vectors_c_ptr, n_iter_c_ptr, &
                                          verbose_c_ptr
        logical(c_bool) :: error_c

        real(rp) :: conv_tol, h_diag(n_param_c)
        real(rp) :: kappa(n_param_c)
        logical :: stable, error, jacobi_davidson
        integer(ip) :: n_param, n_random_trial_vectors, n_iter, verbose
        procedure(hess_x_c_wrapper), pointer :: hess_x
        procedure(precond_c_wrapper), pointer :: precond
        procedure(logger_c_wrapper), pointer :: logger

        ! set optional arguments to default values
        jacobi_davidson = set_default_c_ptr(jacobi_davidson_c_ptr, &
                                            stability_jacobi_davidson_default)
        conv_tol = set_default_c_ptr(conv_tol_c_ptr, stability_conv_tol_default)
        n_random_trial_vectors = set_default_c_ptr(n_random_trial_vectors_c_ptr, &
                                               stability_n_random_trial_vectors_default)
        n_iter = set_default_c_ptr(n_iter_c_ptr, stability_n_iter_default)
        verbose = set_default_c_ptr(verbose_c_ptr, stability_verbose_default)

        ! associate the input C pointer to update_orbs subroutine to a Fortran
        ! procedure pointer
        call c_f_procpointer(cptr=hess_x_c_funptr, fptr=hess_x_before_wrapping)

        ! convert dummy argument to Fortran kind
        n_param = int(n_param_c, kind=ip)

        ! associate procedure pointer to wrapper function
        hess_x => hess_x_c_wrapper

        ! convert dummy argument to Fortran kind
        h_diag = h_diag_c

        ! call stability check with or without preconditioner
        if (c_associated(precond_c_funptr)) then
            call c_f_procpointer(cptr=precond_c_funptr, fptr=precond_before_wrapping)
            precond => precond_c_wrapper
        else
            precond => null()
        end if
        if (c_associated(logger_c_funptr)) then
            call c_f_procpointer(cptr=logger_c_funptr, fptr=logger_before_wrapping)
            logger => logger_c_wrapper
        else
            logger => null()
        end if
        call stability_check(h_diag, hess_x, stable, kappa, error, precond, &
                             jacobi_davidson, conv_tol, n_random_trial_vectors, &
                             n_iter, verbose, logger)

        ! convert return arguments to C kind
        stable_c = stable
        kappa_c = kappa
        error_c = error

    end function stability_check_c_wrapper

    subroutine update_orbs_c_wrapper(kappa, func, grad, h_diag, hess_x, error)
        !
        ! this subroutine wraps the orbital update subroutine to convert C variables to
        ! Fortran variables
        !
        use opentrustregion, only: hess_x_type

        real(rp), intent(in) :: kappa(:)
        real(rp), intent(out) :: func, grad(:), h_diag(:)
        procedure(hess_x_type), intent(out), pointer :: hess_x
        logical, intent(out) :: error

        real(c_double) :: func_c
        real(c_double) :: kappa_c(size(kappa))
        type(c_ptr) :: grad_c_ptr, h_diag_c_ptr
        type(c_funptr) :: hess_x_c_funptr
        real(c_double), pointer :: grad_ptr(:), h_diag_ptr(:)
        logical(c_bool) :: error_c

        ! convert dummy argument to C kind
        kappa_c = kappa

        ! call update_orbs C function
        error_c = update_orbs_before_wrapping(kappa_c, func_c, grad_c_ptr, &
                                              h_diag_c_ptr, hess_x_c_funptr)

        ! associate C pointer to arrays with Fortran pointer
        call c_f_pointer(cptr=grad_c_ptr, fptr=grad_ptr, shape=[size(kappa)])
        call c_f_pointer(cptr=h_diag_c_ptr, fptr=h_diag_ptr, shape=[size(kappa)])

        ! associate the input C pointer to hess_x function to a Fortran procedure
        ! pointer
        call c_f_procpointer(cptr=hess_x_c_funptr, fptr=hess_x_before_wrapping)

        ! associate procedure pointer to wrapper function
        hess_x => hess_x_c_wrapper

        ! convert arguments to Fortran kind
        func = func_c
        grad = grad_ptr
        h_diag = h_diag_ptr
        error = error_c

    end subroutine update_orbs_c_wrapper

    function hess_x_c_wrapper(x, error) result(hess_x)
        !
        ! this function wraps the Hessian linear transformation to convert C variables
        ! to Fortran variables
        !
        real(rp), intent(in) :: x(:)
        logical, intent(out) :: error
        real(rp) :: hess_x(size(x))

        real(c_double) :: x_c(size(x))
        type(c_ptr) :: hess_x_c_ptr
        real(c_double), pointer :: hess_x_ptr(:)
        logical(c_bool) :: error_c

        ! convert trial vector to C kind
        x_c = x

        ! call C function
        error_c = hess_x_before_wrapping(x_c, hess_x_c_ptr)

        ! associate C pointer to arrays with Fortran pointer
        call c_f_pointer(cptr=hess_x_c_ptr, fptr=hess_x_ptr, shape=[size(x)])

        ! convert arguments to Fortran kind
        hess_x = hess_x_ptr
        error = error_c

    end function hess_x_c_wrapper

    function obj_func_c_wrapper(kappa, error) result(obj_func)
        !
        ! this function wraps the objective function to convert C variables to Fortran
        ! variables
        !
        real(rp), intent(in) :: kappa(:)
        logical, intent(out) :: error
        real(rp) :: obj_func

        real(c_double) :: kappa_c(size(kappa)), func_c
        logical(c_bool) :: error_c

        ! convert dummy argument to C kind
        kappa_c = kappa

        ! call obj_func C function
        error_c = obj_func_before_wrapping(kappa_c, func_c)

        ! convert arguments to Fortran kind
        obj_func = func_c
        error = error_c

    end function obj_func_c_wrapper

    function precond_c_wrapper(residual, mu, error) result(precond_residual)
        !
        ! this function wraps the preconditioner function to convert C variables to 
        ! Fortran variables
        !
        real(rp), intent(in) :: residual(:), mu
        logical, intent(out) :: error
        real(rp) :: precond_residual(size(residual))

        real(c_double) :: residual_c(size(residual)), mu_c
        type(c_ptr) :: precond_residual_c_ptr
        real(c_double), pointer :: precond_residual_ptr(:)
        logical(c_bool) :: error_c

        ! convert dummy argument to C kind
        residual_c = residual
        mu_c = mu

        ! call precond C function
        error_c = precond_before_wrapping(residual_c, mu_c, precond_residual_c_ptr)

        ! associate C pointer to arrays with Fortran pointer
        call c_f_pointer(cptr=precond_residual_c_ptr, fptr=precond_residual_ptr, &
                         shape=[size(residual)])

        ! convert arguments to Fortran kind
        precond_residual = precond_residual_ptr
        error = error_c

    end function precond_c_wrapper

    function conv_check_c_wrapper(error) result(converged)
        !
        ! this function wraps the convergence check function to convert C variables to 
        ! Fortran variables
        !
        logical, intent(out) :: error
        logical :: converged

        logical(c_bool) :: error_c, converged_c

        ! call conv_check C function
        error_c = conv_check_before_wrapping(converged_c)

        ! convert arguments to Fortran kind
        converged = converged_c
        error = error_c

    end function conv_check_c_wrapper

    subroutine logger_c_wrapper(message)
        !
        ! this function wraps the logging subroutine to convert C variables to Fortran 
        ! variables
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

    end subroutine logger_c_wrapper

    function set_default_c_ptr_integer(optional_value, default_value) result(variable)
        !
        ! this function sets a default value for reals
        !
        integer(ip), intent(in)  :: default_value
        type(c_ptr), intent(in) :: optional_value

        integer(ip) :: variable

        integer(c_long), pointer :: variable_ptr

        if (c_associated(optional_value)) then
            call c_f_pointer(cptr=optional_value, fptr=variable_ptr)
            variable = int(variable_ptr, kind=ip)
        else
            variable = default_value
        end if

    end function set_default_c_ptr_integer

    function set_default_c_ptr_real(optional_value, default_value) result(variable)
        !
        ! this function sets a default value for reals
        !
        real(rp), intent(in)  :: default_value
        type(c_ptr), intent(in) :: optional_value

        real(rp) :: variable

        real(c_double), pointer :: variable_ptr

        if (c_associated(optional_value)) then
            call c_f_pointer(cptr=optional_value, fptr=variable_ptr)
            variable = variable_ptr
        else
            variable = default_value
        end if

    end function set_default_c_ptr_real

    function set_default_c_ptr_logical(optional_value, default_value) result(variable)
        !
        ! this function sets a default value for logicals
        !
        logical, intent(in)  :: default_value
        type(c_ptr), intent(in) :: optional_value

        logical :: variable

        logical(c_bool), pointer :: variable_ptr

        if (c_associated(optional_value)) then
            call c_f_pointer(cptr=optional_value, fptr=variable_ptr)
            variable = variable_ptr
        else
            variable = default_value
        end if

    end function set_default_c_ptr_logical

end module c_interface
