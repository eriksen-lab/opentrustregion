! Copyright (C) 2025- Jonas Greiner
!
! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.

module opentrustregion_system_tests

    use opentrustregion, only: rp, ip, stderr
    use, intrinsic :: iso_c_binding, only: c_bool, c_ptr, c_char, c_f_pointer, &
                                           c_null_char

    implicit none

    ! dimensions of the Foster-Boys localization of the occupied orbitals of water
    integer(ip), parameter :: n_ao = 13, n_mo = 5, n_param = n_mo * (n_mo - 1) / 2

    ! directory of the test data, set by the test driver
    character(len=:), allocatable :: data_dir

    ! define type for the host context handed to the callback functions of the
    ! Foster-Boys localization, which holds the current orbitals, the integrals and
    ! the intermediates of the Hessian linear transformation
    type :: fb_context_type
        real(rp) :: mo_coeff(n_ao, n_mo) = 0.0_rp, r_ao_ints(3, n_ao, n_ao) = 0.0_rp, &
                    r2_ao_ints(n_ao, n_ao) = 0.0_rp, &
                    r_mo_ints(3, n_mo, n_mo) = 0.0_rp, &
                    rii_rij_rjj_rji(n_mo, n_mo) = 0.0_rp
    end type

    ! error code returned by callback functions reached without their context
    integer(ip), parameter :: missing_context_error = 42

contains

    subroutine set_test_data_path(path) bind(C)
        !
        ! this subroutine sets the path to test data
        !
        type(c_ptr), intent(in), value :: path
        character(kind=c_char), pointer :: c_path(:)
        integer(ip) :: len, i

        ! conda paths need ample space
        call c_f_pointer(path, c_path, [1024])
        len = 0
        do i = 1, size(c_path)
            if (c_path(i) == c_null_char) exit
            len = len + 1
        end do
        allocate(character(len=len) :: data_dir)
        data_dir = transfer(c_path(1:len), data_dir)

    end subroutine set_test_data_path

    function exp_asymm_mat(mat)
        !
        ! this function calculates the matrix exponential of an asymmetric matrix
        !
        real(rp), intent(in) :: mat(:, :)
        real(rp) :: exp_asymm_mat(size(mat, 1), size(mat, 2))

        integer(ip) :: n, lwork, info, i
        real(rp), allocatable :: eigvals(:), rwork(:)
        complex(rp), allocatable :: work(:), eigvecs(:, :), tmp(:, :)

        external :: zheev

        ! size of matrix
        n = size(mat, 1)

        ! convert to Hermitian matrix
        eigvecs = cmplx(0.0_rp, mat, kind=rp)

        ! query optimal workspace size
        lwork = -1
        allocate(eigvals(n), work(1), rwork(3 * n - 2))
        call zheev("V", "U", n, eigvecs, n, eigvals, work, lwork, rwork, info)
        lwork = int(work(1))
        deallocate(work)
        allocate(work(lwork))

        ! perform eigendecomposition
        call zheev("V", "U", n, eigvecs, n, eigvals, work, lwork, rwork, info)
        deallocate(work, rwork)

        ! compute matrix exponential under assumption that eigenvalues are purely
        ! imaginary
        allocate(tmp(n, n))
        do i = 1, n
            tmp(:, i) = eigvecs(:, i) * cmplx(cos(eigvals(i)), sin(eigvals(i)), kind=rp)
        end do
        exp_asymm_mat = real(transpose(matmul(tmp, conjg(transpose(eigvecs)))))
        deallocate(eigvecs, eigvals, tmp)

    end function exp_asymm_mat

    function resolve_fb_context(context, error) result(state)
        !
        ! this function returns the Foster-Boys context handed to a callback function
        ! and an error if no such context was handed over
        !
        class(*), intent(in), pointer :: context
        integer(ip), intent(out) :: error
        type(fb_context_type), pointer :: state

        ! initialize error flag
        error = 0

        ! get context
        state => null()
        if (associated(context)) then
            select type (context)
            type is (fb_context_type)
                state => context
            end select
        end if

        ! report missing context
        if (.not. associated(state)) error = missing_context_error

    end function resolve_fb_context

    subroutine hess_x_fun(x, hess_x, error, context)
        !
        ! this subroutine performs the Hessian linear transformation for Foster-Boys
        ! orbital localization, it cannot be defined within update_orbs as it would
        ! otherwise go out of scope when that subroutine returns
        !
        real(rp), intent(in), target :: x(:)
        real(rp), intent(out), target :: hess_x(:)
        integer(ip), intent(out) :: error
        class(*), intent(in), pointer :: context

        type(fb_context_type), pointer :: state
        real(rp), allocatable :: x_full(:, :), hess_x_full(:, :), tmp2(:, :), &
                                 tmp3(:, :, :)
        integer(ip) :: xyz1, i1, j1, idx1

        ! get Foster-Boys state
        state => resolve_fb_context(context, error)
        if (error /= 0) return

        ! unpack trial vector
        allocate(x_full(n_mo, n_mo))
        x_full = 0.0_rp
        idx1 = 1
        do i1 = 2, n_mo
            do j1 = 1, i1 - 1
                x_full(i1, j1) = x(idx1)
                x_full(j1, i1) = -x(idx1)
                idx1 = idx1 + 1
            end do
        end do

        ! construct intermediates
        allocate(tmp2(3, n_mo), tmp3(3, n_mo, n_mo))
        do xyz1 = 1, 3
            do i1 = 1, n_mo
                tmp2(xyz1, i1) = sum(x_full(:, i1) * state%r_mo_ints(xyz1, i1, :))
                do j1 = 1, n_mo
                    tmp3(xyz1, i1, j1) = &
                        sum(x_full(:, j1) * state%r_mo_ints(xyz1, i1, :))
                end do
            end do
        end do

        ! construct Hessian linear transformation
        allocate(hess_x_full(n_mo, n_mo))
        hess_x_full = matmul(transpose(x_full), transpose(state%rii_rij_rjj_rji))
        do i1 = 1, n_mo
            do j1 = 1, n_mo
                hess_x_full(i1, j1) = hess_x_full(i1, j1) + 2 * sum( &
                    state%r_mo_ints(:, j1, i1) * tmp2(:, j1) - &
                    state%r_mo_ints(:, i1, i1) * tmp3(:, j1, i1) - &
                    state%r_mo_ints(:, j1, i1) * tmp2(:, i1))
            end do
        end do
        deallocate(x_full, tmp2, tmp3)

        ! extract lower triagonal
        idx1 = 1
        do i1 = 2, n_mo
            do j1 = 1, i1 - 1
                hess_x(idx1) = -(hess_x_full(i1, j1) - hess_x_full(j1, i1))
                idx1 = idx1 + 1
            end do
        end do
        deallocate(hess_x_full)

    end subroutine hess_x_fun

    real(rp) function obj_func(kappa, error, context)
        !
        ! this function calculates the Foster-Boys orbital localization objective
        ! function
        !
        real(rp), intent(in), target :: kappa(:)
        integer(ip), intent(out) :: error
        class(*), intent(in), pointer :: context

        type(fb_context_type), pointer :: state
        real(rp), allocatable :: kappa_full(:, :), mo_coeff_tmp(:, :)
        integer(ip) :: xyz, i, j, idx

        ! initialize result in case the Foster-Boys state is missing
        obj_func = 0.0_rp

        ! get Foster-Boys state
        state => resolve_fb_context(context, error)
        if (error /= 0) return

        ! unpack orbital rotation
        allocate(kappa_full(n_mo, n_mo))
        kappa_full = 0.0_rp
        idx = 1
        do i = 2, n_mo
            do j = 1, i - 1
                kappa_full(i, j) = kappa(idx)
                kappa_full(j, i) = -kappa(idx)
                idx = idx + 1
            end do
        end do

        ! rotate orbitals
        mo_coeff_tmp = matmul(state%mo_coeff, exp_asymm_mat(kappa_full))
        deallocate(kappa_full)

        ! compute cost function
        obj_func = 0.0_rp
        do i = 1, n_mo
            obj_func = obj_func + dot_product( &
                mo_coeff_tmp(:, i), matmul(state%r2_ao_ints, mo_coeff_tmp(:, i)))
            do xyz = 1, 3
                obj_func = obj_func - dot_product(mo_coeff_tmp(:, i), matmul( &
                    state%r_ao_ints(xyz, :, :), mo_coeff_tmp(:, i)))**2
            end do
        end do
        deallocate(mo_coeff_tmp)

    end function obj_func

    subroutine update_orbs(kappa, func, grad, h_diag, hess_x_funptr, error, context)
        !
        ! this subroutine updates the orbitals for Foster-Boys orbital localization
        !
        use opentrustregion, only: hess_x_type

        real(rp), intent(in), target :: kappa(:)
        real(rp), intent(out) :: func
        real(rp), intent(out), target :: grad(:), h_diag(:)
        procedure(hess_x_type), intent(inout), pointer :: hess_x_funptr
        integer(ip), intent(out) :: error
        class(*), intent(in), pointer :: context

        type(fb_context_type), pointer :: state
        integer(ip) :: xyz, i, j, idx
        real(rp), allocatable :: kappa_full(:, :), h_diag_tmp(:, :), tmp1(:, :)

        ! get Foster-Boys state
        state => resolve_fb_context(context, error)
        if (error /= 0) return

        ! unpack orbital rotation
        allocate(kappa_full(n_mo, n_mo))
        kappa_full = 0.0_rp
        idx = 1
        do i = 2, n_mo
            do j = 1, i - 1
                kappa_full(i, j) = kappa(idx)
                kappa_full(j, i) = -kappa(idx)
                idx = idx + 1
            end do
        end do

        ! rotate orbitals
        state%mo_coeff = matmul(state%mo_coeff, exp_asymm_mat(kappa_full))
        deallocate(kappa_full)

        ! transform integrals to MO basis
        do xyz = 1, 3
            state%r_mo_ints(xyz, :, :) = matmul(matmul( &
                transpose(state%mo_coeff), state%r_ao_ints(xyz, :, :)), state%mo_coeff)
        end do

        ! compute cost function
        func = 0.0_rp
        do i = 1, n_mo
            func = func + dot_product(state%mo_coeff(:, i), matmul( &
                state%r2_ao_ints, state%mo_coeff(:, i))) - &
                   sum(state%r_mo_ints(:, i, i)**2)
        end do

        ! construct temporary intermediate
        allocate(tmp1(n_mo, n_mo))
        do i = 1, n_mo
            do j = 1, n_mo
                tmp1(i, j) = sum(state%r_mo_ints(:, j, j) * state%r_mo_ints(:, j, i))
            end do
        end do

        ! construct gradient and extract lower triagonal
        idx = 1
        do i = 2, n_mo
            do j = 1, i - 1
                grad(idx) = -2 * (tmp1(i, j) - tmp1(j, i))
                idx = idx + 1
            end do
        end do

        ! construct Hessian diagonal
        allocate(h_diag_tmp(n_mo, n_mo))
        do i = 1, n_mo
            do j = 1, n_mo
                h_diag_tmp(i, j) = &
                    2 * sum(state%r_mo_ints(:, j, j) * state%r_mo_ints(:, i, i) + &
                            state%r_mo_ints(:, j, i) * state%r_mo_ints(:, j, i) + &
                            state%r_mo_ints(:, j, i) * state%r_mo_ints(:, i, j)) - &
                    tmp1(i, i) - tmp1(j, j)
            end do
        end do

        ! extract lower triagonal
        idx = 1
        do i = 2, n_mo
            do j = 1, i - 1
                h_diag(idx) = -2 * (h_diag_tmp(i, j))
                idx = idx + 1
            end do
        end do
        deallocate(h_diag_tmp)

        ! save for Hessian linear transformation
        state%rii_rij_rjj_rji = tmp1 + transpose(tmp1)
        deallocate(tmp1)

        ! get function pointer to Hessian linear transformation
        hess_x_funptr => hess_x_fun

    end subroutine update_orbs

    subroutine setup_fb_context(context, start_file)
        !
        ! this subroutine reads the starting orbitals and the integrals of the
        ! Foster-Boys localization of water into a context
        !
        type(fb_context_type), intent(out) :: context
        character(len=*), intent(in) :: start_file

        integer :: ios

        ! check if test data variable is set
        if (.not. allocated(data_dir)) error stop "Test data directory not set "// &
            "through set_test_data_path subroutine before calling system test."

        ! read raw binary data
        open(unit=10, file=data_dir//"/"//start_file, form="unformatted", &
             access="stream", status="old", action="read", iostat=ios)
        if (ios /= 0) error stop "Error opening file"
        read(10) context%mo_coeff
        close(10)
        open(unit=10, file=data_dir//"/h2o_r_ints.bin", form="unformatted", &
             access="stream", status="old", action="read", iostat=ios)
        if (ios /= 0) error stop "Error opening file"
        read(10) context%r_ao_ints
        close(10)
        open(unit=10, file=data_dir//"/h2o_r2_ints.bin", form="unformatted", &
             access="stream", status="old", action="read", iostat=ios)
        if (ios /= 0) error stop "Error opening file"
        read(10) context%r2_ao_ints
        close(10)

    end subroutine setup_fb_context

    logical function check_h2o_fb_solver(test_name, option)
        !
        ! this function performs the Foster-Boys localization of the occupied orbitals
        ! of water from a guess and from a saddle point with a solver option switched
        ! on and checks that the localized orbitals are a stable minimum with the
        ! reference objective function value
        !
        use opentrustregion, only: update_orbs_type, obj_func_type, &
                                   solver_settings_type, solver, hess_x_type, &
                                   stability_settings_type, stability_check, &
                                   subsystem_solver_options

        character(len=*), intent(in) :: test_name, option

        type(fb_context_type), target :: context
        procedure(update_orbs_type), pointer :: update_orbs_funptr
        procedure(obj_func_type), pointer :: obj_func_funptr
        procedure(hess_x_type), pointer :: hess_x_funptr
        type(solver_settings_type) :: solver_settings
        type(stability_settings_type) :: stability_settings
        class(*), pointer :: context_ptr
        integer(ip) :: error, i_start
        real(rp) :: kappa(n_param), grad(n_param), h_diag(n_param), func
        logical :: stable
        character(len=*), parameter :: start_names(2) = &
            [character(len=12) :: "guess", "saddle point"]
        character(len=*), parameter :: start_files(2) = &
            ["h2o_atomic_mo_coeff.bin", "h2o_saddle_mo_coeff.bin"]
        real(rp), parameter :: ref_func = 6.890557872085_rp, func_tol = 1e-8_rp

        ! assume test passes
        check_h2o_fb_solver = .true.

        ! set function pointers
        update_orbs_funptr => update_orbs
        obj_func_funptr => obj_func

        do i_start = 1, size(start_names)
            ! read starting orbitals and integrals
            call setup_fb_context(context, start_files(i_start))

            ! initialize settings and switch on the option, which is either a
            ! subsystem solver option or a logical setting, the Jacobi-Davidson method
            ! is started immediately since the Davidson method converges before
            ! switching on a problem of this size
            call solver_settings%init(error)
            solver_settings%context => context
            select case (option)
            case ("default")
            case ("line search")
                solver_settings%line_search = .true.
            case ("stability")
                solver_settings%stability = .true.
            case default
                if (.not. any(option == subsystem_solver_options)) then
                    write(stderr, *) test_name//" failed: Unknown solver option."
                    check_h2o_fb_solver = .false.
                    return
                end if
                solver_settings%subsystem_solver = option
                solver_settings%jacobi_davidson_start = 0
            end select

            ! run solver
            call solver(update_orbs_funptr, obj_func_funptr, n_param, error, &
                        solver_settings)
            if (error /= 0) then
                write(stderr, *) test_name//" failed: Solver produced error from "// &
                    "the "//trim(start_names(i_start))//"."
                check_h2o_fb_solver = .false.
                cycle
            end if

            ! get objective function, gradient, Hessian diagonal and Hessian linear
            ! transformation at the localized orbitals
            kappa = 0.0_rp
            context_ptr => context
            call update_orbs(kappa, func, grad, h_diag, hess_x_funptr, error, &
                             context_ptr)
            if (norm2(grad) / sqrt(real(n_param, kind=rp)) > solver_settings%conv_tol) &
                then
                write(stderr, *) test_name//" failed: Gradient not converged from "// &
                    "the "//trim(start_names(i_start))//"."
                check_h2o_fb_solver = .false.
            end if
            if (abs(func - ref_func) > func_tol) then
                write(stderr, *) test_name//" failed: Objective function does not "// &
                    "reach the minimum from the "//trim(start_names(i_start))//"."
                check_h2o_fb_solver = .false.
            end if

            ! check that the localized orbitals are a minimum
            call stability_settings%init(error)
            stability_settings%context => context
            call stability_check(h_diag, hess_x_funptr, stable, error, &
                                 stability_settings)
            if (error /= 0) then
                write(stderr, *) test_name//" failed: Stability check of the "// &
                    "localized orbitals produced error from the "// &
                    trim(start_names(i_start))//"."
                check_h2o_fb_solver = .false.
            end if
            if (.not. stable) then
                write(stderr, *) test_name//" failed: Localized orbitals are not a "// &
                    "stable minimum from the "//trim(start_names(i_start))//"."
                check_h2o_fb_solver = .false.
            end if
        end do

    end function check_h2o_fb_solver

    logical function check_h2o_fb_stability_check(test_name, diag_solver)
        !
        ! this function performs the stability check for the Foster-Boys localization
        ! of the occupied orbitals of water at the minimum and at a saddle point with
        ! a diagonalization solver and checks that the minimum is found to be stable
        ! and the saddle point to be unstable with a returned direction along the
        ! eigenvector of the lowest eigenvalue of the Hessian
        !
        use opentrustregion, only: hess_x_type, stability_settings_type, stability_check

        character(len=*), intent(in) :: test_name, diag_solver

        type(fb_context_type), target :: context
        procedure(hess_x_type), pointer :: hess_x_funptr
        type(stability_settings_type) :: settings
        class(*), pointer :: context_ptr
        integer(ip) :: error, i_point, i, info
        real(rp) :: kappa(n_param), grad(n_param), h_diag(n_param), func, &
                    unit_vector(n_param), hess(n_param, n_param), eigvals(n_param), &
                    work(3 * n_param)
        logical :: stable
        real(rp), parameter :: direction_tol = 1e-6_rp
        character(len=*), parameter :: point_names(2) = &
            [character(len=12) :: "minimum", "saddle point"]
        character(len=*), parameter :: point_files(2) = &
            [character(len=24) :: "h2o_minimum_mo_coeff.bin", "h2o_saddle_mo_coeff.bin"]

        external :: dsyev

        ! assume test passes
        check_h2o_fb_stability_check = .true.

        do i_point = 1, size(point_names)
            ! read orbitals and integrals and get objective function, Hessian diagonal
            ! and Hessian linear transformation at the orbitals
            call setup_fb_context(context, trim(point_files(i_point)))
            context_ptr => context
            kappa = 0.0_rp
            call update_orbs(kappa, func, grad, h_diag, hess_x_funptr, error, &
                             context_ptr)

            ! initialize settings and select the diagonalization solver, the
            ! Jacobi-Davidson method is started immediately since the Davidson method
            ! converges before switching on a problem of this size
            call settings%init(error)
            settings%context => context
            settings%diag_solver = diag_solver
            if (diag_solver == "jacobi-davidson") settings%jacobi_davidson_start = 0

            ! perform stability check
            call stability_check(h_diag, hess_x_funptr, stable, error, settings, kappa)
            if (error /= 0) then
                write(stderr, *) test_name//" failed: Stability check produced "// &
                    "error at the "//trim(point_names(i_point))//"."
                check_h2o_fb_stability_check = .false.
                cycle
            end if

            ! the minimum has to be found to be stable
            if (i_point == 1) then
                if (.not. stable) then
                    write(stderr, *) test_name// &
                        " failed: Minimum not found to be stable."
                    check_h2o_fb_stability_check = .false.
                end if
                cycle
            end if

            ! the saddle point has to be found to be unstable
            if (stable) then
                write(stderr, *) test_name// &
                    " failed: Saddle point not found to be unstable."
                check_h2o_fb_stability_check = .false.
            end if

            ! the returned direction has to be the eigenvector of the lowest eigenvalue
            ! of the Hessian, which is assembled from its linear transformations of the
            ! unit vectors and diagonalized, since the saddle point has several
            ! unstable modes and the objective function decreases along almost any
            ! direction from it
            do i = 1, n_param
                unit_vector = 0.0_rp
                unit_vector(i) = 1.0_rp
                call hess_x_funptr(unit_vector, hess(:, i), error, context_ptr)
                if (error /= 0) exit
            end do
            if (error /= 0) then
                write(stderr, *) test_name//" failed: Hessian linear "// &
                    "transformation produced error at the saddle point."
                check_h2o_fb_stability_check = .false.
                cycle
            end if
            call dsyev("V", "U", n_param, hess, n_param, eigvals, work, &
                       size(work, kind=ip), info)
            if (info /= 0) then
                write(stderr, *) test_name//" failed: Diagonalization of the "// &
                    "Hessian at the saddle point failed."
                check_h2o_fb_stability_check = .false.
                cycle
            end if
            if (abs(abs(dot_product(kappa, hess(:, 1))) - 1.0_rp) > direction_tol) then
                write(stderr, *) test_name//" failed: Returned direction at the "// &
                    "saddle point is not the lowest eigenvector of the Hessian."
                check_h2o_fb_stability_check = .false.
            end if
        end do

    end function check_h2o_fb_stability_check

    logical(c_bool) function test_h2o_fb_solver_default() bind(C)
        !
        ! this function tests the solver for the Foster-Boys localization of the
        ! occupied orbitals of water with the default settings
        !
        test_h2o_fb_solver_default = logical( &
            check_h2o_fb_solver("test_h2o_fb_solver_default", "default"), kind=c_bool)

    end function test_h2o_fb_solver_default

    logical(c_bool) function test_h2o_fb_solver_jacobi_davidson() bind(C)
        !
        ! this function tests the solver for the Foster-Boys localization of the
        ! occupied orbitals of water with the Jacobi-Davidson subsystem solver
        !
        test_h2o_fb_solver_jacobi_davidson = logical(check_h2o_fb_solver( &
            "test_h2o_fb_solver_jacobi_davidson", "jacobi-davidson"), kind=c_bool)

    end function test_h2o_fb_solver_jacobi_davidson

    logical(c_bool) function test_h2o_fb_solver_tcg() bind(C)
        !
        ! this function tests the solver for the Foster-Boys localization of the
        ! occupied orbitals of water with the truncated conjugate gradient subsystem
        ! solver
        !
        test_h2o_fb_solver_tcg = &
            logical(check_h2o_fb_solver("test_h2o_fb_solver_tcg", "tcg"), kind=c_bool)

    end function test_h2o_fb_solver_tcg

    logical(c_bool) function test_h2o_fb_solver_line_search() bind(C)
        !
        ! this function tests the solver for the Foster-Boys localization of the
        ! occupied orbitals of water with line search
        !
        test_h2o_fb_solver_line_search = logical(check_h2o_fb_solver( &
            "test_h2o_fb_solver_line_search", "line search"), kind=c_bool)

    end function test_h2o_fb_solver_line_search

    logical(c_bool) function test_h2o_fb_solver_stability() bind(C)
        !
        ! this function tests the solver for the Foster-Boys localization of the
        ! occupied orbitals of water with the stability check at convergence
        !
        test_h2o_fb_solver_stability = logical(check_h2o_fb_solver( &
            "test_h2o_fb_solver_stability", "stability"), kind=c_bool)

    end function test_h2o_fb_solver_stability

    logical(c_bool) function test_h2o_fb_stability_check_default() bind(C)
        !
        ! this function tests the stability check for the Foster-Boys localization of
        ! the occupied orbitals of water with the default settings
        !
        test_h2o_fb_stability_check_default = logical(check_h2o_fb_stability_check( &
            "test_h2o_fb_stability_check_default", "davidson"), kind=c_bool)

    end function test_h2o_fb_stability_check_default

    logical(c_bool) function test_h2o_fb_stability_check_jacobi_davidson() bind(C)
        !
        ! this function tests the stability check for the Foster-Boys localization of
        ! the occupied orbitals of water with the Jacobi-Davidson diagonalization
        ! solver
        !
        test_h2o_fb_stability_check_jacobi_davidson = logical( &
            check_h2o_fb_stability_check( &
                "test_h2o_fb_stability_check_jacobi_davidson", "jacobi-davidson"), &
            kind=c_bool)

    end function test_h2o_fb_stability_check_jacobi_davidson

end module opentrustregion_system_tests
