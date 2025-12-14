! Copyright (C) 2025- Jonas Greiner
!
! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.

module opentrustregion_system_tests

    use opentrustregion, only: rp, ip, stderr
    use, intrinsic :: iso_c_binding, only: c_bool

    implicit none

    integer(ip) :: n_ao, n_mo, n_param
    real(rp), allocatable :: r_ao_ints(:, :, :), r2_ao_ints(:, :), mo_coeff(:, :), &
                             r_mo_ints(:, :, :), rii_rij_rjj_rji(:, :)
    character(:), allocatable :: data_dir

contains

    subroutine set_test_data_path(path) bind(C)
        !
        ! this subroutine sets the path to test data
        !
        use, intrinsic :: iso_c_binding, only: c_ptr, c_char, c_f_pointer, c_null_char

        type(c_ptr), intent(in), value :: path 
        character(c_char), pointer :: c_path(:)
        integer(ip) :: len, i

        ! conda paths need ample space
        call c_f_pointer(path, c_path, [1024])
        len = 0
        do i = 1, size(c_path)
            if (c_path(i) == c_null_char) exit
            len = len + 1
        end do
        allocate(character(len) :: data_dir)
        data_dir = transfer(c_path(1:len), data_dir)

    end subroutine set_test_data_path

    logical(c_bool) function test_h2o_atomic_fb() bind(C)
        !
        ! this function tests the Foster-Boys localization on water
        !
        use opentrustregion, only: update_orbs_type, obj_func_type, &
                                   solver_settings_type, solver,  hess_x_type, &
                                   stability_settings_type, stability_check

        procedure(update_orbs_type), pointer :: update_orbs_funptr
        procedure(obj_func_type), pointer :: obj_func_funptr
        procedure(hess_x_type), pointer :: hess_x_funptr
        type(solver_settings_type) :: solver_settings
        type(stability_settings_type) :: stability_settings
        integer(ip) :: ios, error
        real(rp), allocatable :: kappa(:), grad(:), h_diag(:)
        real(rp) :: func
        logical :: stable

        ! assume tests pass
        test_h2o_atomic_fb = .true.

        ! set size parameters for system
        n_ao = 13
        n_mo = 5
        n_param = n_mo * (n_mo - 1) / 2

        ! allocate arrays
        allocate(mo_coeff(n_ao, n_mo), r_ao_ints(3, n_ao, n_ao), &
                 r2_ao_ints(n_ao, n_ao), r_mo_ints(3, n_mo, n_mo), &
                 rii_rij_rjj_rji(n_mo, n_mo), kappa(n_param), grad(n_param), &
                 h_diag(n_param))

        ! check if test data variable is set
        if (.not. allocated(data_dir)) error stop "Test data directory not set "// &
            "through set_test_data_path subroutine before calling system test."

        ! read raw binary data
        open (unit=10, file=data_dir//"/h2o_atomic_mo_coeff.bin", form="unformatted", &
              access="stream", status="old", action="read", iostat=ios)
        if (ios /= 0) error stop "Error opening file"
        read (10) mo_coeff
        close (10)
        open (unit=10, file=data_dir//"/h2o_r_ints.bin", form="unformatted", &
              access="stream", status="old", action="read", iostat=ios)
        if (ios /= 0) error stop "Error opening file"
        read (10) r_ao_ints
        close (10)
        open (unit=10, file=data_dir//"/h2o_r2_ints.bin", form="unformatted", &
              access="stream", status="old", action="read", iostat=ios)
        if (ios /= 0) error stop "Error opening file"
        read (10) r2_ao_ints
        close (10)

        ! set function pointers
        update_orbs_funptr => update_orbs
        obj_func_funptr => obj_func

        ! initialize settings
        call solver_settings%init(error)

        ! call solver
        call solver(update_orbs_funptr, obj_func_funptr, n_param, error, &
                    solver_settings)

        ! check if error has occured
        if (error /= 0) then
            write (stderr, *) "test_h2o_atomic_fb failed: Solver subroutine "// &
                "produced error."
            test_h2o_atomic_fb = .false.
        end if

        ! get gradient, Hessian diagonal and Hessian linear transformation function
        ! pointer
        kappa = 0.0_rp
        call update_orbs(kappa, func, grad, h_diag, hess_x_funptr, error)

        ! check if error has occured
        if (error /= 0) then
            write (stderr, *) "test_h2o_atomic_fb failed: Orbital updating "// &
            "subroutine produced error."
            test_h2o_atomic_fb = .false.
        end if

        ! initialize settings
        call stability_settings%init(error)

        ! perform stability check
        call stability_check(h_diag, hess_x_funptr, stable, error, stability_settings)

        ! check if error has occured
        if (error /= 0) then
            write (stderr, *) "test_h2o_atomic_fb failed: Stability check "// &
            "subroutine produced error."
            test_h2o_atomic_fb = .false.
        end if

        ! test if solution is stable
        if (.not. stable) then
            write (stderr, *) "test_h2o_atomic_fb failed: Solver did not converge "// &
                "to minimum."
            test_h2o_atomic_fb = .false.
        end if

        ! deallocate arrays
        deallocate(mo_coeff, r_ao_ints, r2_ao_ints, r_mo_ints, rii_rij_rjj_rji, kappa, &
                   grad, h_diag)

    end function test_h2o_atomic_fb

    logical(c_bool) function test_h2o_saddle_fb() bind(C)
        !
        ! this function tests the Foster-Boys localization on water starting from a
        ! saddle point
        !
        use opentrustregion, only: update_orbs_type, obj_func_type, &
                                   solver_settings_type, solver, hess_x_type, &
                                   stability_settings_type, stability_check

        procedure(update_orbs_type), pointer :: update_orbs_funptr
        procedure(obj_func_type), pointer :: obj_func_funptr
        procedure(hess_x_type), pointer :: hess_x_funptr
        type(solver_settings_type) :: solver_settings
        type(stability_settings_type) :: stability_settings
        integer(ip) :: ios, error
        real(rp), allocatable :: kappa(:), grad(:), h_diag(:)
        real(rp) :: func
        logical :: stable

        ! assume tests pass
        test_h2o_saddle_fb = .true.

        ! set size parameters for system
        n_ao = 13
        n_mo = 5
        n_param = n_mo*(n_mo - 1)/2

        ! allocate arrays
        allocate(mo_coeff(n_ao, n_mo), r_ao_ints(3, n_ao, n_ao), &
                 r2_ao_ints(n_ao, n_ao), r_mo_ints(3, n_mo, n_mo), &
                 rii_rij_rjj_rji(n_mo, n_mo), kappa(n_param), grad(n_param), &
                 h_diag(n_param))

        ! check if test data variable is set
        if (.not. allocated(data_dir)) error stop "Test data directory not set "// &
            "through set_test_data_path subroutine before calling system test."

        ! read raw binary data
        open (unit=10, file=data_dir//"/h2o_saddle_mo_coeff.bin", form="unformatted", &
              access="stream", status="old", action="read", iostat=ios)
        if (ios /= 0) error stop "Error opening file"
        read (10) mo_coeff
        close (10)
        open (unit=10, file=data_dir//"/h2o_r_ints.bin", form="unformatted", &
              access="stream", status="old", action="read", iostat=ios)
        if (ios /= 0) error stop "Error opening file"
        read (10) r_ao_ints
        close (10)
        open (unit=10, file=data_dir//"/h2o_r2_ints.bin", form="unformatted", &
              access="stream", status="old", action="read", iostat=ios)
        if (ios /= 0) error stop "Error opening file"
        read (10) r2_ao_ints
        close (10)

        ! set function pointers
        update_orbs_funptr => update_orbs
        obj_func_funptr => obj_func

        ! initialize settings
        call solver_settings%init(error)

        ! call solver
        call solver(update_orbs_funptr, obj_func_funptr, n_param, error, &
                    solver_settings)

        ! check if error has occured
        if (error /= 0) then
            write (stderr, *) "test_h2o_saddle_fb failed: Solver subroutine "// &
                "produced error."
            test_h2o_saddle_fb = .false.
        end if

        ! get gradient, Hessian diagonal and Hessian linear transformation function
        ! pointer
        kappa = 0.0_rp
        call update_orbs(kappa, func, grad, h_diag, hess_x_funptr, error)

        ! check if error has occured
        if (error /= 0) then
            write (stderr, *) "test_h2o_saddle_fb failed: Orbital update "// &
            "subroutine produced error."
            test_h2o_saddle_fb = .false.
        end if

        ! initialize settings
        call stability_settings%init(error)

        ! perform stability check
        call stability_check(h_diag, hess_x_funptr, stable, error, stability_settings)

        ! check if error has occured
        if (error /= 0) then
            write (stderr, *) "test_h2o_saddle_fb failed: Stability check "// &
            "subroutine produced error."
            test_h2o_saddle_fb = .false.
        end if

        ! test if solution is stable
        if (.not. stable) then
            write (stderr, *) "test_h2o_saddle_fb failed: Solver did not converge "// &
                "to minimum."
            test_h2o_saddle_fb = .false.
        end if

        ! deallocate arrays
        deallocate(mo_coeff, r_ao_ints, r2_ao_ints, r_mo_ints, rii_rij_rjj_rji, kappa, &
                   grad, h_diag)

    end function test_h2o_saddle_fb

    real(rp) function obj_func(kappa, error)
        !
        ! this function calculates the Foster-Boys orbital localization objective 
        ! function
        !
        real(rp), intent(in), target :: kappa(:)
        integer(ip), intent(out) :: error

        real(rp), allocatable :: kappa_full(:, :), mo_coeff_tmp(:, :)
        integer(ip) :: xyz, i, j, idx

        ! initialize error flag
        error = 0

        ! unpack orbital rotation
        allocate(kappa_full(n_mo, n_mo))
        kappa_full = 0.0
        idx = 1
        do i = 2, n_mo
            do j = 1, i - 1
                kappa_full(i, j) = kappa(idx)
                kappa_full(j, i) = -kappa(idx)
                idx = idx + 1
            end do
        end do

        ! rotate orbitals
        mo_coeff_tmp = matmul(mo_coeff, exp_asymm_mat(kappa_full))
        deallocate(kappa_full)

        ! compute cost function
        obj_func = 0.0_rp
        do i = 1, n_mo
            obj_func = obj_func &
                       + dot_product(mo_coeff_tmp(:, i), &
                                     matmul(r2_ao_ints, mo_coeff_tmp(:, i)))
            do xyz = 1, 3
                obj_func = obj_func &
                           - dot_product(mo_coeff_tmp(:, i), &
                                         matmul(r_ao_ints(xyz, :, :), &
                                                mo_coeff_tmp(:, i)))**2
            end do
        end do
        deallocate(mo_coeff_tmp)

    end function obj_func

    subroutine update_orbs(kappa, func, grad, h_diag, hess_x_funptr, error)
        !
        ! this function updates the orbitals for Foster-Boys orbital localization
        !
        use opentrustregion, only: hess_x_type

        real(rp), intent(in), target :: kappa(:)
        real(rp), intent(out) :: func
        real(rp), intent(out), target :: grad(:), h_diag(:)
        procedure(hess_x_type), intent(out), pointer :: hess_x_funptr
        integer(ip), intent(out) :: error

        integer(ip) :: xyz, i, j, idx
        real(rp), allocatable :: kappa_full(:, :), h_diag_tmp(:, :), tmp1(:, :)

        ! initialize error flag
        error = 0

        ! unpack orbital rotation
        allocate(kappa_full(n_mo, n_mo))
        kappa_full = 0.0
        idx = 1
        do i = 2, n_mo
            do j = 1, i - 1
                kappa_full(i, j) = kappa(idx)
                kappa_full(j, i) = -kappa(idx)
                idx = idx + 1
            end do
        end do

        ! rotate orbitals
        mo_coeff = matmul(mo_coeff, exp_asymm_mat(kappa_full))
        deallocate(kappa_full)

        ! transform integrals to MO basis
        do xyz = 1, 3
            r_mo_ints(xyz, :, :) = matmul(matmul(transpose(mo_coeff), &
                                                 r_ao_ints(xyz, :, :)), mo_coeff)
        end do

        ! compute cost function
        func = 0.0_rp
        do i = 1, n_mo
            func = func &
                   + dot_product(mo_coeff(:, i), matmul(r2_ao_ints, mo_coeff(:, i))) &
                   - sum(r_mo_ints(:, i, i)**2)
        end do

        ! construct temporary intermediate
        allocate(tmp1(n_mo, n_mo))
        do i = 1, n_mo
            do j = 1, n_mo
                tmp1(i, j) = sum(r_mo_ints(:, j, j)*r_mo_ints(:, j, i))
            end do
        end do

        ! construct gradient and extract lower triagonal
        idx = 1
        do i = 2, n_mo
            do j = 1, i - 1
                grad(idx) = -2*(tmp1(i, j) - tmp1(j, i))
                idx = idx + 1
            end do
        end do

        ! construct Hessian diagonal
        allocate(h_diag_tmp(n_mo, n_mo))
        do i = 1, n_mo
            do j = 1, n_mo
                h_diag_tmp(i, j) = 2*sum(r_mo_ints(:, j, j)*r_mo_ints(:, i, i) + &
                                         r_mo_ints(:, j, i)*r_mo_ints(:, j, i) + &
                                         r_mo_ints(:, j, i)*r_mo_ints(:, i, j)) - &
                                   tmp1(i, i) - tmp1(j, j)
            end do
        end do

        ! extract lower triagonal
        idx = 1
        do i = 2, n_mo
            do j = 1, i - 1
                h_diag(idx) = -2*(h_diag_tmp(i, j))
                idx = idx + 1
            end do
        end do
        deallocate(h_diag_tmp)

        ! save for Hessian linear transformation
        rii_rij_rjj_rji = tmp1 + transpose(tmp1)
        deallocate(tmp1)

        ! get function pointer to Hessian linear transformation
        hess_x_funptr => hess_x_fun

    end subroutine update_orbs

    subroutine hess_x_fun(x, hess_x, error)
        !
        ! this function performs the Hessian linear transformation for Foster-Boys 
        ! orbital localization, it cannot be defined within update_orbs as it would 
        ! otherwise go out of scope when that subroutine returns
        !
        real(rp), intent(in), target :: x(:)
        real(rp), intent(out), target :: hess_x(:)
        integer(ip), intent(out) :: error

        real(rp), allocatable :: x_full(:, :), hess_x_full(:, :), tmp2(:, :), &
                                 tmp3(:, :, :)
        integer(ip) :: xyz1, i1, j1, idx1

        ! initialize error flag
        error = 0

        ! unpack trial vector
        allocate(x_full(n_mo, n_mo))
        x_full = 0.0
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
                tmp2(xyz1, i1) = sum(x_full(:, i1)*r_mo_ints(xyz1, i1, :))
                do j1 = 1, n_mo
                    tmp3(xyz1, i1, j1) = sum(x_full(:, j1)*r_mo_ints(xyz1, i1, :))
                end do
            end do
        end do

        ! construct Hessian linear transformation
        allocate(hess_x_full(n_mo, n_mo))
        hess_x_full = matmul(transpose(x_full), transpose(rii_rij_rjj_rji))
        do i1 = 1, n_mo
            do j1 = 1, n_mo
                hess_x_full(i1, j1) = hess_x_full(i1, j1) + &
                                      2*sum(r_mo_ints(:, j1, i1)*tmp2(:, j1) - &
                                            r_mo_ints(:, i1, i1)*tmp3(:, j1, i1) - &
                                            r_mo_ints(:, j1, i1)*tmp2(:, i1))
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
            tmp(:, i) = eigvecs(:, i)*cmplx(cos(eigvals(i)), sin(eigvals(i)), kind=rp)
        end do
        exp_asymm_mat = real(transpose(matmul(tmp, conjg(transpose(eigvecs)))))
        deallocate(eigvecs, eigvals, tmp)

    end function exp_asymm_mat

end module opentrustregion_system_tests
