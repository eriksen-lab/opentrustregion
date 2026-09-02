! Copyright (C) 2025- Susi Lehtola
!
! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.
!
! Internal BLAS/LAPACK shim. Provides the six routines the OpenTrustRegion
! core actually calls -- DDOT, DNRM2, DGEMV, DGEMM, DSYEV, DSYSV -- in the
! library's working precision. Selected at build time via
! OPENTRUSTREGION_USE_INTERNAL_BLAS=ON; CMakeLists.txt arranges the textual
! remap `ddot=otr_ddot`, etc., through Fortran_PREPROCESS, so the core source
! is identical between builds that use external vendor BLAS and builds that
! use this shim.
!
! Two reasons to enable this:
! 1. quad-precision builds (no external binary128 BLAS exists)
! 2. portable / minimal builds where pulling in an external BLAS dependency
!    is undesirable
!
! The actual BLAS/LAPACK entry points below are deliberately free-standing
! subprograms (not module procedures) so they take the same linker name shape
! as the BLAS/LAPACK routines they replace. They share kind parameters via
! the small `otr_blas_shim_kinds` module to keep `rp`/`ip` definitions in one
! place.
!
! Two distinct workload regimes:
!   - DSYEV and DSYSV operate on the *reduced-space* matrices, which are
!     tiny (n typically <= 30, always <= ~150). Textbook O(n^3) cyclic
!     Jacobi and partial-pivoting LU are adequate here; vendor LAPACK's
!     blocked algorithms would buy nothing at this size.
!   - DDOT, DNRM2, DGEMV, DGEMM operate on quantities of length / leading
!     dimension n_param, which is the host's orbital-rotation parameter
!     count. n_param can be very large (10^6 - 10^9 in real chemistry use)
!     and these calls do real work even though the calls are tall-skinny
!     (one large dim, one small).
!
! For the BLAS-2/3 routines we rely on gfortran/ifort autovectorization of
! the inner kernels: loops are written column-major over contiguous strides,
! with unit-stride fast paths in the BLAS-1 routines, kept simple enough
! that the compiler reliably issues SIMD code. We do not implement
! cache-blocked GEMM, hand-vectorized intrinsics, or threaded loops -- the
! library does not use OpenMP (the host program owns the threading model).
! Users who want top-tier double-precision BLAS-2/3 performance should
! build with OPENTRUSTREGION_USE_INTERNAL_BLAS=OFF and link a vendor BLAS.
! Quad has no vendor option, so its BLAS-2/3 performance is whatever the
! autovectorized loop achieves -- adequate but not competitive with what
! a hand-blocked binary128 GEMM could do.

module otr_blas_shim_kinds
    use, intrinsic :: iso_fortran_env, only: real64, real128, int32, int64
    implicit none
    private

#ifdef USE_QUAD
    integer, parameter, public :: rp = real128
#else
    integer, parameter, public :: rp = real64
#endif
#ifdef USE_ILP64
    integer, parameter, public :: ip = int64
#else
    integer, parameter, public :: ip = int32
#endif

    real(rp), parameter, public :: zero = 0._rp, one = 1._rp, two = 2._rp
end module otr_blas_shim_kinds

! ----  BLAS-1: DDOT  ----------------------------------------------------------
real(rp) function otr_ddot(n, x, incx, y, incy)
    use otr_blas_shim_kinds, only: rp, ip, zero
    implicit none
    integer(ip), intent(in) :: n, incx, incy
    real(rp), intent(in) :: x(*), y(*)
    integer(ip) :: i, ix, iy
    real(rp) :: acc

    acc = zero
    if (n <= 0) then
        otr_ddot = acc
        return
    end if
    if (incx == 1 .and. incy == 1) then
        do i = 1, n
            acc = acc + x(i) * y(i)
        end do
    else
        ix = 1; if (incx < 0) ix = 1 - (n - 1) * incx
        iy = 1; if (incy < 0) iy = 1 - (n - 1) * incy
        do i = 1, n
            acc = acc + x(ix) * y(iy)
            ix = ix + incx
            iy = iy + incy
        end do
    end if
    otr_ddot = acc
end function otr_ddot

! ----  BLAS-1: DNRM2  ---------------------------------------------------------
! Scaled accumulation for over/underflow safety: nrm2 = scale * sqrt(sumsq).
real(rp) function otr_dnrm2(n, x, incx)
    use otr_blas_shim_kinds, only: rp, ip, zero, one
    implicit none
    integer(ip), intent(in) :: n, incx
    real(rp), intent(in) :: x(*)
    integer(ip) :: i, ix
    real(rp) :: scale, sumsq, ax

    if (n <= 0) then
        otr_dnrm2 = zero
        return
    end if

    scale = zero
    sumsq = one
    ix = 1; if (incx < 0) ix = 1 - (n - 1) * incx
    do i = 1, n
        if (x(ix) /= zero) then
            ax = abs(x(ix))
            if (scale < ax) then
                sumsq = one + sumsq * (scale / ax) ** 2
                scale = ax
            else
                sumsq = sumsq + (ax / scale) ** 2
            end if
        end if
        ix = ix + incx
    end do
    otr_dnrm2 = scale * sqrt(sumsq)
end function otr_dnrm2

! ----  BLAS-2: DGEMV  ---------------------------------------------------------
subroutine otr_dgemv(trans, m, n, alpha, a, lda, x, incx, beta, y, incy)
    use otr_blas_shim_kinds, only: rp, ip, zero, one
    implicit none
    character, intent(in) :: trans
    integer(ip), intent(in) :: m, n, lda, incx, incy
    real(rp), intent(in) :: alpha, beta, a(lda, *), x(*)
    real(rp), intent(inout) :: y(*)
    integer(ip) :: i, j, ix, iy, leny
    real(rp) :: acc, t

    if (trans == 'N' .or. trans == 'n') then
        leny = m
    else
        leny = n
    end if
    if (m == 0 .or. n == 0 .or. (alpha == zero .and. beta == one)) return

    ! y := beta * y
    iy = 1; if (incy < 0) iy = 1 - (leny - 1) * incy
    if (beta == zero) then
        do i = 1, leny
            y(iy) = zero
            iy = iy + incy
        end do
    else if (beta /= one) then
        do i = 1, leny
            y(iy) = beta * y(iy)
            iy = iy + incy
        end do
    end if
    if (alpha == zero) return

    if (trans == 'N' .or. trans == 'n') then
        ! y := alpha * A * x + y
        if (incx == 1 .and. incy == 1) then
            do j = 1, n
                t = alpha * x(j)
                do i = 1, m
                    y(i) = y(i) + t * a(i, j)
                end do
            end do
        else
            ix = 1; if (incx < 0) ix = 1 - (n - 1) * incx
            do j = 1, n
                t = alpha * x(ix)
                iy = 1; if (incy < 0) iy = 1 - (m - 1) * incy
                do i = 1, m
                    y(iy) = y(iy) + t * a(i, j)
                    iy = iy + incy
                end do
                ix = ix + incx
            end do
        end if
    else
        ! y := alpha * A^T * x + y
        if (incx == 1 .and. incy == 1) then
            do j = 1, n
                acc = zero
                do i = 1, m
                    acc = acc + a(i, j) * x(i)
                end do
                y(j) = y(j) + alpha * acc
            end do
        else
            iy = 1; if (incy < 0) iy = 1 - (n - 1) * incy
            do j = 1, n
                acc = zero
                ix = 1; if (incx < 0) ix = 1 - (m - 1) * incx
                do i = 1, m
                    acc = acc + a(i, j) * x(ix)
                    ix = ix + incx
                end do
                y(iy) = y(iy) + alpha * acc
                iy = iy + incy
            end do
        end if
    end if
end subroutine otr_dgemv

! ----  BLAS-3: DGEMM  ---------------------------------------------------------
subroutine otr_dgemm(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc)
    use otr_blas_shim_kinds, only: rp, ip, zero, one
    implicit none
    character, intent(in) :: transa, transb
    integer(ip), intent(in) :: m, n, k, lda, ldb, ldc
    real(rp), intent(in) :: alpha, beta, a(lda, *), b(ldb, *)
    real(rp), intent(inout) :: c(ldc, *)
    integer(ip) :: i, j, l
    real(rp) :: acc
    logical :: nota, notb

    nota = (transa == 'N' .or. transa == 'n')
    notb = (transb == 'N' .or. transb == 'n')

    if (m == 0 .or. n == 0) return

    ! C := beta * C
    if (beta == zero) then
        do j = 1, n
            do i = 1, m
                c(i, j) = zero
            end do
        end do
    else if (beta /= one) then
        do j = 1, n
            do i = 1, m
                c(i, j) = beta * c(i, j)
            end do
        end do
    end if
    if (alpha == zero .or. k == 0) return

    if (nota .and. notb) then
        ! C := alpha * A * B + C
        do j = 1, n
            do l = 1, k
                acc = alpha * b(l, j)
                do i = 1, m
                    c(i, j) = c(i, j) + acc * a(i, l)
                end do
            end do
        end do
    else if ((.not. nota) .and. notb) then
        ! C := alpha * A^T * B + C
        do j = 1, n
            do i = 1, m
                acc = zero
                do l = 1, k
                    acc = acc + a(l, i) * b(l, j)
                end do
                c(i, j) = c(i, j) + alpha * acc
            end do
        end do
    else if (nota .and. (.not. notb)) then
        ! C := alpha * A * B^T + C
        do j = 1, n
            do l = 1, k
                acc = alpha * b(j, l)
                do i = 1, m
                    c(i, j) = c(i, j) + acc * a(i, l)
                end do
            end do
        end do
    else
        ! C := alpha * A^T * B^T + C
        do j = 1, n
            do i = 1, m
                acc = zero
                do l = 1, k
                    acc = acc + a(l, i) * b(j, l)
                end do
                c(i, j) = c(i, j) + alpha * acc
            end do
        end do
    end if
end subroutine otr_dgemm

! ----  LAPACK: DSYEV via cyclic Jacobi  ---------------------------------------
! Computes all eigenvalues and (optionally) eigenvectors of a real symmetric
! matrix using cyclic Jacobi rotations. For the small matrices this routine
! sees (n typically <= 30, always <= ~150), Jacobi converges quadratically and
! is overwhelmingly the simplest correct algorithm.
subroutine otr_dsyev(jobz, uplo, n, a, lda, w, work, lwork, info)
    use otr_blas_shim_kinds, only: rp, ip, zero, one, two
    implicit none
    character, intent(in) :: jobz, uplo
    integer(ip), intent(in) :: n, lda, lwork
    integer(ip), intent(out) :: info
    real(rp), intent(inout) :: a(lda, *)
    real(rp), intent(out) :: w(*), work(*)
    integer(ip) :: i, j, k, sweep, max_sweeps
    real(rp) :: off, off_init, tol, theta, t, c, s, app, aqq, apq, tmp
    logical :: want_vec

    info = 0
    if (jobz /= 'N' .and. jobz /= 'n' .and. jobz /= 'V' .and. jobz /= 'v') info = -1
    if (uplo /= 'U' .and. uplo /= 'u' .and. uplo /= 'L' .and. uplo /= 'l') info = -2
    if (n < 0) info = -3
    if (lda < max(1_ip, n)) info = -5
    if (lwork < 1 .and. lwork /= -1) info = -8
    if (info /= 0) return

    ! workspace query: Jacobi needs no scratch space
    if (lwork == -1) then
        work(1) = one
        return
    end if
    if (n == 0) return

    want_vec = (jobz == 'V' .or. jobz == 'v')

    ! mirror the stored triangle into the full matrix so the rotations below
    ! can read either (i,j) or (j,i)
    if (uplo == 'U' .or. uplo == 'u') then
        do j = 1, n
            do i = 1, j - 1
                a(j, i) = a(i, j)
            end do
        end do
    else
        do j = 1, n
            do i = 1, j - 1
                a(i, j) = a(j, i)
            end do
        end do
    end if

    block
        real(rp), allocatable :: m(:, :), v(:, :)
        allocate(m(n, n))
        do j = 1, n
            do i = 1, n
                m(i, j) = a(i, j)
            end do
        end do
        if (want_vec) then
            allocate(v(n, n))
            v = zero
            do i = 1, n
                v(i, i) = one
            end do
        end if

        max_sweeps = 100_ip
        ! tolerance scales with the Frobenius norm of the matrix
        off_init = zero
        do j = 1, n
            do i = 1, n
                off_init = off_init + m(i, j) * m(i, j)
            end do
        end do
        tol = sqrt(off_init) * epsilon(zero)

        do sweep = 1, max_sweeps
            ! off-diagonal Frobenius norm (squared)
            off = zero
            do j = 2, n
                do i = 1, j - 1
                    off = off + m(i, j) * m(i, j)
                end do
            end do
            if (sqrt(two * off) <= tol) exit

            ! cyclic sweep over upper-triangular off-diagonal elements
            do j = 2, n
                do i = 1, j - 1
                    apq = m(i, j)
                    if (abs(apq) <= tol / real(n, rp)) cycle
                    app = m(i, i)
                    aqq = m(j, j)
                    ! rotation per Rutishauser
                    theta = (aqq - app) / (two * apq)
                    if (theta >= zero) then
                        t = one / (theta + sqrt(one + theta * theta))
                    else
                        t = one / (theta - sqrt(one + theta * theta))
                    end if
                    c = one / sqrt(one + t * t)
                    s = t * c
                    m(i, i) = app - t * apq
                    m(j, j) = aqq + t * apq
                    m(i, j) = zero
                    m(j, i) = zero
                    do k = 1, n
                        if (k == i .or. k == j) cycle
                        tmp = m(k, i)
                        m(k, i) = c * tmp - s * m(k, j)
                        m(k, j) = s * tmp + c * m(k, j)
                        m(i, k) = m(k, i)
                        m(j, k) = m(k, j)
                    end do
                    if (want_vec) then
                        do k = 1, n
                            tmp = v(k, i)
                            v(k, i) = c * tmp - s * v(k, j)
                            v(k, j) = s * tmp + c * v(k, j)
                        end do
                    end if
                end do
            end do
        end do

        if (sweep > max_sweeps) then
            info = 1   ! did not converge -- mirrors LAPACK's `info > 0` convention
        end if

        ! extract eigenvalues and sort ascending; permute eigenvectors to match
        do i = 1, n
            w(i) = m(i, i)
        end do
        do i = 1, n - 1
            k = i
            do j = i + 1, n
                if (w(j) < w(k)) k = j
            end do
            if (k /= i) then
                tmp = w(i); w(i) = w(k); w(k) = tmp
                if (want_vec) then
                    do j = 1, n
                        tmp = v(j, i); v(j, i) = v(j, k); v(j, k) = tmp
                    end do
                end if
            end if
        end do

        if (want_vec) then
            do j = 1, n
                do i = 1, n
                    a(i, j) = v(i, j)
                end do
            end do
            deallocate(v)
        end if
        deallocate(m)
    end block

    work(1) = one
end subroutine otr_dsyev

! ----  LAPACK: DSYSV via plain LU with partial pivoting  ----------------------
! The core only ever feeds DSYSV very small matrices (n <= ~150, typically
! <= 30) with one right-hand side. Bunch-Kaufman is unnecessary at that size;
! plain LU-with-partial-pivoting is correct for indefinite symmetric systems
! and adds only O(n^3 / 3) extra work compared to DSYTRF.
subroutine otr_dsysv(uplo, n, nrhs, a, lda, ipiv, b, ldb, work, lwork, info)
    use otr_blas_shim_kinds, only: rp, ip, zero, one
    implicit none
    character, intent(in) :: uplo
    integer(ip), intent(in) :: n, nrhs, lda, ldb, lwork
    integer(ip), intent(out) :: info, ipiv(*)
    real(rp), intent(inout) :: a(lda, *), b(ldb, *)
    real(rp), intent(out) :: work(*)
    integer(ip) :: i, j, k, p
    real(rp) :: amax, mult, tmp

    info = 0
    if (uplo /= 'U' .and. uplo /= 'u' .and. uplo /= 'L' .and. uplo /= 'l') info = -1
    if (n < 0) info = -2
    if (nrhs < 0) info = -3
    if (lda < max(1_ip, n)) info = -5
    if (ldb < max(1_ip, n)) info = -8
    if (lwork < 1 .and. lwork /= -1) info = -10
    if (info /= 0) return

    if (lwork == -1) then
        work(1) = one
        return
    end if
    if (n == 0 .or. nrhs == 0) return

    ! mirror stored triangle to full so we can do non-symmetric LU
    if (uplo == 'U' .or. uplo == 'u') then
        do j = 1, n
            do i = 1, j - 1
                a(j, i) = a(i, j)
            end do
        end do
    else
        do j = 1, n
            do i = 1, j - 1
                a(i, j) = a(j, i)
            end do
        end do
    end if

    ! LU factorization with partial (row) pivoting
    do k = 1, n
        p = k
        amax = abs(a(k, k))
        do i = k + 1, n
            if (abs(a(i, k)) > amax) then
                amax = abs(a(i, k))
                p = i
            end if
        end do
        ipiv(k) = p
        if (amax == zero) then
            info = k
            return
        end if
        if (p /= k) then
            do j = 1, n
                tmp = a(k, j); a(k, j) = a(p, j); a(p, j) = tmp
            end do
            do j = 1, nrhs
                tmp = b(k, j); b(k, j) = b(p, j); b(p, j) = tmp
            end do
        end if
        do i = k + 1, n
            mult = a(i, k) / a(k, k)
            a(i, k) = mult
            do j = k + 1, n
                a(i, j) = a(i, j) - mult * a(k, j)
            end do
            do j = 1, nrhs
                b(i, j) = b(i, j) - mult * b(k, j)
            end do
        end do
    end do

    ! back-substitution
    do j = 1, nrhs
        do i = n, 1, -1
            tmp = b(i, j)
            do k = i + 1, n
                tmp = tmp - a(i, k) * b(k, j)
            end do
            b(i, j) = tmp / a(i, i)
        end do
    end do

    work(1) = one
end subroutine otr_dsysv
