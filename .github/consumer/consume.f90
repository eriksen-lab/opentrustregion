! Copyright (C) 2025- Jonas Greiner
!
! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.

module consume_callbacks

    use opentrustregion, only: rp, ip

    implicit none

contains

    subroutine hess_x(x, hess_x_x, error, context)
        !
        ! this subroutine describes the Hessian linear transformation of the diagonal
        ! Hessian diag(1, 2, 3)
        !
        real(rp), intent(in), target :: x(:)
        real(rp), intent(out), target :: hess_x_x(:)
        integer(ip), intent(out) :: error
        class(*), intent(in), pointer :: context

        hess_x_x = [1.0_rp, 2.0_rp, 3.0_rp] * x
        error = 0

    end subroutine hess_x

end module consume_callbacks

program consume
    !
    ! this program is built against an installed OpenTrustRegion, runs a stability
    ! check on a diagonal Hessian with positive eigenvalues and fails unless it finds
    ! the point stable without an error
    !
    use opentrustregion, only: rp, ip, hess_x_type, stability_check
    use consume_callbacks, only: hess_x

    implicit none

    procedure(hess_x_type), pointer :: hess_x_funptr
    logical :: stable
    integer(ip) :: error

    ! run stability check
    hess_x_funptr => hess_x
    call stability_check([1.0_rp, 2.0_rp, 3.0_rp], hess_x_funptr, stable, error)
    print *, "Fortran program: error", error, "stable", stable
    if (error /= 0 .or. .not. stable) error stop 1

end program consume
