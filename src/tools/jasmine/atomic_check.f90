!-------------------------------------------------------------------------
! project : jasmine
! program : atomic_check_hmat_real
! source  : atomic_check.f90
! type    : subroutines
! author  : yilin wang (email: qhwyl2006@126.com)
! history : 07/09/2014 by yilin wang
! purpose : do some check
! input   :
! output  :
! status  : unstable
! comment :
!-------------------------------------------------------------------------

!>>> check whether Hamiltonian is real
subroutine atomic_check_hmat_real(ndim, hmat, lreal)
    use constants, only: dp, eps6

    implicit none

    ! external variables
    ! dimension of the Hamiltonian
    integer, intent(in) :: ndim
    ! the Hamiltonian matrix
    complex(dp), intent(in) :: hmat(ndim, ndim)
    ! whether Hamiltonian is real
    logical, intent(out) :: lreal

    ! local variables
    integer :: i, j

    do i=1, ndim
        do j=1, ndim
            if ( aimag(hmat(j,i)) > eps6 ) then
                lreal = .false. 
                return
            endif
        enddo
    enddo

    lreal = .true.

    return
end subroutine atomic_check_hmat_real
