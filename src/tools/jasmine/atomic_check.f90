!>>> check whether Hamiltonian is real
subroutine atomic_check_hmat_real(ndim, hmat, lreal)
    use constants

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
