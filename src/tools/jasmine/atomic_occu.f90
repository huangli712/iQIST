!-------------------------------------------------------------------------
! project : jasmine
! program : atomic_make_occumat_fullspace
! source  : atomic_occu.f90
! type    : subroutines
! author  : yilin wang (email: qhwyl2006@126.com)
! history : 07/09/2014 by yilin wang
! purpose : make occupancy of eigensates of atomic Hamiltonian 
! input   :
! output  :
! status  : unstable
! comment :
!-------------------------------------------------------------------------

!>>> make occupancy for atomic eigenstates, fullspace case
subroutine atomic_make_occumat_fullspace()
    use constants,         only: zero, one
    use control,           only: ncfgs, norbs
    use m_basis_fullspace, only: bin_basis
    use m_glob_fullspace,  only: occu_mat

    implicit none

    ! local variables
    ! loop index over orbits
    integer :: iorb

    ! loop index over configurations
    integer :: ibas

    occu_mat = zero
    do ibas=1,ncfgs
        do iorb=1,norbs
            if (bin_basis(iorb, ibas) .eq. 1) then
                occu_mat(ibas, ibas) = occu_mat(ibas, ibas) + one
            endif
        enddo 
    enddo 

    return
end subroutine atomic_make_occumat_fullspace
