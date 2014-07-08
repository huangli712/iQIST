!>>> make occupancy for atomic eigenstates, fullspace case
subroutine atomic_make_occumat_fullspace()
    use constants
    use control
    use m_basis_fullspace
    use m_glob_fullspace

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
                occu_mat(ibas, ibas) = occu_nmat(ibas, ibas) + one
            endif
        enddo 
    enddo 

    return
end subroutine atomic_make_occumat_fullspace
