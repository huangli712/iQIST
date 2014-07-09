!-------------------------------------------------------------------------
! project : jasmine
! program : atomic_driver_fullspace
!         : atomic_driver_n
!         : atomic_driver_nszps
!         : atomic_driver_njz
!         : atomic_solve_sectors
! source  : atomic_driver.f90
! type    : subroutines
! author  : yilin wang (email: qhwyl2006@126.com)
! history : 07/09/2014 by yilin wang
! purpose : this driver solve atomic problem for different CTQMC trace algorithm
! input   :
! output  :
! status  : unstable
! comment :
!-------------------------------------------------------------------------

!>>> for CTQMC trace algorithm, use full local Hilbert space 
subroutine atomic_driver_fullspace()
    use constants, only: dp
    use control, only: ncfgs
    use m_glob_fullspace, only: hmat, hmat_eigval, hmat_eigvec, &
                                alloc_m_glob_fullspace, dealloc_m_glob_fullspace
                       

    implicit none

    ! local variables
    real(dp) :: tmp_mat(ncfgs, ncfgs)
    logical :: lreal 

    ! allocate memory 
    call alloc_m_glob_fullspace()

    ! build atomic many particle Hamiltonian matrix
    call atomic_mkhmat_fullspace()

    ! check whether the many particle Hamiltonian is real 
    call atomic_check_hmat_real(hmat, lreal)
    if (lreal .eqv. .false.) then
        call atomic_print_error('atomic_driver_fullspace', 'hmat is not real !')
    endif

    ! diagonalize mp_hmat
    tmp_mat = real(hmat)
    call dmat_dsyev(ncfgs, ncfgs, tmp_mat, hmat_eigval, hmat_eigvec)

    ! build fmat
    ! first, build fmat of annihilation operators in Fock basis
    ! then, transform them to the eigen basis
    call atomic_make_annifmat_fullspace()

    ! build occupancy number
    call atomic_make_occumat_fullspace()    

    ! write eigenvalue of hmat to file 'atom.eigval.dat'
    call atomic_write_eigval_fullspace()

    ! write eigenvector of hmat to file 'atom.eigvec.dat'
    call atomic_write_eigvec_fullspace()
 
    ! write eigenvalue of hmat, occupany number of eigenstates and 
    ! fmat of annihilation fermion operators to file "atom.cix"
    ! for begonia, lavender package
    call atomic_write_atomcix_fullspace()

    ! deallocate memory
    call dealloc_m_glob_fullspace()

    return
end subroutine atomic_driver_fullspace
  
!>>> for CTQMC trace algorithm: use good quantum number, total electrons N
subroutine atomic_driver_n()
    use m_glob_sectors, only: dealloc_m_glob_sectors

    implicit none

    ! make all the sectors, allocate m_glob_sectors memory inside
    call atomic_mksectors_n()

    ! solve the atomic problem for good quantum numbers algorithm
    call atomic_solve_sectors()

    ! free memory
    call dealloc_m_glob_sectors()

    return
end subroutine atomic_driver_n

!>>> for CTQMC trace algorithm: use good quantum number, total electrons N, Sz, PS
subroutine atomic_driver_nszps()
    use m_glob_sectors, only: dealloc_m_glob_sectors

    implicit none

    ! make all the sectors, allocate m_glob_sectors memory inside
    call atomic_mksectors_nszps()

    ! solve the atomic problem for good quantum numbers algorithm
    call atomic_solve_sectors()

    ! free memory
    call dealloc_m_glob_sectors()

    return
end subroutine atomic_driver_nszps

!>>> for CTQMC trace algorithm: use good quantum number, total electrons N, Jz
subroutine atomic_driver_njz()
    use m_glob_sectors, only: dealloc_m_glob_sectors

    implicit none

    ! make all the sectors, allocate m_glob_sectors memory inside
    call atomic_mksectors_njz()

    ! solve the atomic problem for good quantum numbers algorithm
    call atomic_solve_sectors()

    ! free memory
    call dealloc_m_glob_sectors()

    return
end subroutine atomic_driver_njz

subroutine atomic_solve_sectors()
    implicit none

    ! make atomic Hamiltonian
    call atomic_mkhmat_sectors()

    ! diagonalize sector one by one
    call atomic_diag_hmat_sectors()

    ! make fmat of both creation and annihilation operators
    call atomic_make_fmat_sectors()

    ! write eigenvalue to file 'atom.eigval.dat'
    call atomic_write_eigval_sectors()

    ! write eigenvector to file 'atom.eigvec.dat'
    call atomic_write_eigvec_sectors()

    ! write information of sectors to file 'atom.cix'
    call atomic_write_atomcix_sectors()

    return
end subroutine atomic_solve_sectors
