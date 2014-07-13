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
    use constants,        only: dp, mystd
    use control,          only: ncfgs
    use m_glob_fullspace, only: hmat, hmat_eigval, hmat_eigvec, &
                                alloc_m_glob_fullspace, dealloc_m_glob_fullspace

    implicit none

    ! local variables
    real(dp) :: tmp_mat(ncfgs, ncfgs)
    logical :: lreal 

    ! allocate memory 
    write(mystd, "(2X,a)") "jasmine >>> allocate global memory ..."
    write(mystd,*)
    call alloc_m_glob_fullspace()

    ! build atomic many particle Hamiltonian matrix
    write(mystd, "(2X,a)") "jasmine >>> make atomic many particle Hamiltonian ..."
    write(mystd,*)
    call atomic_mkhmat_fullspace()

    ! check whether the many particle Hamiltonian is real 
    write(mystd, "(2X,a)") "jasmine >>> check whether Hamiltonian is real or not ..."
    write(mystd,*)
    call atomic_check_hmat_real(ncfgs, hmat, lreal)
    if (lreal .eqv. .false.) then
        call atomic_print_error('atomic_driver_fullspace', 'hmat is not real !')
    else
        write(mystd, "(2X,a)") "jasmine >>> the Hamiltonian is real"
        write(mystd,*)
    endif

    ! diagonalize hmat
    write(mystd, "(2X,a)") "jasmine >>> diagonalize atomic Hamiltonian ..."
    write(mystd,*)
    tmp_mat = real(hmat)
    call dmat_dsyev(ncfgs, ncfgs, tmp_mat, hmat_eigval, hmat_eigvec)

    ! build fmat
    ! first, build fmat of annihilation operators in Fock basis
    ! then, transform them to the eigen basis
    write(mystd, "(2X,a)") "jasmine >>> make fmat for annihilation fermion operators ... "
    write(mystd,*)
    call atomic_make_annifmat_fullspace()

    write(mystd, "(2X,a)") "jasmine >>> make occupancy number of atomic eigenstates ... "
    write(mystd,*)
    ! build occupancy number
    call atomic_make_occumat_fullspace()    

    write(mystd, "(2X,a)") "jasmine >>> write eigenvalue, eigenvector, and atom.cix to files ..."
    write(mystd,*)
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
    use constants, only: mystd
    use m_glob_sectors, only: dealloc_m_glob_sectors

    implicit none

    ! make all the sectors, allocate m_glob_sectors memory inside
    write(mystd, "(2X,a)") "jasmine >>> determine sectors by good quantum numbers N ... "
    write(mystd,*)
    call atomic_mksectors_n()

    ! solve the atomic problem for good quantum numbers algorithm
    call atomic_solve_sectors()

    ! free memory
    call dealloc_m_glob_sectors()

    return
end subroutine atomic_driver_n

!>>> for CTQMC trace algorithm: use good quantum number, total electrons N
subroutine atomic_driver_nsz()
    use constants, only: mystd
    use m_glob_sectors, only: dealloc_m_glob_sectors

    implicit none

    ! make all the sectors, allocate m_glob_sectors memory inside
    write(mystd, "(2X,a)") "jasmine >>> determine sectors by good quantum numbers N ... "
    write(mystd,*)
    call atomic_mksectors_nsz()

    ! solve the atomic problem for good quantum numbers algorithm
    call atomic_solve_sectors()

    ! free memory
    call dealloc_m_glob_sectors()

    return
end subroutine atomic_driver_nsz


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
    use constants, only: mystd
    use m_glob_sectors, only: dealloc_m_glob_sectors

    implicit none

    ! make all the sectors, allocate m_glob_sectors memory inside
    write(mystd, "(2X,a)") "jasmine >>> determine sectors by good quantum numbers N, Jz ... "
    write(mystd,*)
    call atomic_mksectors_njz()

    ! solve the atomic problem for good quantum numbers algorithm
    call atomic_solve_sectors()

    ! free memory
    call dealloc_m_glob_sectors()

    return
end subroutine atomic_driver_njz

subroutine atomic_solve_sectors()
    use constants, only: mystd
    use m_glob_sectors

    implicit none
    ! local variables
    integer :: i
    logical :: lreal

    ! make atomic Hamiltonian
    write(mystd, "(2X,a)") "jasmine >>> make atomic Hamiltonian for each sector ... "
    write(mystd,*)
    call atomic_mkhmat_sectors()

    ! check whether the many particle Hamiltonian is real 
    write(mystd, "(2X,a)") "jasmine >>> check whether Hamiltonian is real or not ..."
    write(mystd,*)
    do i=1, nsectors 
        call atomic_check_hmat_real(sectors(i)%ndim, sectors(i)%myham, lreal)
        if (lreal .eqv. .false.) then
            call atomic_print_error('atomic_driver_fullspace', 'hmat is not real !')
        endif
    enddo
    write(mystd, "(2X,a)") "jasmine >>> the Hamiltonian is real"
    write(mystd,*)

    ! diagonalize sector one by one
    write(mystd, "(2X,a)") "jasmine >>> diagonalize atomic Hamiltonian for each sector ... "
    write(mystd,*)
    call atomic_diag_hmat_sectors()

    ! make fmat of both creation and annihilation operators
    write(mystd, "(2X,a)") "jasmine >>> make fmat for each sector ..."
    write(mystd,*)
    call atomic_make_fmat_sectors()

    write(mystd, "(2X,a)") "jasmine >>> write eigenvalue, eigenvector, and atom.cix to files ... "
    write(mystd,*)
    ! write eigenvalue to file 'atom.eigval.dat'
    call atomic_write_eigval_sectors()

    ! write eigenvector to file 'atom.eigvec.dat'
    call atomic_write_eigvec_sectors()

    ! write information of sectors to file 'atom.cix'
    call atomic_write_atomcix_sectors()

    return
end subroutine atomic_solve_sectors
