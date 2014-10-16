!!!-------------------------------------------------------------------------
!!! project : jasmine
!!! program : atomic_f_driver
!!!           atomic_s_driver
!!! source  : atomic_driver.f90
!!! type    : subroutines
!!! author  : yilin wang (email: qhwyl2006@126.com)
!!! history : 07/09/2014 by yilin wang
!!!           08/13/2014 by yilin wang
!!!           08/22/2014 by yilin wang
!!! purpose : solve atomic problem for different CTQMC trace algorithms
!!! status  : unstable
!!! comment :
!!!-------------------------------------------------------------------------

!!>>> atomic_f_driver: CTQMC direct matrices multiplications 
!!>>> trace algorithm, use full Hilbert space 
  subroutine atomic_f_driver()
     use constants, only : dp, mystd
     use control, only : ncfgs

     use m_glob_fullspace, only : hmat, hmat_eigval, hmat_eigvec
     use m_glob_fullspace, only : alloc_m_glob_fullspace, dealloc_m_glob_fullspace
  
     implicit none
  
! local variables
! a temp matrix
     real(dp) :: tmp_mat(ncfgs, ncfgs)

! whether the Hamiltonian is real ? 
     logical :: lreal 
  
! allocate memory 
     write(mystd, "(2X,a)") "jasmine >>> allocate memory of global variables for fullspace case ..."
     write(mystd,*)
     call alloc_m_glob_fullspace()
  
! build atomic many particle Hamiltonian matrix
     write(mystd, "(2X,a)") "jasmine >>> make atomic many particle Hamiltonian ..."
     write(mystd,*)
     call atomic_make_fhmat()
  
! check whether the many particle Hamiltonian is real 
     write(mystd, "(2X,a)") "jasmine >>> check whether Hamiltonian is real or not ..."
     write(mystd,*)
     call atomic_check_realmat(ncfgs, hmat, lreal)
     if (lreal .eqv. .false.) then
         call s_print_error('atomic_f_driver', 'hmat is not real !')
     else
         write(mystd, "(2X,a)") "jasmine >>> the atomic Hamiltonian is real"
         write(mystd,*)
     endif
  
! diagonalize hmat
     write(mystd, "(2X,a)") "jasmine >>> diagonalize the atomic Hamiltonian ..."
     write(mystd,*)
     tmp_mat = real(hmat)
     call s_eig_sy(ncfgs, ncfgs, tmp_mat, hmat_eigval, hmat_eigvec)
  
! build fmat
! first, build fmat of annihilation operators in Fock basis
! then, transform them to the eigen basis
     write(mystd, "(2X,a)") "jasmine >>> make fmat for annihilation fermion operators ... "
     write(mystd,*)
     call atomic_make_ffmat()

! build occupancy number
     write(mystd, "(2X,a)") "jasmine >>> make occupancy number of atomic eigenstates ... "
     write(mystd,*)
     call atomic_make_foccu()    

! write eigenvalues of hmat to file 'atom.eigval.dat'
     write(mystd, "(2X,a)") "jasmine >>> write eigenvalue, eigenvector, and atom.cix to files ..."
     write(mystd,*)
     call atomic_dump_feigval()
  
! write eigenvectors of hmat to file 'atom.eigvec.dat'
     call atomic_dump_feigvec()
   
! write eigenvalue of hmat, occupany number of eigenstates and 
! fmat of annihilation fermion operators to file "atom.cix"
! this is for begonia, lavender codes of iQIST package
     call atomic_dump_fcix()
  
! deallocate memory
     write(mystd, "(2X,a)") "jasmine >>> free memory of global variables for fullspace case ... "
     write(mystd,*)
     call dealloc_m_glob_fullspace()
  
     return
  end subroutine atomic_f_driver
  
!!>>> atomic_s_driver: CTQMC trace algorithm: use good quantum numbers (GQNs)
  subroutine atomic_s_driver()
     use constants, only : mystd
     use m_glob_sectors, only : nsectors, sectors, dealloc_m_glob_sectors
  
     implicit none

! local variables
! loop index
     integer :: i

! whether the Hamiltonian is real ?
     logical :: lreal

! make all the sectors, allocate m_glob_sectors memory inside
     write(mystd, "(2X,a)") "jasmine >>> determine sectors by good quantum numbers (GQNs)... "
     write(mystd,*)
     call atomic_mksectors()
 
! make atomic Hamiltonian
     write(mystd, "(2X,a)") "jasmine >>> make atomic Hamiltonian for each sector ... "
     write(mystd,*)
     call atomic_make_shmat()
  
! check whether the many particle Hamiltonian is real 
     write(mystd, "(2X,a)") "jasmine >>> check whether Hamiltonian is real or not ..."
     write(mystd,*)
     do i=1, nsectors 
         call atomic_check_realmat(sectors(i)%ndim, sectors(i)%myham, lreal)
         if (lreal .eqv. .false.) then
             call s_print_error('atomic_solve_sectors', 'hmat is not real !')
         endif
     enddo
     write(mystd, "(2X,a)") "jasmine >>> the Hamiltonian is real"
     write(mystd,*)
  
! diagonalize Hamiltonian of each sector one by one
     write(mystd, "(2X,a)") "jasmine >>> diagonalize atomic Hamiltonian for each sector ... "
     write(mystd,*)
     call atomic_diaghmat_sectors()
  
! make fmat of both creation and annihilation operators for each sector
     write(mystd, "(2X,a)") "jasmine >>> make fmat for each sector ..."
     write(mystd,*)
     call atomic_make_sfmat()
  
     write(mystd, "(2X,a)") "jasmine >>> write eigenvalue, eigenvector, and atom.cix to files ... "
     write(mystd,*)
! write eigenvalues to file 'atom.eigval.dat'
     call atomic_dump_seigval()
  
! write eigenvectors to file 'atom.eigvec.dat'
     call atomic_dump_seigvec()
  
! write information of sectors to file 'atom.cix'
     call atomic_dump_scix()

! free memory
     write(mystd, "(2X,a)") "jasmine >>> free memory for sectors case ..."
     write(mystd,*)
     call dealloc_m_glob_sectors()
 
     return
  end subroutine atomic_s_driver
