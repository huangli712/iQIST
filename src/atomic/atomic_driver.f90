!!!-----------------------------------------------------------------------
!!! project : jasmine
!!! program : atomic_f_driver
!!!           atomic_s_driver
!!! source  : atomic_driver.f90
!!! type    : subroutines
!!! author  : yilin wang (email:qhwyl2006@126.com)
!!! history : 07/09/2014 by yilin wang (created)
!!!           08/17/2015 by li huang (last modified)
!!! purpose : kernel drivers for atomic eigenvalue problem solver
!!! status  : unstable
!!! comment :
!!!-----------------------------------------------------------------------

!!>>> atomic_f_driver: solve the atomic eigenvalue problem using full
!!>>> Hilbert space diagonalization
!!>>> note: the output files are only compatible with BEGONIA and LAVENDER
!!>>> commponents in the iQIST software package
  subroutine atomic_f_driver()
     use constants, only : dp, eps6, mystd

     use m_cntr, only : ncfgs
     use m_full, only : hmat, eval, evec
     use m_full, only : alloc_m_full, dealloc_m_full

     implicit none

! local variables
! starting time
     real(dp) :: time_begin

! ending time
     real(dp) :: time_end

! allocate memory for global variables
     write(mystd,'(2X,a)') 'allocate memory for global variables in full Hilbert space'
     call cpu_time(time_begin) ! record starting time
     call alloc_m_full()
     call cpu_time(time_end)   ! record ending   time
     write(mystd,'(2X,a,f10.3,a)') 'time:', time_end - time_begin, 's'
     write(mystd,*)

! build atomic many particle Hamiltonian matrix
     write(mystd,'(2X,a)') 'assemble atomic many particle Hamiltonian'
     call cpu_time(time_begin) ! record starting time
     call atomic_make_fhmat()
     call cpu_time(time_end)   ! record ending   time
     write(mystd,'(2X,a,f10.3,a)') 'time:', time_end - time_begin, 's'
     write(mystd,*)

! check whether the many particle Hamiltonian is real
     write(mystd,"(2X,a)") 'check whether the atomic Hamiltonian is real or not'
     call cpu_time(time_begin) ! record starting time
     if ( any( abs( aimag(hmat) ) > eps6 ) ) then
         call s_print_error('atomic_f_driver','hmat is not real!')
     endif ! back if ( any( abs( aimag(hmat) ) > eps6 ) ) block
     call cpu_time(time_end)   ! record ending   time
     write(mystd,'(2X,a,f10.3,a)') 'time:', time_end - time_begin, 's'
     write(mystd,*)

! diagonalize hmat
     write(mystd,'(2X,a)') 'diagonalize the atomic Hamiltonian'
     call cpu_time(time_begin) ! record starting time
     call s_eig_sy(ncfgs, ncfgs, real(hmat), eval, evec)
     call cpu_time(time_end)   ! record ending   time
     write(mystd,'(2X,a,f10.3,a)') 'time:', time_end - time_begin, 's'
     write(mystd,*)

! build F-matrix
! first, build fmat of annihilation operators in Fock basis
! then, transform them to the eigen basis
     write(mystd,'(2X,a)') 'build F-matrix for annihilation fermion operators'
     call cpu_time(time_begin) ! record starting time
     call atomic_make_ffmat()
     call cpu_time(time_end)   ! record ending   time
     write(mystd,'(2X,a,f10.3,a)') 'time:', time_end - time_begin, 's'
     write(mystd,*)

! build occupancy number
     write(mystd,'(2X,a)') 'compute occupancy number of atomic eigenstates'
     call cpu_time(time_begin) ! record starting time
     call atomic_make_foccu()
     call cpu_time(time_end)   ! record ending   time
     write(mystd,'(2X,a,f10.3,a)') 'time:', time_end - time_begin, 's'
     write(mystd,*)

! build Sz
     write(mystd,'(2X,a)') 'compute magnetic moment of atomic eigenstates'
     call cpu_time(time_begin) ! record starting time
     call atomic_make_fspin()
     call cpu_time(time_end)   ! record ending   time
     write(mystd,'(2X,a,f10.3,a)') 'time:', time_end - time_begin, 's'
     write(mystd,*)

! write eigenvalues of hmat to file 'atom.eigval.dat'
! write eigenvectors of hmat to file 'atom.eigvec.dat'
! write eigenvalue of hmat, occupany number of eigenstates and fmat of
! annihilation fermion operators to file 'atom.cix'
     write(mystd,'(2X,a)') 'write eigenvalue, eigenvector, and F-matrix to files'
     call cpu_time(time_begin) ! record starting time
     call atomic_dump_feigval()
     call atomic_dump_feigvec()
     call atomic_dump_fcix()
     call cpu_time(time_end)   ! record ending   time
     write(mystd,'(2X,a,f10.3,a)') 'time:', time_end - time_begin, 's'
     write(mystd,*)

! deallocate memory
     write(mystd,'(2X,a)') 'deallocate memory for global variables in full Hilbert space'
     call cpu_time(time_begin) ! record starting time
     call dealloc_m_full()
     call cpu_time(time_end)   ! record ending   time
     write(mystd,'(2X,a,f10.3,a)') 'time:', time_end - time_begin, 's'
     write(mystd,*)

     return
  end subroutine atomic_f_driver

!!>>> atomic_s_driver: solve the atomic eigenvalue problem using good
!!>>> quantum numbers (GQNs) algorithm and sector-by-sector diagonalization
!!>>> note: the output files are only compatible with PANSY and MANJUSHAKA
!!>>> commponents in the iQIST software package
  subroutine atomic_s_driver()
     use constants, only : dp, eps6, mystd

     use m_sector, only : nsectors
     use m_sector, only : sectors
     use m_sector, only : dealloc_m_sector

     implicit none

! local variables
! loop index
     integer :: i

! starting time
     real(dp) :: time_begin

! ending time
     real(dp) :: time_end

! make all the sectors, allocate sectors memory inside
     write(mystd,'(2X,a)') 'determine sectors using good quantum numbers'
     call cpu_time(time_begin) ! record starting time
     call atomic_make_sectors()
     call cpu_time(time_end)   ! record ending   time
     write(mystd,'(2X,a,f10.3,a)') 'time:', time_end - time_begin, 's'
     write(mystd,*)

! make atomic Hamiltonian
     write(mystd,'(2X,a)') 'assemble atomic Hamiltonian for all sectors'
     call cpu_time(time_begin) ! record starting time
     call atomic_make_shmat()
     call cpu_time(time_end)   ! record ending   time
     write(mystd,'(2X,a,f10.3,a)') 'time:', time_end - time_begin, 's'
     write(mystd,*)

! check whether the many particle Hamiltonian is real
     write(mystd,'(2X,a)') 'check whether the atomic Hamiltonian is real or not'
     call cpu_time(time_begin) ! record starting time
     do i=1,nsectors
         if ( any( abs( aimag(sectors(i)%hmat) ) > eps6 ) ) then
             call s_print_error('atomic_s_driver','hmat is not real!')
         endif ! back if ( any( abs( aimag(sectors(i)%hmat) ) > eps6 ) ) block
     enddo ! over i={1,nsectors} loop
     call cpu_time(time_end)   ! record ending   time
     write(mystd,'(2X,a,f10.3,a)') 'time:', time_end - time_begin, 's'
     write(mystd,*)

! diagonalize Hamiltonian of each sector one by one
     write(mystd,'(2X,a)') 'diagonalize the atomic Hamiltonian for all sectors'
     call cpu_time(time_begin) ! record starting time
     call atomic_diag_shmat()
     call cpu_time(time_end)   ! record ending   time
     write(mystd,'(2X,a,f10.3,a)') 'time:', time_end - time_begin, 's'
     write(mystd,*)

! make F-matrix of both creation and annihilation operators for each sector
     write(mystd,'(2X,a)') 'build F-matrix for all sectors'
     call cpu_time(time_begin) ! record starting time
     call atomic_make_sfmat()
     call cpu_time(time_end)   ! record ending   time
     write(mystd,'(2X,a,f10.3,a)') 'time:', time_end - time_begin, 's'
     write(mystd,*)

! write eigenvalues to file 'atom.eigval.dat'
! write eigenvectors to file 'atom.eigvec.dat'
! write information of sectors to file 'atom.cix'
     write(mystd,'(2X,a)') 'write eigenvalue, eigenvector, and F-matrix to files'
     call cpu_time(time_begin) ! record starting time
     call atomic_dump_seigval()
     call atomic_dump_seigvec()
     call atomic_dump_scix()
     call cpu_time(time_end)   ! record ending   time
     write(mystd,'(2X,a,f10.3,a)') 'time:', time_end - time_begin, 's'
     write(mystd,*)

! deallocate memory
     write(mystd,'(2X,a)') 'deallocate memory for global variables in sectors'
     call cpu_time(time_begin) ! record starting time
     call dealloc_m_sector()
     call cpu_time(time_end)   ! record ending   time
     write(mystd,'(2X,a,f10.3,a)') 'time:', time_end - time_begin, 's'
     write(mystd,*)

     return
  end subroutine atomic_s_driver
