!!!-----------------------------------------------------------------------
!!! project : iqist @ jasmine
!!! program : atomic_dispatcher
!!!           atomic_f_driver
!!!           atomic_s_driver
!!! source  : atomic_driver.f90
!!! type    : subroutines
!!! author  : yilin wang (email:qhwyl2006@126.com)
!!! history : 07/09/2014 by yilin wang (created)
!!!           06/11/2024 by li huang (last modified)
!!! purpose : try to launch various computational tasks for the atomic
!!!           eigenvalue problem solver.
!!! status  : unstable
!!! comment :
!!!-----------------------------------------------------------------------

!!
!! @sub atomic_dispatcher
!!
!! launch various computational tasks. it is the core subroutine of the
!! atomic eigenvalue problem solver
!!
  subroutine atomic_dispatcher()
     use constants, only : mystd

     use control, only : ictqmc

     implicit none

!! local variables
     ! used to draw a dashed line
     character (len=1) :: dash(54)

!! [body

     ! setup dash
     dash = '-'

     write(mystd,'(2X,a)') 'start initialization'
     write(mystd,'(2X,54a1)') dash ! dashed line
     write(mystd,*)

     ! make Fock basis for the full many particle Hiblert space
     write(mystd,'(2X,a)') 'make Fock basis'
     call atomic_build_fock()
     write(mystd,*)

     ! make single particle matrix
     ! such as spin-orbit coupling and crystal field splitting
     write(mystd,'(2X,a)') 'make single particle matrix'
     call atomic_build_spmat()
     write(mystd,*)

     ! make natural eigenbasis. the onsite energy (CFS + SOC) matrix of
     ! impurity will become diagonal in this basis.
     write(mystd,'(2X,a)') 'make natural eigenbasis'
     call atomic_build_natural()
     write(mystd,*)

     ! launch the computational kernel
     select case (ictqmc)

         ! diagonalize the atomic Hamiltonian directly
         ! good quantum numbers -> N/A
         case (1)
             write(mystd,'(2X,a)') 'start exact diagonalization'
             write(mystd,'(2X,54a1)') dash ! dashed line
             call atomic_f_driver()

         ! subspace diagonalization by using good quantum numbers
         ! good quantum numbers -> N
         !
         ! crystal field splitting: no limits
         ! spin-orbit coupling: no limits
         ! Coulomb interaction: no limits
         case (2)
             write(mystd,'(2X,a)') 'start subspace diagonalization (N)'
             write(mystd,'(2X,54a1)') dash ! dashed line
             call atomic_s_driver()

         ! subspace diagonalization by using good quantum numbers
         ! good quantum numbers -> N, Sz
         !
         ! crystal field splitting: no limits
         ! spin-orbit coupling: no
         ! Coulomb interaction: no limits
         case (3)
             write(mystd,'(2X,a)') 'start subspace diagonalization (N, Sz)'
             write(mystd,'(2X,54a1)') dash ! dashed line
             call atomic_s_driver()

         ! subspace diagonalization by using good quantum numbers
         ! good quantum numbers -> N, Sz, PS
         !
         ! crystal field splitting: no limits
         ! spin-orbit coupling: no
         ! Coulomb interaction: Kanamori type
         case (4)
             write(mystd,'(2X,a)') 'start subspace diagonalization (N, Sz, PS)'
             write(mystd,'(2X,54a1)') dash ! dashed line
             call atomic_s_driver()

         ! subspace diagonalization by using good quantum numbers
         ! good quantum numbers -> N, Jz
         !
         ! crystal field splitting: no
         ! spin-orbit coupling: yes
         ! Coulomb interaction: no limits
         case (5)
             write(mystd,'(2X,a)') 'start subspace diagonalization (N, Jz)'
             write(mystd,'(2X,54a1)') dash ! dashed line
             call atomic_s_driver()

         ! subspace diagonalization by using good quantum numbers
         ! good quantum numbers -> N + Sz (or Jz) + AP
         !
         ! crystal field splitting: no limits
         ! spin-orbit coupling: no limits
         ! Coulomb interaction: no limits
         case (6)
             write(mystd,'(2X,a)') 'start subspace diagonalization (N, AP)'
             write(mystd,'(2X,54a1)') dash ! dashed line
             call atomic_s_driver()

         case default
             call s_print_error('atomic_dispatcher','this task is not supported')

     end select

!! body]

     return
  end subroutine atomic_dispatcher

!!
!! @sub atomic_f_driver
!!
!! solve the atomic eigenvalue problem by direct diagonalization in the
!! full Hilbert space. this subroutine is called when ictqmc == 1
!!
  subroutine atomic_f_driver()
     use constants, only : dp
     use constants, only : mystd

     use m_fock, only : cat_alloc_fock_eigen
     use m_fock, only : cat_free_fock_eigen

     implicit none

!! local variables
     ! starting time
     real(dp) :: time_begin

     ! ending time
     real(dp) :: time_end

!! [body

     ! allocate memory for atomic eigenstates
     write(mystd,*)
     write(mystd,'(2X,a)') 'allocate memory for atomic eigenstates'
     !
     call cpu_time(time_begin) ! record starting time
     call cat_alloc_fock_eigen()
     call cpu_time(time_end)   ! record ending   time
     !
     write(mystd,'(2X,a,f10.3,a)') 'time:', time_end - time_begin, 's'
     write(mystd,*)

     ! build the atomic Hamiltonian
     write(mystd,'(2X,a)') 'assemble atomic Hamiltonian'
     !
     call cpu_time(time_begin) ! record starting time
     call atomic_make_fhmat()
     call cpu_time(time_end)   ! record ending   time
     !
     write(mystd,'(2X,a,f10.3,a)') 'time:', time_end - time_begin, 's'
     write(mystd,*)

     ! check whether the atomic Hamiltonian is real
     write(mystd,"(2X,a)") 'check atomic Hamiltonian'
     !
     call cpu_time(time_begin) ! record starting time
     call atomic_check_fhmat()
     call cpu_time(time_end)   ! record ending   time
     !
     write(mystd,'(2X,a,f10.3,a)') 'time:', time_end - time_begin, 's'
     write(mystd,*)

     ! diagonalize the atomic Hamiltonian
     write(mystd,'(2X,a)') 'diagonalize atomic Hamiltonian'
     !
     call cpu_time(time_begin) ! record starting time
     call atomic_diag_fhmat()
     call cpu_time(time_end)   ! record ending   time
     !
     write(mystd,'(2X,a,f10.3,a)') 'time:', time_end - time_begin, 's'
     write(mystd,*)

     ! calculate annihilation operator matrix
     !
     ! at first, build annihilation operator matrix in the Fock basis
     ! and then, it should be transformed to the atomic eigenbasis
     write(mystd,'(2X,a)') 'compute f operator in atomic eigenbasis'
     !
     call cpu_time(time_begin) ! record starting time
     call atomic_make_ffmat()
     call cpu_time(time_end)   ! record ending   time
     !
     write(mystd,'(2X,a,f10.3,a)') 'time:', time_end - time_begin, 's'
     write(mystd,*)

     ! calculate occupancy matrix in atomic eigenbasis
     !
     ! the occupancy matrix is built in the Fock basis, and then is
     ! rotated to the atomic eigenbasis
     write(mystd,'(2X,a)') 'compute density matrix in atomic eigenbasis'
     !
     call cpu_time(time_begin) ! record starting time
     call atomic_make_foccu()
     call cpu_time(time_end)   ! record ending   time
     !
     write(mystd,'(2X,a,f10.3,a)') 'time:', time_end - time_begin, 's'
     write(mystd,*)

     ! calculate Sz in atomic eigenbasis
     !
     ! the Sz matrix is built in the Fock basis, and then is rotated to
     ! the atomic eigenbasis
     write(mystd,'(2X,a)') 'compute magnetic moment in atomic eigenbasis'
     !
     call cpu_time(time_begin) ! record starting time
     call atomic_make_fspin()
     call cpu_time(time_end)   ! record ending   time
     !
     write(mystd,'(2X,a,f10.3,a)') 'time:', time_end - time_begin, 's'
     write(mystd,*)

     ! write essential data to external files for reference
     !
     ! write eigenvalues of hmat to file 'atom.eigval.dat'
     ! write eigenvectors of hmat to file 'atom.eigvec.dat'
     ! write f operator to file 'atom.cix'
     write(mystd,'(2X,a)') 'write atomic eigenstates'
     !
     call cpu_time(time_begin) ! record starting time
     call atomic_dump_feigval()
     call atomic_dump_feigvec()
     call atomic_dump_fcix()
     call cpu_time(time_end)   ! record ending   time
     !
     write(mystd,'(2X,a,f10.3,a)') 'time:', time_end - time_begin, 's'
     write(mystd,*)

     ! deallocate memory for atomic eigenstates
     write(mystd,'(2X,a)') 'deallocate memory for atomic eigenstates'
     !
     call cpu_time(time_begin) ! record starting time
     call cat_free_fock_eigen()
     call cpu_time(time_end)   ! record ending   time
     !
     write(mystd,'(2X,a,f10.3,a)') 'time:', time_end - time_begin, 's'
     write(mystd,*)

!! body]

     return
  end subroutine atomic_f_driver

!!
!! @sub atomic_s_driver
!!
!! solve the atomic eigenvalue problem using good quantum numbers (GQNs)
!! algorithm and subspace diagonalization. this subroutine is called
!! when 2 <= ictqmc <= 6
!!
  subroutine atomic_s_driver()
     use constants, only : dp
     use constants, only : mystd

     use control, only : ictqmc

     use m_fock, only : cat_alloc_hmat_only
     use m_fock, only : cat_free_hmat_only

     use m_sector, only : cat_free_sectors

     implicit none

!! local variables
     ! starting time
     real(dp) :: time_begin

     ! ending time
     real(dp) :: time_end

!! [body

     ! for automatic partition algorithm only
     if ( ictqmc == 6 ) then
         ! allocate memory for atomic Hamiltonian
         write(mystd,*)
         write(mystd,'(2X,a)') 'allocate memory for atomic Hamiltonian'
         !
         call cpu_time(time_begin) ! record starting time
         call cat_alloc_hmat_only()
         call cpu_time(time_end)   ! record ending   time
         !
         write(mystd,'(2X,a,f10.3,a)') 'time:', time_end - time_begin, 's'
         write(mystd,*)

         ! build the atomic Hamiltonian
         write(mystd,'(2X,a)') 'assemble atomic Hamiltonian'
         !
         call cpu_time(time_begin) ! record starting time
         call atomic_make_fhmat()
         call cpu_time(time_end)   ! record ending   time
         !
         write(mystd,'(2X,a,f10.3,a)') 'time:', time_end - time_begin, 's'
         write(mystd,*)

         ! check whether the atomic Hamiltonian is real
         write(mystd,"(2X,a)") 'check atomic Hamiltonian'
         !
         call cpu_time(time_begin) ! record starting time
         call atomic_check_fhmat()
         call cpu_time(time_end)   ! record ending   time
         !
         write(mystd,'(2X,a,f10.3,a)') 'time:', time_end - time_begin, 's'
     endif ! back if ( ictqmc == 6 ) block

     ! make the subspaces (sectors)
     ! the memory will be allocated automatically
     write(mystd,*)
     write(mystd,'(2X,a)') 'construct subspaces for atomic eigenstates'
     !
     call cpu_time(time_begin) ! record starting time
     if ( ictqmc /= 6 ) then
         call atomic_make_sectors()
     else
         call automatic_partition()
     endif ! back if ( ictqmc /= 6 ) block
     call cpu_time(time_end)   ! record ending   time
     !
     write(mystd,'(2X,a,f10.3,a)') 'time:', time_end - time_begin, 's'
     write(mystd,*)

     ! build the atomic Hamiltonian
     write(mystd,'(2X,a)') 'assemble atomic Hamiltonian'
     !
     call cpu_time(time_begin) ! record starting time
     call atomic_make_shmat()
     call cpu_time(time_end)   ! record ending   time
     !
     write(mystd,'(2X,a,f10.3,a)') 'time:', time_end - time_begin, 's'
     write(mystd,*)

     ! check whether the atomic Hamiltonian is real
     write(mystd,'(2X,a)') 'check atomic Hamiltonian'
     !
     call cpu_time(time_begin) ! record starting time
     call atomic_check_shmat()
     call cpu_time(time_end)   ! record ending   time
     !
     write(mystd,'(2X,a,f10.3,a)') 'time:', time_end - time_begin, 's'
     write(mystd,*)

     ! diagonalize the atomic Hamiltonian
     write(mystd,'(2X,a)') 'diagonalize atomic Hamiltonian'
     !
     call cpu_time(time_begin) ! record starting time
     call atomic_diag_shmat()
     call cpu_time(time_end)   ! record ending   time
     !
     write(mystd,'(2X,a,f10.3,a)') 'time:', time_end - time_begin, 's'
     write(mystd,*)

     ! calculate both creation and annihilation operators matrices
     write(mystd,'(2X,a)') 'compute f^+ and f operators in atomic eigenbasis'
     !
     call cpu_time(time_begin) ! record starting time
     call atomic_make_sfmat()
     call cpu_time(time_end)   ! record ending   time
     !
     write(mystd,'(2X,a,f10.3,a)') 'time:', time_end - time_begin, 's'
     write(mystd,*)

     ! write essential data to external files
     !
     ! write eigenvalues of hmat to file 'atom.eigval.dat'
     ! write eigenvectors of hmat to file 'atom.eigvec.dat'
     ! write f^+ and f operators to file 'atom.cix'
     write(mystd,'(2X,a)') 'write atomic eigenstates'
     !
     call cpu_time(time_begin) ! record starting time
     call atomic_dump_seigval()
     call atomic_dump_seigvec()
     call atomic_dump_scix()
     call cpu_time(time_end)   ! record ending   time
     !
     write(mystd,'(2X,a,f10.3,a)') 'time:', time_end - time_begin, 's'
     write(mystd,*)

     ! deallocate memory for atomic eigenstates
     write(mystd,'(2X,a)') 'deallocate memory for atomic eigenstates'
     !
     call cpu_time(time_begin) ! record starting time
     call cat_free_sectors()
     call cpu_time(time_end)   ! record ending   time
     !
     write(mystd,'(2X,a,f10.3,a)') 'time:', time_end - time_begin, 's'
     write(mystd,*)

     ! for automatic partition algorithm only
     if ( ictqmc == 6 ) then
         ! deallocate memory for atomic Hamiltonian
         write(mystd,'(2X,a)') 'deallocate memory for atomic Hamiltonian'
         !
         call cpu_time(time_begin) ! record starting time
         call cat_free_hmat_only()
         call cpu_time(time_end)   ! record ending   time
         !
         write(mystd,'(2X,a,f10.3,a)') 'time:', time_end - time_begin, 's'
         write(mystd,*)
     endif ! back if ( ictqmc == 6 ) block

!! body]

     return
  end subroutine atomic_s_driver
