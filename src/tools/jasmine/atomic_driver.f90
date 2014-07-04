!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
! project : jasmine
! program : atomic_driver.f90
! history : 07/05/2014
! authors : yilin wang (qhwyl2006@126.com)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
  subroutine atomic_driver_fullspace()
     use constants
     use control

     implicit none
  
  end subroutine atomic_driver_fullspace
  
  subroutine atomic_driver_n()
     use constants
     use control

     implicit none

     write(mystd, "(2X, a)") ">>> hahaha, atomic_driver_n has not been implemented!"   

     return
  end subroutine atomic_driver_n

  subroutine atomic_driver_nszps()
     use constants
     use control

     implicit none

     write(mystd, "(2X, a)") ">>> hahaha, atomic_driver_nszps has not been implemented!"   

     return
  end subroutine atomic_driver_nszps

  subroutine atomic_driver_njz()
      use constants
      use control
      use mod_global
  
      use mod_subspace
      
      implicit none
  
      integer :: iorb
      integer :: ibas, jbas
      integer :: i, j
  
      ! allocate memory for global arrays and initialize them
      call atomic_setup_array()
  
      ! setup key matrix in atomic problems
      call atomic_jovian(nstat, cemat, somat, cumat)
  
      ! construct Fock basis sheet (just follow nhtong's note)
      call atomic_make_basis(norbs, ncfgs, ntots, nstat, basis, invcd, invsn)
  
      call atomic_make_good()
      call atomic_make_subspaces()
      call atomic_make_towhich()
      call atomic_trunk_space()
  
      ! setup transformation matrix from real orbitals to |j2,jz> single particle basis
      call atomic_make_amtrx(norbs, amtrx)
  
      ! construct atomic Hamiltonian matrix in natural many body basis
      call atomic_tran_cumat(norbs, amtrx, cumat, cumat_t)
      cumat = cumat_t
  
      ! build Hamiltonian matrix for each subspace
      call atomic_make_hmtrx()
  
      ! diagonalize each subspace
      call atomic_diag_hmtrx()
  
      ! build fmat
      call atomic_build_cfmat()
  
      call atomic_write_eigval()
      call atomic_write_eigvec()
      call atomic_write_fmat()
      call atomic_write_ctqmc()
  
      ! deallocate memory and finalize them
      call atomic_final_array()
  
      ! print the ending information for atomic problems
      call atomic_print_footer()
  
      return
  end subroutine atomic_driver_njz


