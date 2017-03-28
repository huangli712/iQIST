!!!-----------------------------------------------------------------------
!!! project : jasmine
!!! program : cat_init_atomic
!!!           cat_exec_atomic
!!!           cat_stop_atomic
!!! source  : atomic_open.f90
!!! type    : subroutines
!!! author  : li huang (email:lihuang.dmft@gmail.com)
!!! history : 08/12/2015 by li huang (created)
!!!           03/29/2017 by li huang (last modified)
!!! purpose : to provide necessary application programming interface for
!!!           the atomic eigenvalue problem solver
!!! status  : unstable
!!! comment :
!!!-----------------------------------------------------------------------

# if !defined (PYAPI)

!!>>> cat_init_atomic: initialize the atomic eigenvalue problem solver
!!>>> fortran version
  subroutine cat_init_atomic(I_solver)
     use capi, only : T_jasmine

     use m_cntr ! ALL
     use m_full, only : alloc_m_full_basis
     use m_spmat, only : alloc_m_spmat

     implicit none

! external arguments
! type structure of generic atomic eigenvalue problem solver
     type (T_jasmine), intent(in) :: I_solver

! print the running header
     call atomic_print_header()

! setup I_solver: integer parameters
     ibasis = I_solver%ibasis
     ictqmc = I_solver%ictqmc
     icu    = I_solver%icu
     icf    = I_solver%icf
     isoc   = I_solver%isoc
     nband  = I_solver%nband
     nspin  = I_solver%nspin
     norbs  = I_solver%norbs
     ncfgs  = I_solver%ncfgs
     nmini  = I_solver%nmini
     nmaxi  = I_solver%nmaxi

! setup I_solver: real parameters
     Uc     = I_solver%Uc
     Uv     = I_solver%Uv
     Jz     = I_solver%Jz
     Js     = I_solver%Js
     Jp     = I_solver%Jp
     Ud     = I_solver%Ud
     Jh     = I_solver%Jh
     mune   = I_solver%mune
     lambda = I_solver%lambda

! check validity of control parameters
     call atomic_check_config()

! print the summary of control parameters
     call atomic_print_summary()

! allocate memory for basis-related matrices
     call alloc_m_full_basis()

! allocate memory for single particle matrices
     call alloc_m_spmat()

     return
  end subroutine cat_init_atomic

# else   /* PYAPI */

!!>>> cat_init_atomic: initialize the atomic eigenvalue problem solver
!!>>> python version
  subroutine cat_init_atomic()
     use m_full, only : alloc_m_full_basis
     use m_spmat, only : alloc_m_spmat

     implicit none

! print the running header
     call atomic_print_header()

! set control parameters
     call atomic_config()

! check validity of control parameters
     call atomic_check_config()

! print the summary of control parameters
     call atomic_print_summary()

! allocate memory for basis-related matrices
     call alloc_m_full_basis()

! allocate memory for single particle matrices
     call alloc_m_spmat()

     return
  end subroutine cat_init_atomic

# endif  /* PYAPI */

!!>>> cat_exec_atomic: execute the atomic eigenvalue problem solver
  subroutine cat_exec_atomic()
     use constants, only : mystd

     use m_cntr, only : ictqmc

     implicit none

! make Fock basis for the full many particle Hiblert space
     write(mystd,'(2X,a)') 'make Fock basis'
     call atomic_make_fock()
     write(mystd,*)

! make single particle related matrices, including crystal field (CF),
! spin-orbit coupling (SOC), and Coulomb interaction U.
! when writing these matrices, we should define a single particle basis,
! there are four basis we will use (take 5-orbitals system for example)
! (1) real orbital basis
!     for example, |dz2,up>, |dz2,dn>,
!                  |dxz,up>, |dxz,dn>,
!                  |dyz,up>, |dyz,dn>,
!                  |dx2-y2,up>, |dx2-y2,dn>,
!                  |dxy,up>, |dxy,dn>
! (2) |lz,sz> complex orbital basis (the complex spherical functions)
!     for example, |2,-2,up>, |2,-2,dn>,
!                  |2,-1,up>, |2,-1,dn>,
!                  |2, 0,up>, |2, 0,dn>,
!                  |2, 1,up>, |2, 1,dn>,
!                  |2, 2,up>, |2, 2,dn>
! (3) |j2,jz> orbital basis (eigenstates of j2, jz)
!     for example, |3/2,-3/2>, |3/2,3/2>,
!                  |3/2,-1/2>, |3/2,1/2>,
!                  |5/2,-5/2>, |5/2,5/2>,
!                  |5/2,-3/2>, |5/2,3/2>,
!                  |5/2,-1/2>, |5/2,1/2>,
! (4) the so-called natural basis, on which the onsite energy of impurity
!     is diagonal. we have to diagonalize CF + SOC to obtain natural basis
! Note that the CF is always defined in real orbital basis, SOC is always
! defined in complex orbital basis, and Coulomb interaction U is defined
! in real orbital basis or complex orbital basis which depends on the form
! of Coulomb interaction, so we often need to transform them between two
! different basis sets
     write(mystd,'(2X,a)') 'make single particle matrices'
     call atomic_make_spmat()
     write(mystd,*)

! make natural basis
     write(mystd,'(2X,a)') 'make natural eigenbasis'
     call atomic_make_natural()
     write(mystd,*)

! call the drivers to perform different tasks
     select case(ictqmc)

! task 0: diagonalize the atomic Hamiltonian in full Hilbert space
         case (0)
             write(mystd,'(2X,a)') 'start full diagonalization'
             call atomic_f_driver()

! task 1: diagonalize the atomic Hamiltonian in full Hilbert space
         case (1)
             write(mystd,'(2X,a)') 'start full diagonalization'
             call atomic_f_driver()

! task 2: use good quantum numbers
! total number of electrons: N
! for the case of crystal field (CF) plus spin-orbital coupling (SOC)
         case (2)
             write(mystd,'(2X,a)') 'start sector-by-sector diagonalization (N)'
             call atomic_s_driver()

! task 3: use good quantum numbers
! total number of electrons: N
! z component of spin: Sz
! for the case without SOC and Slater parameterized Coulomb interaction
         case (3)
             write(mystd,'(2X,a)') 'start sector-by-sector diagonalization (N, Sz)'
             call atomic_s_driver()

! task 4: use good quantum numbers
! total number of electrons: N
! z component of spin: Sz
! PS number
! for the case without SOC and Kanamori parametrized Coulomb interaction
         case (4)
             write(mystd,'(2X,a)') 'start sector-by-sector diagonalization (N, Sz, PS)'
             call atomic_s_driver()

! task 5: use good quantum numbers
! total number of electrons: N
! z component of spin-orbit momentum: Jz
! for the case with SOC, and no CF
         case (5)
             write(mystd,'(2X,a)') 'start sector-by-sector diagonalization (N, Jz)'
             call atomic_s_driver()

         case default
             call s_print_error('atomic_main','this computational mode is not supported')

     end select

     return
  end subroutine cat_exec_atomic

!!>>> cat_stop_atomic: stop the atomic eigenvalue problem solver
  subroutine cat_stop_atomic()
     use m_full, only : dealloc_m_full_basis
     use m_spmat, only : dealloc_m_spmat

     implicit none

! deallocate memory
     call dealloc_m_spmat()
     call dealloc_m_full_basis()

! print footer
     call atomic_print_footer()

     return
  end subroutine cat_stop_atomic
