!!!=========+=========+=========+=========+=========+=========+=========+!
!!! JASMINE @ iQIST                                                      !
!!!                                                                      !
!!! An atomic eigenvalue problem solver which is used to generate input  !
!!! file (atom.cix) for the hybridization expansion version continuous   !
!!! time quantum Monte Carlo (CTQMC) quantum impurity solver             !
!!! author  : Yilin Wang (at IOP/CAS)                                    !
!!!           Li Huang (at IOP/CAS & SPCLab/CAEP & UNIFR)                !
!!! version : v2014.10.20T                                               !
!!! status  : WARNING: IN TESTING STAGE, USE IT IN YOUR RISK             !
!!! comment : the atomic solver is based on Dr. Liang Du's rambutan code !
!!!           any question, please contact with huangli712@gmail.com     !
!!!=========+=========+=========+=========+=========+=========+=========+!

!!
!!
!! WARNING
!! =======
!!
!! If you want to obtain an executable program, please go to src/build/,
!! open make.sys and comment out the API flag. On the contrary, if you
!! want to compile jasmine as a library, please activate the API flag.
!!
!! Introduction
!! ============
!!
!! Usage
!! =====
!!
!! # ./atomic or bin/jasmine.x
!!
!! Input
!! =====
!!
!! Output
!! ======
!!
!! Documents
!! =========
!!
!! For more details, please go to iqist/doc/manual directory.
!!
!!
  program atomic_main
     use constants, only : mystd

     use control, only : ictqmc
     use m_full, only : alloc_m_full_basis
     use m_full, only : dealloc_m_full_basis
     use m_spmat, only : alloc_m_spmat
     use m_spmat, only : dealloc_m_spmat

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

! make Single Particle related MATrix
! including crystal field (CF), spin-orbit coupling (SOC), Coulomb interaction U,
! when writing these matrices, we should define a single particle basis,
! there are four basis we will use (take 5-orbitals system for example)
! (1) real orbital basis
!     for example, |dz2,up>, |dz2,dn>, |dxz,up>, |dxz,dn>, |dyz,up>, |dyz,dn>,
!                  |dx2-y2,up>, |dx2-y2,dn>, |dxy,up>, |dxy,dn>
! (2) |lz,sz> complex orbital basis (the complex spherical functions)
!     for example, |2,-2,up>, |2,-2,dn>, |2,-1,up>, |2,-1,dn>, |2,0,up>, |2,0,dn>,
!                  |2, 1,up>, |2, 1,dn>, |2, 2,up>, |2, 2,dn>
! (3) |j2,jz> orbital basis (eigenstates of j^2, jz)
!     for example, |3/2,-3/2>, |3/2,-1/2>, |3/2, 1/2>, |3/2,3/2>,
!                  |5/2,-5/2>, |5/2,-3/2>, |5/2,-1/2>, |5/2,1/2>, |5/2,3/2>, |5/2,5/2>
! (4) the so-called natural basis, on which the on-site energy of impurity
!     is diagonal. we have to diagonalize CF + SOC to obtain natural basis
     write(mystd,'(2X,a)') 'JASMINE >>> prepare basis set and single particle matrix'
     write(mystd,'(2X,a)') 'make crystal field, spin-orbital coupling, and Coulomb interaction U'
     call atomic_make_spmat()
     write(mystd,*)

! make Fock basis for the full many particle Hiblert space
     write(mystd,'(2X,a)') 'make Fock basis'
     call atomic_make_fock()
     write(mystd,*)

! make natural basis
     write(mystd,'(2X,a)') "make natural basis"
     call atomic_make_natural()
     write(mystd,*)

! call the drivers to perform different tasks
     select case(ictqmc)

! task 1: diagonalize the atomic Hamiltonian in full Hilbert space
         case (1)
             write(mystd,"(2X,a)") "JASMINE >>> using direct diagonalization"
             call atomic_f_driver()

! task 2: use good quantum numbers
! total number of electrons: N
! for the case of crystal field (CF) plus spin-orbital coupling (SOC)
         case (2)
             write(mystd,"(2X,a)") "JASMINE >>> using good quantum number N"
             call atomic_s_driver()

! task 3: use good quantum numbers
! total number of electrons: N
! z component of spin: Sz
! for the case without SOC and Slater parameterized Coulomb interaction
         case (3)
             write(mystd,"(2X,a)") "JASMINE >>> using good quantum numbers N, Sz"
             call atomic_s_driver()

! task 4: use good quantum numbers
! total number of electrons: N
! z component of spin: Sz
! PS number
! for the case without SOC and Kanamori parametrized Coulomb interaction
         case (4)
             write(mystd,"(2X,a)") "JASMINE >>>> using good quantum numbers N, Sz, PS"
             call atomic_s_driver()

! task 5: use good quantum numbers
! total number of electrons: N
! z component of spin-orbit momentum: Jz
! for the case with SOC, and no CF
         case (5)
             write(mystd,"(2X,a)") "JASMINE >>> using good quantum numbers N, Jz"
             call atomic_s_driver()

         case default
             call s_print_error('atomic_main','this computational mode is not supported')

     end select

! deallocate memory
     call dealloc_m_spmat()
     call dealloc_m_full_basis()

! print footer
     call atomic_print_footer()

  end program atomic_main
