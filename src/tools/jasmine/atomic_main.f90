!!!-------------------------------------------------------------------------
!!! project : jasmine
!!! program : main
!!! source  : atomic_main.f90
!!! type    : program
!!! author  : yilin wang (email: qhwyl2006@126.com)
!!! history : 07/09/2014 by yilin wang
!!! purpose : main program of jasmine
!!! input   :
!!! output  :
!!! status  : unstable
!!! comment :
!!!-------------------------------------------------------------------------

  program main
     use constants,         only : mystd
     use control,           only : ictqmc

     use m_full, only : alloc_m_basis_fullspace, dealloc_m_basis_fullspace
     use m_spmat,           only : alloc_m_spmat, dealloc_m_spmat
  
     implicit none
  
! print the running header 
     call atomic_print_header()
  
! set control parameters
     write(mystd, "(2X,a)") "jasmine >>> setting control parameters ..."
     write(mystd,*)
     call atomic_config()

! print the summary of control parameters 
     call atomic_print_summary()

! check validity of control parameters 
     write(mystd, "(2X,a)") "jasmine >>> checking validity of control parameters ..."
     write(mystd,*)
     call atomic_check_config()
     write(mystd,*)

! allocate memory for basis-related matrices
     call alloc_m_basis_fullspace()

! allocate memory for single particle matrices
     call alloc_m_spmat()
  
! make Single Particle related MATrix
! including crystal field (CF), spin-orbit coupling (SOC), Coulomb interaction U,
! when writing these matrices, we should define a single particle basis, 
! there are four basis we will use (take 5-orbitals system for example)
! (1) real orbital basis
!     for example, |dz2,up>, |dz2,dn>, |dxz,up>, |dxz,dn>, |dyz,up>, |dyz,dn>, 
!                  |dx2-y2,up>, |dx2-y2,dn>, |dxy,up>, |dxy,dn>
! (2) complex orbital (the complex spherical functions) basis
!     for example, |2,-2,up>, |2,-2,dn>, |2,-1,up>, |2,-1,dn>, |2,0,up>, |2,0,dn>, 
!                  |2, 1,up>, |2, 1,dn>, |2, 2,up>, |2, 2,dn>
! (3) |j2,jz> (eigenstates of j^2, jz) orbital basis
!     for example, |3/2,-3/2>, |3/2,-1/2>, |3/2, 1/2>, |3/2,3/2>, 
!                  |5/2,-5/2>, |5/2,-3/2>, |5/2,-1/2>, |5/2,1/2>, |5/2,3/2>, |5/2,5/2>
! (4) the so-called natural basis, on which the on-site energy of impurity is diagonal
!     we diagonalize CF + SOC to obtain natural basis
     write(mystd, "(2X,a)") "jasmine >>> make crystal field, spin-orbital coupling, &
                                                       and Coulomb interaction U ..."
     write(mystd,*)
     call atomic_make_spmat()
  
! make natural basis
     write(mystd, "(2X,a)") "jasmine >>> make natural basis, the Fock basis will be defined on it ..."
     write(mystd,*)
     call atomic_make_natural()
  
! make Fock basis for the full many particle Hiblert space
     write(mystd, "(2X,a)") "jasmine >>> make Fock basis ..."
     write(mystd,*)
     call atomic_make_fock()
  
! call the drivers for different CTQMC algorithm
     select case(ictqmc)
! itask 1: diagonalize the full Hilbert space
         case(1) 
             write(mystd, "(2X,a)") "jasmine >>> CTQMC trace algorithm: full space &
                                                           matrices multiplications."
             write(mystd, *)
             call atomic_f_driver()
  
! itask 2: use good quantum numbers
! total number of electrons: N
! for the case of crystal field (CF) plus spin-orbital coupling (SOC)
         case(2) 
             write(mystd, "(2X,a)") "jasmine >>> CTQMC trace algorithm: use good quantum number N."
             write(mystd, *)
             call atomic_s_driver()
  
! itask 3: use good quantum numbers
! total number of electrons: N 
! z component of spin: Sz
! for the case without SOC and Slater parameterized Coulomb interaction
         case(3) 
             write(mystd, "(2X,a)") "jasmine >>> CTQMC trace algorithm: use good quantum numbers N, Sz."
             write(mystd, *)
             call atomic_s_driver()
  
! itask 4: use good quantum numbers
! total number of electrons: N 
! z component of spin: Sz
! PS number
! for the case without SOC and Kanamori parametrized Coulomb interaction
         case(4) 
             write(mystd, "(2X,a)") "jasmine >>> CTQMC trace algorithm: use good quantum numbers N, Sz, PS."
             write(mystd, *)
             call atomic_s_driver()
  
! itask 5: use good quantum numbers
! total number of electrons: N
! z component of spin-orbit momentum: Jz
! for the case with SOC, and no CF
         case(5) 
             write(mystd, "(2X,a)") "jasmine >>> CTQMC trace algorithm: use good quantum numbers N, Jz."
             write(mystd, *)
             call atomic_s_driver()
  
     end select 
  
! deallocate memory
     call dealloc_m_spmat()
     call dealloc_m_basis_fullspace()
  
! print footer
     call atomic_print_footer()
  
  end program main
