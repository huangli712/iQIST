!!!=========+=========+=========+=========+=========+=========+=========+!
!!! iQIST @ JASMINE                                                      !
!!!                                                                      !
!!! An atomic eigenvalue problem solver which is used to generate input  !
!!! file (atom.cix) for the hybridization expansion version continuous   !
!!! time quantum Monte Carlo (CTQMC) quantum impurity solver             !
!!!                                                                      !
!!! author  : Yilin Wang (University of Science and Technology of China) !
!!!           Li Huang (China Academy of Engineering Physics)            !
!!! status  : WARNING: IN TESTING STAGE, USE IT IN YOUR RISK             !
!!! comment : the atomic solver is based on Dr. Liang Du's rambutan code !
!!!           any question, please contact with huangli@caep.cn          !
!!!=========+=========+=========+=========+=========+=========+=========+!

  program atomic_main
     use constants, only : mystd

     use control, only : ictqmc
     use m_fock, only : cat_alloc_fock_basis
     use m_fock, only : cat_free_fock_basis
     use m_spmat, only : cat_alloc_spmat
     use m_spmat, only : cat_free_spmat

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
     call cat_alloc_fock_basis()

! allocate memory for single particle matrices
     call cat_alloc_spmat()

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

! deallocate memory
     call cat_free_spmat()
     call cat_free_fock_basis()

! print footer
     call atomic_print_footer()

  end program atomic_main
