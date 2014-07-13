!-------------------------------------------------------------------------
! project : jasmine
! program : main
! source  : atomic_main.f90
! type    : program
! author  : yilin wang (email: qhwyl2006@126.com)
! history : 07/09/2014 by yilin wang
! purpose : the main program of jasmine
! input   :
! output  :
! status  : unstable
! comment :
!-------------------------------------------------------------------------
program main
    use constants,         only: mystd
    use control,           only: ictqmc
    use m_basis_fullspace, only: dealloc_m_basis_fullspace
    use m_spmat,           only: dealloc_m_spmat

    implicit none

    ! print the running header 
    call atomic_print_header()

    ! set control parameters
    write(mystd, "(2X,a)") "jasmine >>> set control parameters ..."
    write(mystd,*)
    call atomic_config()

    ! check validity of control parameters 
    write(mystd, "(2X,a)") "jasmine >>> checking validity of control parameters ..."
    write(mystd,*)
    call atomic_check_config()

    ! print the summary of control parameters 
    call atomic_print_summary()

    ! make Single Particle related MATrix
    ! including crystal field (CF), spin-orbital coupling (SOC), Coulomb interaction U
    ! when writing these matrices, we should define a single particle 
    ! orbital basis, there are four basis we will use
    ! (1) real orbital basis
    !     for example, |dz2,up>, |dz2,dn>, |dxz,up>, |dxz,dn>, |dyz,up>, |dyz,dn>, 
    !                  |dx2-y2,up>, |dx2-y2,dn>, |dxy,up>, |dxy,dn>
    ! (2) complex orbital (the complex spherical functions) basis
    !     for example, |2,-2,up>, |2,-2,dn>, |2,-1,up>, |2,-1,dn>, |2,0,up>, |2,0,dn>, 
    !                  |2,1,up>, |2,1,dn>, |2,2,up>, |2,2,dn>
    ! (3) |j2,jz> (eigenstates of j^2, jz) orbital basis
    !     for example, |3/2,-3/2>, |3/2,-1/2>, |3/2,1/2>, |3/2,3/2>, 
    !                  |5/2,-5/2>, |5/2,-3/2>,|5/2,-1/2>,|5/2,1/2>,|5/2,3/2>,|5/2,5/2>
    ! (4) the so-called natural basis, on which the on-site energy of impurity is diagonal
    !     we diagonalize CF + SOC to obtain natural basis
    write(mystd, "(2X,a)") "jasmine >>> make crystal field, spin-orbital coupling, and Coulomb interaction U ..."
    write(mystd,*)
    call atomic_make_spmat()

    ! make natural basis
    write(mystd, "(2X,a)") "jasmine >>> make natural basis, the Fock basis will be defined on it ..."
    write(mystd,*)
    call atomic_make_natural()

    ! make Fock basis for the full many particle Hiblert space
    write(mystd, "(2X,a)") "jasmine >>> make Fock basis ..."
    write(mystd,*)
    call atomic_make_basis_fullspace()

    ! call the drivers for different CTQMC algorithm
    select case(ictqmc)
        ! itask 1: diagonalize the full Hilbert space
        case(1) 
            write(mystd, "(2X,a)") "jasmine >>> CTQMC trace algorithm: full space matrices multiplication."
            write(mystd, *)
            call atomic_driver_fullspace()

        ! itask 2: use good quantum numbers
        ! total number of electrons: N
        ! for the case of crystal field (CF) plus spin-orbital coupling (SOC)
        case(2) 
            write(mystd, "(2X,a)") "jasmine >>> CTQMC trace algorithm: use good quantum number N."
            write(mystd, *)
            call atomic_driver_n()

        ! itask 3: use good quantum numbers
        ! total number of electrons: N 
        ! spin: Sz
        ! for the case without SOC
        case(3) 
            write(mystd, "(2X,a)") "jasmine >>> CTQMC trace algorithm: use good quantum numbers N, Sz."
            write(mystd, *)
            call atomic_driver_nsz()

        ! itask 4: use good quantum numbers
        ! total number of electrons: N 
        ! spin: Sz
        ! PS number
        ! for the case without SOC
        case(4) 
            write(mystd, "(2X,a)") "jasmine >>> CTQMC trace algorithm: use good quantum numbers N, Sz, PS."
            write(mystd, *)
            call atomic_driver_nszps()

        ! itask 5: use good quantum numbers
        ! total number of electrons: N
        ! Jz
        ! for the case with SOC, and no CF
        case(5) 
            write(mystd, "(2X,a)") "jasmine >>> CTQMC trace algorithm: use good quantum numbers N, Jz."
            write(mystd, *)
            call atomic_driver_njz() 

    end select 

    ! deallocate memory
    call dealloc_m_spmat()
    call dealloc_m_basis_fullspace()

    ! print footer
    call atomic_print_footer()

end program main
