!-------------------------------------------------------------------------
! project : jasmine
! program : main
! source  : atomic_main.f90
! type    : program
! author  : yilin wang (email: qhwyl2006@126.com)
! history : 07/06/2014 by yilin wang
! purpose : the main program of atomic problem
! input   :
! output  :
! status  : unstable
! comment :
!-------------------------------------------------------------------------
program main
    use control, only ictqmc
    use m_basis_fullspace
    use m_spmat

    implicit none

    ! print the running header for atomic problems
    call atomic_print_header()

    ! setup important parameters for atomic problems
    call atomic_config()

    ! print the running summary for atomic problems
    call atomic_print_summary()

    ! make Single Particle related MATrix
    ! cystal field, spin-orbital coupling, Coulomb interaction U
    call atomic_make_spmat()

    ! make natural basis
    ! natural basis is the basis on which the on-site impurity 
    ! energy matrix is diagonal 
    call atomic_make_natural()

    ! make Fock BASIS for the FULL many particle Hiblert  space
    call atomic_make_basis_full()

    ! call the driver
    select case(ictqmc)
        ! itask 1: diagonalize the full Hilbert space
        case(1) call atomic_driver_fullspace()

        ! itask 2: use good quantum numbers
        ! total number of electrons: N
        ! for the case of crystal field (CF) plus spin-orbital coupling (SOC)
        case(2) call atomic_driver_n()

        ! itask 3: use good quantum numbers
        ! total number of electrons: N 
        ! spin: Sz
        ! PS number
        ! for the case without SOC
        case(3) call atomic_driver_nszps()

        ! itask 4: use good quantum numbers
        ! total number of electrons: N
        ! Jz
        ! for the case with SOC, and no CF
        case(4) call atomic_driver_njz() 

    end select 

    ! deallocate memory
    call dealloc_m_spmat()
    call dealloc_m_basis_fullspace()

    ! print footer
    call atomic_print_footer()

end program main
