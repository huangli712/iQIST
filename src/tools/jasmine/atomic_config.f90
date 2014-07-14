!-------------------------------------------------------------------------
! project : jasmine
! program : atomic_config
! source  : atomic_config.f90
! type    : subroutines
! author  : yilin wang (email: qhwyl2006@126.com)
! history : 07/09/2014 by yilin wang
! purpose : set config parameters 
! input   :
! output  :
! status  : unstable
! comment :
!-------------------------------------------------------------------------

!>>> read config parameters from file 'atom.config.in'
subroutine atomic_config()
    use constants,  only: dp, mytmp
    use control

    implicit none

    ! local variables
    logical :: exists
    
    !----------------------------------------------------------------
    itask  = 1           ! type of task
    ictqmc = 1           ! type of CTQMC algorithm
    icf    = 0           ! type of crystal field
    isoc   = 0           ! type of spin-orbital coupling (SOC)
    icu    = 1           ! type of Coulomb interaction
    
    !---------------------------------------------------------------- 
    nband = 1            ! number of bands
    nspin = 2            ! number of spins
    norbs = nband*nspin  ! number of orbits
    ncfgs = 2**norbs     ! number of many-body configurations
 
    !----------------------------------------------------------------
    Uc = 4.00_dp         ! intraorbital Coulomb interaction
    Uv = 2.00_dp         ! interorbital Coulomb interaction
    Jz = 1.00_dp         ! Hund's exchange interaction
    Js = 1.00_dp         ! spin-flip interaction
    Jp = 1.00_dp         ! pair-hopping interaction

    !----------------------------------------------------------------
    F0 = 4.00_dp         ! F0
    F2 = 6.00_dp         ! F2
    F4 = 3.75_dp         ! F4
    F6 = 0.00_dp         ! F6

    !----------------------------------------------------------------
    lambda = 0.00_dp     ! spin-orbit coupling parameter
    mune   = 0.00_dp     ! chemical potential

    ! read from input file if it exists
    exists = .false.

    ! inquire the input file status
    inquire( file="atom.config.in", exist=exists )

    ! read parameters from dft.atom.in
    if ( exists .eqv. .true. ) then
        open( mytmp, file="atom.config.in")
        read(mytmp, *) ! skip header
        !----------------------------------------------------------------
        read(mytmp, *) ! skip header
        read(mytmp, *)  nband
        norbs = nband * nspin
        ncfgs = 2 ** norbs 

        read(mytmp, *)  itask
        read(mytmp, *)  ictqmc
        read(mytmp, *)  icf
        read(mytmp, *)  isoc
        read(mytmp, *)  icu
        !----------------------------------------------------------------
        read(mytmp, *) ! skip header
        read(mytmp, *)  Uc
        read(mytmp, *)  Uv
        read(mytmp, *)  Jz
        read(mytmp, *)  Js
        read(mytmp, *)  Jp
        !----------------------------------------------------------------
        read(mytmp, *) ! skip header
        read(mytmp, *)  Ud
        read(mytmp, *)  JH

        F0 = Ud
        if (nband == 5) then
            F2 = JH * 14.0_dp / 1.625_dp 
            F4 = 0.625_dp * F2
        elseif(nband == 7) then
            F2 = JH * 6435.0_dp / (286.0_dp + (195.0_dp * 451.0_dp / 675.0_dp) + (250.0_dp * 1001.0_dp / 2025.0_dp))
            F4 = 451.0_dp / 675.0_dp * F2
            F6 = 1001.0_dp / 2025.0_dp * F2
        endif
        !----------------------------------------------------------------
        read(mytmp, *) ! skip header
        read(mytmp, *)  lambda
        read(mytmp, *)  mune
 
        close(mytmp)
    else
        call atomic_print_error('atomic_config', 'no file atom.config.in !')
    endif

    return
end subroutine atomic_config
