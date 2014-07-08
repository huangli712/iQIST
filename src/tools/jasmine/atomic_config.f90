!-------------------------------------------------------------------------
! project : jasmine
! program : atomic_config
! source  : atomic_config.f90
! type    : subroutine
! author  : yilin wang (email: qhwyl2006@126.com)
! history : 07/06/2014 by yilin wang
! purpose : set config parameters 
! input   :
! output  :
! status  : unstable
! comment :
!-------------------------------------------------------------------------

!>>> read config parameters from file 'atom.config.in'
subroutine atomic_config()
    use constants
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
    norbs = 2            ! number of orbits
    ncfgs = 4            ! number of many-body configurations
 
    !----------------------------------------------------------------
    Uc = 4.00_dp         ! intraorbital Coulomb interaction
    Uv = 2.00_dp         ! interorbital Coulomb interaction
    Jz = 1.00_dp         ! Hund's exchange interaction
    Js = 1.00_dp         ! spin-flip interaction
    Jp = 1.00_dp         ! pair-hopping interaction

    !----------------------------------------------------------------
    F0 = 4.00_dp         ! F0
    F2 = 0.00_dp         ! F2
    F4 = 0.00_dp         ! F4
    F6 = 0.00_dp         ! F6

    !----------------------------------------------------------------
    lambda = 0.00_dp     ! spin-orbit coupling parameter

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
        read(mytmp, *)  itask
        read(mytmp, *)  ictqmc
        read(mytmp, *)  icf
        read(mytmp, *)  isoc
        read(mytmp, *)  icu

        !----------------------------------------------------------------
        read(mytmp, *) ! skip header
        read(mytmp, *)  nband
        read(mytmp, *)  nspin
        read(mytmp, *)  norbs
        read(mytmp, *)  ncfgs
 
        !----------------------------------------------------------------
        read(mytmp, *) ! skip header
        read(mytmp, *)  Uc
        read(mytmp, *)  Uv
        read(mytmp, *)  Jz
        read(mytmp, *)  Js
        read(mytmp, *)  Jp

        !----------------------------------------------------------------
        read(mytmp, *) ! skip header
        read(mytmp, *)  F0
        read(mytmp, *)  F2
        read(mytmp, *)  F4
        read(mytmp, *)  F6

        !----------------------------------------------------------------
        read(mytmp, *) ! skip header
        read(mytmp, *)  lambda
 
        close(mytmp)
    else
        call atomic_print_error('atomic_config', 'no file atom.config.in !')
    endif

    return
end subroutine atomic_config


