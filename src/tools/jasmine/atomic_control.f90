!!!-----------------------------------------------------------------------
!!! project : jasmine
!!! program : control    module
!!! source  : atomic_control.f90
!!! type    : module
!!! author  : yilin wang (email:qhwyl2006@126.com)
!!! history : 07/09/2014 by yilin wang
!!!           08/22/2014 by yilin wang
!!!           10/16/2014 by li huang
!!! purpose : define global control parameters for the atomic eigenvalue
!!!           problem solver
!!! status  : unstable
!!! comment :
!!!-----------------------------------------------------------------------

  module control
     use constants, only : dp

     implicit none

!!========================================================================
!!>>> integer variables                                                <<<
!!========================================================================

! control flag: type of task
! 1: make natural basis inside of this program
! 2: make natural basis outside of this program
     integer, public, save :: itask  = 1

! control flag: type of CTQMC trace algorithm
! 1: general matrices multiplication
! 2: good quantum numbers: N
! 3: good quantum numbers: N, Sz
! 4: good quantum numbers: N, Sz, PS
! 5: good quantum numbers: N, Jz
     integer, public, save :: ictqmc = 1

! control flag: type of Coulomb interaction U
! 1: Kanamori parameters (Uc, Uv, Jz, Js, Jp)
! 2: Slater-Cordon parameters (Ud, Jh => F0, F2, F4, F6)
     integer, public, save :: icu    = 1

! control flag: type of crystal field
! 0: no crystal field
! 1: diagonal crystal field
! 2: non-diagonal crystal field
     integer, public, save :: icf    = 0

! control flag: type of spin-orbit coupling (SOC)
! 0: no SOC
! 1: on-site atomic SOC, H_soc = \lambda * L*S
     integer, public, save :: isoc   = 0

!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

! number of bands
     integer, public, save :: nband  = 1

! number of spins
     integer, public, save :: nspin  = 2

! number of orbitals
     integer, public, save :: norbs  = 2

! number of many-body configurations, the dimension of Hilbert space
     integer, public, save :: ncfgs  = 4

!!========================================================================
!!>>> real variables                                                   <<<
!!========================================================================

! the following are useful when icu = 1
! intraorbital Coulomb interaction
     real(dp), public, save :: Uc    = 2.0_dp

! interorbital Coulomb interaction
     real(dp), public, save :: Uv    = 2.0_dp

! Hund's exchange interaction
     real(dp), public, save :: Jz    = 0.0_dp

! spin-flip interaction
     real(dp), public, save :: Js    = 0.0_dp

! pair-hopping interaction
     real(dp), public, save :: Jp    = 0.0_dp

!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

! the following are useful when icu = 2
! Coulomb parameters
     real(dp), public, save :: Ud    = 2.0_dp

! Hund's exchange parameters
     real(dp), public, save :: Jh    = 0.0_dp

!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

! chemical potential, used to shift energy level
! only useful for model calculation
     real(dp), public, save :: mune  = 0.0_dp

! SOC strength
     real(dp), public, save :: lambda= 0.0_dp

  end module control
