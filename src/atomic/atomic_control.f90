!!!-----------------------------------------------------------------------
!!! project : jasmine
!!! program : m_cntr     module
!!! source  : atomic_control.f90
!!! type    : module
!!! author  : yilin wang (email:qhwyl2006@126.com)
!!! history : 07/09/2014 by yilin wang (created)
!!!           08/17/2015 by li huang (last modified)
!!! purpose : define global control parameters for the atomic eigenvalue
!!!           problem solver
!!! status  : unstable
!!! comment :
!!!-----------------------------------------------------------------------

  module m_cntr
     use constants, only : dp

     implicit none

!!========================================================================
!!>>> integer variables                                                <<<
!!========================================================================

! control flag: where is the source for natural basis (the eigenstate of
! crystal field + spin-orbital coupling)
! 1: make natural basis inside of this program
! 2: make natural basis outside of this program
     integer, public, save :: ibasis = 1

! control flag: type of atomic Hamiltonian matrix diagonalization
! 0: direct diagonalization in full Hilbert space (for camellia code)
! 1: direct diagonalization in full Hilbert space (for begonia and lavender codes)
! 2: good quantum numbers: N
! 3: good quantum numbers: N, Sz
! 4: good quantum numbers: N, Sz, PS
! 5: good quantum numbers: N, Jz
     integer, public, save :: ictqmc = 1

! control flag: type of Coulomb interaction U
! 1: Kanamori parameters (Uc, Uv, Jz, Js, Jp), isotropic Hund's rule coupling
! 2: Slater-Cordon parameters (Ud, Jh => F0, F2, F4, F6)
! 3: Kanamori parameters (Uc, Uv, Jz, Js, Jp), anisotropic Hund's rule coupling
     integer, public, save :: icu    = 1

! control flag: type of crystal field (CF)
! 0: no crystal field
! 1: diagonal crystal field
! 2: non-diagonal crystal field
     integer, public, save :: icf    = 0

! control flag: type of spin-orbit coupling (SOC)
! 0: no SOC
! 1: onsite atomic SOC, H_soc = \lambda * L*S
     integer, public, save :: isoc   = 0

!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

! number of bands
     integer, public, save :: nband  = 1

! number of spins, it should not be changed
     integer, public, save :: nspin  = 2

! number of orbitals
     integer, public, save :: norbs  = 2

! number of many-body configurations, the dimension of Hilbert space
     integer, public, save :: ncfgs  = 4

!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

! the minimal total occupancy N will be kept
     integer, public, save :: nmini = 0

! the maximal total occupancy N will be kept
     integer, public, save :: nmaxi = 2

!!========================================================================
!!>>> real variables                                                   <<<
!!========================================================================

! the following parameters are useful when icu = 1 or icu = 3
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

! the following parameters are useful when icu = 2. they are used to
! calculate the F0, F2, F4, and F6.
! Coulomb parameter
     real(dp), public, save :: Ud    = 2.0_dp

! Hund's exchange parameter
     real(dp), public, save :: Jh    = 0.0_dp

!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

! chemical potential, used to shift energy level
! only useful for model calculation
     real(dp), public, save :: mune  = 0.0_dp

! SOC strength
     real(dp), public, save :: lambda= 0.0_dp

  end module m_cntr
