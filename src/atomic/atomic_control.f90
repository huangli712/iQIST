!!!-----------------------------------------------------------------------
!!! project : iqist @ jasmine
!!! program : control module
!!!           version module
!!! source  : atomic_control.f90
!!! type    : module
!!! author  : yilin wang (email:qhwyl2006@126.com)
!!! history : 07/09/2014 by yilin wang (created)
!!!           01/03/2024 by li huang (last modified)
!!! purpose : define global control parameters for the atomic eigenvalue
!!!           problem solver
!!! status  : unstable
!!! comment :
!!!-----------------------------------------------------------------------

!!========================================================================
!!>>> module control                                                   <<<
!!========================================================================

!!
!! @mod control
!!
!! define the control parameters and dimensional parameters.
!!
  module control
     use constants, only : dp

     implicit none

!!========================================================================
!!>>> character variables                                              <<<
!!========================================================================

!!
!! @var cname
!!
!! code name of the current atomic eigenvalue problem solver
!!
     character(len = 07), public, save :: cname = 'JASMINE'

!!========================================================================
!!>>> integer variables                                                <<<
!!========================================================================

!!
!! @var ibasis
!!
!! control flag, where is the source for natural basis (the eigenstate of
!! crystal field + spin-orbital coupling)
!!
!! if ibasis == 1:
!!     make natural basis inside of this program
!!
!! if ibasis == 2:
!!     make natural basis outside of this program
!!
     integer, public, save :: ibasis = 1

!!
!! @var ictqmc
!!
!! control flag, type of atomic Hamiltonian matrix diagonalization
!!
!! if ictqmc == 0:
!!     direct diagonalization in full Hilbert space (for camellia code)
!!
!! if ictqmc == 1:
!!     direct diagonalization in full Hilbert space (for begonia and lavender codes)
!!
!! if ictqmc == 2:
!!     good quantum numbers: N
!!
!! if ictqmc == 3:
!!     good quantum numbers: N, Sz
!!
!! if ictqmc == 4:
!!     good quantum numbers: N, Sz, PS
!!
!! if ictqmc == 5:
!!     good quantum numbers: N, Jz
!!
     integer, public, save :: ictqmc = 1

!!
!! @var icu
!!
!! control flag, type of Coulomb interaction matrix
!!
!! if icu == 1:
!!     Kanamori parameters (Uc, Uv, Jz, Js, Jp), isotropic Hund's rule coupling
!!
!! if icu == 2:
!!     Slater-Cordon parameters (Ud, Jh => F0, F2, F4, F6)
!!
!! if icu == 3:
!!     Kanamori parameters (Uc, Uv, Jz, Js, Jp), anisotropic Hund's rule coupling
!!
     integer, public, save :: icu    = 1

!!
!! @var icf
!!
!! control flag, type of crystal field (CF)
!!
!! if icf == 0:
!!     no crystal field
!!
!! if icf == 1:
!!     diagonal crystal field
!!
!! if icf == 2:
!!     non-diagonal crystal field
!!
     integer, public, save :: icf    = 0

!!
!! @var isoc
!!
!! control flag, type of spin-orbit coupling (SOC)
!!
!! if isoc == 0:
!!     no SOC
!!
!! if isoc == 1:
!!     onsite atomic SOC, H_soc = \lambda * L*S
!!
     integer, public, save :: isoc   = 0

!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

!!
!! @var nband
!!
!! number of correlated bands
!!
     integer, public, save :: nband  = 1

!!
!! @var nspin
!!
!! number of spin projections, it should not be changed
!!
     integer, public, save :: nspin  = 2

!!
!! @var norbs
!!
!! number of correlated orbitals (= nband * nspin)
!!
     integer, public, save :: norbs  = 2

!!
!! @var ncfgs
!!
!! number of many-body configurations, the dimension of Hilbert space
!!
     integer, public, save :: ncfgs  = 4

!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

!!
!! @var nmini
!!
!! the lower boundary of occupancy N
!!
     integer, public, save :: nmini = 0

!!
!! @var nmaxi
!!
!! the upper boundary of occupancy N
!!
     integer, public, save :: nmaxi = 2

!!========================================================================
!!>>> real variables                                                   <<<
!!========================================================================

!!
!! the following parameters are useful when icu = 1 or icu = 3
!!

!!
!! @var Uc
!!
!! intra-orbital Coulomb interaction
!!
     real(dp), public, save :: Uc    = 2.0_dp

!!
!! @var Uv
!!
!! inter-orbital Coulomb interaction
!!
     real(dp), public, save :: Uv    = 2.0_dp

!!
!! @var Jz
!!
!! Hund's exchange interaction in z axis
!!
     real(dp), public, save :: Jz    = 0.0_dp

!!
!! @var Js
!!
!! spin-flip interaction
!!
     real(dp), public, save :: Js    = 0.0_dp

!!
!! @var Jp
!!
!! pair-hopping interaction
!!
     real(dp), public, save :: Jp    = 0.0_dp

!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

!!
!! the following parameters are useful when icu = 2. they are used to
!! calculate the F0, F2, F4, and F6.
!!

!!
!! @var Ud
!!
!! Coulomb interaction parameter
!!
     real(dp), public, save :: Ud    = 2.0_dp

!!
!! @var Jh
!!
!! Hund's exchange parameter
!!
     real(dp), public, save :: Jh    = 0.0_dp

!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

!!
!! @var mune
!!
!! chemical potential, used to shift energy level. it is only useful for
!! model calculation
!!
     real(dp), public, save :: mune  = 0.0_dp

!!
!! @var lambda
!!
!! SOC strength
!!
     real(dp), public, save :: lambda= 0.0_dp

  end module control

!!========================================================================
!!>>> module version                                                   <<<
!!========================================================================

!!
!! @mod version
!!
!! define the semantic version string.
!!
  module version
     implicit none

!!
!! @var V_FULL
!!
!! version string, version number + date info. + status info.
!!
     character(len=20), public, parameter :: V_FULL = 'v0.8.2 @ 2023.07.06D'

!!
!! @var V_CURR
!!
!! version string, only version number
!!
     character(len=06), public, parameter :: V_CURR = 'v0.8.2'

!!
!! @var V_DATE
!!
!! version string, only date info.
!!
     character(len=11), public, parameter :: V_DATE = '2023.07.06'

!!
!! @var V_STAT
!!
!! version string, only status info., D means devel, T testing, R released.
!!
     character(len=01), public, parameter :: V_STAT = 'D'

!!
!! @var V_AUTH
!!
!! version string, author info.
!!
     character(len=11), public, parameter :: V_AUTH = 'by li huang'

!!
!! @var V_INST
!!
!! version string, affiliation info.
!!
     character(len=36), public, parameter :: V_INST = 'China Academy of Engineering Physics'

!!
!! @var V_MAIL
!!
!! version string, email info.
!!
     character(len=22), public, parameter :: V_MAIL = 'huangli@caep.cn'

!!
!! @var V_GPL3
!!
!! version string, license info.
!!
     character(len=36), public, parameter :: V_GPL3 = 'GNU General Public License version 3'

  end module version
