!!!-----------------------------------------------------------------------
!!! project : iqist @ jasmine
!!! program : control module
!!!           version module
!!! source  : atomic_control.f90
!!! type    : module
!!! author  : yilin wang (email:qhwyl2006@126.com)
!!! history : 07/09/2014 by yilin wang (created)
!!!           01/31/2024 by li huang (last modified)
!!! purpose : define global control parameters for the atomic eigenvalue
!!!           problem solver.
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
!! code name of the atomic eigenvalue problem solver
!!
     character(len = 07), public, save :: cname = 'JASMINE'

!!========================================================================
!!>>> integer variables                                                <<<
!!========================================================================

!!
!! @var ibasis
!!
!! control flag. how to build the natural eigenbasis (eigenstates of
!! crystal field splitting + spin-orbit coupling)
!!
!! if ibasis == 1:
!!     make natural eigenbasis inside of this program. the crystal field
!!     splitting and spin-orbit coupling are built separately
!!
!! if ibasis == 2:
!!     make natural eigenbasis outside of this program. the crystal field
!!     splitting and spin-orbit coupling are built as a whole
!!
     integer, public, save :: ibasis = 1

!!
!! @var ictqmc
!!
!! control flag. how to diagonalize the atomic Hamiltonian matrix
!!
!! if ictqmc == 1:
!!     direct diagonalization in full Hilbert space
!!
!! if ictqmc == 2:
!!     subspace diagonalization using good quantum numbers (N)
!!
!! if ictqmc == 3:
!!     subspace diagonalization using good quantum numbers (N and Sz)
!!
!! if ictqmc == 4:
!!     subspace diagonalization using good quantum numbers (N, Sz, and PS)
!!
!! if ictqmc == 5:
!!     subspace diagonalization using good quantum numbers (N and Jz)
!!
!! we note that the format of the atom.cix file exactly depends on the
!! ictqmc parameter. if ictqmc == 1, the generated atom.cix file is
!! only suitable for the lavender code. if ictqmc > 1, the generated
!! atom.cix file is only suitable for the manjushaka code. the two
!! atom.cix files are not compatible with each other. the lavender code
!! has already been deprecated, so we retain the option (ictqmc == 1)
!! only for internal reference
!!
     integer, public, save :: ictqmc = 1

!!
!! @var icu
!!
!! control flag. type of Coulomb interaction matrix
!!
!! if icu == 1:
!!     Kanamori type interaction. the Hund's rule coupling is isotropic.
!!     it needs the Uc, Uv, Jz, Js, Jp parameters to build the Coulomb
!!     interaction matrix
!!
!! if icu == 2:
!!     Slater-Cordon type interaction. it needs the Ud, Jh parameters to
!!     evaluate the Slater integrals (F0, F2, F4, F6). and they are used
!!     build the Coulomb interaction matrix
!!
!! if icu == 3:
!!     Kanamori type interaction. the Hund's rule coupling is anisotropic.
!!     it needs the Uc, Uv, Jz, Js, Jp parameters to build the Coulomb
!!     interaction matrix
!!
!! we note that if icu == 3, perhaps users need to modify the
!!     atomic_util.f90/atomic_make_hund()
!! subroutine to customize the Hund's rule coupling matrix
!!
     integer, public, save :: icu    = 1

!!
!! @var icf
!!
!! control flag. type of crystal field splitting (CFS)
!!
!! if icf == 0:
!!     without crystal field splitting
!!
!! if icf == 1:
!!     diagonal crystal field splitting
!!
!! if icf == 2:
!!     non-diagonal crystal field splitting
!!
     integer, public, save :: icf    = 0

!!
!! @var isoc
!!
!! control flag, type of spin-orbit coupling (SOC)
!!
!! if isoc == 0:
!!     without spin-orbit coupling
!!
!! if isoc == 1:
!!     onsite atomic spin-orbit coupling, H_{soc} = \lambda * L \cdot S
!!
     integer, public, save :: isoc   = 0

!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

     integer, public, save :: nband  = 1

     integer, public, save :: nspin  = 2

     integer, public, save :: norbs  = 2

     integer, public, save :: ncfgs  = 4

!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

     integer, public, save :: nmini = 0

     integer, public, save :: nmaxi = 2

!!========================================================================
!!>>> real variables                                                   <<<
!!========================================================================

     real(dp), public, save :: Uc    = 2.0_dp

     real(dp), public, save :: Uv    = 2.0_dp

     real(dp), public, save :: Jz    = 0.0_dp

     real(dp), public, save :: Js    = 0.0_dp

     real(dp), public, save :: Jp    = 0.0_dp

!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

     real(dp), public, save :: Ud    = 2.0_dp

     real(dp), public, save :: Jh    = 0.0_dp

!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

     real(dp), public, save :: mune  = 0.0_dp

!!
!! @var lambda
!!
!! strength of spin-orbit coupling
!!
     real(dp), public, save :: lambda= 0.0_dp

  end module control

!!========================================================================
!!>>> module version                                                   <<<
!!========================================================================

!!
!! @mod version
!!
!! define the semantic version string for the jasmine code.
!!
  module version
     implicit none

!!
!! @var V_FULL
!!
!! version string, version number + date info. + status info.
!!
     character(len=20), public, parameter :: V_FULL = 'v0.8.5 @ 2024.01.30D'

!!
!! @var V_CURR
!!
!! version string, only version number
!!
     character(len=06), public, parameter :: V_CURR = 'v0.8.5'

!!
!! @var V_DATE
!!
!! version string, only date info.
!!
     character(len=11), public, parameter :: V_DATE = '2024.01.30'

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
