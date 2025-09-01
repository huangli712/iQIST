!!!-----------------------------------------------------------------------
!!! project : iqist @ narcissus
!!! program : control module
!!!           version module
!!! source  : ctqmc_control.f90
!!! type    : module
!!! author  : li huang (email:huangli@caep.cn)
!!! history : 09/15/2009 by li huang (created)
!!!           05/30/2025 by li huang (last modified)
!!! purpose : define global control parameters for hybridization expansion
!!!           version continuous time quantum Monte Carlo (CTQMC) quantum
!!!           impurity solver and dynamical mean field theory (DMFT) self-
!!!           consistent engine.
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
!! code name of the current quantum impurity solver
!!
     character(len = 09), public, save :: cname = 'NARCISSUS'

!!========================================================================
!!>>> integer variables                                                <<<
!!========================================================================

!!
!! @var isscf
!!
!! control flag, define the running scheme of the code
!!
!! if isscf == 1:
!!     one-shot non-self-consistent scheme, usually used in the density
!!     functional theory plus dynamical mean field theory case or used
!!     to solve the quantum impurity model separately
!!
!! if isscf == 2:
!!     self-consistent scheme, used in the dynamical mean field theory
!!     case. the code implements a typical dynamical mean field theory
!!     self-consistent loop for solving the Hubbard model in the bethe
!!     lattice (semicircular density of state)
!!
     integer, public, save :: isscf  = 1

!!
!! @var isscr
!!
!! control flag, define whether the Coulomb interaction U is dynamic
!!
!! if isscr == 1:
!!     static interaction
!!
!! if isscr == 2:
!!     dynamic interaction, for plasmon pole model
!!
!! if isscr == 3:
!!     dynamic interaction, for ohmic model
!!
!! if isscr == 4:
!!     dynamic interaction, for realistic materials
!!
     integer, public, save :: isscr  = 1

!!
!! @var isbnd
!!
!! control flag, define symmetry of the impurity model (band part)
!!
!! if isbnd == 1:
!!     the bands are not symmetrized
!!
!! if isbnd == 2:
!!     the bands are symmetrized according to symmetry matrix
!!
     integer, public, save :: isbnd  = 1

!!
!! @var isspn
!!
!! control flag, define symmetry of the impurity model (spin part)
!!
!! if isspn == 1:
!!     let spin up and spin down states evolve independently
!!
!! if isspn == 2:
!!     enforce spin up = spin down
!!
     integer, public, save :: isspn  = 1

!!
!! @var iswor
!!
!! control flag, define which algorithm will be used to do the measurement
!!
!! if iswor == 1:
!!     without worm algorithm, fast but unreliable
!!
!! if iswor == 2:
!!     with worm algorithm, slow but more reliable. note that only some
!!     selected physical observables can be measured by this algorithm
!!
!! note:
!!
!!     this feature has not been implemented so far
!!
     integer, public, save :: iswor  = 1

!!
!! @var isort
!!
!! control flag, define which basis will be used to do the measurement
!!
!! if isort == 1:
!!     using standard representation
!!
!! if isort == 2:
!!     using legendre orthogonal polynomial representation
!!
!! if isort == 3:
!!     using intermediate (singular value decomposition) representation
!!
     integer, public, save :: isort  = 1

!!
!! @var isobs
!!
!! control flag, it is used to tell the code which physical observables
!! should be measured. we just use the following rules to judge:
!!
!! rule 1:
!!     isobs is firstly converted to a binary representation. for example,
!!     10_10 is converted to 1010_2, 15_10 is converted to 1111_2, etc
!!
!! rule 2:
!!     then we examine the bits one by one. if it is 1, then we try to do
!!     the calculation. if it is 0, then we ignore the calculation. for
!!     example, we just use the second bit (from right side to left side)
!!     to represent the calculation of kinetic energy fluctuation. so, if
!!     isobs is 10_10 (1010_2), we will try to compute the kinetic energy
!!     fluctuation. if isobs is 13_10 (1101_2), we will not calculate it
!!     since the second bit is 0
!!
!! the following are the definitions of bit representation:
!!
!! if p == 1:
!!     do nothing
!!
!! if p == 2:
!!     calculate kinetic energy fluctuation < k^2 > - < k >^2
!!
!! if p == 3:
!!     calculate fidelity susceptibility < k_l k_r > - < k_l > < k_r >
!!
!! if p == 4:
!!     calculate powers of local magnetization < S^n_z >
!!
!! if p >= 5:
!!     reserved
!!
!! example:
!!   ( 1 1 1 0 1 0 1 0 1)_2
!! p = 9 8 7 6 5 4 3 2 1
!!
     integer, public, save :: isobs  = 1

!!
!! @var issus
!!
!! control flag, it is used to tell the code whether we should calculate
!! the charge or spin susceptibility. we just use the following rules to
!! make a judgement:
!!
!! rule 1:
!!     issus is firstly converted to a binary representation. for example,
!!     10_10 is converted to 1010_2, 15_10 is converted to 1111_2, etc
!!
!! rule 2:
!!     then we examine the bits one by one. if it is 1, then we try to do
!!     the calculation. if it is 0, then we ignore the calculation. for
!!     example, we just use the second bit (from right side to left side)
!!     to represent the calculation of spin-spin correlation function. so,
!!     if issus is 10_10 (1010_2), we will try to calculate the spin-spin
!!     correlation function. if issus is 13_10 (1101_2), since the second
!!     bit is 0 we will not calculate it
!!
!! the following are the definitions of bit representation:
!!
!! if p == 1:
!!     do nothing
!!
!! if p == 2:
!!     calculate spin-spin correlation function (imaginary time)
!!
!! if p == 3:
!!     calculate charge-charge correlation function (imaginary time)
!!
!! if p == 4:
!!     calculate spin-spin correlation function (matsubara frequency)
!!
!! if p == 5:
!!     calculate charge-charge correlation function (matsubara frequency)
!!
!! if p >= 6:
!!     reserved
!!
!! example:
!!   ( 1 1 1 0 1 0 1 0 1)_2
!! p = 9 8 7 6 5 4 3 2 1
!!
     integer, public, save :: issus  = 1

!!
!! @var isvrt
!!
!! control flag, it is used to tell the code whether we should measure
!! the two-particle green's functions. we just use the following rules
!! to judge:
!!
!! rule 1:
!!     isvrt is firstly converted to a binary representation. for example,
!!     10_10 is converted to 1010_2, 15_10 is converted to 1111_2, etc
!!
!! rule 2:
!!     then we examine the bits one by one. if it is 1, then we try to do
!!     the calculation. if it is 0, then we ignore the calculation. for
!!     example, we just use the second bit (from right side to left side)
!!     to represent the calculation of two-particle green's function for
!!     the particle-hole channel in AABB block structure. so, if isvrt is
!!     10_10 (1010_2), we will try to compute the specified two-particle
!!     green's function. if isvrt is 13_10 (1101_2), we will not calculate
!!     it since the second bit is 0
!!
!! the following are the definitions of bit representation:
!!
!! if p == 1:
!!     do nothing
!!
!! if p == 2:
!!     calculate two-particle green's function
!!     block: AABB
!!     channel: particle-hole
!!
!! if p == 3:
!!     calculate two-particle green's function
!!     block: ABBA
!!     channel: particle-hole
!!
!! if p == 4:
!!     calculate two-particle green's function
!!     block: AABB
!!     channel: particle-particle
!!
!! if p == 5:
!!     calculate two-particle green's function
!!     block: ABBA
!!     channel: particle-particle
!!
!! if p >= 6:
!!     reserved
!!
!! note 1:
!!
!!     p = 2 and p = 3 could not be enabled at the same time
!!
!! note 2:
!!
!!     p = 4 and p = 5 could not be enabled at the same time
!!
!! example:
!!   ( 1 1 1 0 1 0 1 0 1)_2
!! p = 9 8 7 6 5 4 3 2 1
!!
     integer, public, save :: isvrt  = 1

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
!! number of spin projections
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
!! number of atomic eigenstates (= 2**norbs)
!!
     integer, public, save :: ncfgs  = 4

!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

!!
!! @var niter
!!
!! number of self-consistent iterations for dynamical mean field theory
!! simulation combined with continuous time quantum Monte Carlo quantum
!! impurity solver
!!
     integer, public, save :: niter  = 20

!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

!!
!! @var lemax
!!
!! maximum expansion order for legendre orthogonal polynomial
!!
     integer, public, save :: lemax  = 32

!!
!! @var legrd
!!
!! number of mesh points for legendre orthogonal polynomial in [-1,1] range
!!
     integer, public, save :: legrd  = 20001

!!
!! @var svmax
!!
!! maximum expansion order for svd orthogonal polynomial
!!
     integer, public, save :: svmax  = 32

!!
!! @var svgrd
!!
!! number of mesh points for svd orthogonal polynomial in [-1,1] range
!!
     integer, public, save :: svgrd  = 2001

!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

!!
!! @var mkink
!!
!! maximum perturbation expansion order
!!
     integer, public, save :: mkink  = 1024

!!
!! @var mfreq
!!
!! maximum number of matsubara frequency points
!!
     integer, public, save :: mfreq  = 8193

!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

!!
!! @var nffrq
!!
!! number of fermionic frequencies for the two-particle green's function
!!
     integer, public, save :: nffrq  = 32

!!
!! @var nbfrq
!!
!! number of bosonic frequncies for the two-particle green's function
!!
     integer, public, save :: nbfrq  = 8

!!
!! @var nfreq
!!
!! number of matsubara frequencies sampled by continuous time quantum
!! Monte Carlo quantum impurity solver directly. the values for the other
!! points should be evaluated by using the other tricks
!!
     integer, public, save :: nfreq  = 128

!!
!! @var ntime
!!
!! number of imaginary time slices sampled by continuous time quantum
!! Monte Carlo quantum impurity solver
!!
     integer, public, save :: ntime  = 1024

!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

!!
!! @var nflip
!!
!! flip period for spin up and spin down states. some care must be taken
!! to prevent the system from being trapped in a state which breaks a
!! symmetry of local hamiltonian when it should not be. to avoid this
!! unphysical trapping, we introduce "flip" moves, which exchange the
!! operators corresponding, for example, to up and down spins in a given
!! orbital. note that nflip could be negative. if nflip /= 0, then the
!! absolute value of nflip is the flip period
!!
!! in this code, nowadays the following flip schemes are supported:
!!
!! if nflip == 0:
!!     means infinite long period to do flip. do not flip the spins
!!
!! if nflip >  0:
!!     flip intra-orbital spins one by one (90%) and globally (10%)
!!
!! if nflip <  0:
!!     flip intra-orbital spins globally (90%) and one by one (10%)
!!
     integer, public, save :: nflip  = 20000

!!
!! @var ntherm
!!
!! number of thermalization steps
!!
     integer, public, save :: ntherm = 200000

!!
!! @var nsweep
!!
!! number of Monte Carlo sampling steps
!!
     integer, public, save :: nsweep = 20000000

!!
!! @var nwrite
!!
!! output period for quantum impurity solver
!!
     integer, public, save :: nwrite = 2000000

!!
!! @var nclean
!!
!! clean update period for quantum impurity solver
!!
     integer, public, save :: nclean = 100000

!!
!! @var nmonte
!!
!! how often to sample the physical observables. it would be adjusted in
!! the ctqmc_try_warming() subroutine automatically via rough estimation
!! of the integrated autocorrelation time for the total occupation number
!!
     integer, public, save :: nmonte = 10

!!
!! @var ncarlo
!!
!! how often to sample the physical observables
!!
!! note:
!!
!!     it is reserved for the future
!!
     integer, public, save :: ncarlo = 10

!!========================================================================
!!>>> real variables                                                   <<<
!!========================================================================

!!
!! @var Uc
!!
!! intra-orbital Coulomb interaction
!!
     real(dp), public, save :: Uc    = 4.00_dp

!!
!! @var Jz
!!
!! Hund's exchange interaction in z axis
!!
     real(dp), public, save :: Jz    = 0.00_dp

!!
!! @var lc
!!
!! strength of dynamical screening effect. its meaning depends on isscr
!!
!! if isscr == 1:
!!     lc is ignored
!!
!! if isscr == 2:
!!     lc just means the model parameter \lambda
!!
!! if isscr == 3:
!!     lc just means the model parameter \alpha
!!
!! if isscr == 4:
!!     lc is ignored
!!
     real(dp), public, save :: lc    = 1.00_dp

!!
!! @var wc
!!
!! screening frequency. its meaning depends on isscr
!!
!! if isscr == 1:
!!     wc is ignored
!!
!! if isscr == 2:
!!     wc just means the model parameter \omega^{'}
!!
!! if isscr == 3:
!!     wc just means the model parameter \omega_{c}
!!
!! if isscr == 4:
!!     wc is ignored
!!
     real(dp), public, save :: wc    = 1.00_dp

!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

!!
!! @var mune
!!
!! chemical potential or fermi level
!!
     real(dp), public, save :: mune  = 2.00_dp

!!
!! @var beta
!!
!! inversion of temperature
!!
     real(dp), public, save :: beta  = 8.00_dp

!!
!! @var part
!!
!! hopping parameter t for Hubbard model
!!
     real(dp), public, save :: part  = 0.50_dp

!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

!!
!! @var alpha
!!
!! mixing factor for dynamical mean field theory self-consistent engine
!!
     real(dp), public, save :: alpha = 0.70_dp

!!========================================================================
!!>>> MPI related common variables                                     <<<
!!========================================================================

!!
!! @var nprocs
!!
!! number of processors: default value 1
!!
     integer, public, save :: nprocs = 1

!!
!! @var myid
!!
!! the id of current process: default value 0
!!
     integer, public, save :: myid   = 0

!!
!! @var master
!!
!! denote as the controller process: default value 0
!!
     integer, public, save :: master = 0

!!
!! @var cid
!!
!! the id of current process in cartesian topology (cid == myid)
!!
     integer, public, save :: cid    = 0

!!
!! @var cx
!!
!! the x coordinates of current process in cartesian topology
!!
     integer, public, save :: cx     = 0

!!
!! @var cy
!!
!! the y coordinates of current process in cartesian topology
!!
     integer, public, save :: cy     = 0

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
     character(len=20), public, parameter :: V_FULL = 'v0.8.6 @ 2025.05.30D'

!!
!! @var V_CURR
!!
!! version string, only version number
!!
     character(len=06), public, parameter :: V_CURR = 'v0.8.6'

!!
!! @var V_DATE
!!
!! version string, only date info.
!!
     character(len=11), public, parameter :: V_DATE = '2025.05.30'

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
