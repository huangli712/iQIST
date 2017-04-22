!!!-----------------------------------------------------------------------
!!! project : narcissus
!!! program : control    module
!!! source  : ctqmc_control.f90
!!! type    : module
!!! author  : li huang (email:lihuang.dmft@gmail.com)
!!! history : 09/15/2009 by li huang (created)
!!!           04/22/2017 by li huang (last modified)
!!! purpose : define global control parameters for hybridization expansion
!!!           version continuous time quantum Monte Carlo (CTQMC) quantum
!!!           impurity solver and dynamical mean field theory (DMFT) self-
!!!           consistent engine
!!! status  : unstable
!!! comment :
!!!-----------------------------------------------------------------------

  module control
     use constants, only : dp

     implicit none

!!========================================================================
!!>>> character variables                                              <<<
!!========================================================================

!!
!! @var cname
!!
!! the code name of the current quantum impurity solver
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
!!     to solve the quantum impurity model
!!
!! if isscf == 2:
!!     self-consistent scheme, used in the dynamical mean field theory
!!     case. the code implements a dynamical mean field self-consistent
!!     loop for solving the hubbard model in the bethe lattice
!!
     integer, public, save :: isscf  = 2

!!
!! @var isscr
!!
!! control flag, define whether the Coulomb interaction U is dynamic
!!
!! if isscr == 1:
!!     static interaction
!!
!! if isscr == 2:
!!     dynamic interaction, for palsmon pole model
!!
!! if isscr == 3
!!     dynamic interaction, for ohmic model
!!
!! if isscr ==99
!!     dynamic interaction, for realistic materials
!!
     integer, public, save :: isscr  = 1

!!
!! @var isbnd
!!
!! control flag, define symmetry of the model (band part)
!!
!! if isbnd == 1:
!!     the bands are not symmetrized
!!
!! if isbnd == 2:
!!     the bands are symmetrized according to symmetry matrix
!!
     integer, public, save :: isbnd  = 2

!!
!! @var isspn
!!
!! control flag, define symmetry of the model (spin part)
!!
!! if isspn == 1:
!!     enforce spin up = spin down
!!
!! if isspn == 2:
!!     let spin up and spin down states evolve independently
!!
     integer, public, save :: isspn  = 1

!!
!! @var iswor
!!
!! control flag, define which algorithm should be used to do measurement
!!
!! if iswor == 1:
!!     without worm algorithm, fast but unreliable
!!
!! if iswor == 2
!!     with worm algorithm, slow but reliable
!!
     integer, public, save :: iswor  = 1

!!
!! @var isbin
!!
!! control flag, define how to accumulate the imaginary-time impurity
!! green's function G(\tau) data
!!
!! if isbin == 1:
!!     without data binning mode
!!
!! if isbin == 2:
!!     with data binning mode
!!
     integer, public, save :: isbin  = 2

!!
!! @var isort
!!
!! control flag, define which basis should be used to do measurement
!!
!! if isort == 1:
!!     use standard representation
!!
!! if isort == 2
!!     use legendre orthogonal polynomial representation
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
!!     calculate fidelity susceptibility
!!
!! if p == 4:
!!     calculate < S^n_z >, powers of local magnetization
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
!!     calculate spin-spin correlation function (time space)
!!
!! if p == 3:
!!     calculate orbital-orbital correlation function (time space)
!!
!! if p == 4:
!!     calculate spin-spin correlation function (frequency space)
!!
!! if p == 5:
!!     calculate orbital-orbital correlation function (frequency space)
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
!!     to represent the calculation of two-particle green's function. so,
!!     if isvrt is 10_10 (1010_2), we will try to compute the two-particle
!!     green's function. if isvrt is 13_10 (1101_2), we will not calculate
!!     it since the second bit is 0
!!
!! the following are the definitions of bit representation:
!!
!! if p == 1:
!!     do nothing
!!
!! if p == 2:
!!     calculate two-particle green's function and vertex function
!!
!! if p == 3:
!!     calculate particle-particle pair susceptibility
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
!! number of spin projection
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
!! number of atomic states (= 2**norbs)
!!
     integer, public, save :: ncfgs  = 4

!!
!! @var niter
!!
!! maximum number of self-consistent iterations for the continuous time
!! quantum Monte Carlo quantum impurity solver plus dynamical mean field
!! theory calculation
!!
     integer, public, save :: niter  = 20

!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

! maximum order for legendre polynomial
     integer, public, save :: lemax  = 32

! number of mesh points for legendre polynomial in [-1,1] range
     integer, public, save :: legrd  = 20001

!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

! maximum perturbation expansion order
     integer, public, save :: mkink  = 1024

! maximum number of matsubara frequency point
     integer, public, save :: mfreq  = 8193

!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

! number of matsubara frequency for the two-particle green's function
     integer, public, save :: nffrq  = 32

! number of bosonic frequncy for the two-particle green's function
     integer, public, save :: nbfrq  = 8

! number of matsubara frequency sampling by continuous time quantum Monte
! Carlo quantum impurity solver
!
! note: the rest (mfreq - nfreq + 1 points) values are evaluated by using
! Hubbard-I approximation
     integer, public, save :: nfreq  = 128

! number of imaginary time slice sampling by continuous time quantum Monte
! Carlo quantum impurity solver
     integer, public, save :: ntime  = 1024

! flip period for spin up and spin down states
!
! note: care must be taken to prevent the system from being trapped in a
! state which breaks a symmetry of local hamiltonian when it should not
! be. to avoid unphysical trapping, we introduce "flip" moves, which
! exchange the operators corresponding, for example, to up and down spins
! in a given orbital.
!
! note: in this code, nowadays the following flip schemes are supported
!     if cflip = 1, flip inter-orbital spins randomly;
!     if cflip = 2, flip intra-orbital spins one by one;
!     if cflip = 3, flip intra-orbital spins globally.
! here cflip is an internal variable.
!
! note: we use the sign of nflip to control flip schemes
!     if nflip = 0, means infinite long period to do flip
!     if nflip > 0, combine cflip = 2 (80%) and cflip = 3 (20%)
!     if nflip < 0, combine cflip = 1 (80%) and cflip = 3 (20%)
!
! note: if nflip /= 0, the absolute value of nflip is the flip period
!
! note: when cflip = 1, the symmetry of all orbitals must be taken into
! consideration, otherwise the code may be trapped by a deadlock.
     integer, public, save :: nflip  = 20000

! maximum number of thermalization steps
     integer, public, save :: ntherm = 200000

! maximum number of quantum Monte Carlo sampling steps
     integer, public, save :: nsweep = 20000000

! output period for quantum impurity solver
     integer, public, save :: nwrite = 2000000

! clean update period for quantum impurity solver
     integer, public, save :: nclean = 100000

! how often to sampling the physical observables
     integer, public, save :: nmonte = 10

! how often to sampling the physical observables
     integer, public, save :: ncarlo = 10

!!========================================================================
!!>>> real variables                                                   <<<
!!========================================================================

! average Coulomb interaction
     real(dp), public, save :: U     = 4.00_dp

! intraorbital Coulomb interaction
     real(dp), public, save :: Uc    = 4.00_dp

! interorbital Coulomb interaction, Uv = Uc - 2 * Jz for t2g system
     real(dp), public, save :: Uv    = 4.00_dp

! Hund's exchange interaction in z axis (Jz = Js = Jp = J)
     real(dp), public, save :: Jz    = 0.00_dp

! spin-flip term
     real(dp), public, save :: Js    = 0.00_dp

! pair-hopping term
     real(dp), public, save :: Jp    = 0.00_dp

! strength of dynamical screening effect
! when isscr == 01, lc is ignored
! when isscr == 02, lc just means the control parameter \lambda 
! when isscr == 03, lc just means the control parameter \alpha
! when isscr == 99, lc just means the shift for interaction matrix
     real(dp), public, save :: lc    = 1.00_dp

! screening frequency
! when isscr == 01, wc is ignored. 
! when isscr == 02, wc just means the control parameter \omega^{'}
! when isscr == 03, wc just means the control parameter \omega_{c}
! when isscr == 99, wc just means the shift for chemical potential
     real(dp), public, save :: wc    = 1.00_dp

!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

! chemical potential or fermi level
!
! note: it should/can be replaced with eimp
     real(dp), public, save :: mune  = 2.00_dp

! inversion of temperature
     real(dp), public, save :: beta  = 8.00_dp

! coupling parameter t for Hubbard model
     real(dp), public, save :: part  = 0.50_dp

! mixing parameter for dynamical mean field theory self-consistent engine
     real(dp), public, save :: alpha = 0.70_dp

!!========================================================================
!!>>> MPI related common variables                                     <<<
!!========================================================================

! number of processors: default value 1
     integer, public, save :: nprocs = 1

! the id of current process: default value 0
     integer, public, save :: myid   = 0

! denote as the controller process: default value 0
     integer, public, save :: master = 0

! the id of current process in cartesian topology (cid == myid)
     integer, public, save :: cid    = 0

! the x coordinates of current process in cartesian topology
     integer, public, save :: cx     = 0

! the y coordinates of current process in cartesian topology
     integer, public, save :: cy     = 0

  end module control
