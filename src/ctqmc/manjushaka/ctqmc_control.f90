!!!-----------------------------------------------------------------------
!!! project : manjushaka
!!! program : control    module
!!! source  : ctqmc_control.f90
!!! type    : module
!!! author  : li huang (email:lihuang.dmft@gmail.com)
!!!           yilin wang (email:qhwyl2006@126.com)
!!! history : 09/15/2009 by li huang (created)
!!!           05/10/2017 by li huang (last modified)
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
!! code name of the current quantum impurity solver
!!
     character(len = 10), public, save :: cname = 'MANJUSHAKA'

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
!!     to solve the quantum impurity model once
!!
!! if isscf == 2:
!!     self-consistent scheme, used in the dynamical mean field theory
!!     case. the code implements a typical dynamical mean field theory
!!     self-consistent loop for solving the Hubbard model in the bethe
!!     lattice (semicircular density of state)
!!
     integer, public, save :: isscf  = 1

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
!! @var isbin
!!
!! control flag, define how to accumulate data for imaginary time
!! impurity green's function G(\tau)
!!
!! if isbin == 1:
!!     without data binning mode
!!
!! if isbin == 2:
!!     with data binning mode
!!
     integer, public, save :: isbin  = 1

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
!!     selected physical observables support this algorithm
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
!!     calculate charge-charge correlation function (time space)
!!
!! if p == 4:
!!     calculate spin-spin correlation function (frequency space)
!!
!! if p == 5:
!!     calculate charge-charge correlation function (frequency space)
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
!!     calculate two-particle green's function
!!
!! if p == 3:
!!     calculate particle-particle pairing susceptibility
!!
!! example:
!!   ( 1 1 1 0 1 0 1 0 1)_2
!! p = 9 8 7 6 5 4 3 2 1
!!
     integer, public, save :: isvrt  = 1

! control flag: the efficient algorithm for calculate the trace
! if ifast == 1, use divide-and-conquer algorithm (see npart as well)
! if ifast == 2, use classic time evolution algorithm, not implemented
! if ifast == 3, use skip listing algorithm, not implemented
     integer, public, save :: ifast  = 1

! control flag: the mode how to truncate the Hilbert space
! if itrun == 1, don't truncate it
! if itrun == 2, truncate high energy states
     integer, public, save :: itrun  = 1

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
!! number of self-consistent iterations for the continuous time quantum
!! Monte Carlo quantum impurity solver plus dynamical mean field theory
!! simulation
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

! number of parts that the imaginary time axis is split
!
! note: all operators in the imaginary time axis are grouped into npart
! parts according to their time values, in each Monte Carlo steps, only
! those changed parts are carefully dealt with, not all the parts.
!
! note: 2\sqrt{3 <k> nband} ~ 4\sqrt{3 <k> nband} may be the optimal value
! for npart to achieve maximum performance.
     integer, public, save :: npart  = 4

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

! how often to sampling the gmat and paux (nmat and nnmat)
!
! note: the measure periods for schi, sschi, ochi, oochi, g2_re, g2_im,
! h2_re, h2_im, ps_re, and ps_im are also controlled by nmonte parameter.
     integer, public, save :: nmonte = 10

! how often to sampling the gtau and prob
!
! note: the measure period for ftau is also controlled by ncarlo parameter.
     integer, public, save :: ncarlo = 10

!!========================================================================
!!>>> real variables                                                   <<<
!!========================================================================

! note: U, Uc, Uv, Jz, Js, and Jp are not used by this quantum impurity
! solver actually. we keep them here is just for reference.

! average Coulomb interaction
     real(dp), public, save :: U     = 4.00_dp

! intraorbital Coulomb interaction
     real(dp), public, save :: Uc    = 4.00_dp

! Hund's exchange interaction in z axis (Jz = Js = Jp = J)
     real(dp), public, save :: Jz    = 0.00_dp

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
