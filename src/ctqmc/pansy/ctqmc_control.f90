!!!-----------------------------------------------------------------------
!!! project : pansy
!!! program : control    module
!!! source  : ctqmc_control.f90
!!! type    : module
!!! author  : li huang (email:lihuang.dmft@gmail.com)
!!!           yilin wang (email:qhwyl2006@126.com)
!!! history : 09/15/2009 by li huang
!!!           02/23/2010 by li huang
!!!           11/11/2014 by yilin wang
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
!!>>> integer variables                                                <<<
!!========================================================================

! control flag: running mode
! if isscf == 1, one-shot non-self-consistent scheme, used in local density
! approximation plus dynamical mean field theory case
! if isscf == 2, self-consistent scheme, used in normal model hamiltonian
! plus dynamical mean field theory case
     integer, public, save :: isscf  = 2

! control flag: symmetry of bands
! if issun == 1, the bands are not symmetrized
! if issun == 2, the bands are symmetrized according to symmetry matrix
     integer, public, save :: issun  = 2

! control flag: symmetry of spin orientation
! if isspn == 1, enforce spin up = spin down
! if isspn == 2, let spin up and spin down states evolve independently
     integer, public, save :: isspn  = 1

! control flag: impurity green's function binning mode
! if isbin == 1, without binning mode
! if isbin == 2, with binning mode
     integer, public, save :: isbin  = 2

!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

! number of correlated bands
     integer, public, save :: nband  = 1

! number of spin projection
     integer, public, save :: nspin  = 2

! number of correlated orbitals (= nband * nspin)
     integer, public, save :: norbs  = 2

! number of atomic states (= 2**norbs)
     integer, public, save :: ncfgs  = 4

! maximum number of continuous time quantum Monte Carlo quantum impurity
! solver plus dynamical mean field theory self-consistent iterations
     integer, public, save :: niter  = 20

!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

! maximum perturbation expansion order
     integer, public, save :: mkink  = 1024

! maximum number of matsubara frequency point
     integer, public, save :: mfreq  = 8193

!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

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
     integer, public, save :: nmonte = 10

! how often to sampling the gtau and prob
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

! interorbital Coulomb interaction, Uv = Uc - 2 * Jz for t2g system
     real(dp), public, save :: Uv    = 4.00_dp

! Hund's exchange interaction in z axis (Jz = Js = Jp = J)
     real(dp), public, save :: Jz    = 0.00_dp

! spin-flip term
     real(dp), public, save :: Js    = 0.00_dp

! pair-hopping term
     real(dp), public, save :: Jp    = 0.00_dp

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
