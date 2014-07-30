!-------------------------------------------------------------------------
! project : manjushaka
! program : control    module
! source  : ctqmc_control.f90
! type    : module
! author  : li huang (email:huangli712@yahoo.com.cn)
! history : 09/15/2009 by li huang
!           09/20/2009 by li huang
!           11/01/2009 by li huang
!           12/01/2009 by li huang
!           02/23/2010 by li huang
! purpose : define global control parameters for hybridization expansion
!           version continuous time quantum Monte Carlo (CTQMC) quantum
!           impurity solver and dynamical mean field theory (DMFT) self-
!           consistent engine
! input   :
! output  :
! status  : unstable
! comment :
!-------------------------------------------------------------------------

  module control
     use constants, only : dp

     implicit none

!=========================================================================
!>>> integer variables                                                 <<<
!=========================================================================

! control flag: running mode
! if isscf == 1, one-shot non-self-consistent scheme, used in local density
! approximation plus dynamical mean field theory case
! if isscf == 2, self-consistent scheme, used in normal model hamiltonian
! plus dynamical mean field theory case
     integer, public, save :: isscf  = 1

! control flag: symmetry of bands
! if issun == 1, the bands are not symmetrized
! if issun == 2, the bands are symmetrized according to symmetry matrix
     integer, public, save :: issun  = 1

! control flag: symmetry of spin orientation
! if isspn == 1, enforce spin up = spin down
! if isspn == 2, let spin up and spin down states evolve independently
     integer, public, save :: isspn  = 1

! control flag: impurity green's function binning mode
! if isbin == 1, without binning mode
! if isbin == 2, with binning mode
     integer, public, save :: isbin  = 1

! control flag: apply orthogonal polynomial representation to perform measurement
! if isort == 1, use normal representation to measure G(\tau)
! if isort == 2, use legendre polynomial to measure G(\tau)
! if isort == 3, use chebyshev polynomial (the second kind) to measure G(\tau)
! if isort == 4, use normal representation to measure G(\tau) and F(\tau)
! if isort == 5, use legendre polynomial to measure G(\tau) and F(\tau)
! if isort == 6, use chebyshev polynomial (the second kind) to measure G(\tau) and F(\tau)
! note: if isort \in [1,3], we use ctqmc_make_hub1() to calculate the self
! energy function, or else we use ctqmc_make_hub2().
! note: as for the kernel polynomial representation, the default dirichlet
! kernel is applied automatically. if you want to choose the other kernel,
! please check the ctqmc_make_gtau() subroutine.
! note: isort == 4, 5, and 6 are not implemeted so far.
     integer, public, save :: isort  = 1

! control flag: whether we measure the high order correlation function
! if isvrt == 1, do nothing
! if isvrt == 2, calculate spin-spin correlation function
! if isvrt == 3, calculate orbital-orbital correlation function
! if isvrt == 4, calculate both two-particle green's function and vertex function
! if isvrt == 5, calculate both two-particle green's function and vertex function
! note: when isvrt == 4 and isvrt == 5, both the two-particle green's and
! vertex functions are computed by using two different algorithms.
! note: when isvrt == 4, both the solver.twop.dat and solver.vrtx.dat files
! are written, but the data contained in solver.vrtx.dat file is not right.
! note: when isvrt == 5, both the solver.twop.dat and solver.vrtx.dat files
! are written, the data contained in solver.twop.dat are less accurate than
! those in solver.vrtx.dat, more specifically, the irreducible part of two-
! particle green's function and vertex function.
! note: isvrt == 2, 3, and 5 are not implemented so far.
     integer, public, save :: isvrt  = 1

! control flag: which trace algorithm to be used
! if iskip == 0, use npart algorithm
! if iskip == 1, use skip lists algorithm
     integer, public, save :: iskip

!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

! number of correlated bands
     integer, public, save :: nband  = 1

! number of spin projection
     integer, public, save :: nspin  = 2

! number of correlated orbitals (= nband * nspin)
     integer, public, save :: norbs  = 2

! number of atomic states (= 2**norbs)
     integer, public, save :: ncfgs  = 4

! maximum allowed number of non-zero elements in F-matrix
     integer, public, save :: nzero  = 128

! maximum number of continuous time quantum Monte Carlo quantum impurity
! solver plus dynamical mean field theory self-consistent iterations
     integer, public, save :: niter  = 20

!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

! maximum order for legendre polynomial
     integer, public, save :: lemax  = 32

! number of mesh points for legendre polynomial in [-1,1] range
     integer, public, save :: legrd  = 20001

! maximum order for chebyshev polynomial
     integer, public, save :: chmax  = 32

! number of mesh points for chebyshev polynomial in [-1,1] range
     integer, public, save :: chgrd  = 20001

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
! note: the rest (mfreq - nfreq + 1 points) values are evaluated by using
! Hubbard-I approximation
     integer, public, save :: nfreq  = 128

! number of imaginary time slice sampling by continuous time quantum Monte
! Carlo quantum impurity solver
     integer, public, save :: ntime  = 1024

! number of parts that the imaginary time axis is split
! note: all operators in the imaginary time axis are grouped into npart
! parts according to their time values, in each Monte Carlo steps, only
! those changed parts are carefully dealt with, not all the parts.
! note: 2\sqrt{3 <k> nband} ~ 4\sqrt{3 <k> nband} may be the optimal value
! for npart to achieve maximum performance
     integer, public, save :: npart  = 16

! maximum level for skip lists
     integer, public, save :: mlevl  = 8

! flip period for spin up and spin down states
! note: care must be taken to prevent the system from being trapped in a
! state which breaks a symmetry of local hamiltonian when it should not
! be. to avoid unphysical trapping, we introduce "flip" moves, which
! exchange the operators corresponding, for example, to up and down spins
! in a given orbital.
! note: in this code, nowadays the following flip schemes are supported
!     if cflip = 1, flip inter-orbital spins randomly;
!     if cflip = 2, flip intra-orbital spins one by one;
!     if cflip = 3, flip intra-orbital spins globally.
! note: we use the sign of nflip to control flip schemes
!     if nflip = 0, means infinite long period to do flip
!     if nflip > 0, combine cflip = 2 (80%) and cflip = 3 (20%)
!     if nflip < 0, combine cflip = 1 (80%) and cflip = 3 (20%)
! note: if nflip /= 0, the absolute value of nflip is the flip period
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
! note: the measure periods for schi, sschi, ochi, oochi, g2_re, g2_im,
! h2_re, and h2_im are also controlled by nmonte parameter.
     integer, public, save :: nmonte = 10

! how often to sampling the gtau and prob
! note: the measure period for ftau is also controlled by ncarlo parameter.
     integer, public, save :: ncarlo = 10

!=========================================================================
!>>> real variables                                                    <<<
!=========================================================================

! note: U, Uc, Uv, Jz, Js, and Jp are not used by this quantum impurity
! solver actually. we keep them here is just for reference
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

! chemical potential or fermi level
! note: it should be replaced with eimp
     real(dp), public, save :: mune  = 2.00_dp

! inversion of temperature
     real(dp), public, save :: beta  = 8.00_dp

! coupling parameter t for Hubbard model
     real(dp), public, save :: part  = 0.50_dp

! mixing parameter for dynamical mean field theory self-consistent engine
     real(dp), public, save :: alpha = 0.70_dp

!=========================================================================
!>>> MPI related common variables                                      <<<
!=========================================================================

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
