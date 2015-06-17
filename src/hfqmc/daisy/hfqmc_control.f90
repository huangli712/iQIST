!!!-----------------------------------------------------------------------
!!! project : daisy
!!! program : control    module
!!! source  : hfqmc_control.f90
!!! type    : module
!!! author  : li huang (email:lihuang.dmft@gmail.com)
!!! history : 10/24/2008 by li huang
!!!           03/25/2010 by li huang
!!!           12/04/2014 by li huang
!!! purpose : define global control parameters for Hirsch-Fye quantum
!!!           Monte Carlo (HFQMC) quantum impurity solver and dynamical
!!!           mean field theory (DMFT) self-consistent engine
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

! maximum number of Hirsch-Fye quantum Monte Carlo quantum impurity solver
! plus dynamical mean field theory self-consistent iterations
     integer, public, save :: niter  = 20

!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

! maximum number of delayed update steps
! if mstep == 1, using traditional update algorithm
! if mstep >= 1, using delayed update algorithm to improve the efficiency
     integer, public, save :: mstep  = 16

! maximum number of matsubara frequency point
     integer, public, save :: mfreq  = 8193

!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

! number of auxiliary ising-like fields (= norbs * (norbs - 1) / 2)
     integer, public, save :: nsing  = 1

! number of imaginary time slice, 64, 128 or 256 in common
     integer, public, save :: ntime  = 64

! maximum number of thermalization steps
     integer, public, save :: ntherm = 100

! maximum number of quantum Monte Carlo sampling steps
     integer, public, save :: nsweep = 240000

! clean update period for quantum impurity solver
     integer, public, save :: nclean = 100

! how many QMC steps between two measurements
! note: in order to relieve the data correlation between two successive
! QMC measurements, we make an effective measurement after ncarlo accepted
! ising spin flips. noted by li huang, 2007/03/15
     integer, public, save :: ncarlo = 10

!!========================================================================
!!>>> real variables                                                   <<<
!!========================================================================

! intraorbital Coulomb interaction
     real(dp), public, save :: Uc    = 4.00_dp

! Hund's exchange interaction in z axis (Jz = Js = Jp = J)
     real(dp), public, save :: Jz    = 0.00_dp

!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

! chemical potential or fermi level
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
