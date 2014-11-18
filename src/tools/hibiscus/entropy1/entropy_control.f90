!!!-----------------------------------------------------------------------
!!! project : hibiscus/entropy1
!!! program : control    module
!!! source  : entropy_control.f90
!!! type    : module
!!! author  : li huang (email:huangli712@gmail.com)
!!! history : 01/08/2011 by li huang
!!!           01/26/2011 by li huang
!!!           11/17/2014 by li huang
!!! purpose : define global control parameters for classic maximum entropy
!!!           method code
!!! status  : unstable
!!! comment :
!!!-----------------------------------------------------------------------

  module control
     use constants, only : dp

     implicit none

!!========================================================================
!!>>> integer variables                                                <<<
!!========================================================================

! number of imaginary time slices sampling by continuous time or hirsh-fye
! quantum Monte Carlo quantum impurity solver
     integer, public, save :: ntime = 129

! number of frequency points on half axis, energy range can be expressed by
! [ -wstep * nwmax, wstep * nwmax ]
     integer, public, save :: nwmax = 200

! maximum number of cycles for classic maximum entropy method
     integer, public, save :: niter = 20

! number of smooth runs for classic maximum entropy method
     integer, public, save :: ntune = 20

! number of annealing steps per classic maximum entropy method cycle
     integer, public, save :: nstep = 4000

! number of bands
     integer, public, save :: nband = 1

! number of orbitals
     integer, public, save :: norbs = 2

! the way the default model function is build
! if ntype == 0, gaussian model
! if ntype == 1, flat model
     integer, public, save :: ntype = 1

!!========================================================================
!!>>> real variables                                                   <<<
!!========================================================================

! initial alpha parameter
     real(dp), public, save :: ainit = 1200._dp

! it is the deviation from the average green's function
     real(dp), public, save :: devia = 0.001_dp

! \beta, inversion of real temperature
     real(dp), public, save :: beta  = 10.00_dp

! gauss broadening parameter, used to build the default model
     real(dp), public, save :: sigma = 1.600_dp

! delta frequency, step of real frequency grid
     real(dp), public, save :: wstep = 0.025_dp

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
