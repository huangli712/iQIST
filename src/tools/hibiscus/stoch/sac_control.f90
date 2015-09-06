!!!-----------------------------------------------------------------------
!!! project : hibiscus/stoch
!!! program : control    module
!!! source  : sac_control.f90
!!! type    : module
!!! author  : li huang (email:lihuang.dmft@gmail.com)
!!! history : 01/08/2011 by li huang (created)
!!!           08/17/2015 by li huang (last modified)
!!! purpose : define global control parameters for stochastic analytic
!!!           continuation code
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
     integer, public, save :: ntime = 1024

! number of frequency points on half axis, energy range can be expressed by
! [ -wstep * nwmax, wstep * nwmax ]
     integer, public, save :: nwmax = 128

! number of slices of x in [0,1]
     integer, public, save :: ngrid = 10001

! number of configurations, dimension for r_{\gamma} and a_{\gamma}
     integer, public, save :: ngamm = 1024

! number of alpha parameters used in parallel tempering
! note: it must be an even number, since we need to exchange configurations
! between different alpha channel
     integer, public, save :: nalph = 10

! maximum number of thermalization steps
     integer, public, save :: nwarm = 4000

! maximum number of quantum Monte Carlo sampling steps
     integer, public, save :: nstep = 4000000

! output period for stochastic analytic continuation code
     integer, public, save :: ndump = 40000

! measurement scheme
! if ltype == 1, normal measurement
! if ltype == 2, using legendre polynomial representation
     integer, public, save :: ltype = 1

! maximum order for legendre polynomial
     integer, public, save :: lemax = 64

! number of mesh points for legendre polynomial in [-1,1] range
     integer, public, save :: legrd = 20001

!!========================================================================
!!>>> real variables                                                   <<<
!!========================================================================

! initial alpha parameter
     real(dp), public, save :: ainit = 1.00_dp

! \alpha_(p+1) / \alpha_p = R, used to build alpha parameter list
     real(dp), public, save :: ratio = 2.00_dp

! \beta, inversion of real temperature
     real(dp), public, save :: beta  = 10.0_dp

! lorentz broadening parameter \eta, used to represent delta function
     real(dp), public, save :: eta1  = 0.02_dp
     real(dp), public, save :: eta2  = 4E-4_dp

! gauss broadening parameter, used to build the default model
     real(dp), public, save :: sigma = 1.00_dp

! frequency step, used to build the frequency mesh
     real(dp), public, save :: wstep = 0.05_dp

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
