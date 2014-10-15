!!!---------------------------------------------------------------
!!! project : maxent
!!! program : control  module
!!! source  : maxent_control.f90
!!! type    : module
!!! author  : yilin wang (email: qhwyl2006@126.com)
!!! history : 05/29/2013 by yilin wang
!!!           10/13/2014 by yilin wang
!!! purpose : define some control variables for maxent
!!! status  : unstable
!!! comment :
!!!---------------------------------------------------------------

  module control
     use constants, only : dp, one

     implicit none

!!==============================================================
!!>>> integer control variables
!!==============================================================

! default model type
! 0 for flat model, 1 for Gauss model, 2 for reading from file
     integer, public, save :: imode = 0

! whether to diagonalize the covariance matrix type
! 0 for no, 1 for yes
     integer, public, save :: icov = 0

! number of Monte Carlo data bins
     integer, public, save :: nbins = 1000

! number of time slice in each Monte Carlo data bins
     integer, public, save :: ntime = 129

! number of frequency mesh points of half axis
     integer, public, save :: nwhf = 1000

! number of frequency mesh points of the total axis
     integer, public, save :: nw = 2001

! dimension of alpha mesh
     integer, public, save :: nalpha = 100

! slice for histogram of data bins
     integer, public, save :: slice = 100

!!==============================================================
!!>>> real control variables
!!==============================================================

! inversion of temperatue
     real(dp), public, save :: beta = 10.0_dp

! step size for frequency mesh, for fixed step size, the 
! range of frequency is [-nw/2 * step : nw/2 * step]
     real(dp), public, save :: step = 0.02_dp

! if imode == 1, Gauss default model is used, then sigma is
! the standard deviation of Gaussian distribution
     real(dp), public, save :: sigma = 5.0_dp

! max value of parameter \alpha
     real(dp), public, save :: max_alpha = 100.0_dp

! low value of parameter \alpha
     real(dp), public, save :: min_alpha = one

  end module control
