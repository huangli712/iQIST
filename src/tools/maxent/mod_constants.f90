!---------------------------------------------------------------
! project : maxent
! program : constants  module
! source  : mod_constants.f90
! type    : module
! author  : yilin wang (email: qhwyl2006@126.com)
! history : 05/29/2013 by yilin wang
! purpose : define some constants for maxent
! input   :
! output  :
! status  : unstable
! comment :
!---------------------------------------------------------------

  module constants
      implicit none
!===============================================================
! integer constants
!===============================================================
! single precision of float numbers
      integer, public, parameter :: sp = kind(1.0)

! double precision of float numbers
      integer, public, parameter :: dp = kind(1.0d0)

! file unit for standard output
      integer, public, parameter :: mystd = 6

! file unit for ordinary files
      integer, public, parameter :: mytmp = 100

! file unit for log file
      integer, public, parameter :: mylog = 101

! master node in mpi case
      integer, public, parameter :: master = 0

!===============================================================
! real constants
!===============================================================
! zero in double precision
      real(dp), public, parameter :: zero = 0.0_dp

! one in double precision
      real(dp), public, parameter :: one  = 1.0_dp

! two in double precision
      real(dp), public, parameter :: two  = 2.0_dp

! three in double precision
      real(dp), public, parameter :: three  = 3.0_dp

! half in double precision
      real(dp), public, parameter :: half = 0.5_dp

! one third in double precision
      real(dp), public, parameter :: one_third = one / three

! one quarter in double precision
      real(dp), public, parameter :: one_quart = 0.25_dp

! pi in double precision
      real(dp), public, parameter :: pi = 3.14159265358979323846_dp

! two*pi in double precision
      real(dp), public, parameter :: twopi = two * pi

! some tolerance in double precision form
      real(dp), public, parameter :: eps1  = 1E-1 
      real(dp), public, parameter :: eps2  = 1E-2 
      real(dp), public, parameter :: eps3  = 1E-3 
      real(dp), public, parameter :: eps4  = 1E-4 
      real(dp), public, parameter :: eps5  = 1E-5 
      real(dp), public, parameter :: eps6  = 1E-6 
      real(dp), public, parameter :: eps7  = 1E-7 
      real(dp), public, parameter :: eps8  = 1E-8 
      real(dp), public, parameter :: eps9  = 1E-9 
      real(dp), public, parameter :: eps10 = 1E-10 
      real(dp), public, parameter :: eps11 = 1E-11 
      real(dp), public, parameter :: eps12 = 1E-12

!===============================================================
! complex constants
!===============================================================
! zero in complex form
      complex(dp), public, parameter ::  czero = dcmplx(zero, zero) 

! real one in complex form
      complex(dp), public, parameter :: cone = dcmplx(one, zero)

! imaginary one 
      complex(dp), public, parameter :: ci = dcmplx(zero, one) 

  end module constants

