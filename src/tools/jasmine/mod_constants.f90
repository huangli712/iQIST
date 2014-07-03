!=========================================================================!
! project : rambutan
! program : constants
! history : Apr 27, 2011
! authors : duliang (duleung@gmail.com)
! purpose : constants for general use
! comment : copied from lihuang's codes
!=========================================================================!
  module constants
     implicit none

! define numerical precision consatnts
     integer, public, parameter :: sp = kind(1.0)
     integer, public, parameter :: dp = kind(1.0d0)

! define real(dp) numerical constants
     real(dp), public, parameter :: pi   = dacos(-1.0d0)
     real(dp), public, parameter :: zero = 0.0_dp
     real(dp), public, parameter :: half = 0.5_dp
     real(dp), public, parameter :: one  = 1.0_dp
     real(dp), public, parameter :: two  = 2.0_dp
     real(dp), public, parameter :: eps1  = 1.0D-1
     real(dp), public, parameter :: eps2  = 1.0D-2
     real(dp), public, parameter :: eps3  = 1.0D-3
     real(dp), public, parameter :: eps4  = 1.0D-4
     real(dp), public, parameter :: eps5  = 1.0D-5
     real(dp), public, parameter :: eps6  = 1.0D-6
     real(dp), public, parameter :: eps9  = 1.0D-9

! define complex(dp) numerical constants
     complex(dp), public, parameter :: ci    = dcmplx(0.0_dp, 1.0_dp)
     complex(dp), public, parameter :: cone  = dcmplx(1.0_dp, 0.0_dp)
     complex(dp), public, parameter :: czero = dcmplx(0.0_dp, 0.0_dp)

! unit transformation factor
     real(dp), public, parameter :: ev2cm = 8065.50000_dp
     real(dp), public, parameter :: ha2ev = 27.2113845_dp

! define file handler unit
     integer, public, parameter :: mystd = 6
     integer, public, parameter :: myout = 77
     integer, public, parameter :: mytmp = 88
     integer, public, parameter :: mytst = 99

     ! the max number of subspaces, 1024 is enough for (N,Jz)
     integer, public, parameter :: maxsubs = 1024
     ! the max dimension of each subspaces, 1024 is enough
     integer, public, parameter :: maxdim  = 1024

  end module constants
