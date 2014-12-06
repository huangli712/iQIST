!!!-----------------------------------------------------------------------
!!! project : lilac
!!! program : dapi
!!!           dapi@T_daisy
!!!           dapi@init_hfqmc
!!!           dapi@exec_hfqmc
!!!           dapi@stop_hfqmc
!!! source  : hfqmc_api.f90
!!! type    : module
!!! author  : li huang (email:huangli712@gmail.com)
!!! history : 12/06/2014 by li huang
!!! purpose : the purpose of this module is to define a generic and robust
!!!           application programming interface (API) for the Hirsch-Fye
!!!           quantum Monte Carlo impurity solver
!!! status  : unstable
!!! comment :
!!!-----------------------------------------------------------------------

  module dapi
     implicit none

!!========================================================================
!!>>> declare global parameters                                        <<<
!!========================================================================

! dp: number precision, double precision for reals
     integer, private, parameter :: dp = kind(1.0d0)

!!========================================================================
!!>>> declare global data structure                                    <<<
!!========================================================================

! note: now f2py does not support derived types, so we have to comment
! out them when f2py is used.

# if !defined (F2PY)

! define type T_jasmine, which is used to describe the control parameters
! for the jasmine code
     public :: T_jasmine
     type :: T_jasmine
         integer :: ibasis
         integer :: ictqmc
         integer :: icu
         integer :: icf
         integer :: isoc

         integer :: nband
         integer :: nspin
         integer :: norbs
         integer :: ncfgs

         integer :: nmini
         integer :: nmaxi

         real(dp) :: Uc
         real(dp) :: Uv
         real(dp) :: Jz
         real(dp) :: Js
         real(dp) :: Jp
         real(dp) :: Ud
         real(dp) :: Jh
         real(dp) :: mune
         real(dp) :: lambda
     end type T_jasmine

# endif  /* F2PY */

!!========================================================================
!!>>> declare accessibility for module routines                        <<<
!!========================================================================

     public  :: init_hfqmc
     public  :: exec_hfqmc
     public  :: stop_hfqmc

  contains ! encapsulated functionality

  subroutine init_hfqmc()
  end subroutine init_hfqmc

  subroutine exec_hfqmc()
  end subroutine exec_hfqmc

  subroutine stop_hfqmc()
  end subroutine stop_hfqmc

  end module dapi
