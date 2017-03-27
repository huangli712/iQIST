!!!-----------------------------------------------------------------------
!!! project : CAPI (Common Application Programming Interface)
!!! program : japi
!!! source  : i_atomic.f90
!!! type    : module
!!! author  : li huang (email:lihuang.dmft@gmail.com)
!!! history : 10/29/2014 by li huang (created)
!!!           03/24/2017 by li huang (last modified)
!!! purpose : the purpose of this module is to define a generic and robust
!!!           application programming interface (API) for atomic eigenvalue
!!!           problem solver.
!!! status  : unstable
!!! comment :
!!!-----------------------------------------------------------------------

  module japi
     implicit none

!!========================================================================
!!>>> declare private parameters                                       <<<
!!========================================================================

! dp: number precision, double precision for reals
     integer, private, parameter :: dp = kind(1.0d0)

!!========================================================================
!!>>> declare global data structure                                    <<<
!!========================================================================

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

  end module japi
