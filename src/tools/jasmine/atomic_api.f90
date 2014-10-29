!!!-----------------------------------------------------------------------
!!! project : lilac
!!! program : api
!!!           api@T_jasmine
!!!           api@init_atomic
!!!           api@exec_atomic
!!!           api@stop_atomic
!!! source  : atomic_api.f90
!!! type    : module
!!! author  : li huang (email:huangli712@gmail.com)
!!! history : 10/29/2014 by li huang
!!! purpose : the purpose of this module is to define a generic and robust
!!!           application programming interface (API) for atomic eigenvalue
!!!           problem solver
!!! status  : unstable
!!! comment :
!!!-----------------------------------------------------------------------

  module api
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
! for jasmine code
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

     public  :: init_atomic
     public  :: exec_atomic
     public  :: stop_atomic

  contains ! encapsulated functionality

# if !defined (F2PY)

!!>>> init_atomic: initialize the atomic eigenvalue problem solver
!!>>> fortran version
  subroutine init_atomic(I_jasmine)
     implicit none

! external arguments
! type structure of generic solver
     class(*), intent(in) :: I_jasmine

     call cat_init_atomic(I_jasmine)

     return
  end subroutine init_atomic

# else   /* F2PY */

!!>>> init_atomic: initialize the atomic eigenvalue problem solver
!!>>> python version
  subroutine init_atomic()
     implicit none

     call cat_init_atomic()

     return
  end subroutine init_atomic

# endif  /* F2PY */

  end module api
