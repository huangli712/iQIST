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
  end module api
