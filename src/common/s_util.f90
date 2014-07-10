!!!-----------------------------------------------------------------------
!!! project : CSSL (Common Service Subroutines Library)
!!! program : s_assert
!!! source  : s_util.f90
!!! type    : subroutines
!!! author  : li huang (email:huangli712@gmail.com)
!!! history : 07/10/2014 by li huang
!!!           07/10/2014 by li huang
!!! purpose : these subroutines are used to
!!! status  : unstable
!!! comment :
!!!-----------------------------------------------------------------------

!!>>> s_assert: fortran version of assert
  subroutine s_assert(condition)
     implicit none

! external arguments
     logical, intent(in) :: condition

! if condition == .false., it aborts the program.
     if ( .not. condition ) then
         call s_print_error("s_assert", "assert failed.")
     endif

     return
  end subroutine s_assert
