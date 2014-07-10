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

pure function count_chars(str,chars) result(n)
  character(len=*), intent(in) :: str, chars
  integer :: n
  integer :: pos,inc
  n = 0
  pos = 0
  do
     inc = scan(str(pos + 1:),chars)
     if (inc == 0) exit
     pos = pos + inc
     n = n + 1
  end do
end function count_chars


function real2str(num) result(str)
  integer, parameter, private :: max_real_len = 16
  real, intent(in) :: num
  character(len=:) ,allocatable :: str
  character(len=max_real_len) :: base_str
  write(base_str,*) num
  str = trim(adjustl(base_str))
end function real2str

function dble2str(num) result(str)
  integer, parameter, private :: max_dble_len = 28
  double precision, intent(in) :: num
  character(len=:) ,allocatable :: str
  character(len=max_dble_len) :: base_str
  write(base_str,*) num
  str = trim(adjustl(base_str))
end function dble2str

function int2str(num) result(str)
  integer, parameter, private :: max_int_len = 12
  integer, intent(in) :: num
  character(len=:) ,allocatable :: str
  character(len=max_int_len) :: base_str
  write(base_str,*) num
  str = trim(adjustl(base_str))
end function int2str
