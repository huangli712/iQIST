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

! string_count_substring --
!     Return the number of times a substring occurs
! Arguments:
!     string     String to be examined
!     substr     Substring in question
! Result:
!     Number of occurrences
! Note:
!     Trailing blanks _are_ taken into account.
!     Possible overlaps are ignored:
!     string = 'ababab' and substr = 'abab'
!     will give the answer 1, not 2
!
function string_count_substring( string, substr ) result( count )
    character(len=*)           :: string
    character(len=*)           :: substr

    integer                    :: start
    integer                    :: count

    count  = 0
    start  = 0
    do
        start = index( string(start+1:), substr )
        count = count + 1
        if ( start == 0 ) exit
    enddo
end function string_count_substring

function real2str(num) result(str)
  integer, parameter :: max_real_len = 16
  real, intent(in) :: num
  character(len=:) ,allocatable :: str
  character(len=max_real_len) :: base_str
  write(base_str,*) num
  str = trim(adjustl(base_str))
end function real2str

function dble2str(num) result(str)
  integer, parameter :: max_dble_len = 28
  double precision, intent(in) :: num
  character(len=:) ,allocatable :: str
  character(len=max_dble_len) :: base_str
  write(base_str,*) num
  str = trim(adjustl(base_str))
end function dble2str

function int2str(num) result(str)
  integer, parameter :: max_int_len = 12
  integer, intent(in) :: num
  character(len=:) ,allocatable :: str
  character(len=max_int_len) :: base_str
  write(base_str,*) num
  str = trim(adjustl(base_str))
end function int2str

function upcase(s) result(t)
! Returns string 's' in uppercase
character(*), intent(in) :: s
character(len(s)) :: t
integer :: i, diff
t = s; diff = ichar('A')-ichar('a')
do i = 1, len(t)
    if (ichar(t(i:i)) >= ichar('a') .and. ichar(t(i:i)) <= ichar('z')) then
        ! if lowercase, make uppercase
        t(i:i) = char(ichar(t(i:i)) + diff)
    end if
end do
end function

function lowcase(s) result(t)
! Returns string 's' in lowercase
character(*), intent(in) :: s
character(len(s)) :: t
integer :: i, diff
t = s; diff = ichar('A')-ichar('a')
do i = 1, len(t)
    if (ichar(t(i:i)) >= ichar('A') .and. ichar(t(i:i)) <= ichar('Z')) then
        ! if uppercase, make lowercase
        t(i:i) = char(ichar(t(i:i)) - diff)
    end if
end do
end function
