!!!-----------------------------------------------------------------------
!!! project : CSSL (Common Service Subroutines Library)
!!! program : s_assert
!!!           s_str_upcase
!!!           s_str_lowcase
!!! source  : s_util.f90
!!! type    : subroutines
!!! author  : li huang (email:huangli712@gmail.com)
!!! history : 07/10/2014 by li huang
!!! purpose : these subroutines are used to
!!! status  : unstable
!!! comment :
!!!-----------------------------------------------------------------------

!!========================================================================
!!>>> assertion checker                                                <<<
!!========================================================================

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

!!========================================================================
!!>>> string manipulation                                              <<<
!!========================================================================

!!>>> s_str_upcase: returns string 's' in uppercase
  function s_str_upcase(s) result(t)
     implicit none

! external arguments
! input string
     character(*), intent(in) :: s

! output string
     character(len(s)) :: t

! local variables
! loop index
     integer :: i

! difference between 'A' and 'a'
     integer :: diff

     t = s; diff = ichar('A')-ichar('a')

! if lowercase, make uppercase
     do i=1,len(t)
         if ( ichar(t(i:i)) >= ichar('a') .and. ichar(t(i:i)) <= ichar('z') ) then
             t(i:i) = char(ichar(t(i:i)) + diff)
         endif ! back if block
     enddo ! over i={1, len(t)} loop

     return
  end function s_str_upcase

!!>>> s_str_lowcase: returns string 's' in lowercase
  function s_str_lowcase(s) result(t)
     implicit none

! external arguments
! input string
     character(*), intent(in) :: s

! output string
     character(len(s)) :: t

! local variables
! loop index
     integer :: i

! difference between 'A' and 'a'
     integer :: diff

     t = s; diff = ichar('A')-ichar('a')

! if uppercase, make lowercase
     do i=1,len(t)
         if ( ichar(t(i:i)) >= ichar('A') .and. ichar(t(i:i)) <= ichar('Z') ) then
             t(i:i) = char(ichar(t(i:i)) - diff)
         endif ! back if block
     enddo ! over i={1,len(t)} loop

     return
  end function s_str_lowcase

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


