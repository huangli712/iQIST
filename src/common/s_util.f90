!!!-----------------------------------------------------------------------
!!! project : CSSL (Common Service Subroutines Library)
!!! program : s_assert
!!!           s_str_upcase
!!!           s_str_lowcase
!!!           s_str_count
!!!           s_str_double
!!!           s_str_integer
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

!!>>> s_str_count: return the number of times a substring occurs
  function s_str_count(string, substr) result( count )
     implicit none

! external arguments
! string to be examined
     character(len=*), intent(in) :: string

! substring in question
     character(len=*), intent(in) :: substr

! return value, number of occurrences
     integer :: count

! local variables
! position to start the match
     integer :: start

     count = 0
     start = 0
     do
         start = index( string(start+1:), substr )
         count = count + 1
         if ( start == 0 ) EXIT
     enddo ! over do loop

     return
  end function s_str_count

!!>>> s_str_double: convert a real number to a string
  function s_str_double(num) result(str)
     implicit none

! external arguments
! input double precision real number
     double precision, intent(in) :: num

! return value: a string
     character(len=:) ,allocatable :: str

! local variables
! auxiliary string
     character(len=28) :: base_str

     write(base_str,*) num
     str = trim( adjustl(base_str) )

     return
  end function s_str_double

!!>>> s_str_integer: convert a integer number to a string
  function s_str_integer(num) result(str)
     implicit none

! external arguments
! input integer number
     integer, intent(in) :: num

! return value: a string
     character(len=:) ,allocatable :: str

! local variables
! auxiliary string
     character(len=12) :: base_str

     write(base_str,*) num
     str = trim( adjustl(base_str) )

     return
  end function s_str_integer
