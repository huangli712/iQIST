!!!-----------------------------------------------------------------------
!!! project : CSSL (Common Service Subroutines Library)
!!! program : s_assert
!!!           s_str_upcase
!!!           s_str_lowcase
!!!           s_str_count
!!!           s_str_double
!!!           s_str_integer
!!!           s_str_compress
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

!!>>> s_str_compress: return a copy of an input string with all whitespace
!!>>> (spaces and tabs) removed.
  function s_str_compress(input_string) result (output_string)
     implicit none

! external arguments
! character string to be compressed.
     character( * ), intent(in) :: input_string

! return values
! input string with all whitespace removed before the first non-whitespace
! character, and from in-between non-whitespace characters.
     character( len( input_string ) ) :: output_string

! local parameters
     integer, parameter :: IACHAR_SPACE = 32
     integer, parameter :: IACHAR_TAB   = 9

! local variables
     integer :: i
     integer :: j
     integer :: curr_char

!
! Definitions of a space and a tab character are made for the ASCII collating
! sequence. Each single character of the input string is checked against
! these definitions using the IACHAR() intrinsic. If the input string
! character DOESNOT correspond to a space or tab, it is not copied to
! the output string.
!
! Note that for input that ONLY has spaces or tabs BEFORE the first useful
! character, the output of this function is the same as the ADJUSTL() instrinsic.
!

! Initialise output string
     output_string = ' '

! initialise output string "useful" length counter
     j = 0

! loop over string elements
     do i=1,len(input_string)
! convert the current character to its position in the ASCII collating sequence
         curr_char = iachar( input_string(i:i) )
! if the character is NOT a space ' ' or a tab '->|', copy it to the output string.
         if ( curr_char /= IACHAR_SPACE .and. curr_char /= IACHAR_TAB ) then
             j = j + 1
             output_string(j:j) = input_string(i:i)
         endif ! back if block
     enddo ! over i={1,len(input_string)} loop

     return
  end function s_str_compress
