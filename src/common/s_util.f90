!!!-----------------------------------------------------------------------
!!! project : CSSL (Common Service Subroutines Library)
!!! program : s_assert
!!!           s_str_upcase
!!!           s_str_lowcase
!!!           s_str_count
!!!           s_str_compress
!!! source  : s_util.f90
!!! type    : subroutines
!!! author  : li huang (email:huangli712@gmail.com)
!!! history : 07/10/2014 by li huang
!!!           07/14/2014 by li huang
!!! purpose : these subroutines are used to provide some useful facilities
!!!           including string manipulation, date time information, etc.
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
  subroutine s_str_upcase(s)
     implicit none

! external arguments
! input/output string
     character(*), intent(inout) :: s

! local variables
! loop index
     integer :: i

! difference between 'A' and 'a'
     integer :: diff

     diff = ichar('A') - ichar('a')

! if lowercase, make uppercase
     do i=1,len(s)
         if ( ichar(s(i:i)) >= ichar('a') .and. ichar(s(i:i)) <= ichar('z') ) then
             s(i:i) = char(ichar(s(i:i)) + diff)
         endif ! back if block
     enddo ! over i={1,len(s)} loop

     return
  end subroutine s_str_upcase

!!>>> s_str_lowcase: returns string 's' in lowercase
  subroutine s_str_lowcase(s)
     implicit none

! external arguments
! input/output string
     character(*), intent(inout) :: s

! local variables
! loop index
     integer :: i

! difference between 'A' and 'a'
     integer :: diff

     diff = ichar('A') - ichar('a')

! if uppercase, make lowercase
     do i=1,len(s)
         if ( ichar(s(i:i)) >= ichar('A') .and. ichar(s(i:i)) <= ichar('Z') ) then
             s(i:i) = char(ichar(s(i:i)) - diff)
         endif ! back if block
     enddo ! over i={1,len(s)} loop

     return
  end subroutine s_str_lowcase

!!>>> s_str_count: return the number of times a substring occurs
  subroutine s_str_count(string, substr, count)
     implicit none

! external arguments
! string to be examined
     character(len=*), intent(in) :: string

! substring in question
     character(len=*), intent(in) :: substr

! return value, number of occurrences
     integer, intent(out) :: count

! local variables
! position to start the match
     integer :: start
     integer :: offset

     count = 0
     start = 0
     do
         offset = index( string(start+1:), substr )
         if ( offset == 0 ) EXIT
         start = start + offset
         count = count + 1
     enddo ! over do loop

     return
  end subroutine s_str_count

!!>>> s_str_compress: return a copy of an input string with all whitespace
!!>>> (spaces and tabs) removed.
  subroutine s_str_compress(string)
     implicit none

! external arguments
! character string to be compressed.
     character( * ), intent(inout) :: string

! local parameters
! ASCII number for tab space ' ' and tab 
     integer, parameter :: SPACE = 32
     integer, parameter :: TAB   = 9
     integer, parameter :: NUL   = 0

! local variables
! loop index
     integer :: i
     integer :: j

! ASCII number for current character
     integer :: curr_char

! return values
! input string with all whitespace removed before the first non-whitespace
! character, and from in-between non-whitespace characters.
     character( len( string ) ) :: output

!
! definitions of a space and a tab character are made for the ASCII collating
! sequence. Each single character of the input string is checked against
! these definitions using the IACHAR() intrinsic. If the input string
! character DOESNOT correspond to a space or tab, it is not copied to
! the output string.
!
! Note that for input that ONLY has spaces or tabs BEFORE the first useful
! character, the output of this function is the same as the ADJUSTL() instrinsic.
!

! Initialise output string
     output = ' '

! initialise output string "useful" length counter
     j = 0

! loop over string elements
     do i=1,len(string)
! convert the current character to its position in the ASCII collating sequence
         curr_char = iachar( string(i:i) )
! if the character is NOT a space ' ' or a tab '->|', copy it to the output string.
         if ( curr_char /= SPACE .and. curr_char /= TAB .and. curr_char /= NUL ) then
             j = j + 1
             output(j:j) = string(i:i)
         endif ! back if block
     enddo ! over i={1,len(input_string)} loop

! copy output string to input string
     string = output

     return
  end subroutine s_str_compress
