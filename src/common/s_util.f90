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
!!>>> assertion checker                                                <<<
!!========================================================================

!!>>> s_sorter: using bubble sort algorithm to sort a real dataset, the slowest algorithm
  subroutine s_sorter(nsize, list)
     use constants, only : dp

     implicit none

! external arguments
! grab the number of values from the calling code
     integer, intent(in) :: nsize

! dataset to be sorted
     real(dp), intent(inout) :: list(nsize)

! local variables
! dataset index
     integer  :: i = 0
     integer  :: j = 0

! dummy variables
     real(dp) :: swap

! basically we just loop through every element to compare it against
! every other element
! this loop increments i which is our starting point for the comparison
     sort_loop1: do i=nsize,1,-1
! this loop increments j which is the ending point for the comparison
         sort_loop2: do j=1,i-1
! swap the two elements here
             exchange: if ( list(j) > list(j+1) ) then
                 swap = list(j)
                 list(j) = list(j+1)
                 list(j+1) = swap
             endif exchange
         enddo sort_loop2 ! over j={1,i-1} loop
     enddo sort_loop1 ! over i={nsize,1,-1} loop

     return
  end subroutine s_sorter

!!>>> s_qsorter: sets up for the quick sort recursive method
  subroutine s_qsorter(nsize, list)
     use constants, only : dp

     implicit none

! external arguments
! grab the number of values from the calling code
     integer, intent(in) :: nsize

! dataset to be sorted
     real(dp), intent(inout) :: list(nsize)

! kicks off the recursive process
     call s_qscorer(1, nsize, nsize, list)

     return
  end subroutine s_qsorter

!!>>> s_qscorer: this is the actually recursive portion of the quicksort algorithm
!!>>> note: do not call it directly, please use s_qsorter() insteadly
  recursive &
  subroutine s_qscorer(pstart, pend, nsize, list)
     use constants, only : dp

     implicit none

! external arguments
! start point
     integer, intent(in) :: pstart

! end point
     integer, intent(in) :: pend

! size of array
     integer, intent(in) :: nsize

! dataset to be sorted
     real(dp), intent(inout) :: list(nsize)

! local variables
! used to find out list(left) > kaux and list(right) < kaux
     integer  :: left, right

! used to record list(pstart)
     real(dp) :: kaux

! used to swap data
     real(dp) :: taux

! setup left and right
     left = pstart
     right = pend + 1

! only in right > left, the data is to be sorted
     if ( right > left ) then

! record list(pstart) at first
         kaux = list(pstart)

         do while ( .true. )

! find out where list(left) < kaux
             do while ( .true. )
                 left = left + 1
                 if ( list(left)  > kaux .or. left  >= pend   ) EXIT
             enddo ! over do while loop

! find out where list(right) > kaux
             do while ( .true. )
                 right = right - 1
                 if ( list(right) < kaux .or. right <= pstart ) EXIT
             enddo ! over do while loop

! we should ensure right is larger than left
             if ( right <= left ) EXIT

! exchange data between list(left) and list(right)
             taux = list(left)
             list(left) = list(right)
             list(right) = taux

         enddo ! over do while loop

! exchange data between list(pstart) and list(right)
        list(pstart) = list(right)
        list(right) = kaux

! sort data from pstart to right-1
        call s_qscorer(pstart, right-1, nsize, list)

! sort data from right+1 to pend
        call s_qscorer(right+1, pend, nsize, list)

     endif ! back if ( right > left ) block

     return
  end subroutine s_qscorer

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
