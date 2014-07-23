!!!-----------------------------------------------------------------------
!!! project : CSSL (Common Service Subroutines Library)
!!! program : s_print_error
!!!           s_print_exception
!!!           s_print_message
!!! source  : s_error.f90
!!! type    : subroutines
!!! author  : li huang (email:huangli712@gmail.com)
!!! history : 09/15/2009 by li huang
!!!           07/08/2014 by li huang
!!!           07/23/2014 by li huang
!!! purpose : these subroutines are used to display the (error/exception/
!!!           normal) messages in the console, and then STOP orCONTINUE
!!!           the code according to the error level.
!!! status  : unstable
!!! comment :
!!!-----------------------------------------------------------------------

!!>>> s_print_error: print the error information and STOP the program
  subroutine s_print_error(sub, msg)
     implicit none

! external arguments
! subroutine name
     character(len=*), intent(in) :: sub

! error message
     character(len=*), intent(in) :: msg

! print error information
     write(*,'(2X,4a)') 'fatal error occurred in ', sub, ': ', msg

! TERMINATE THE PROGRAM
!-------------------------------------------------------------------------
     STOP
!-------------------------------------------------------------------------

     return
  end subroutine s_print_error

!!>>> s_print_exception: print normal runtime exceptional information, and continue
  subroutine s_print_exception(sub, msg)
     implicit none

! external arguments
! subroutine name
     character(len=*), intent(in) :: sub

! exception message
     character(len=*), intent(in) :: msg

! print error information
     write(*,'(2X,4a)') 'runtime exception occurred in ', sub, ': ', msg

! CONTINUE/PAUSE THE PROGRAM
!-------------------------------------------------------------------------
     CONTINUE ! OR PAUSE
!-------------------------------------------------------------------------

     return
  end subroutine s_print_exception

!!>>> s_print_message: print normal runtime message to the console
  subroutine s_print_message(sub, msg)
     implicit none

! external arguments
! subroutine name
     character(len=*), intent(in) :: sub

! runtime message
     character(len=*), intent(in) :: msg

! print error information
     write(*,'(2X,4a)') 'instant message from ', sub, ': ', msg

     return
  end subroutine s_print_message
