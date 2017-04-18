!!!-----------------------------------------------------------------------
!!! project : CSSL (Common Service Subroutines Library)
!!! program : s_print_error
!!!           s_write_error
!!!           s_print_exception
!!!           s_write_exception
!!!           s_print_message
!!!           s_write_message
!!! source  : s_error.f90
!!! type    : subroutines
!!! author  : li huang (email:lihuang.dmft@gmail.com)
!!! history : 09/15/2009 by li huang (created)
!!!           04/19/2017 by li huang (last modified)
!!! purpose : these subroutines are used to display the (error/exception/
!!!           normal) messages in the console, and then STOP or CONTINUE
!!!           the code according to the error level.
!!! status  : unstable
!!! comment :
!!!-----------------------------------------------------------------------

!!
!!
!! Introduction
!! ============
!!
!! 1. display the error, exception, or message on the terminal
!! -----------------------------------------------------------
!!
!! subroutine s_print_error(...)
!! subroutine s_print_exception(...)
!! subroutine s_print_message(...)
!!
!! 2. write the error, exception, or message to external file
!! ----------------------------------------------------------
!!
!! subroutine s_write_error(...)
!! subroutine s_write_exception(...)
!! subroutine s_write_message(...)
!!
!!

!!
!! @sub s_print_error
!!
!! print the error information and STOP the program
!!
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

!!
!! @sub s_write_error
!!
!! write the error information and STOP the program
!!
  subroutine s_write_error(sub, msg, file_unit)
     implicit none

! external arguments
! subroutine name
     character(len=*), intent(in) :: sub

! error message
     character(len=*), intent(in) :: msg

! file handler
     integer, intent(in)          :: file_unit

! print error information
     write(file_unit,'(2X,4a)') 'fatal error occurred in ', sub, ': ', msg

! TERMINATE THE PROGRAM
!-------------------------------------------------------------------------
     STOP
!-------------------------------------------------------------------------

     return
  end subroutine s_write_error

!!
!! @sub s_print_exception
!!
!! print normal runtime exceptional information, and continue
!!
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

!!
!! @sub s_write_exception
!!
!! write normal runtime exceptional information, and continue
!!
  subroutine s_write_exception(sub, msg, file_unit)
     implicit none

! external arguments
! subroutine name
     character(len=*), intent(in) :: sub

! exception message
     character(len=*), intent(in) :: msg

! file handler
     integer, intent(in)          :: file_unit

! print error information
     write(file_unit,'(2X,4a)') 'runtime exception occurred in ', sub, ': ', msg

! CONTINUE/PAUSE THE PROGRAM
!-------------------------------------------------------------------------
     CONTINUE ! OR PAUSE
!-------------------------------------------------------------------------

     return
  end subroutine s_write_exception

!!
!! @sub s_print_message
!!
!! print normal runtime message to the console
!!
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

!!
!! @sub s_write_message
!!
!! write normal runtime message to the console
!!
  subroutine s_write_message(sub, msg, file_unit)
     implicit none

! external arguments
! subroutine name
     character(len=*), intent(in) :: sub

! runtime message
     character(len=*), intent(in) :: msg

! file handler
     integer, intent(in)          :: file_unit

! print error information
     write(file_unit,'(2X,4a)') 'instant message from ', sub, ': ', msg

     return
  end subroutine s_write_message
