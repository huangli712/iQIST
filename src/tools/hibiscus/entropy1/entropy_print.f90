!-------------------------------------------------------------------------
! project : hibiscus
! program : entropy_print_header
!           entropy_print_footer
!           entropy_print_summary
!           entropy_print_error
!           entropy_print_exception
! source  : entropy_print.f90
! type    : subroutines
! author  : li huang (email:huangli712@yahoo.com.cn)
! history : 01/08/2011 by li huang
!           01/10/2011 by li huang
!           01/20/2011 by li huang
!           01/26/2011 by li huang
! purpose : provide printing infrastructure for classic maximum entropy
!           method code
! input   :
! output  :
! status  : very unstable
! comment :
!-------------------------------------------------------------------------

!>>> print the startup information for classic maximum entropy method code
  subroutine entropy_print_header()
     use constants
     use control, only : nprocs

     implicit none

     write(mystd,'(2X,a)') 'HIBISCUS'
     write(mystd,'(2X,a)') '>>> A Classic Maximum Entropy Method Code for Imaginary Time Data'
     write(mystd,*)

     write(mystd,'(2X,a)') 'version: 2011.08.18T            '
     write(mystd,'(2X,a)') 'develop: by li huang, CAEP & IOP'
     write(mystd,'(2X,a)') 'support: huangli712@yahoo.com.cn'
     write(mystd,'(2X,a)') 'license: GPL2 and later versions'
     write(mystd,*)

     write(mystd,'(2X,a)') 'HIBISCUS >>> running'

# if defined (MPI)

     write(mystd,'(2X,a,i4)') 'HIBISCUS >>> parallelism: Yes >>> processors:', nprocs

# else   /* MPI */

     write(mystd,'(2X,a,i4)') 'HIBISCUS >>> parallelism: No  >>> processors:', 1

# endif  /* MPI */

     write(mystd,*)

     return
  end subroutine entropy_print_header

!>>> print the ending information for classic maximum entropy method code
  subroutine entropy_print_footer()
     use constants

     implicit none

! used to record the time information
     real(dp) :: tot_time

! obtain time information
     call cpu_time(tot_time)

     write(mystd,'(2X,a,f10.2,a)') 'HIBISCUS >>> total time spent:', tot_time, 's'
     write(mystd,*)

     write(mystd,'(2X,a)') 'HIBISCUS >>> I am tired and want to go to bed. Bye!'
     write(mystd,'(2X,a)') 'HIBISCUS >>> ending'

     return
  end subroutine entropy_print_footer

!>>> print the running parameters, only for reference
  subroutine entropy_print_summary()
     use constants
     use control

     implicit none

     write(mystd,'(2X,a)') 'HIBISCUS >>> parameters list:'

     write(mystd,'(2(4X,a,i10)  )') 'ntime:', ntime, 'niter:', niter
     write(mystd,'(2(4X,a,i10)  )') 'nwmax:', nwmax, 'ntune:', ntune
     write(mystd,'(2(4X,a,i10)  )') 'ntype:', ntype, 'nstep:', nstep
     write(mystd,'(2(4X,a,i10)  )') 'nband:', nband, 'norbs:', norbs

     write(mystd,'(2(4X,a,f10.5))') 'ainit:', ainit, 'devia:', devia
     write(mystd,'(2(4X,a,f10.5))') 'wstep:', wstep, 'sigma:', sigma
     write(mystd,'(2(4X,a,f10.5))') 'beta :', beta , 'temp :', ev2k / beta

     write(mystd,*)

     return
  end subroutine entropy_print_summary

!>>> print the error information and STOP the program
  subroutine entropy_print_error(sub, msg)
     use constants

     implicit none

! external arguments
! subroutine name
     character(len=*), intent(in) :: sub

! error message
     character(len=*), intent(in) :: msg

! print error information
     write(mystd,'(2X,4a)') 'fatal error occurred in ', sub, ': ', msg

! TERMINATE THE PROGRAM
!-------------------------------------------------------------------------
     STOP
!-------------------------------------------------------------------------

     return
  end subroutine entropy_print_error

!>>> print normal runtime exceptional information, and continue
  subroutine entropy_print_exception(sub, msg)
     use constants

     implicit none

! external arguments
! subroutine name
     character(len=*), intent(in) :: sub

! exception message
     character(len=*), intent(in) :: msg

! print error information
     write(mystd,'(2X,4a)') 'runtime exception occurred in ', sub, ': ', msg

! CONTINUE/PAUSE THE PROGRAM
!-------------------------------------------------------------------------
     CONTINUE ! OR PAUSE
!-------------------------------------------------------------------------

     return
  end subroutine entropy_print_exception
