!!!-----------------------------------------------------------------------
!!! project : hibiscus/entropy
!!! program : entropy_print_header
!!!           entropy_print_footer
!!!           entropy_print_summary
!!! source  : entropy_print.f90
!!! type    : subroutines
!!! author  : li huang (email:lihuang.dmft@gmail.com)
!!! history : 01/08/2011 by li huang
!!!           01/26/2011 by li huang
!!!           11/17/2014 by li huang
!!! purpose : provide printing infrastructure for classic maximum entropy
!!!           method code
!!! status  : very unstable
!!! comment :
!!!-----------------------------------------------------------------------

!!>>> entropy_print_header: print the startup information for classic
!!>>> maximum entropy method code
  subroutine entropy_print_header()
     use constants, only : mystd

     use control, only : nprocs

     implicit none

! string for current date and time
     character (len = 20) :: date_time_string

! obtain current date and time
     call s_time_builder(date_time_string)

     write(mystd,'(2X,a)') 'HIBISCUS/entropy'
     write(mystd,'(2X,a)') '>>> A Classic Maximum Entropy Method Code for Imaginary Time Data'
     write(mystd,*)

     write(mystd,'(2X,a)') 'Version: 2015.01.06T '//'(built at '//__TIME__//" "//__DATE__//')'
     write(mystd,'(2X,a)') 'Develop: by li huang (at IOP/CAS & SPCLab/CAEP & UNIFR)'
     write(mystd,'(2X,a)') 'Support: lihuang.dmft@gmail.com'
     write(mystd,'(2X,a)') 'License: GNU General Public License version 3'
     write(mystd,*)

     write(mystd,'(2X,a)') 'HIBISCUS/entropy >>> start running at '//date_time_string

# if defined (MPI)

     write(mystd,'(2X,a,i4)') 'HIBISCUS/entropy >>> parallelism: Yes >>> processors:', nprocs

# else   /* MPI */

     write(mystd,'(2X,a,i4)') 'HIBISCUS/entropy >>> parallelism: No  >>> processors:', 1

# endif  /* MPI */

     write(mystd,*)

     return
  end subroutine entropy_print_header

!!>>> entropy_print_footer: print the ending information for classic
!!>>> maximum entropy method code
  subroutine entropy_print_footer()
     use constants, only : dp, mystd

     implicit none

! string for current date and time
     character (len = 20) :: date_time_string

! used to record the time usage information
     real(dp) :: tot_time

! obtain time usage information
     call cpu_time(tot_time)

! obtain current date and time
     call s_time_builder(date_time_string)

     write(mystd,'(2X,a,f10.2,a)') 'HIBISCUS/entropy >>> total time spent:', tot_time, 's'
     write(mystd,*)

     write(mystd,'(2X,a)') 'HIBISCUS/entropy >>> I am tired and want to go to bed. Bye!'
     write(mystd,'(2X,a)') 'HIBISCUS/entropy >>> happy ending at '//date_time_string

     return
  end subroutine entropy_print_footer

!!>>> entropy_print_summary: print the running parameters,
!!>>> only for reference
  subroutine entropy_print_summary()
     use constants, only : mystd, ev2k

     use control ! ALL

     implicit none

     write(mystd,'(2X,a)') 'HIBISCUS/entropy >>> parameters list:'

     write(mystd,'(2(4X,a,i10)  )') 'ntime:', ntime, 'niter:', niter
     write(mystd,'(2(4X,a,i10)  )') 'nwmax:', nwmax, 'ntune:', ntune
     write(mystd,'(2(4X,a,i10)  )') 'ntype:', ntype, 'nstep:', nstep
     write(mystd,'(2(4X,a,i10)  )') 'nband:', nband, 'norbs:', norbs

     write(mystd,'(2(4X,a,f10.5))') 'ainit:', ainit, 'devia:', devia
     write(mystd,'(2(4X,a,f10.5))') 'wstep:', wstep, 'sigma:', sigma
     write(mystd,'(2(4X,a,f10.5))') 'beta :', beta , 'temp :', ev2k/beta

     write(mystd,*)

     return
  end subroutine entropy_print_summary
