!!!-----------------------------------------------------------------------
!!! project : azalea
!!! program : dt_print_header
!!!           dt_print_footer
!!!           dt_print_summary
!!!           dt_print_control
!!!           dt_print_runtime
!!!           dt_print_it_info
!!! source  : dt_print.f90
!!! type    : subroutines
!!! author  : li huang (email:lihuang.dmft@gmail.com)
!!! history : 09/15/2009 by li huang (created)
!!!           01/04/2018 by li huang (last modified)
!!! purpose : provide printing infrastructure for dual fermion engine.
!!! status  : unstable
!!! comment :
!!!-----------------------------------------------------------------------

!!
!! @sub dt_print_header
!!
!! print the startup information for dual fermion engine 
!!
  subroutine dt_print_header()
     use constants, only : mystd

     use version, only : V_FULL
     use version, only : V_AUTH
     use version, only : V_MAIL
     use version, only : V_GPL3

     use control, only : cname
     use control, only : nprocs

     implicit none

! local variables
! string for current date and time
     character (len = 20) :: date_time_string

! obtain current date and time
     call s_time_builder(date_time_string)

# if defined (MPI)

     write(mystd,'(2X,a)') cname//' (parallelized edition)'

# else   /* MPI */

     write(mystd,'(2X,a)') cname//' (sequential edition)'

# endif  /* MPI */

     write(mystd,'(2X,a)') 'A Modern Dual Fermion Framework For Quantum Lattice Models'
     write(mystd,*)

     write(mystd,'(2X,a)') 'Version: '//V_FULL//' (built at '//__TIME__//" "//__DATE__//')'
     write(mystd,'(2X,a)') 'Develop: '//V_AUTH
     write(mystd,'(2X,a)') 'Support: '//V_MAIL
     write(mystd,'(2X,a)') 'License: '//V_GPL3
     write(mystd,*)

     write(mystd,'(2X,a)') 'start running at '//date_time_string

# if defined (MPI)

     write(mystd,'(2X,a,i4)') 'currently using cpu cores:', nprocs

# else   /* MPI */

     write(mystd,'(2X,a,i4)') 'currently using cpu cores:', 1

# endif  /* MPI */

     return
  end subroutine dt_print_header

!!
!! @sub dt_print_footer
!!
!! print the ending information for dual fermion engine
!!
  subroutine dt_print_footer()
     use constants, only : dp
     use constants, only : mystd

     use control, only : cname

     implicit none

! local variables
! string for current date and time
     character (len = 20) :: date_time_string

! used to record the time usage information
     real(dp) :: tot_time

! obtain time usage information
     call cpu_time(tot_time)

! obtain current date and time
     call s_time_builder(date_time_string)

     write(mystd,'(2X,a,f10.2,a)') cname//' >>> total time spent:', tot_time, 's'
     write(mystd,*)

     write(mystd,'(2X,a)') cname//' >>> I am tired and want to go to bed. Bye!'
     write(mystd,'(2X,a)') cname//' >>> happy ending at '//date_time_string

     return
  end subroutine dt_print_footer

!!
!! @sub dt_print_summary
!!
!! print the running parameters, only for reference
!!
  subroutine dt_print_summary()
     use constants, only : mystd, ev2k

     use control ! ALL

     implicit none

     write(mystd,'(2X,a)') cname//' >>> parameters list:'

     write(mystd,*)

     return
  end subroutine dt_print_summary

  subroutine dt_print_runtime()
     implicit none

     return
  end subroutine dt_print_runtime

  subroutine dt_print_it_info()
     implicit none

     return
  end subroutine dt_print_it_info
