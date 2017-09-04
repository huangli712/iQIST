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
     use constants, only : mystd

     use control ! ALL

     implicit none

     write(mystd,'(2X,a)') 'configuration parameters -> lattice model'
     write(mystd,'(2X,a)') '----------------------------------------------------'
     write(mystd,'(4X,a8,i10,  2X,a8)') 'nband  /', nband , 'type / i'
     write(mystd,'(4X,a8,i10,  2X,a8)') 'nspin  /', nspin , 'type / i'
     write(mystd,'(4X,a8,i10,  2X,a8)') 'norbs  /', norbs , 'type / i'
     write(mystd,'(4X,a8,i10,  2X,a8)') 'nkpts  /', nkpts , 'type / i'
     write(mystd,'(4X,a8,i10,  2X,a8)') 'nkp_x  /', nkp_x , 'type / i'
     write(mystd,'(4X,a8,i10,  2X,a8)') 'nkp_y  /', nkp_y , 'type / i'
     write(mystd,'(4X,a8,i10,  2X,a8)') 'nkp_z  /', nkp_z , 'type / i'
     write(mystd,'(4X,a8,f10.5,2X,a8)') 'mune   /', mune  , 'type / d'
     write(mystd,'(4X,a8,f10.5,2X,a8)') 'beta   /', beta  , 'type / d'
     write(mystd,'(4X,a8,f10.5,2X,a8)') 'part   /', part  , 'type / d'

     write(mystd,'(2X,a)') 'configuration parameters -> dual fermion engine'
     write(mystd,'(2X,a)') '----------------------------------------------------'
     write(mystd,'(4X,a8,i10,  2X,a8)') 'nffrq  /', nffrq , 'type / i'
     write(mystd,'(4X,a8,i10,  2X,a8)') 'nbfrq  /', nbfrq , 'type / i'
     write(mystd,'(4X,a8,i10,  2X,a8)') 'ndfit  /', ndfit , 'type / i'
     write(mystd,'(4X,a8,i10,  2X,a8)') 'nbsit  /', nbsit , 'type / i'
     write(mystd,'(4X,a8,f10.5,2X,a8)') 'dfmix  /', dfmix , 'type / d'
     write(mystd,'(4X,a8,f10.5,2X,a8)') 'bsmix  /', bsmix , 'type / d'

     write(mystd,*)

     return
  end subroutine dt_print_summary

!!
!! @sub dt_print_control
!!
!! print the control parameters, only for reference
!!
  subroutine dt_print_control()
     implicit none

     return
  end subroutine dt_print_control

!!
!! @sub dt_print_runtime
!!
!! print the runtime information, including some physical observables and
!! statistic data, only for reference
!!
  subroutine dt_print_runtime()
     implicit none

     return
  end subroutine dt_print_runtime

!!
!! @sub dt_print_it_info
!!
!! print the iteration information to the screen
!!
  subroutine dt_print_it_info()
     implicit none

     return
  end subroutine dt_print_it_info
