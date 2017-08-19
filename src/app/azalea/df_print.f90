

  subroutine df_print_header()
     use constants, only : mystd
     use version, only : FULL_VER

     use df_control, only : cname
     use df_control, only : nprocs

     implicit none

! string for current date and time
     character (len = 20) :: date_time_string

! obtain current date and time
     call s_time_builder(date_time_string)

     write(mystd,'(2X,a)') cname
     write(mystd,'(2X,a)') '>>> A Modern Dual Fermion Framework For Quantum Lattice Models'
     write(mystd,*)

     write(mystd,'(2X,a)') 'Version: '//FULL_VER//' (built at '//__TIME__//" "//__DATE__//')'
     write(mystd,'(2X,a)') 'Develop: by li huang (at SPCLab/CAEP)'
     write(mystd,'(2X,a)') 'Support: lihuang.dmft@gmail.com'
     write(mystd,'(2X,a)') 'License: GNU General Public License version 3'
     write(mystd,*)

     write(mystd,'(2X,a)') cname//' >>> start running at '//date_time_string

# if defined (MPI)

     write(mystd,'(2X,a,i4)') cname//' >>> parallelism: Yes >>> processors:', nprocs

# else   /* MPI */

     write(mystd,'(2X,a,i4)') cname//' >>> parallelism: No  >>> processors:', 1

# endif  /* MPI */

     write(mystd,*)

     return
  end subroutine df_print_header

  subroutine df_print_footer()
     use constants, only : dp, mystd

     use df_control, only : cname

     implicit none

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
  end subroutine df_print_footer

  subroutine df_print_summary()
     use constants, only : mystd, ev2k

     use df_control ! ALL

     implicit none

     write(mystd,'(2X,a)') cname//' >>> parameters list:'

     write(mystd,*)

     return
  end subroutine df_print_summary

  subroutine df_print_runtime()
     implicit none

     return
  end subroutine df_print_runtime

  subroutine df_print_it_info()
     implicit none

     return
  end subroutine df_print_it_info
