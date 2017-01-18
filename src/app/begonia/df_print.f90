

  subroutine df_print_header()
     implicit none

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
     implicit none

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
