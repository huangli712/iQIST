!!!-----------------------------------------------------------------------
!!! project : iqist @ jasmine
!!! program : atomic_print_header
!!!           atomic_print_footer
!!!           atomic_print_summary
!!! source  : atomic_print.f90
!!! type    : subroutines
!!! author  : yilin wang (email:qhwyl2006@126.com)
!!! history : 07/09/2014 by yilin wang (created)
!!!           06/11/2024 by li huang (last modified)
!!! purpose : provide printing infrastructure for the atomic eigenvalue
!!!           problem solver. the subroutines will print many runtime
!!!           information about the solver to the terminal.
!!! status  : unstable
!!! comment :
!!!-----------------------------------------------------------------------

!!
!! @sub atomic_print_header
!!
!! print the runtime information (start) to the terminal
!!
  subroutine atomic_print_header()
     use constants, only : mystd

     use version, only : V_FULL
     use version, only : V_AUTH
     use version, only : V_INST
     use version, only : V_MAIL
     use version, only : V_GPL3

     use control, only : cname

!! local variables
     ! string for current date and time
     character (len = 20) :: date_time_string

!! [body

     ! obtain current date and time
     call s_time_builder(date_time_string)

     write(mystd,'(2X,a)') cname//' (Sequential Edition)'
     write(mystd,'(2X,a)') 'A Modern Atomic Eigenvalue Problem Solver'
     write(mystd,*)

     write(mystd,'(2X,a)') 'Version: '//V_FULL//' (built at '//__TIME__//' '//__DATE__//')'
     write(mystd,'(2X,a)') 'Develop: '//V_AUTH//' ('//V_INST//')'
     write(mystd,'(2X,a)') 'Support: '//V_MAIL
     write(mystd,'(2X,a)') 'License: '//V_GPL3
     write(mystd,*)

     write(mystd,'(2X,a)') 'start running at '//date_time_string
     write(mystd,'(2X,a,i4)') 'currently using cpu cores:', 1

!! body]

     return
  end subroutine atomic_print_header

!!
!! @sub atomic_print_footer
!!
!! print the runtime information (end) to the terminal
!!
  subroutine atomic_print_footer()
     use constants, only : dp
     use constants, only : mystd

     use control, only : cname

     implicit none

!! local variables
     ! string for current date and time
     character (len = 20) :: date_time_string

     ! used to record the time usage information
     real(dp) :: tot_time

!! [body

     ! obtain time usage information
     call cpu_time(tot_time)

     ! obtain current date and time
     call s_time_builder(date_time_string)

     write(mystd,'(2X,a,f10.2,a)') cname//' >>> total time spent:', tot_time, 's'
     write(mystd,*)

     write(mystd,'(2X,a)') cname//' >>> I am tired and want to go to bed. Bye!'
     write(mystd,'(2X,a)') cname//' >>> happy ending at '//date_time_string

!! body]

     return
  end subroutine atomic_print_footer

!!
!! @sub atomic_print_summary
!!
!! print the control parameters to the terminal
!!
  subroutine atomic_print_summary()
     use constants, only : mystd

     use control ! ALL

     implicit none

!! [body

     write(mystd,'(2X,a)') '[configuration parameters] -> core control'
     write(mystd,'(2X,a)') '-----------------------------------------------------'
     write(mystd,'(4X,a16,i10,  2X,a8)') 'ibasis / value :', ibasis , 'type : i'
     write(mystd,'(4X,a16,i10,  2X,a8)') 'ictqmc / value :', ictqmc , 'type : i'
     write(mystd,'(4X,a16,i10,  2X,a8)') 'icu    / value :', icu    , 'type : i'
     write(mystd,'(4X,a16,i10,  2X,a8)') 'icf    / value :', icf    , 'type : i'
     write(mystd,'(4X,a16,i10,  2X,a8)') 'isoc   / value :', isoc   , 'type : i'
     !
     write(mystd,'(2X,a)') '[configuration parameters] -> atomic Hamiltonian'
     write(mystd,'(2X,a)') '-----------------------------------------------------'
     write(mystd,'(4X,a16,i10,  2X,a8)') 'nband  / value :', nband  , 'type : i'
     write(mystd,'(4X,a16,i10,  2X,a8)') 'nspin  / value :', nspin  , 'type : i'
     write(mystd,'(4X,a16,i10,  2X,a8)') 'norbs  / value :', norbs  , 'type : i'
     write(mystd,'(4X,a16,i10,  2X,a8)') 'ncfgs  / value :', ncfgs  , 'type : i'
     write(mystd,'(4X,a16,i10,  2X,a8)') 'nmini  / value :', nmini  , 'type : i'
     write(mystd,'(4X,a16,i10,  2X,a8)') 'nmaxi  / value :', nmaxi  , 'type : i'
     write(mystd,'(4X,a16,f10.5,2X,a8)') 'Uc     / value :', Uc     , 'type : d'
     write(mystd,'(4X,a16,f10.5,2X,a8)') 'Uv     / value :', Uv     , 'type : d'
     write(mystd,'(4X,a16,f10.5,2X,a8)') 'Jz     / value :', Jz     , 'type : d'
     write(mystd,'(4X,a16,f10.5,2X,a8)') 'Js     / value :', Js     , 'type : d'
     write(mystd,'(4X,a16,f10.5,2X,a8)') 'Jp     / value :', Jp     , 'type : d'
     write(mystd,'(4X,a16,f10.5,2X,a8)') 'Ud     / value :', Ud     , 'type : d'
     write(mystd,'(4X,a16,f10.5,2X,a8)') 'Jh     / value :', Jh     , 'type : d'
     write(mystd,'(4X,a16,f10.5,2X,a8)') 'mune   / value :', mune   , 'type : d'
     write(mystd,'(4X,a16,f10.5,2X,a8)') 'lambda / value :', lambda , 'type : d'
     !
     write(mystd,*)

!! body]

     return
  end subroutine atomic_print_summary
