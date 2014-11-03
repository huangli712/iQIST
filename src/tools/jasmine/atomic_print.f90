!!!-----------------------------------------------------------------------
!!! project : jasmine
!!! program : atomic_print_header
!!!           atomic_print_footer
!!!           atomic_print_summary
!!! source  : atomic_print.f90
!!! type    : subroutines
!!! author  : yilin wang (email: qhwyl2006@126.com)
!!! history : 07/09/2014 by yilin wang
!!!           08/22/2014 by yilin wang
!!!           10/20/2014 by li huang
!!! purpose : provide printing infrastructure for the atomic eigenvalue
!!!           problem solver
!!! status  : unstable
!!! comment :
!!!-----------------------------------------------------------------------

!!>>> atomic_print_header: print the running header
  subroutine atomic_print_header()
     use constants, only : mystd

! string for current date and time
     character (len = 20) :: date_time_string

! obtain current date and time
     call s_time_builder(date_time_string)

     write(mystd,'(2X,a)') 'JASMINE'
     write(mystd,'(2X,a)') '>>> An Atomic Eigenvalue Problem Solver'
     write(mystd,*)

     write(mystd,'(2X,a)') 'Version: 2014.10.29T '//'(built at '//__TIME__//" "//__DATE__//')'
     write(mystd,'(2X,a)') 'Develop: by yilin wang (at IOP/CAS)'
     write(mystd,'(2X,a)') '         by li huang (at IOP/CAS & SPCLab/CAEP & UNIFR)'
     write(mystd,'(2X,a)') 'Support: qhwyl2006@126.com'
     write(mystd,'(2X,a)') 'License: GNU General Public License version 3'
     write(mystd,*)

     write(mystd,'(2X,a)') 'JASMINE >>> start running at '//date_time_string
     write(mystd,*)

     return
  end subroutine atomic_print_header

!!>>> atomic_print_footer: print running footer
  subroutine atomic_print_footer()
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

     write(mystd,'(2X,a,f10.2,a)') 'JASMINE >>> total time spent:', tot_time, 's'
     write(mystd,*)

     write(mystd,'(2X,a)') 'JASMINE >>> I am tired and want to go to bed. Bye!'
     write(mystd,'(2X,a)') 'JASMINE >>> happy ending at '//date_time_string

     return
  end subroutine atomic_print_footer

!!>>> atomic_print_summary: print summary of parameters list
  subroutine atomic_print_summary()
     use constants, only : mystd

     use control ! ALL

     implicit none

     write(mystd,'(2X,a)') 'JASMINE >>> parameters list:'

     write(mystd,'(2(4X,a,i10))')   'ibasis :', ibasis , 'ictqmc :', ictqmc
     write(mystd,'(2(4X,a,i10))')   'icu    :', icu    , 'icf    :', icf
     write(mystd,'(1(4X,a,i10))')   'isoc   :', isoc

     write(mystd,'(2(4X,a,i10))')   'nband  :', nband  , 'nspin  :', nspin
     write(mystd,'(2(4X,a,i10))')   'norbs  :', norbs  , 'ncfgs  :', ncfgs
     write(mystd,'(2(4X,a,i10))')   'nmini  :', nmini  , 'nmaxi  :', nmaxi 

     write(mystd,'(2(4X,a,f10.5))') 'Uc     :', Uc     , 'Uv     :', Uv
     write(mystd,'(2(4X,a,f10.5))') 'Jz     :', Jz     , 'Js     :', Js
     write(mystd,'(1(4X,a,f10.5))') 'Jp     :', Jp
     write(mystd,'(2(4X,a,f10.5))') 'Ud     :', Ud     , 'Jh     :', Jh

     write(mystd,'(2(4X,a,f10.5))') 'mune   :', mune   , 'lambda :', lambda

     write(mystd,*)

     return
  end subroutine atomic_print_summary
