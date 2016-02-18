!!!-----------------------------------------------------------------------
!!! project : hibiscus/stoch
!!! program : sac_print_header
!!!           sac_print_footer
!!!           sac_print_summary
!!!           sac_print_runtime
!!! source  : sac_print.f90
!!! type    : subroutines
!!! author  : li huang (email:lihuang.dmft@gmail.com)
!!! history : 01/08/2011 by li huang (created)
!!!           08/17/2015 by li huang (last modified)
!!! purpose : provide printing infrastructure for stochastic analytic
!!!           continuation code
!!! status  : very unstable
!!! comment :
!!!-----------------------------------------------------------------------

!!>>> sac_print_header: print the startup information for stochastic
!!>>> analytic continuation code
  subroutine sac_print_header()
     use constants, only : mystd

     use control, only : nprocs

     implicit none

! string for current date and time
     character (len = 20) :: date_time_string

! obtain current date and time
     call s_time_builder(date_time_string)

     write(mystd,'(2X,a)') 'HIBISCUS/stoch'
     write(mystd,'(2X,a)') '>>> A Stochastic Analytic Continuation Code for Imaginary Time Data'
     write(mystd,*)

     write(mystd,'(2X,a)') 'Version: 2016.02.13T '//'(built at '//__TIME__//" "//__DATE__//')'
     write(mystd,'(2X,a)') 'Develop: by li huang (at IOP/CAS & SPCLab/CAEP & UNIFR)'
     write(mystd,'(2X,a)') 'Support: lihuang.dmft@gmail.com'
     write(mystd,'(2X,a)') 'License: GNU General Public License version 3'
     write(mystd,*)

     write(mystd,'(2X,a)') 'HIBISCUS/stoch >>> start running at '//date_time_string

# if defined (MPI)

     write(mystd,'(2X,a,i4)') 'HIBISCUS/stoch >>> parallelism: Yes >>> processors:', nprocs

# else   /* MPI */

     write(mystd,'(2X,a,i4)') 'HIBISCUS/stoch >>> parallelism: No  >>> processors:', 1

# endif  /* MPI */

     write(mystd,*)

     return
  end subroutine sac_print_header

!!>>> sac_print_footer: print the ending information for stochastic
!!>>> analytic continuation code
  subroutine sac_print_footer()
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

     write(mystd,'(2X,a,f10.2,a)') 'HIBISCUS/stoch >>> total time spent:', tot_time, 's'
     write(mystd,*)

     write(mystd,'(2X,a)') 'HIBISCUS/stoch >>> I am tired and want to go to bed. Bye!'
     write(mystd,'(2X,a)') 'HIBISCUS/stoch >>> happy ending at '//date_time_string

     return
  end subroutine sac_print_footer

!!>>> sac_print_summary: print the running parameters, only for reference
  subroutine sac_print_summary()
     use constants, only : mystd, ev2k

     use control ! ALL

     implicit none

     write(mystd,'(2X,a)') 'HIBISCUS/stoch >>> parameters list:'

     write(mystd,'(2(4X,a,i10)  )') 'ntime:', ntime, 'nalph:', nalph
     write(mystd,'(2(4X,a,i10)  )') 'nwmax:', nwmax, 'nwarm:', nwarm
     write(mystd,'(2(4X,a,i10)  )') 'ngrid:', ngrid, 'nstep:', nstep
     write(mystd,'(2(4X,a,i10)  )') 'ngamm:', ngamm, 'ndump:', ndump
     write(mystd,'(2(4X,a,i10)  )') 'ltype:', ltype, 'lemax:', lemax

     write(mystd,'(2(4X,a,f10.5))') 'ainit:', ainit, 'wstep:', wstep
     write(mystd,'(2(4X,a,f10.5))') 'ratio:', ratio, 'sigma:', sigma
     write(mystd,'(2(4X,a,f10.5))') 'eta1 :', eta1 , 'eta2 :', eta2
     write(mystd,'(2(4X,a,f10.5))') 'beta :', beta , 'temp :', ev2k/beta

     write(mystd,*)

     return
  end subroutine sac_print_summary

!!>>> sac_print_runtime: print the runtime information, including the
!!>>> statistic data, only for reference
  subroutine sac_print_runtime(step, time_start, time_end)
     use constants, only : dp, zero, mystd

     use control, only : nalph, nstep
     use context, only : move_accept, move_tcount
     use context, only : swap_accept, swap_tcount

     implicit none

! external arguments
! current monte carlo sweep number
     real(dp), intent(in) :: step

! start time point
     real(dp), intent(in) :: time_start

! end time point
     real(dp), intent(in) :: time_end

! local variables
! loop index
     integer :: i

     write(mystd,'(4X,a,i8,2X,a,i8)')   '>>> iter:', int(step), ' of ', nstep
     do i=1,nalph
         if ( nalph == 1 ) then
             write(mystd,'(2(8X,a,f10.4))') 'move:', move_accept(i) / move_tcount(i), &
                                            'swap:', zero
         else
             write(mystd,'(2(8X,a,f10.4))') 'move:', move_accept(i) / move_tcount(i), &
                                            'swap:', swap_accept(i) / swap_tcount(i)
         endif ! back if ( nalph == 1 ) block
     enddo ! over i={1,nalph} loop
     write(mystd,'(4X,a,f10.4,a)') '>>> time used in this interval:', time_end - time_start, 's'
     write(mystd,*)

     return
  end subroutine sac_print_runtime
