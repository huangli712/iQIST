!!!-----------------------------------------------------------------------
!!! project : hibiscus/stoch
!!! program : sac_print_header
!!!           sac_print_footer
!!!           sac_print_summary
!!!           sac_print_runtime
!!! source  : sac_print.f90
!!! type    : subroutines
!!! author  : li huang (email:huangli712@gmail.com)
!!! history : 01/08/2011 by li huang
!!!           01/10/2011 by li huang
!!!           11/18/2014 by li huang
!!! purpose : provide printing infrastructure for stochastic analytic
!!!           continuation code
!!! status  : very unstable
!!! comment :
!!!-----------------------------------------------------------------------

!>>> print the startup information for stochastic analytic continuation code
  subroutine sac_print_header()
     use constants
     use control, only : nprocs

     implicit none

     write(mystd,'(2X,a)') 'HIBISCUS'
     write(mystd,'(2X,a)') '>>> A Stochastic Analytic Continuation Code for Imaginary Time Data'
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
  end subroutine sac_print_header

!>>> print the ending information for stochastic analytic continuation code
  subroutine sac_print_footer()
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
  end subroutine sac_print_footer

!>>> print the running parameters, only for reference
  subroutine sac_print_summary()
     use constants
     use control

     implicit none

     write(mystd,'(2X,a)') 'HIBISCUS >>> parameters list:'

     write(mystd,'(2(4X,a,i10)  )') 'ntime:', ntime, 'nalph:', nalph
     write(mystd,'(2(4X,a,i10)  )') 'nwmax:', nwmax, 'nwarm:', nwarm
     write(mystd,'(2(4X,a,i10)  )') 'ngrid:', ngrid, 'nstep:', nstep
     write(mystd,'(2(4X,a,i10)  )') 'ngamm:', ngamm, 'ndump:', ndump
     write(mystd,'(2(4X,a,i10)  )') 'ltype:', ltype, 'lemax:', lemax

     write(mystd,'(2(4X,a,f10.5))') 'ainit:', ainit, 'eta1 :', eta1
     write(mystd,'(2(4X,a,f10.5))') 'ratio:', ratio, 'sigma:', sigma
     write(mystd,'(2(4X,a,f10.5))') 'beta :', beta , 'wstep:', wstep

     write(mystd,*)

     return
  end subroutine sac_print_summary

!>>> print the runtime information, including statistic data, only for reference
  subroutine sac_print_runtime(step, time_start, time_end)
     use constants
     use control
     use context

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
         endif
     enddo ! over i={1,nalph} loop
     write(mystd,'(4X,a,f10.4,a)') '>>> time used in this interval:', time_end - time_start, 's'
     write(mystd,*)

     return
  end subroutine sac_print_runtime
