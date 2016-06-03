!!!-----------------------------------------------------------------------
!!! project : daisy
!!! program : hfqmc_print_header
!!!           hfqmc_print_footer
!!!           hfqmc_print_summary
!!!           hfqmc_print_runtime
!!!           hfqmc_print_it_info
!!!           hfqmc_print_diagram
!!! source  : hfqmc_print.f90
!!! type    : subroutines
!!! author  : li huang (email:lihuang.dmft@gmail.com)
!!! history : 10/24/2008 by li huang (created)
!!!           08/17/2015 by li huang (last modified)
!!! purpose : provide printing infrastructure for Hirsch-Fye quantum Monte
!!!           Carlo (HFQMC) quantum impurity solver
!!! status  : unstable
!!! comment :
!!!-----------------------------------------------------------------------

!!>>> hfqmc_print_header: print the startup information for Hirsch-Fye
!!>>> quantum Monte Carlo quantum impurity solver plus dynamical mean
!!>>> field theory self-consistent engine
  subroutine hfqmc_print_header()
     use constants, only : mystd

     use control, only : cname
     use control, only : nprocs

     implicit none

! string for current date and time
     character (len = 20) :: date_time_string

! obtain current date and time
     call s_time_builder(date_time_string)

     write(mystd,'(2X,a)') cname
     write(mystd,'(2X,a)') '>>> A DMFT Engine With Hirsch-Fye Quantum Monte Carlo Impurity Solver'
     write(mystd,*)

     write(mystd,'(2X,a)') 'Version: 2016.02.13T '//'(built at '//__TIME__//" "//__DATE__//')'
     write(mystd,'(2X,a)') 'Develop: by li huang (at IOP/CAS & SPCLab/CAEP & UNIFR)'
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
  end subroutine hfqmc_print_header

!!>>> hfqmc_print_footer: print the ending information for Hirsch-Fye
!!>>> quantum Monte Carlo quantum impurity solver plus dynamical mean
!!>>> field theory self-consistent engine
  subroutine hfqmc_print_footer()
     use constants, only : dp, mystd

     use control, only : cname

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
  end subroutine hfqmc_print_footer

!!>>> hfqmc_print_summary: print the running parameters, only for reference
  subroutine hfqmc_print_summary()
     use constants, only : mystd, ev2k

     use control ! ALL

     implicit none

     write(mystd,'(2X,a)') cname//' >>> parameters list:'

     write(mystd,'(2(4X,a,i10))')   'isscf :', isscf  , 'isbin :', isbin
     write(mystd,'(2(4X,a,i10))')   'issun :', issun  , 'isspn :', isspn

     write(mystd,'(2(4X,a,i10))')   'mstep :', mstep  , 'mfreq :', mfreq
     write(mystd,'(2(4X,a,i10))')   'nband :', nband  , 'nspin :', nspin
     write(mystd,'(2(4X,a,i10))')   'norbs :', norbs  , 'niter :', niter

     write(mystd,'(2(4X,a,i10))')   'nsing :', nsing  , 'ntime :', ntime
     write(mystd,'(2(4X,a,i10))')   'ntherm:', ntherm , 'nsweep:', nsweep
     write(mystd,'(2(4X,a,i10))')   'nclean:', nclean , 'ncarlo:', ncarlo

     write(mystd,'(2(4X,a,f10.5))') 'Uc    :', Uc     , 'Jz    :', Jz
     write(mystd,'(2(4X,a,f10.5))') 'mune  :', mune   , 'beta  :', beta
     write(mystd,'(2(4X,a,f10.5))') 'part  :', part   , 'temp  :', ev2k/beta

     write(mystd,*)

     return
  end subroutine hfqmc_print_summary

!!>>> hfqmc_print_runtime: print the runtime information, including the
!!>>> iteration messages and statistic data, only for reference
  subroutine hfqmc_print_runtime(iter, nstep, accept, reject, tcount)
     use constants, only : dp, one, mystd

     use control, only : cname
     use control, only : norbs
     use control, only : nsweep
     use context, only : ktep

     implicit none

! external arguments
! current iteration number
     integer, intent(in)  :: iter

! current QMC effective sweep count
     integer, intent(in)  :: nstep

! accepted QMC flip count
     real(dp), intent(in) :: accept

! rejected QMC flip count
     real(dp), intent(in) :: reject

! total QMC flip count
     real(dp), intent(in) :: tcount

! local variables
! loop index
     integer :: i

! about iteration number
     write(mystd,'(2X,a,i3,2(a,i10))') cname//' >>> iter:', iter, ' sweep:', nstep, ' of ', nsweep

! about update action
     write(mystd,'(4X,a)')        'hfqmc sampling statistics:'
     write(mystd,'(4X,a,3i12)')   'count:', int(tcount), int(accept), int(reject)
     write(mystd,'(4X,a,3f12.5)') 'ratio:', one, accept / tcount, reject / tcount

     write(mystd,'(4X,a)')        'delayed update statistics:'
     write(mystd,'(4X,a,10i6)')   'count:', (ktep(i), i=1,norbs)

     return
  end subroutine hfqmc_print_runtime

!!>>> hfqmc_print_it_info: print the iteration information to the screen
  subroutine hfqmc_print_it_info(iter)
     use constants, only : mystd

     use control, only : cname

     implicit none

! external arguments
! current iteration number
     integer, intent(in) :: iter

! according to the value of iter, we can judge whether the impurity solver
! is in the binning mode.
     if ( iter /= 999 ) then
         write(mystd,'(2X,a,i3,a)') cname//' >>> DMFT iter:', iter, ' <<< SELFING'
     else
         write(mystd,'(2X,a,i3,a)') cname//' >>> DMFT iter:', iter, ' <<< BINNING'
     endif ! back if ( iter /= 999 ) block

     return
  end subroutine hfqmc_print_it_info

!!>>> hfqmc_diagram_plotting: write out a snapshot for the current diagram
!!>>> configuration, the results can be used to make a dynamical video.
  subroutine hfqmc_print_diagram(iter, nstep, nfast)
     use constants, only : mystd, mytmp, zero, one

     use control, only : niter
     use control, only : nsing, ntime, nsweep
     use context, only : imat

     implicit none

! external arguments
! current iteration number
     integer, intent(in) :: iter

! current QMC effective sweep count
     integer, intent(in) :: nstep

! internal fast cycle numbers, used to control the size of data bins
     integer, intent(in) :: nfast

! local variables
! loop index
     integer :: i
     integer :: j

! write the snapshot
! open data file: solver.diag.dat
     open(mytmp, file='solver.diag.dat', form='formatted', status='unknown', position='append')

! write diagram info
     write(mytmp,'(2(a,i4))') '>> cur_iter:', iter, ' tot_iter:', niter
     write(mytmp,'(2(a,i4))') '>> cur_diag:', nstep/nfast, ' tot_diag:', nsweep/nfast

! write the position of operators
     do i=1,ntime
         write(mytmp,'(i4)', advance = 'no') i
         do j=1,nsing
             if ( imat(i,j) > zero ) then
                 write(mytmp,'(f5.1)', advance = 'no') +one
             else
                 write(mytmp,'(f5.1)', advance = 'no') -one
             endif ! back if ( imat(i,j) > zero ) block
         enddo ! over j={1,nsing} loop
         write(mytmp,*)
     enddo ! over i={1,ntime} loop

! write two blank lines
     write(mytmp,*)
     write(mytmp,*)

! close data file
     close(mytmp)

! write the message to the terminal
     write(mystd,'(4X,a)') '>>> quantum impurity solver config: saving'

     return
  end subroutine hfqmc_print_diagram
