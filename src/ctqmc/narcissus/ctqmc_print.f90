!!!-----------------------------------------------------------------------
!!! project : narcissus
!!! program : ctqmc_print_header
!!!           ctqmc_print_footer
!!!           ctqmc_print_summary
!!!           ctqmc_print_runtime
!!!           ctqmc_print_it_info
!!! source  : ctqmc_print.f90
!!! type    : subroutines
!!! author  : li huang (email:lihuang.dmft@gmail.com)
!!! history : 09/15/2009 by li huang (created)
!!!           04/29/2017 by li huang (last modified)
!!! purpose : provide printing infrastructure for hybridization expansion
!!!           version continuous time quantum Monte Carlo (CTQMC) quantum
!!!           impurity solver and dynamical mean field theory (DMFT) self
!!!           -consistent engine
!!! status  : unstable
!!! comment :
!!!-----------------------------------------------------------------------

!!
!! @sub ctqmc_print_header
!!
!! print the startup information for continuous time quantum Monte Carlo
!! quantum impurity solver plus dynamical mean field theory
!! self-consistent engine
!!
  subroutine ctqmc_print_header()
     use constants, only : mystd

     use version, only : FULL_VER
     use version, only : AUTH_VER
     use version, only : MAIL_VER
     use version, only : GPL3_VER

     use control, only : cname
     use control, only : nprocs

     implicit none

! local variables
! string for current date and time
     character (len = 20) :: date_time_string

! obtain current date and time
     call s_time_builder(date_time_string)

     write(mystd,'(2X,a)') cname
     write(mystd,'(2X,a)') '>>> A Modern Continuous Time Quantum Monte Carlo Impurity Solver'
     write(mystd,*)

     write(mystd,'(2X,a)') 'Version: '//FULL_VER//' (built at '//__TIME__//" "//__DATE__//')'
     write(mystd,'(2X,a)') 'Develop: '//AUTH_VER
     write(mystd,'(2X,a)') 'Support: '//MAIL_VER
     write(mystd,'(2X,a)') 'License: '//GPL3_VER
     write(mystd,*)

     write(mystd,'(2X,a)') cname//' >>> start running at '//date_time_string

# if defined (MPI)

     write(mystd,'(2X,a,i4)') cname//' >>> parallelism: Yes >>> processors:', nprocs

# else   /* MPI */

     write(mystd,'(2X,a,i4)') cname//' >>> parallelism: No  >>> processors:', 1

# endif  /* MPI */

     write(mystd,*)

     return
  end subroutine ctqmc_print_header

!!
!! @sub ctqmc_print_footer
!!
!! print the ending information for continuous time quantum Monte Carlo
!! quantum impurity solver plus dynamical mean field theory
!! self-consistent engine
!!
  subroutine ctqmc_print_footer()
     use constants, only : dp, mystd

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
  end subroutine ctqmc_print_footer

!!
!! @sub ctqmc_print_summary
!!
!! print the running parameters, only for reference
!!
  subroutine ctqmc_print_summary()
     use constants, only : mystd, ev2k

     use control ! ALL

     implicit none

     write(mystd,'(2X,a)') cname//' >>> parameters list:'

     write(mystd,'(2(4X,a,i10))')   'isscf :', isscf  , 'isscr :', isscr
     write(mystd,'(2(4X,a,i10))')   'isbnd :', isbnd  , 'isspn :', isspn
     write(mystd,'(2(4X,a,i10))')   'isbin :', isbin  , 'iswor :', iswor
     write(mystd,'(2(4X,a,i10))')   'isort :', isort  , 'isobs :', isobs
     write(mystd,'(2(4X,a,i10))')   'issus :', issus  , 'isvrt :', isvrt

     write(mystd,'(2(4X,a,i10))')   'lemax :', lemax  , 'legrd :', legrd
     write(mystd,'(2(4X,a,i10))')   'mkink :', mkink  , 'mfreq :', mfreq
     write(mystd,'(2(4X,a,i10))')   'nband :', nband  , 'nspin :', nspin
     write(mystd,'(2(4X,a,i10))')   'norbs :', norbs  , 'ncfgs :', ncfgs
     write(mystd,'(2(4X,a,i10))')   'niter :', niter  , 'nfreq :', nfreq
     write(mystd,'(2(4X,a,i10))')   'nffrq :', nffrq  , 'nbfrq :', nbfrq
     write(mystd,'(2(4X,a,i10))')   'ntime :', ntime  , 'nflip :', nflip
     write(mystd,'(2(4X,a,i10))')   'ntherm:', ntherm , 'nsweep:', nsweep
     write(mystd,'(2(4X,a,i10))')   'nclean:', nclean , 'nwrite:', nwrite
     write(mystd,'(2(4X,a,i10))')   'nmonte:', nmonte , 'ncarlo:', ncarlo

     write(mystd,'(2(4X,a,f10.5))') 'U     :', U      , 'Uc    :', Uc
     write(mystd,'(2(4X,a,f10.5))') 'Js    :', Js     , 'Uv    :', Uv
     write(mystd,'(2(4X,a,f10.5))') 'Jp    :', Jp     , 'Jz    :', Jz
     write(mystd,'(2(4X,a,f10.5))') 'lc    :', lc     , 'wc    :', wc
     write(mystd,'(2(4X,a,f10.5))') 'mune  :', mune   , 'beta  :', beta
     write(mystd,'(2(4X,a,f10.5))') 'part  :', part   , 'alpha :', alpha
     write(mystd,'(1(4X,a,f10.5))') 'system temperature:', ev2k / beta

     write(mystd,*)

     return
  end subroutine ctqmc_print_summary

!!
!! @sub ctqmc_print_runtime
!!
!! print the runtime information, including physical observables and
!! statistic data, only for reference
!!
  subroutine ctqmc_print_runtime(iter, cstep)
     use constants, only : one, half, mystd

     use control, only : cname
     use control, only : nsweep, nmonte
     use context, only : ins_t, ins_a, ins_r
     use context, only : rmv_t, rmv_a, rmv_r
     use context, only : lsh_t, lsh_a, lsh_r
     use context, only : rsh_t, rsh_a, rsh_r
     use context, only : rfl_t, rfl_a, rfl_r
     use context, only : paux

     implicit none

! external arguments
! current self-consistent iteration number
     integer, intent(in) :: iter

! current QMC sweeping steps
     integer, intent(in) :: cstep

! local variables
! integer dummy variables
     integer :: istat

! about iteration number
     write(mystd,'(2X,a,i3,2(a,i10))') cname//' >>> iter:', iter, ' sweep:', cstep, ' of ', nsweep

! about auxiliary physical observables
     istat = cstep / nmonte
     write(mystd,'(4X,a)')        'auxiliary system observables:'
     write(mystd,'(2(4X,a,f10.5))') 'etot :', paux(1) / istat, 'epot :', paux(2) / istat
     write(mystd,'(2(4X,a,f10.5))') 'ekin :', paux(3) / istat, '<Sz> :', paux(4) / istat
     write(mystd,'(2(4X,a,f10.5))') '<N1> :', paux(5) / istat, '<N2> :', paux(6) / istat
     write(mystd,'(2(4X,a,e10.5))') '<K2> :', paux(7) / istat, '<K3> :', paux(8) / istat
     write(mystd,'(1(4X,a,e10.5))') '<K4> :', paux(9) / istat

! about insert action
     if ( ins_t <= half ) ins_t = -one ! if insert is disable
     write(mystd,'(4X,a)')        'insert kink statistics:'
     write(mystd,'(4X,a,3i10)')   'count:', int( ins_t ), int( ins_a ), int( ins_r )
     write(mystd,'(4X,a,3f10.5)') 'ratio:', one, ins_a / ins_t, ins_r / ins_t

! about remove action
     if ( rmv_t <= half ) rmv_t = -one ! if remove is disable
     write(mystd,'(4X,a)')        'remove kink statistics:'
     write(mystd,'(4X,a,3i10)')   'count:', int( rmv_t ), int( rmv_a ), int( rmv_r )
     write(mystd,'(4X,a,3f10.5)') 'ratio:', one, rmv_a / rmv_t, rmv_r / rmv_t

! about lshift action
     if ( lsh_t <= half ) lsh_t = -one ! if lshift is disable
     write(mystd,'(4X,a)')        'lshift kink statistics:'
     write(mystd,'(4X,a,3i10)')   'count:', int( lsh_t ), int( lshift_accept ), int( lshift_reject )
     write(mystd,'(4X,a,3f10.5)') 'ratio:', one, lshift_accept / lsh_t, lshift_reject / lsh_t

! about rshift action
     if ( rsh_t <= half ) rsh_t = -one ! if rshift is disable
     write(mystd,'(4X,a)')        'rshift kink statistics:'
     write(mystd,'(4X,a,3i10)')   'count:', int( rsh_t ), int( rshift_accept ), int( rshift_reject )
     write(mystd,'(4X,a,3f10.5)') 'ratio:', one, rshift_accept / rsh_t, rshift_reject / rsh_t

! about reflip action
     if ( reflip_tcount <= half ) reflip_tcount = -one ! if reflip is disable
     write(mystd,'(4X,a)')        'global flip statistics:'
     write(mystd,'(4X,a,3i10)')   'count:', int( reflip_tcount ), int( reflip_accept ), int( reflip_reject )
     write(mystd,'(4X,a,3f10.5)') 'ratio:', one, reflip_accept / reflip_tcount, reflip_reject / reflip_tcount

     return
  end subroutine ctqmc_print_runtime

!!
!! @sub ctqmc_print_it_info
!!
!! print the iteration information to the screen
!!
  subroutine ctqmc_print_it_info(iter)
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
  end subroutine ctqmc_print_it_info
