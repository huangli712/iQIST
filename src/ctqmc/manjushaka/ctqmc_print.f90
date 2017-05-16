!!!-----------------------------------------------------------------------
!!! project : manjushaka
!!! program : ctqmc_print_header
!!!           ctqmc_print_footer
!!!           ctqmc_print_summary
!!!           ctqmc_print_control
!!!           ctqmc_print_runtime
!!!           ctqmc_print_it_info
!!! source  : ctqmc_print.f90
!!! type    : subroutines
!!! author  : li huang (email:lihuang.dmft@gmail.com)
!!! history : 09/15/2009 by li huang (created)
!!!           05/16/2017 by li huang (last modified)
!!! purpose : provide printing infrastructure for hybridization expansion
!!!           version continuous time quantum Monte Carlo (CTQMC) quantum
!!!           impurity solver and dynamical mean field theory (DMFT) self
!!!           -consistent engine.
!!! status  : unstable
!!! comment :
!!!-----------------------------------------------------------------------

!!
!! @sub ctqmc_print_header
!!
!! print the startup information for continuous time quantum Monte Carlo
!! quantum impurity solver plus dynamical mean field theory engine
!!
  subroutine ctqmc_print_header()
     use constants, only : mystd

     use control, only : cname
     use control, only : nprocs

     implicit none

! string for current date and time
     character (len = 20) :: date_time_string

! obtain current date and time
     call s_time_builder(date_time_string)

     write(mystd,'(2X,a)') cname
     write(mystd,'(2X,a)') '>>> A Modern Continuous Time Quantum Monte Carlo Impurity Solver'
     write(mystd,*)

     write(mystd,'(2X,a)') 'Version: 2016.02.13T '//'(built at '//__TIME__//" "//__DATE__//')'
     write(mystd,'(2X,a)') 'Develop: by li huang (at IOP/CAS & SPCLab/CAEP & UNIFR) and yilin wang (at IOP/CAS)'
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
  end subroutine ctqmc_print_header

!!>>> ctqmc_print_footer: print the ending information for continuous time
!!>>> quantum Monte Carlo quantum impurity solver plus dynamical mean field
!!>>> theory self-consistent engine
  subroutine ctqmc_print_footer()
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
  end subroutine ctqmc_print_footer

!!>>> ctqmc_print_summary: print the running parameters, only for reference
  subroutine ctqmc_print_summary()
     use constants, only : mystd, ev2k

     use control ! ALL

     implicit none

     write(mystd,'(2X,a)') cname//' >>> parameters list:'

     write(mystd,'(2(4X,a,i10))')   'isscf :', isscf  , 'isbin :', isbin
     write(mystd,'(2(4X,a,i10))')   'issun :', issun  , 'isspn :', isspn
     write(mystd,'(2(4X,a,i10))')   'isort :', isort  , 'issus :', issus
     write(mystd,'(1(4X,a,i10))')   'isvrt :', isvrt
     write(mystd,'(2(4X,a,i10))')   'ifast :', ifast  , 'itrun :', itrun 

     write(mystd,'(2(4X,a,i10))')   'lemax :', lemax  , 'legrd :', legrd
     write(mystd,'(2(4X,a,i10))')   'chmax :', chmax  , 'chgrd :', chgrd
     write(mystd,'(2(4X,a,i10))')   'mkink :', mkink  , 'mfreq :', mfreq
     write(mystd,'(2(4X,a,i10))')   'nband :', nband  , 'nspin :', nspin
     write(mystd,'(2(4X,a,i10))')   'norbs :', norbs  , 'ncfgs :', ncfgs
     write(mystd,'(1(4X,a,i10))')   'niter :', niter
     write(mystd,'(2(4X,a,i10))')   'nffrq :', nffrq  , 'nbfrq :', nbfrq
     write(mystd,'(2(4X,a,i10))')   'nfreq :', nfreq  , 'ntime :', ntime
     write(mystd,'(2(4X,a,i10))')   'npart :', npart  , 'nflip :', nflip

     write(mystd,'(2(4X,a,i10))')   'ntherm:', ntherm , 'nsweep:', nsweep
     write(mystd,'(2(4X,a,i10))')   'nclean:', nclean , 'nwrite:', nwrite
     write(mystd,'(2(4X,a,i10))')   'nmonte:', nmonte , 'ncarlo:', ncarlo

     write(mystd,'(2(4X,a,f10.5))') 'U     :', U      , 'Uc    :', Uc
     write(mystd,'(2(4X,a,f10.5))') 'Js    :', Js     , 'Uv    :', Uv
     write(mystd,'(2(4X,a,f10.5))') 'Jp    :', Jp     , 'Jz    :', Jz
     write(mystd,'(2(4X,a,f10.5))') 'mune  :', mune   , 'beta  :', beta
     write(mystd,'(2(4X,a,f10.5))') 'part  :', part   , 'temp  :', ev2k/beta

     write(mystd,*)

     return
  end subroutine ctqmc_print_summary

!!>>> ctqmc_print_runtime: print the runtime information, including physical
!!>>> observables and statistic data, only for reference
  subroutine ctqmc_print_runtime(iter, cstep)
     use constants, only : dp, one, half, mystd

     use control, only : cname
     use control, only : nsweep, nmonte
     use context, only : cnegs, caves
     use context, only : insert_tcount, insert_accept, insert_reject
     use context, only : remove_tcount, remove_accept, remove_reject
     use context, only : lshift_tcount, lshift_accept, lshift_reject
     use context, only : rshift_tcount, rshift_accept, rshift_reject
     use context, only : reflip_tcount, reflip_accept, reflip_reject
     use context, only : paux

     implicit none

! external arguments
! current self-consistent iteration number
     integer, intent(in) :: iter

! current QMC sweeping steps
     integer, intent(in) :: cstep

! local variables
! real(dp) dummy variables
     real(dp) :: raux

! about iteration number
     write(mystd,'(2X,a,i3,2(a,i10))') cname//' >>> iter:', iter, ' sweep:', cstep, ' of ', nsweep

! about auxiliary physical observables
     raux = real(caves) / nmonte
     write(mystd,'(4X,a)')        'auxiliary system observables:'
     write(mystd,'(2(4X,a,f10.5))') 'etot :', paux(1) / raux, 'epot :', paux(2) / raux
     write(mystd,'(2(4X,a,f10.5))') 'ekin :', paux(3) / raux, '<Sz> :', paux(4) / raux
     write(mystd,'(2(4X,a,f10.5))') '<N1> :', paux(5) / raux, '<N2> :', paux(6) / raux
     write(mystd,'(2(4X,a,e10.5))') '<K2> :', paux(7) / raux, '<K3> :', paux(8) / raux
     write(mystd,'(1(4X,a,e10.5))') '<K4> :', paux(9) / raux

! about insert action
     if ( insert_tcount <= half ) insert_tcount = -one ! if insert is disable
     write(mystd,'(4X,a)')        'insert kink statistics:'
     write(mystd,'(4X,a,3i10)')   'count:', int( insert_tcount ), int( insert_accept ), int( insert_reject )
     write(mystd,'(4X,a,3f10.5)') 'ratio:', one, insert_accept / insert_tcount, insert_reject / insert_tcount

! about remove action
     if ( remove_tcount <= half ) remove_tcount = -one ! if remove is disable
     write(mystd,'(4X,a)')        'remove kink statistics:'
     write(mystd,'(4X,a,3i10)')   'count:', int( remove_tcount ), int( remove_accept ), int( remove_reject )
     write(mystd,'(4X,a,3f10.5)') 'ratio:', one, remove_accept / remove_tcount, remove_reject / remove_tcount

! about lshift action
     if ( lshift_tcount <= half ) lshift_tcount = -one ! if lshift is disable
     write(mystd,'(4X,a)')        'lshift kink statistics:'
     write(mystd,'(4X,a,3i10)')   'count:', int( lshift_tcount ), int( lshift_accept ), int( lshift_reject )
     write(mystd,'(4X,a,3f10.5)') 'ratio:', one, lshift_accept / lshift_tcount, lshift_reject / lshift_tcount

! about rshift action
     if ( rshift_tcount <= half ) rshift_tcount = -one ! if rshift is disable
     write(mystd,'(4X,a)')        'rshift kink statistics:'
     write(mystd,'(4X,a,3i10)')   'count:', int( rshift_tcount ), int( rshift_accept ), int( rshift_reject )
     write(mystd,'(4X,a,3f10.5)') 'ratio:', one, rshift_accept / rshift_tcount, rshift_reject / rshift_tcount

! about reflip action
     if ( reflip_tcount <= half ) reflip_tcount = -one ! if reflip is disable
     write(mystd,'(4X,a)')        'global flip statistics:'
     write(mystd,'(4X,a,3i10)')   'count:', int( reflip_tcount ), int( reflip_accept ), int( reflip_reject )
     write(mystd,'(4X,a,3f10.5)') 'ratio:', one, reflip_accept / reflip_tcount, reflip_reject / reflip_tcount

! about negative sign
     write(mystd,'(4X,a,i10)')    'negative sign counter:', cnegs
     write(mystd,'(4X,a,f10.5)')  'averaged sign sampler:', caves / real(cstep)

     return
  end subroutine ctqmc_print_runtime

!!>>> ctqmc_print_it_info: print the iteration information to the screen
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
