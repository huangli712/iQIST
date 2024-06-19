!!!-----------------------------------------------------------------------
!!! project : iqist @ manjushaka
!!! program : ctqmc_print_header
!!!           ctqmc_print_footer
!!!           ctqmc_print_summary
!!!           ctqmc_print_control
!!!           ctqmc_print_runtime
!!!           ctqmc_print_it_info
!!! source  : ctqmc_print.f90
!!! type    : subroutines
!!! author  : li huang (email:huangli@caep.cn)
!!! history : 09/15/2009 by li huang (created)
!!!           06/20/2017 by li huang (last modified)
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

     use version, only : V_FULL
     use version, only : V_AUTH
     use version, only : V_INST
     use version, only : V_MAIL
     use version, only : V_GPL3

     use control, only : cname
     use control, only : nprocs

     implicit none

!! local variables
     ! string for current date and time
     character (len = 20) :: date_time_string

!! [body

     ! obtain current date and time
     call s_time_builder(date_time_string)

# if defined (MPI)

     write(mystd,'(2X,a)') cname//' (Parallelized Edition)'

# else   /* MPI */

     write(mystd,'(2X,a)') cname//' (Sequential Edition)'

# endif  /* MPI */

     write(mystd,'(2X,a)') 'A Modern Continuous Time Quantum Monte Carlo Impurity Solver'
     write(mystd,*)

     write(mystd,'(2X,a)') 'Version: '//V_FULL//' (built at '//__TIME__//' '//__DATE__//')'
     write(mystd,'(2X,a)') 'Develop: '//V_AUTH//' ('//V_INST//')'
     write(mystd,'(2X,a)') 'Support: '//V_MAIL
     write(mystd,'(2X,a)') 'License: '//V_GPL3
     write(mystd,*)

     write(mystd,'(2X,a)') 'start running at '//date_time_string

# if defined (MPI)

     write(mystd,'(2X,a,i4)') 'currently using cpu cores:', nprocs

# else   /* MPI */

     write(mystd,'(2X,a,i4)') 'currently using cpu cores:', 1

# endif  /* MPI */

!! body]

     return
  end subroutine ctqmc_print_header

!!
!! @sub ctqmc_print_footer
!!
!! print the ending information for continuous time quantum Monte Carlo
!! quantum impurity solver plus dynamical mean field theory engine
!!
  subroutine ctqmc_print_footer()
     use constants, only : dp
     use constants, only : mystd

     use control, only : cname

     implicit none

!! local variables
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

!! body]

     return
  end subroutine ctqmc_print_footer

!!
!! @sub ctqmc_print_summary
!!
!! print the running parameters, only for reference
!!
  subroutine ctqmc_print_summary()
     use constants, only : mystd

     use control ! ALL

     implicit none

     write(mystd,'(2X,a)') 'configuration parameters for global control'
     write(mystd,'(2X,a)') '----------------------------------------------------'
     write(mystd,'(4X,a,i10,  a4,a10)') 'isscf  /', isscf , '/', 'integer'
     write(mystd,'(4X,a,i10,  a4,a10)') 'isbnd  /', isbnd , '/', 'integer'
     write(mystd,'(4X,a,i10,  a4,a10)') 'isspn  /', isspn , '/', 'integer'
     write(mystd,'(4X,a,i10,  a4,a10)') 'isopt  /', isopt , '/', 'integer'
     write(mystd,'(4X,a,i10,  a4,a10)') 'iscut  /', iscut , '/', 'integer'
     !write(mystd,'(4X,a,i10,  a4,a10)') 'isbin  /', isbin , '/', 'integer'
     write(mystd,'(4X,a,i10,  a4,a10)') 'iswor  /', iswor , '/', 'integer'
     write(mystd,'(4X,a,i10,  a4,a10)') 'isort  /', isort , '/', 'integer'
     write(mystd,'(4X,a,i10,  a4,a10)') 'isobs  /', isobs , '/', 'integer'
     write(mystd,'(4X,a,i10,  a4,a10)') 'issus  /', issus , '/', 'integer'
     write(mystd,'(4X,a,i10,  a4,a10)') 'isvrt  /', isvrt , '/', 'integer'

     write(mystd,'(2X,a)') 'configuration parameters for self-consistent engine'
     write(mystd,'(2X,a)') '----------------------------------------------------'
     write(mystd,'(4X,a,i10,  a4,a10)') 'niter  /', niter , '/', 'integer'
     write(mystd,'(4X,a,f10.5,a4,a10)') 'alpha  /', alpha , '/', ' double'

     write(mystd,'(2X,a)') 'configuration parameters for quantum impurity model'
     write(mystd,'(2X,a)') '----------------------------------------------------'
     write(mystd,'(4X,a,i10,  a4,a10)') 'nband  /', nband , '/', 'integer'
     write(mystd,'(4X,a,i10,  a4,a10)') 'nspin  /', nspin , '/', 'integer'
     write(mystd,'(4X,a,i10,  a4,a10)') 'norbs  /', norbs , '/', 'integer'
     write(mystd,'(4X,a,i10,  a4,a10)') 'ncfgs  /', ncfgs , '/', 'integer'
     write(mystd,'(4X,a,f10.5,a4,a10)') 'Uc     /', Uc    , '/', ' double'
     write(mystd,'(4X,a,f10.5,a4,a10)') 'Jz     /', Jz    , '/', ' double'
     write(mystd,'(4X,a,f10.5,a4,a10)') 'mune   /', mune  , '/', ' double'
     write(mystd,'(4X,a,f10.5,a4,a10)') 'beta   /', beta  , '/', ' double'
     write(mystd,'(4X,a,f10.5,a4,a10)') 'part   /', part  , '/', ' double'

     write(mystd,'(2X,a)') 'configuration parameters for quantum impurity solver'
     write(mystd,'(2X,a)') '----------------------------------------------------'
     write(mystd,'(4X,a,i10,  a4,a10)') 'lemax  /', lemax , '/', 'integer'
     write(mystd,'(4X,a,i10,  a4,a10)') 'legrd  /', legrd , '/', 'integer'
     write(mystd,'(4X,a,i10,  a4,a10)') 'mkink  /', mkink , '/', 'integer'
     write(mystd,'(4X,a,i10,  a4,a10)') 'mfreq  /', mfreq , '/', 'integer'
     write(mystd,'(4X,a,i10,  a4,a10)') 'nffrq  /', nffrq , '/', 'integer'
     write(mystd,'(4X,a,i10,  a4,a10)') 'nbfrq  /', nbfrq , '/', 'integer'
     write(mystd,'(4X,a,i10,  a4,a10)') 'nfreq  /', nfreq , '/', 'integer'
     write(mystd,'(4X,a,i10,  a4,a10)') 'ntime  /', ntime , '/', 'integer'
     write(mystd,'(4X,a,i10,  a4,a10)') 'npart  /', npart , '/', 'integer'
     write(mystd,'(4X,a,i10,  a4,a10)') 'nflip  /', nflip , '/', 'integer'
     write(mystd,'(4X,a,i10,  a4,a10)') 'ntherm /', ntherm, '/', 'integer'
     write(mystd,'(4X,a,i10,  a4,a10)') 'nsweep /', nsweep, '/', 'integer'
     write(mystd,'(4X,a,i10,  a4,a10)') 'nwrite /', nwrite, '/', 'integer'
     write(mystd,'(4X,a,i10,  a4,a10)') 'nclean /', nclean, '/', 'integer'
     write(mystd,'(4X,a,i10,  a4,a10)') 'nmonte /', nmonte, '/', 'integer'
     write(mystd,'(4X,a,i10,  a4,a10)') 'ncarlo /', ncarlo, '/', 'integer'

     write(mystd,*)

!! body]

     return
  end subroutine ctqmc_print_summary

!!
!! @sub ctqmc_print_control
!!
!! print the control parameters, only for reference
!!
  subroutine ctqmc_print_control()
     use constants, only : mystd

     use control, only : cname               ! code name
                                             !
     use control, only : isscf               ! control running scheme
     use control, only : isbnd, isspn        ! control symmetry
     use control, only : isopt, iscut        ! control optimized algorithms
     !use control, only : isbin, iswor, isort ! control measurement tricks
     use control, only : iswor, isort ! control measurement tricks
     use control, only : isobs, issus, isvrt ! control physical observables

     implicit none

! local variables
! loop index
     integer :: i

! predefined strings for control parameters
     character (len = 4) :: scf(2) = ['nscf', 'scf']
     character (len = 3) :: bnd(2) = ['no', 'yes']
     character (len = 3) :: spn(2) = ['no', 'yes']
     character (len = 4) :: opt(2) = ['lazy', 'tbd']
     character (len = 3) :: cut(2) = ['no', 'yes']
     character (len = 3) :: bin(2) = ['no', 'yes']
     character (len = 3) :: wor(2) = ['no', 'yes']
     character (len = 3) :: ort(3) = ['std', 'leg', 'svd']
     character (len = 8) :: obs(4) = ['none', 'kinetic', 'fidelity', 'binder']
     character (len = 4) :: sus(5) = ['none', 'sp_t', 'ch_t', 'sp_w', 'ch_w']
     character (len = 4) :: vrt(3) = ['none', 'twop', 'pair']

! predefined strings for control parameters
     character (len = 99) :: str_obs
     character (len = 99) :: str_sus
     character (len = 99) :: str_vrt

! build str_obs according to isobs
     str_obs = ''
     do i=1,size(obs)
         if ( btest(isobs, i-1) ) then
             str_obs = ( trim( str_obs ) // ' ' // trim( obs(i) ) )
         endif ! back if ( btest(isobs, i-1) ) block
     enddo ! over i={1,size(obs)} loop
     str_obs = adjustl(str_obs)

! build str_sus according to issus
     str_sus = ''
     do i=1,size(sus)
         if ( btest(issus, i-1) ) then
             str_sus = ( trim( str_sus ) // ' ' // trim( sus(i) ) )
         endif ! back if ( btest(issus, i-1) ) block
     enddo ! over i={1,size(sus)} loop
     str_sus = adjustl(str_sus)

! build str_vrt according to isvrt
     str_vrt = ''
     do i=1,size(vrt)
         if ( btest(isvrt, i-1) ) then
             str_vrt = ( trim( str_vrt ) // ' ' // trim( vrt(i) ) )
         endif ! back if ( btest(isvrt, i-1) ) block
     enddo ! over i={1,size(vrt)} loop
     str_vrt = adjustl(str_vrt)

! write control parameters
     write(mystd,'(2X,a)') cname//' >>> CTQMC quantum impurity solver running'

     write(mystd,'(4X,a,i4,3X,2a)') 'self-consistent scheme  /', isscf, '/ ', scf(isscf)
     write(mystd,'(4X,a,i4,3X,2a)') 'symmetry (band part)    /', isbnd, '/ ', bnd(isbnd)
     write(mystd,'(4X,a,i4,3X,2a)') 'symmetry (spin part)    /', isspn, '/ ', spn(isspn)
     write(mystd,'(4X,a,i4,3X,2a)') 'fast trace evaluation   /', isopt, '/ ', opt(isopt)
     write(mystd,'(4X,a,i4,3X,2a)') 'dynamic truncation      /', iscut, '/ ', cut(iscut)
     !write(mystd,'(4X,a,i4,3X,2a)') 'data binning            /', isbin, '/ ', bin(isbin)
     write(mystd,'(4X,a,i4,3X,2a)') 'worm algorithm          /', iswor, '/ ', wor(iswor)
     write(mystd,'(4X,a,i4,3X,2a)') 'advanced basis          /', isort, '/ ', ort(isort)

     write(mystd,'(4X,a,i4,3X,2a)') 'fidelity susceptibility /', isobs, '/ ', trim(str_obs)
     write(mystd,'(4X,a,i4,3X,2a)') 'sp/ch susceptibility    /', issus, '/ ', trim(str_sus)
     write(mystd,'(4X,a,i4,3X,2a)') 'two-particle quantities /', isvrt, '/ ', trim(str_vrt)

     write(mystd,*)

!! body]

     return
  end subroutine ctqmc_print_control

!!
!! @sub ctqmc_print_runtime
!!
!! print the runtime information, including some physical observables and
!! statistic data, only for reference
!!
  subroutine ctqmc_print_runtime(iter, cstep)
     use constants, only : dp, one, half, mystd

     use control, only : nsweep
     use control, only : nmonte

     use context, only : cnegs, caves
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
! real(dp) dummy variables
     real(dp) :: raux

! about iteration number
     write(mystd,'(4X,a,i3,2(a,i10))') '>>> iter:', iter, ' sweep:', cstep, ' of ', nsweep

! about auxiliary physical observables
     raux = real(caves) / nmonte
     write(mystd,'(4X,a)')        'auxiliary system observables:'
     write(mystd,'(2(4X,a,f10.5))') 'etot :', paux(1) / raux, 'epot :', paux(2) / raux
     write(mystd,'(2(4X,a,f10.5))') 'ekin :', paux(3) / raux, '<Sz> :', paux(4) / raux
     write(mystd,'(2(4X,a,f10.5))') '<N1> :', paux(5) / raux, '<N2> :', paux(6) / raux
     write(mystd,'(2(4X,a,e10.3))') '<K2> :', paux(7) / raux, '<K3> :', paux(8) / raux
     write(mystd,'(1(4X,a,e10.3))') '<K4> :', paux(9) / raux

! about insert action
     if ( ins_t <= half ) ins_t = -one ! if insert is disable
     write(mystd,'(4X,a)')        'C_Z SPACE / insert kink statistics:'
     write(mystd,'(4X,a,3i10)')   'count:', int( ins_t ), int( ins_a ), int( ins_r )
     write(mystd,'(4X,a,3f10.5)') 'ratio:', one, ins_a / ins_t, ins_r / ins_t

! about remove action
     if ( rmv_t <= half ) rmv_t = -one ! if remove is disable
     write(mystd,'(4X,a)')        'C_Z SPACE / remove kink statistics:'
     write(mystd,'(4X,a,3i10)')   'count:', int( rmv_t ), int( rmv_a ), int( rmv_r )
     write(mystd,'(4X,a,3f10.5)') 'ratio:', one, rmv_a / rmv_t, rmv_r / rmv_t

! about lshift action
     if ( lsh_t <= half ) lsh_t = -one ! if lshift is disable
     write(mystd,'(4X,a)')        'C_Z SPACE / lshift kink statistics:'
     write(mystd,'(4X,a,3i10)')   'count:', int( lsh_t ), int( lsh_a ), int( lsh_r )
     write(mystd,'(4X,a,3f10.5)') 'ratio:', one, lsh_a / lsh_t, lsh_r / lsh_t

! about rshift action
     if ( rsh_t <= half ) rsh_t = -one ! if rshift is disable
     write(mystd,'(4X,a)')        'C_Z SPACE / rshift kink statistics:'
     write(mystd,'(4X,a,3i10)')   'count:', int( rsh_t ), int( rsh_a ), int( rsh_r )
     write(mystd,'(4X,a,3f10.5)') 'ratio:', one, rsh_a / rsh_t, rsh_r / rsh_t

! about reflip action
     if ( rfl_t <= half ) rfl_t = -one ! if reflip is disable
     write(mystd,'(4X,a)')        'C_Z SPACE / global flip statistics:'
     write(mystd,'(4X,a,3i10)')   'count:', int( rfl_t ), int( rfl_a ), int( rfl_r )
     write(mystd,'(4X,a,3f10.5)') 'ratio:', one, rfl_a / rfl_t, rfl_r / rfl_t

! about negative sign
     write(mystd,'(4X,a,i10)')    'negative sign counter:', cnegs
     write(mystd,'(4X,a,f10.5)')  'averaged sign sampler:', caves / real(cstep)

!! body]

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
     use control, only : isscf

     implicit none

! external arguments
! current iteration number
     integer, intent(in) :: iter

! according to the value of isscf, we can judge the self-consistent scheme
! of the current simulation
     if ( isscf /= 2 ) then
         write(mystd,'(2X,a,i3,a)') cname//' >>> SCF CYCLE:', iter, ' <<< ONESHOT'
     else
         write(mystd,'(2X,a,i3,a)') cname//' >>> SCF CYCLE:', iter, ' <<< LOOPING'
     endif ! back if ( isscf /= 2 ) block

!! body]

     return
  end subroutine ctqmc_print_it_info
