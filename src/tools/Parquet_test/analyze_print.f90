!-------------------------------------------------------------------------
! project : lavender
! program : ctqmc_print_header
!           ctqmc_print_footer
!           analyze_print_summary
!           ctqmc_print_runtime
!           analyze_print_error
!           ctqmc_print_exception
! source  : analyze_print.f90
! type    : subroutines
! author  : li huang (email:huangli712@yahoo.com.cn)
! history : 09/15/2009 by li huang
!           09/20/2009 by li huang
!           12/01/2009 by li huang
!           02/21/2010 by li huang
!	    10/30/2013 by ziyang meng
! purpose : provide printing infrastructure for the analyze code of the 
!	    hybridization expansion version continuous time quantum 
!           Monte Carlo (CTQMC) quantum impurity solver
! input   :
! output  :
! status  : very unstable
! comment :
!-------------------------------------------------------------------------

!>>> print the startup information for continuous time quantum Monte Carlo
! quantum impurity solver plus dynamical mean field theory self-consistent
! engine
  subroutine ctqmc_print_header()
     use constants
     use control, only : nprocs

     implicit none

     write(mystd,'(2X,a)') 'LAVENDER'
     write(mystd,'(2X,a)') '>>> A DMFT Engine With Continuous Time Quantum Monte Carlo Impurity Solver'
     write(mystd,*)

     write(mystd,'(2X,a)') 'version: 2011.08.18T            '
     write(mystd,'(2X,a)') 'develop: by li huang, CAEP & IOP'
     write(mystd,'(2X,a)') 'support: huangli712@yahoo.com.cn'
     write(mystd,'(2X,a)') 'license: GPL2 and later versions'
     write(mystd,*)

     write(mystd,'(2X,a)') 'LAVENDER >>> running'

# if defined (MPI)

     write(mystd,'(2X,a,i4)') 'LAVENDER >>> parallelism: Yes >>> processors:', nprocs

# else   /* MPI */

     write(mystd,'(2X,a,i4)') 'LAVENDER >>> parallelism: No  >>> processors:', 1

# endif  /* MPI */

     write(mystd,*)

     return
  end subroutine ctqmc_print_header

!>>> print the ending information for continuous time quantum Monte Carlo
! quantum impurity solver plus dynamical mean field theory self-consistent
! engine
  subroutine ctqmc_print_footer()
     use constants

     implicit none

! used to record the time information
     real(dp) :: tot_time

! obtain time information
     call cpu_time(tot_time)

     write(mystd,'(2X,a,f10.2,a)') 'LAVENDER >>> total time spent:', tot_time, 's'
     write(mystd,*)

     write(mystd,'(2X,a)') 'LAVENDER >>> I am tired and want to go to bed. Bye!'
     write(mystd,'(2X,a)') 'LAVENDER >>> ending'

     return
  end subroutine ctqmc_print_footer

!>>> print the running parameters, only for reference
  subroutine analyze_print_summary()
     use constants
     use control

     implicit none

     write(mystd,'(2X,a)') 'LAVENDER >>> parameters list:'

     write(mystd,'(2(4X,a,i10))')   'isscf :', isscf  , 'isbin :', isbin
     write(mystd,'(2(4X,a,i10))')   'issun :', issun  , 'isspn :', isspn
     write(mystd,'(2(4X,a,i10))')   'isort :', isort  , 'isvrt :', isvrt

     write(mystd,'(2(4X,a,i10))')   'lemax :', lemax  , 'legrd :', legrd
     write(mystd,'(2(4X,a,i10))')   'chmax :', chmax  , 'chgrd :', chgrd
     write(mystd,'(2(4X,a,i10))')   'mkink :', mkink  , 'mfreq :', mfreq
     write(mystd,'(2(4X,a,i10))')   'nband :', nband  , 'nspin :', nspin
     write(mystd,'(2(4X,a,i10))')   'norbs :', norbs  , 'ncfgs :', ncfgs
     write(mystd,'(2(4X,a,i10))')   'nzero :', nzero  , 'niter :', niter
     write(mystd,'(2(4X,a,i10))')   'nfreq :', nfreq  , 'ntime :', ntime
     write(mystd,'(2(4X,a,i10))')   'nffrq :', nffrq  , 'nbfrq :', nbfrq
     write(mystd,'(2(4X,a,i10))')   'npart :', npart  , 'nflip :', nflip
     write(mystd,'(2(4X,a,i10))')   'ktime :', ktime

     write(mystd,'(2(4X,a,i10))')   'ntherm:', ntherm , 'nsweep:', nsweep
     write(mystd,'(2(4X,a,i10))')   'nclean:', nclean , 'nwrite:', nwrite
     write(mystd,'(2(4X,a,i10))')   'nmonte:', nmonte , 'ncarlo:', ncarlo

     write(mystd,'(2(4X,a,f10.5))') 'U     :', U      , 'Uc    :', Uc
     write(mystd,'(2(4X,a,f10.5))') 'Js    :', Js     , 'Uv    :', Uv
     write(mystd,'(2(4X,a,f10.5))') 'Jp    :', Jp     , 'Jz    :', Jz
     write(mystd,'(2(4X,a,f10.5))') 'mune  :', mune   , 'beta  :', beta
     write(mystd,'(2(4X,a,f10.5),a)') 'part  :', part , 'temp  :', ev2k/beta, ' Kelvin'
     write(mystd,'(2(4X,a,f10.5))') 't1    :', t1     , 'lambda:', lambda

     write(mystd,*)

     return
  end subroutine analyze_print_summary



!>>> print the error information and STOP the program
  subroutine analyze_print_error(sub, msg)
     use constants

     implicit none

! external arguments
! subroutine name
     character(len=*), intent(in) :: sub

! error message
     character(len=*), intent(in) :: msg

! print error information
     write(mystd,'(2X,4a)') 'fatal error occurred in ', sub, ': ', msg

! TERMINATE THE PROGRAM
!-------------------------------------------------------------------------
     STOP
!-------------------------------------------------------------------------

     return
  end subroutine analyze_print_error

!>>> print normal runtime exceptional information, and continue
  subroutine ctqmc_print_exception(sub, msg)
     use constants

     implicit none

! external arguments
! subroutine name
     character(len=*), intent(in) :: sub

! exception message
     character(len=*), intent(in) :: msg

! print error information
     write(mystd,'(2X,4a)') 'runtime exception occurred in ', sub, ': ', msg

! CONTINUE/PAUSE THE PROGRAM
!-------------------------------------------------------------------------
     CONTINUE ! OR PAUSE
!-------------------------------------------------------------------------

     return
  end subroutine ctqmc_print_exception
