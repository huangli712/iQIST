!-------------------------------------------------------------------------
! project : daisy
! program : hfqmc_print_header
!           hfqmc_print_footer
!           hfqmc_print_summary
! source  : hfqmc_print.f90
! type    : subroutines
! author  : li huang (email:huangli712@yahoo.com.cn)
! history : 10/24/2008 by li huang
!           10/28/2008 by li huang
!           11/04/2008 by li huang
!           12/20/2008 by li huang
!           12/27/2008 by li huang
!           12/30/2008 by li huang
!           01/03/2009 by li huang
!           01/06/2009 by li huang
!           03/18/2009 by li huang
!           04/16/2009 by li huang
!           04/19/2009 by li huang
!           08/11/2009 by li huang
!           08/23/2009 by li huang
!           12/23/2009 by li huang
!           02/26/2010 by li huang
!           03/08/2010 by li huang
!           03/25/2010 by li huang
! purpose : provide printing infrastructure for Hirsch-Fye quantum Monte
!           Carlo (HFQMC) quantum impurity solver
! input   :
! output  :
! status  : very unstable
! comment :
!-------------------------------------------------------------------------

!>>> print the startup information for Hirsch-Fye quantum Monte Carlo
! quantum impurity solver plus dynamical mean field theory self-consistent
! engine
  subroutine hfqmc_print_header()
     use constants
     use control, only : nprocs

     implicit none

! string for current date and time
     character (len = 20) :: date_time_string

! obtain current date and time
     call s_time_builder(date_time_string)

     write(mystd,'(2X,a)') 'DAISY'
     write(mystd,'(2X,a)') '>>> A DMFT Engine With Hirsch-Fye Quantum Monte Carlo Impurity Solver'
     write(mystd,*)

     write(mystd,'(2X,a)') 'version: 2012.08.20T '//'(built at '//__TIME__//" "//__DATE__//')'
     write(mystd,'(2X,a)') 'develop: by li huang, CAEP & IOP'
     write(mystd,'(2X,a)') 'support: huangli712@yahoo.com.cn'
     write(mystd,'(2X,a)') 'license: GPL2 and later versions'
     write(mystd,*)

     write(mystd,'(2X,a)') 'DAISY >>> start running at '//date_time_string

# if defined (MPI)

     write(mystd,'(2X,a,i4)') 'DAISY >>> parallelism: Yes >>> processors:', nprocs

# else   /* MPI */

     write(mystd,'(2X,a,i4)') 'DAISY >>> parallelism: No  >>> processors:', 1

# endif  /* MPI */

     write(mystd,*)

     return
  end subroutine hfqmc_print_header

!>>> print the ending information for Hirsch-Fye quantum Monte Carlo
! quantum impurity solver plus dynamical mean field theory self-consistent
! engine
  subroutine hfqmc_print_footer()
     use constants

     implicit none

! string for current date and time
     character (len = 20) :: date_time_string

! used to record the time usage information
     real(dp) :: tot_time

! obtain time usage information
     call cpu_time(tot_time)

! obtain current date and time
     call s_time_builder(date_time_string)

     write(mystd,'(2X,a,f10.2,a)') 'DAISY >>> total time spent:', tot_time, 's'
     write(mystd,*)

     write(mystd,'(2X,a)') 'DAISY >>> I am tired and want to go to bed. Bye!'
     write(mystd,'(2X,a)') 'DAISY >>> happy ending at '//date_time_string

     return
  end subroutine hfqmc_print_footer

!>>> print the running parameters, only for reference
  subroutine hfqmc_print_summary()
     use constants
     use control

     implicit none

     write(mystd,'(2X,a)') 'DAISY >>> parameters list:'

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
