!-------------------------------------------------------------------------!
! project : strawberry
! program : rtgw_print_header
!           rtgw_print_footer
!           rtgw_print_summary
! history : Dec 21, 2010
! author  : xidai and duliang (email: duleung@gmail.com)
! purpose :
! comment :
!-------------------------------------------------------------------------!
!>>> print the heading information
  subroutine atomic_print_header()
     use constants
     use control

     implicit none

     write(mystd,'(2X,A)') ''
     write(mystd,'(2X,A)') '>>> atomic solver for rotational invariant gutzwiller method'
     write(mystd,*)

     write(mystd,'(2X,A)') 'version: '
     write(mystd,'(2X,A)') 'develop: '
     write(mystd,'(2X,A)') 'support: '
     write(mystd,*)

     write(mystd,'(2X,A)') 'strawberry >>> running'

# if defined (MPI)

     write(mystd,'(2X,A,I3)') 'strawberry >>> parallel: Y >>> nprocs:', nprocs

# else   /* MPI */

     write(mystd,'(2X,A,I3)') 'strawberry >>> parallel: N >>> nprocs:', 1

# endif  /* MPI */

     write(mystd,*)

     return
  end subroutine atomic_print_header

!>>> print the ending information
  subroutine atomic_print_footer()
     use constants
     use control

     implicit none

! used to record the time information
     real(dp) :: time

! obtain time information
     call cpu_time(time)

     write(mystd,'(2X,A,F10.2,A)') 'strawberry >>> total time spent:', time, 's'
     write(mystd,*)

     write(mystd,'(2X,A)') 'strawberry >>> hope you good luck. Bye!'
     write(mystd,'(2X,A)') 'strawberry >>> ending'

     return
  end subroutine atomic_print_footer

!>>> print the running parameters
  subroutine atomic_print_summary()
     use constants
     use control

     implicit none

     write(mystd, '(2X, A)') 'clemalis >>> parameters list:'

     write(mystd, '(2(4X, A, I10))')   'nband :', nband , 'nspin :', nspin
     write(mystd, '(2(4X, A, I10))')   'norbs :', norbs , 'ntots :', ntots
     write(mystd, '(2(4X, A, I10))')   'ncfgs :', ncfgs , 'ncfgs :', ncfgs

     write(mystd, '(2(4X, A, F10.5))') 'Uc    :', Uc    , 'Uv    :', Uv
     write(mystd, '(2(4X, A, F10.5))') 'Jz    :', Jz    , 'Js    :', Js
     write(mystd, '(2(4X, A, F10.5))') 'Jp    :', Jp    , 'J     :', Jz

     write(mystd, *)

     return
  end subroutine atomic_print_summary
