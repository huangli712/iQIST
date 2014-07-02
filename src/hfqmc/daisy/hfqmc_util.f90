!-------------------------------------------------------------------------
! project : daisy
! program : hfqmc_dmat_inv
!           hfqmc_zmat_inv
!           hfqmc_time_builder
!           hfqmc_time_analyzer
! source  : hfqmc_util.f90
! type    : functions & subroutines
! author  : li huang (email:huangli712@yahoo.com.cn)
! history : 12/24/2009 by li huang
!           02/26/2010 by li huang
!           03/04/2010 by li huang
!           08/25/2010 by li huang
! purpose : to provide utility functions and subroutines for Hirsch-Fye
!           quantum Monte Carlo (HFQMC) quantum impurity solver
! input   :
! output  :
! status  : unstable
! comment :
!-------------------------------------------------------------------------

# define prefix '>>> used time: '

# define iolst1 prefix , mday, ' d ', mhou, ' h ', mmin, ' m in this iteration.'
# define iolst2 prefix , mhou, ' h ', mmin, ' m in this iteration.'
# define iolst3 prefix , mmin, ' m ', msec, ' s in this iteration.'
# define iolst4 prefix , msec, ' s in this iteration.'

# define iolst5 prefix , nday, ' d ', nhou, ' h ', nmin, ' m in total iteration.'
# define iolst6 prefix , nhou, ' h ', nmin, ' m in total iteration.'
# define iolst7 prefix , nmin, ' m ', nsec, ' s in total iteration.'
# define iolst8 prefix , nsec, ' s in total iteration.'

!>>> invert real(dp) matrix using lapack subroutines
  subroutine hfqmc_dmat_inv(ndim, dmat)
     use constants, only : dp

     implicit none

! external arguments
! dimension of dmat matrix
     integer, intent(in) :: ndim

! object matrix, on entry, it contains the original matrix, on exit,
! it is destroyed and replaced with the inversed matrix
     real(dp), intent(inout) :: dmat(ndim,ndim)

! local variables
! error flag
     integer  :: ierror

! working arrays for lapack subroutines
     integer  :: ipiv(ndim)
     real(dp) :: work(ndim)

! computes the LU factorization of a general m-by-n matrix, need lapack
! package, dgetrf subroutine
     call dgetrf(ndim, ndim, dmat, ndim, ipiv, ierror)
     if ( ierror /= 0 ) then
         call hfqmc_print_error('hfqmc_dmat_inv','error in lapack subroutine dgetrf')
     endif

! computes the inverse of an LU-factored general matrix, need lapack
! package, dgetri subroutine
     call dgetri(ndim, dmat, ndim, ipiv, work, ndim, ierror)
     if ( ierror /= 0 ) then
         call hfqmc_print_error('hfqmc_dmat_inv','error in lapack subroutine dgetri')
     endif

     return
  end subroutine hfqmc_dmat_inv

!>>> invert complex(dp) matrix using lapack subroutines
  subroutine hfqmc_zmat_inv(ndim, zmat)
     use constants, only : dp

     implicit none

! external arguments
! dimension of zmat matrix
     integer, intent(in) :: ndim

! object matrix, on entry, it contains the original matrix, on exit,
! it is destroyed and replaced with the inversed matrix
     complex(dp), intent(inout) :: zmat(ndim,ndim)

! local variables
! error flag
     integer     :: ierror

! working arrays for lapack subroutines
     integer     :: ipiv(ndim)
     complex(dp) :: work(ndim)

! computes the LU factorization of a general m-by-n matrix, need lapack
! package, zgetrf subroutine
     call zgetrf(ndim, ndim, zmat, ndim, ipiv, ierror)
     if ( ierror /= 0 ) then
         call hfqmc_print_error('hfqmc_zmat_inv','error in lapack subroutine zgetrf')
     endif

! computes the inverse of an LU-factored general matrix, need lapack
! package, zgetri subroutine
     call zgetri(ndim, zmat, ndim, ipiv, work, ndim, ierror)
     if ( ierror /= 0 ) then
         call hfqmc_print_error('hfqmc_zmat_inv','error in lapack subroutine zgetri')
     endif

     return
  end subroutine hfqmc_zmat_inv

!>>> returns a string containing date and time in human-readable format
  subroutine hfqmc_time_builder(date_time_string)
     implicit none

! external arguments
! output date and time
     character (len = 20), intent(out) :: date_time_string

! local variables
! used to extract data from a standard fortran call: date_and_time()
     integer :: date_time(8)

! string for current date
     character (len = 12) :: cdate

! string for current time
     character (len = 08) :: ctime

! month array
     character (len = 03) :: months(12)

! init the month array
     months( 1) = 'Jan'; months( 2) = 'Feb'; months( 3) = 'Mar'
     months( 4) = 'Apr'; months( 5) = 'May'; months( 6) = 'Jun'
     months( 7) = 'Jul'; months( 8) = 'Aug'; months( 9) = 'Sep'
     months(10) = 'Oct'; months(11) = 'Nov'; months(12) = 'Dec'

! obtain current date and time
     call date_and_time(values = date_time)

! convert date and time from integer to string
     write (cdate,'(1X,a3,1X,i2,1X,i4)') months(date_time(2)), date_time(3), date_time(1)
     write (ctime, '(i2,":",i2,":",i2)') date_time(5), date_time(6), date_time(7)

! build final output string by concating them
     date_time_string = ctime // cdate

     return
  end subroutine hfqmc_time_builder

!>>> used to print the iteration timing information about Hirsch-Fye
! quantum Monte Carlo quantum impurity solver.
  subroutine hfqmc_time_analyzer(time_iter, time_niter)
     use constants

     implicit none

! external arguments
! time used in this iteration
     real(dp), intent(in) :: time_iter

! time used in total iteration
     real(dp), intent(in) :: time_niter

! local variables
! number of days
     integer  :: mday, nday

! number of hours
     integer  :: mhou, nhou

! number of minutes
     integer  :: mmin, nmin

! number of seconds
     real(dp) :: msec, nsec

     mday = time_iter / 86400
     msec = time_iter - 86400 * mday
     mhou = msec / 3600
     msec = msec - 3600 * mhou
     mmin = msec / 60
     msec = msec - 60 * mmin

     nday = time_niter / 86400
     nsec = time_niter - 86400 * nday
     nhou = nsec / 3600
     nsec = nsec - 3600 * nhou
     nmin = nsec / 60
     nsec = nsec - 60 * nmin

     if      ( mday > 0 ) then
         write(mystd,'(4X,3(a,i2),a)')     iolst1

     else if ( mhou > 0 ) then
         write(mystd,'(4X,2(a,i2),a)')     iolst2

     else if ( mmin > 0 ) then
         write(mystd,'(4X,a,i2,a,f5.2,a)') iolst3

     else
         write(mystd,'(4X,a,f5.2,a)')      iolst4

     endif ! back if ( mday > 0 ) block

     if      ( nday > 0 ) then
         write(mystd,'(4X,3(a,i2),a)')     iolst5

     else if ( nhou > 0 ) then
         write(mystd,'(4X,2(a,i2),a)')     iolst6

     else if ( nmin > 0 ) then
         write(mystd,'(4X,a,i2,a,f5.2,a)') iolst7

     else
         write(mystd,'(4X,a,f5.2,a)')      iolst8

     endif ! back if ( nday > 0 ) block

     return
  end subroutine hfqmc_time_analyzer
