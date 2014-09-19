!-------------------------------------------------------------------------
! project : lavender_analyze
! program : analyze_dmat_inv
!           analyze_zmat_inv
!	    analyze_dmat_prd
!	    analyze_zmat_prd
! source  : analyze_util.f90
! type    : functions & subroutines
! author  : li huang (email:huangli712@yahoo.com.cn)
! history : 10/01/2008 by li huang
!           02/08/2009 by li huang
!           09/23/2009 by li huang
!           09/26/2009 by li huang
!           11/17/2009 by li huang
!           11/21/2009 by li huang
!           12/18/2009 by li huang
!           12/22/2009 by li huang
!           12/29/2009 by li huang
!           01/12/2010 by li huang
!           02/27/2010 by li huang
!           06/08/2010 by li huang
!           06/22/2010 by li huang
!           10/30/2013 by ziyang meng
! purpose : to provide utility functions and subroutines for the analyze code
!           of the hybridization expansion version continuous time quantum Monte 
!           Carlo (CTQMC) quantum impurity solver
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
  subroutine analyze_dmat_inv(ndim, dmat)
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
         call analyze_print_error('analyze_dmat_inv','error in lapack subroutine dgetrf')
     endif

! computes the inverse of an LU-factored general matrix, need lapack
! package, dgetri subroutine
     call dgetri(ndim, dmat, ndim, ipiv, work, ndim, ierror)
     if ( ierror /= 0 ) then
         call analyze_print_error('analyze_dmat_inv','error in lapack subroutine dgetri')
     endif

     return
  end subroutine analyze_dmat_inv

!>>> invert complex(dp) matrix using lapack subroutines
  subroutine analyze_zmat_inv(ndim, zmat)
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
         call analyze_print_error('analyze_zmat_inv','error in lapack subroutine zgetrf')
     endif

! computes the inverse of an LU-factored general matrix, need lapack
! package, zgetri subroutine
     call zgetri(ndim, zmat, ndim, ipiv, work, ndim, ierror)
     if ( ierror /= 0 ) then
         call analyze_print_error('analyze_zmat_inv','error in lapack subroutine zgetri')
     endif

     return
  end subroutine analyze_zmat_inv
  
  !>>> real(dp) matrix multiplication using lapack subroutines
  subroutine analyze_dmat_prd(ndim, dmat1, dmat2, dmat3)
     use constants

     implicit none

! external arguments
! dimension of dmat matrix
     integer, intent(in) :: ndim

! object matrix, it contains the original matrix
     real(dp), intent(in) :: dmat1(ndim,ndim), dmat2(ndim,ndim)
! object matrix, it contains the multiplication result
     real(dp), intent(out) :: dmat3(ndim,ndim)
     
! local variables
! error flag
!     integer  :: ierror

! computes the matrix multiplication, need lapack package, dgemm subroutine
     call dgemm('N', 'N', ndim, ndim, ndim, one, dmat1, ndim, dmat2, ndim, zero, dmat3, ndim)

     return
  end subroutine analyze_dmat_prd
  
  !>>> complex(dp) matrix multiplication using lapack subroutines
  subroutine analyze_zmat_prd(ndim, zmat1, zmat2, zmat3)
     use constants

     implicit none

! external arguments
! dimension of zmat matrix
     integer, intent(in) :: ndim

! object matrix, it contains the original matrix
     real(dp), intent(in) :: zmat1(ndim,ndim), zmat2(ndim,ndim)
! object matrix, it contains the multiplication result
     real(dp), intent(out) :: zmat3(ndim,ndim)
     
! local variables
! error flag
!     integer     :: ierror

! computes the matrix multiplication, need lapack package, zgemm subroutine
     call zgemm('N', 'N', ndim, ndim, ndim, cone, zmat1, ndim, zmat2, ndim, czero, zmat3, ndim)
     
     return
  end subroutine analyze_zmat_prd
  

