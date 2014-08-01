!!!-----------------------------------------------------------------------
!!! project : CSSL (Common Service Subroutines Library)
!!! program : s_zeros_i
!!!           s_zeros_d
!!!           s_zeros_z
!!!           s_ones_i
!!!           s_ones_d
!!!           s_ones_z
!!!           s_any_i
!!!           s_any_d
!!!           s_any_z
!!!           s_eye_i
!!!           s_eye_d
!!!           s_eye_z
!!!           s_identity_i
!!!           s_identity_d
!!!           s_identity_z
!!!           s_diag_i
!!!           s_diag_d
!!!           s_diag_z
!!!           s_inv_d
!!!           s_inv_z
!!!           s_det_d
!!!           s_det_z
!!!           s_trace_d
!!!           s_trace_z
!!!           s_eig_dg
!!!           s_eig_zg
!!!           s_eigvals_dg
!!!           s_eigvals_zg
!!!           s_eig_sy
!!!           s_eig_he
!!!           s_eigvals_sy
!!!           s_eigvals_he
!!! source  : s_matrix.f90
!!! type    : subroutines
!!! author  : li huang (email:huangli712@gmail.com)
!!! history : 07/10/2014 by li huang
!!!           07/26/2014 by li huang
!!!           08/01/2014 by li huang
!!! purpose : these subroutines are used to encapsulate some important and
!!!           frequently used linear algebra operations.
!!! status  : unstable
!!! comment :
!!!-----------------------------------------------------------------------

!!------------------------------------------------------------------------
!!>>> matrix manipulation: build zeros/ones/any matrix                 <<<
!!------------------------------------------------------------------------

!!>>> s_zeros_i: build an integer matrix with all elements are zero
  subroutine s_zeros_i(A)
     implicit none

! external arguments
! input/output matrix
     integer, intent(out) :: A(:,:)

     A = 0

     return
  end subroutine s_zeros_i

!!>>> s_zeros_d: build a real(dp) matrix with all elements are zero
  subroutine s_zeros_d(A)
     use constants, only : dp, zero

     implicit none

! external arguments
! input/output matrix
     real(dp), intent(out) :: A(:,:)

     A = zero

     return
  end subroutine s_zeros_d

!!>>> s_zeros_z: build a complex(dp) matrix with all elements are zero
  subroutine s_zeros_z(A)
     use constants, only : dp, czero

     implicit none

! external arguments
! input/output matrix
     complex(dp), intent(out) :: A(:,:)

     A = czero

     return
  end subroutine s_zeros_z

!!>>> s_ones_i: build an integer matrix with all elements are one
  subroutine s_ones_i(A)
     implicit none

! external arguments
! input/output matrix
     integer, intent(out) :: A(:,:)

     A = 1

     return
  end subroutine s_ones_i

!!>>> s_ones_d: build a real(dp) matrix with all elements are one
  subroutine s_ones_d(A)
     use constants, only : dp, one

     implicit none

! external arguments
! input/output matrix
     real(dp), intent(out) :: A(:,:)

     A = one

     return
  end subroutine s_ones_d

!!>>> s_ones_z: build a complex(dp) matrix with all elements are one
  subroutine s_ones_z(A)
     use constants, only : dp, cone

     implicit none

! external arguments
! input/output matrix
     complex(dp), intent(out) :: A(:,:)

     A = cone

     return
  end subroutine s_ones_z

!!>>> s_any_i: build an integer matrix with all elements are given by i
  subroutine s_any_i(A, i)
     implicit none

! external arguments
! input/output matrix
     integer, intent(out) :: A(:,:)

! value of matrix element
     integer, intent(in)  :: i

     A = i

     return
  end subroutine s_any_i

!!>>> s_any_d: build a real(dp) matrix with all elements are given by d
  subroutine s_any_d(A, d)
     use constants, only : dp

     implicit none

! external arguments
! input/output matrix
     real(dp), intent(out) :: A(:,:)

! value of matrix element
     real(dp), intent(in)  :: d

     A = d

     return
  end subroutine s_any_d

!!>>> s_any_z: build a complex(dp) matrix with all elements are given by z
  subroutine s_any_z(A, z)
     use constants, only : dp

     implicit none

! external arguments
! input/output matrix
     complex(dp), intent(out) :: A(:,:)

! value of matrix element
     complex(dp), intent(in)  :: z

     A = z

     return
  end subroutine s_any_z
!!------------------------------------------------------------------------
!!>>> matrix manipulation: build diagonal matrix                       <<<
!!------------------------------------------------------------------------



!!>>> s_inv_dmat: invert real(dp) matrix using lapack subroutines
  subroutine s_inv_dmat(ndim, dmat)
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
     call DGETRF(ndim, ndim, dmat, ndim, ipiv, ierror)
     if ( ierror /= 0 ) then
         call s_print_error('s_inv_dmat','error in lapack subroutine dgetrf')
     endif

! computes the inverse of an LU-factored general matrix, need lapack
! package, dgetri subroutine
     call DGETRI(ndim, dmat, ndim, ipiv, work, ndim, ierror)
     if ( ierror /= 0 ) then
         call s_print_error('s_inv_dmat','error in lapack subroutine dgetri')
     endif

     return
  end subroutine s_inv_dmat

!!>>> s_inv_zmat: invert complex(dp) matrix using lapack subroutines
  subroutine s_inv_zmat(ndim, zmat)
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
     call ZGETRF(ndim, ndim, zmat, ndim, ipiv, ierror)
     if ( ierror /= 0 ) then
         call s_print_error('s_inv_zmat','error in lapack subroutine zgetrf')
     endif

! computes the inverse of an LU-factored general matrix, need lapack
! package, zgetri subroutine
     call ZGETRI(ndim, zmat, ndim, ipiv, work, ndim, ierror)
     if ( ierror /= 0 ) then
         call s_print_error('s_inv_zmat','error in lapack subroutine zgetri')
     endif

     return
  end subroutine s_inv_zmat

!!>>> s_det_dmat: calculate the determinant of a real(dp) matrix
  subroutine s_det_dmat(ndim, dmat, ddet)
     use constants, only : dp, one, cone

     implicit none

! external arguments
! dimension of dmat matrix
     integer, intent(in) :: ndim

! determinant of dmat matrix
     real(dp), intent(out) :: ddet

! object matrix, on entry, it contains the original matrix, on exit,
! it is destroyed and replaced with the L and U matrix
     real(dp), intent(inout) :: dmat(ndim,ndim)

! local variables
! loop index
     integer  :: i

! error flag
     integer  :: ierror

! size of working array work
     integer  :: lwork

! working arrays for lapack subroutines: dgetrf
     integer  :: ipiv(ndim)

! working arrays for lapack subroutines: dgeev
     real(dp) :: work(4*ndim)

! real and imaginary parts of the computed eigenvalues
     real(dp) :: wi(ndim)
     real(dp) :: wr(ndim)

! left and right eigenvectors
     real(dp) :: vl(ndim,ndim)
     real(dp) :: vr(ndim,ndim)

! dummy arrays, used to save dmat
     real(dp) :: amat(ndim,ndim)

! used to calculate determinant
     complex(dp) :: cres

! setup lwork
     lwork = 4 * ndim

! copy dmat to amat at first
     amat = dmat

!-------------------------------------------------------------------------
! method A: preferred method
!-------------------------------------------------------------------------
! computes the LU factorization of a general m-by-n matrix, need lapack
! package, dgetrf subroutine
     call DGETRF(ndim, ndim, dmat, ndim, ipiv, ierror)
     if ( ierror /= 0 ) then
         call s_print_exception('s_det_dmat','error in lapack subroutine dgetrf')
     endif

! calculate determinant
     ddet = one
     do i=1,ndim
         if ( ipiv(i) == i ) then
             ddet = ddet * ( +dmat(i,i) )
         else
             ddet = ddet * ( -dmat(i,i) )
         endif
     enddo ! over i={1,ndim} loop

! everything is ok!
     if ( ierror == 0 ) RETURN

!-------------------------------------------------------------------------
! method B: as a backup
!-------------------------------------------------------------------------
! diagonalize amat to obtain its eigenvalues: wr and wi
     call DGEEV('N', 'N', ndim, amat, ndim, wr, wi, vl, ndim, vr, ndim, work, lwork, ierror)
     if ( ierror /= 0 ) then
         call s_print_error('s_det_dmat','error in lapack subroutine dgeev')
     endif

! evaluate the final determinant
     cres = cone
     do i=1,ndim
         cres = cres * dcmplx( wr(i), wi(i) )
     enddo ! over i={1,ndim} loop
     ddet = cres

     return
  end subroutine s_det_dmat

!!>>> s_det_zmat: calculate the determinant of a complex(dp) matrix
  subroutine s_det_zmat(ndim, zmat, zdet)
     use constants, only : dp, cone

     implicit none

! external arguments
! dimension of zmat matrix
     integer, intent(in) :: ndim

! determinant of zmat matrix
     real(dp), intent(out) :: zdet

! object matrix, on entry, it contains the original matrix, on exit,
! it is destroyed and replaced with the L and U matrix
     real(dp), intent(inout) :: zmat(ndim,ndim)

! local variables
! loop index
     integer :: i

! error flag
     integer :: ierror

! working arrays for lapack subroutines
     integer :: ipiv(ndim)

! computes the LU factorization of a general m-by-n matrix, need lapack
! package, zgetrf subroutine
     call ZGETRF(ndim, ndim, zmat, ndim, ipiv, ierror)
     if ( ierror /= 0 ) then
         call s_print_error('s_det_zmat','error in lapack subroutine zgetrf')
     endif

! calculate determinant
     zdet = cone
     do i=1,ndim
         if ( ipiv(i) == i ) then
             zdet = zdet * ( +zmat(i,i) )
         else
             zdet = zdet * ( -zmat(i,i) )
         endif
     enddo ! over i={1,ndim} loop

     return
  end subroutine s_det_zmat
