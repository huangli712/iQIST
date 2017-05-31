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
!!!           s_trace_d
!!!           s_trace_z
!!!           s_det_d
!!!           s_det_z
!!!           s_inv_d
!!!           s_inv_z
!!!           s_eig_dg
!!!           s_eig_zg
!!!           s_eigvals_dg
!!!           s_eigvals_zg
!!!           s_eig_sy
!!!           s_eig_he
!!!           s_eigvals_sy
!!!           s_eigvals_he
!!!           s_solve_dg
!!!           s_solve_zg
!!!           s_solve_sy
!!!           s_solve_he
!!!           s_svd_dg
!!!           s_svd_zg
!!! source  : s_matrix.f90
!!! type    : subroutines
!!! author  : li huang (email:lihuang.dmft@gmail.com)
!!! history : 07/10/2014 by li huang (created)
!!!           05/31/2017 by li huang (last modified)
!!! purpose : these subroutines are used to encapsulate some important and
!!!           frequently used linear algebra operations.
!!! status  : unstable
!!! comment :
!!!-----------------------------------------------------------------------

!!
!!
!! Introduction
!! ============
!!
!! 1. build constants (0) matrix
!! -----------------------------
!!
!! subroutine s_zeros_i(...)
!! subroutine s_zeros_d(...)
!! subroutine s_zeros_z(...)
!!
!! 2. build constants (1) matrix
!! -----------------------------
!!
!! subroutine s_ones_i(...)
!! subroutine s_ones_d(...)
!! subroutine s_ones_z(...)
!!
!! 3. build constants (any values) matrix
!! --------------------------------------
!!
!! subroutine s_any_i(...)
!! subroutine s_any_d(...)
!! subroutine s_any_z(...)
!!
!! 4. build diagonal matrix
!! ------------------------
!!
!! subroutine s_eye_i(...)
!! subroutine s_eye_d(...)
!! subroutine s_eye_z(...)
!!
!! 5. build identity matrix
!! ------------------------
!!
!! subroutine s_identity_i(...)
!! subroutine s_identity_d(...)
!! subroutine s_identity_z(...)
!!
!! 6. build diagonal matrix from vector
!! ------------------------------------
!!
!! subroutine s_diag_i(...)
!! subroutine s_diag_d(...)
!! subroutine s_diag_z(...)
!!
!! 7. calculate trace for matrix
!! -----------------------------
!!
!! subroutine s_trace_d(...)
!! subroutine s_trace_z(...)
!!
!! 8. calculate determinant for matrix
!! -----------------------------------
!!
!! subroutine s_det_d(...)
!! subroutine s_det_z(...)
!!
!! 9. calculate matrix inversion
!! -----------------------------
!!
!! subroutine s_inv_d(...)
!! subroutine s_inv_z(...)
!!
!! 10. general eigensystem problem
!! -------------------------------
!!
!! subroutine s_eig_dg(...)
!! subroutine s_eig_zg(...)
!! subroutine s_eigvals_dg(...)
!! subroutine s_eigvals_zg(...)
!!
!! 11. symmetric eigensystem problem
!! ---------------------------------
!!
!! subroutine s_eig_sy(...)
!! subroutine s_eig_he(...)
!! subroutine s_eigvals_sy(...)
!! subroutine s_eigvals_he(...)
!!
!! 12. linear equation solver
!! --------------------------
!!
!! subroutine s_solve_dg(...)
!! subroutine s_solve_zg(...)
!! subroutine s_solve_sy(...)
!! subroutine s_solve_he(...)
!!
!! 13. general singular value decomposition
!! ----------------------------------------
!!
!! subroutine s_svd_dg(...)
!! subroutine s_svd_zg(...)
!!
!! Note: _i means integer version, _d real(dp) version, and _z complex(dp)
!! version. _dg means real(dp) general version, _zg complex(dp) general
!! version, _sy real(dp) symmetric version, _he complex(dp) Hermitian version.
!!
!!

!!========================================================================
!!>>> matrix construction: build zeros/ones/any matrix                 <<<
!!========================================================================

!!
!! @sub s_zeros_i
!!
!! build an integer matrix with all elements are zero
!!
  subroutine s_zeros_i(n, A)
     implicit none

! external arguments
! size of matrix
     integer, intent(in)  :: n

! input/output matrix
     integer, intent(out) :: A(n,n)

     A = 0

     return
  end subroutine s_zeros_i

!!
!! @sub s_zeros_d
!!
!! build a real(dp) matrix with all elements are zero
!!
  subroutine s_zeros_d(n, A)
     use constants, only : dp
     use constants, only : zero

     implicit none

! external arguments
! size of matrix
     integer, intent(in)   :: n

! input/output matrix
     real(dp), intent(out) :: A(n,n)

     A = zero

     return
  end subroutine s_zeros_d

!!
!! @sub s_zeros_z
!!
!! build a complex(dp) matrix with all elements are zero
!!
  subroutine s_zeros_z(n, A)
     use constants, only : dp
     use constants, only : czero

     implicit none

! external arguments
! size of matrix
     integer, intent(in)      :: n

! input/output matrix
     complex(dp), intent(out) :: A(n,n)

     A = czero

     return
  end subroutine s_zeros_z

!!
!! @sub s_ones_i
!!
!! build an integer matrix with all elements are one
!!
  subroutine s_ones_i(n, A)
     implicit none

! external arguments
! size of matrix
     integer, intent(in)  :: n

! input/output matrix
     integer, intent(out) :: A(n,n)

     A = 1

     return
  end subroutine s_ones_i

!!
!! @sub s_ones_d
!!
!! build a real(dp) matrix with all elements are one
!!
  subroutine s_ones_d(n, A)
     use constants, only : dp
     use constants, only : one

     implicit none

! external arguments
! size of matrix
     integer, intent(in)   :: n

! input/output matrix
     real(dp), intent(out) :: A(n,n)

     A = one

     return
  end subroutine s_ones_d

!!
!! @sub s_ones_z
!!
!! build a complex(dp) matrix with all elements are one
!!
  subroutine s_ones_z(n, A)
     use constants, only : dp
     use constants, only : cone

     implicit none

! external arguments
! size of matrix
     integer, intent(in)      :: n

! input/output matrix
     complex(dp), intent(out) :: A(n,n)

     A = cone

     return
  end subroutine s_ones_z

!!
!! @sub s_any_i
!!
!! build an integer matrix with all elements are given by i
!!
  subroutine s_any_i(n, i, A)
     implicit none

! external arguments
! size of matrix
     integer, intent(in)  :: n

! value of matrix element
     integer, intent(in)  :: i

! input/output matrix
     integer, intent(out) :: A(n,n)

     A = i

     return
  end subroutine s_any_i

!!
!! @sub s_any_d
!!
!! build a real(dp) matrix with all elements are given by d
!!
  subroutine s_any_d(n, d, A)
     use constants, only : dp

     implicit none

! external arguments
! size of matrix
     integer, intent(in)   :: n

! value of matrix element
     real(dp), intent(in)  :: d

! input/output matrix
     real(dp), intent(out) :: A(n,n)

     A = d

     return
  end subroutine s_any_d

!!
!! @sub s_any_z
!!
!! build a complex(dp) matrix with all elements are given by z
!!
  subroutine s_any_z(n, z, A)
     use constants, only : dp

     implicit none

! external arguments
! size of matrix
     integer, intent(in)      :: n

! value of matrix element
     complex(dp), intent(in)  :: z

! input/output matrix
     complex(dp), intent(out) :: A(n,n)

     A = z

     return
  end subroutine s_any_z

!!========================================================================
!!>>> matrix construction: build diagonal matrix                       <<<
!!========================================================================

!!
!! @sub s_eye_i
!!
!! build integer matrix with ones on the diagonal and zeros elsewhere
!!
  subroutine s_eye_i(n, k, A)
     implicit none

! external arguments
! size of matrix
     integer, intent(in)  :: n

! index of the diagonal: 0 refers to the main diagonal, a positive value
! refers to an upper diagonal, and a negative value to a lower diagonal.
     integer, intent(in)  :: k

! input/output matrix
     integer, intent(out) :: A(n,n)

! local variables
! loop index
     integer :: i

     A = 0
     do i=1,n
         if ( i - k < 1 .or. i - k > n ) CYCLE
         A(i,i-k) = 1
     enddo ! over i={1,n} loop

     return
  end subroutine s_eye_i

!!
!! @sub s_eye_d
!!
!! build real(dp) matrix with ones on the diagonal and zeros elsewhere
!!
  subroutine s_eye_d(n, k, A)
     use constants, only : dp
     use constants, only : zero, one

     implicit none

! external arguments
! size of matrix
     integer, intent(in)   :: n

! index of the diagonal: 0 refers to the main diagonal, a positive value
! refers to an upper diagonal, and a negative value to a lower diagonal.
     integer, intent(in)   :: k

! input/output matrix
     real(dp), intent(out) :: A(n,n)

! local variables
! loop index
     integer :: i

     A = zero
     do i=1,n
         if ( i - k < 1 .or. i - k > n ) CYCLE
         A(i,i-k) = one
     enddo ! over i={1,n} loop

     return
  end subroutine s_eye_d

!!
!! @sub s_eye_z
!!
!! build complex(dp) matrix with ones on the diagonal and zeros elsewhere
!!
  subroutine s_eye_z(n, k, A)
     use constants, only : dp
     use constants, only : czero, cone

     implicit none

! external arguments
! size of matrix
     integer, intent(in)      :: n

! index of the diagonal: 0 refers to the main diagonal, a positive value
! refers to an upper diagonal, and a negative value to a lower diagonal.
     integer, intent(in)      :: k

! input/output matrix
     complex(dp), intent(out) :: A(n,n)

! local variables
! loop index
     integer :: i

     A = czero
     do i=1,n
         if ( i - k < 1 .or. i - k > n ) CYCLE
         A(i,i-k) = cone
     enddo ! over i={1,n} loop

     return
  end subroutine s_eye_z

!!
!! @sub s_identity_i
!!
!! build integer identity matrix
!!
  subroutine s_identity_i(n, A)
     implicit none

! external arguments
! size of matrix
     integer, intent(in)  :: n

! input/output matrix
     integer, intent(out) :: A(n,n)

! local variables
! loop index
     integer :: i

     A = 0
     do i=1,n
         A(i,i) = 1
     enddo ! over i={1,n} loop

     return
  end subroutine s_identity_i

!!
!! @sub s_identity_d
!!
!! build real(dp) identity matrix
!!
  subroutine s_identity_d(n, A)
     use constants, only : dp
     use constants, only : zero, one

     implicit none

! external arguments
! size of matrix
     integer, intent(in)   :: n

! input/output matrix
     real(dp), intent(out) :: A(n,n)

! local variables
! loop index
     integer :: i

     A = zero
     do i=1,n
         A(i,i) = one
     enddo ! over i={1,n} loop

     return
  end subroutine s_identity_d

!!
!! @sub s_identity_z
!!
!! build complex(dp) identity matrix
!!
  subroutine s_identity_z(n, A)
     use constants, only : dp
     use constants, only : czero, cone

     implicit none

! external arguments
! size of matrix
     integer, intent(in)      :: n

! input/output matrix
     complex(dp), intent(out) :: A(n,n)

! local variables
! loop index
     integer :: i

     A = czero
     do i=1,n
         A(i,i) = cone
     enddo ! over i={1,n} loop

     return
  end subroutine s_identity_z

!!
!! @sub s_diag_i
!!
!! build integer diagonal matrix from a vector
!!
  subroutine s_diag_i(n, v, A)
     implicit none

! external arguments
! size of matrix
     integer, intent(in)  :: n

! input integer vector
     integer, intent(in)  :: v(n)

! output integer diagonal matrix
     integer, intent(out) :: A(n,n)

! local variables
! loop index
     integer :: i

     A = 0
     do i=1,n
         A(i,i) = v(i)
     enddo ! over i={1,n} loop

     return
  end subroutine s_diag_i

!!
!! @sub s_diag_d
!!
!! build real(dp) diagonal matrix from a vector
!!
  subroutine s_diag_d(n, v, A)
     use constants, only : dp
     use constants, only : zero

     implicit none

! external arguments
! size of matrix
     integer, intent(in)   :: n

! input real(dp) vector
     real(dp), intent(in)  :: v(n)

! output real(dp) diagonal matrix
     real(dp), intent(out) :: A(n,n)

! local variables
! loop index
     integer :: i

     A = zero
     do i=1,n
         A(i,i) = v(i)
     enddo ! over i={1,n} loop

     return
  end subroutine s_diag_d

!!
!! @sub s_diag_z
!!
!! build complex(dp) diagonal matrix from a vector
!!
  subroutine s_diag_z(n, v, A)
     use constants, only : dp
     use constants, only : czero

     implicit none

! external arguments
! size of matrix
     integer, intent(in)      :: n

! input complex(dp) vector
     complex(dp), intent(in)  :: v(n)

! output complex(dp) diagonal matrix
     complex(dp), intent(out) :: A(n,n)

! local variables
! loop index
     integer :: i

     A = czero
     do i=1,n
         A(i,i) = v(i)
     enddo ! over i={1,n} loop

     return
  end subroutine s_diag_z

!!========================================================================
!!>>> matrix query: return matrix's trace or determinant               <<<
!!========================================================================

!!
!! @sub s_trace_d
!!
!! return trace for a real(dp) array
!!
  subroutine s_trace_d(n, A, tr)
     use constants, only : dp
     use constants, only : zero

     implicit none

! external arguments
! size of matrix
     integer, intent(in)   :: n

! output matrix's trace
     real(dp), intent(out) :: tr

! input real(dp) matrix
     real(dp), intent(in)  :: A(n,n)

! local variables
! loop index
     integer :: i

     tr = zero
     do i=1,n
         tr = tr + A(i,i)
     enddo ! over i={1,n} loop

     return
  end subroutine s_trace_d

!!
!! @sub s_trace_z
!!
!! return trace for a complex(dp) array
!!
  subroutine s_trace_z(n, A, tr)
     use constants, only : dp
     use constants, only : czero

     implicit none

! external arguments
! size of matrix
     integer, intent(in)      :: n

! output matrix's trace
     complex(dp), intent(out) :: tr

! input complex(dp) matrix
     complex(dp), intent(in)  :: A(n,n)

! local variables
! loop index
     integer :: i

     tr = czero
     do i=1,n
         tr = tr + A(i,i)
     enddo ! over i={1,n} loop

     return
  end subroutine s_trace_z

!!
!! @sub s_det_d
!!
!! calculate the determinant of a real(dp) matrix
!!
  subroutine s_det_d(ndim, dmat, ddet)
     use constants, only : dp
     use constants, only : one, cone

     implicit none

! external arguments
! dimension of dmat matrix
     integer, intent(in)     :: ndim

! determinant of dmat matrix
     real(dp), intent(out)   :: ddet

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

! used to calculate determinant
     complex(dp) :: cres

! working arrays for lapack subroutines: dgetrf
     integer, allocatable  :: ipiv(:)

! working arrays for lapack subroutines: dgeev
     real(dp), allocatable :: work(:)

! real and imaginary parts of the computed eigenvalues
     real(dp), allocatable :: wi(:)
     real(dp), allocatable :: wr(:)

! left and right eigenvectors
     real(dp), allocatable :: vl(:,:)
     real(dp), allocatable :: vr(:,:)

! dummy arrays, used to save dmat
     real(dp), allocatable :: amat(:,:)

! setup lwork
     lwork = 4 * ndim

! allocate memory
     allocate(ipiv(ndim),      stat=ierror)
     allocate(work(lwork),     stat=ierror)
     allocate(wi(ndim),        stat=ierror)
     allocate(wr(ndim),        stat=ierror)
     allocate(vl(ndim,ndim),   stat=ierror)
     allocate(vr(ndim,ndim),   stat=ierror)
     allocate(amat(ndim,ndim), stat=ierror)
     if ( ierror /= 0 ) then
         call s_print_error('s_det_d','can not allocate enough memory')
     endif ! back if ( ierror /= 0 ) block

! copy dmat to amat at first
     amat = dmat

!-------------------------------------------------------------------------
! method A: preferred method
!-------------------------------------------------------------------------
! computes the LU factorization of a general m-by-n matrix, need lapack
! package, dgetrf subroutine
     call DGETRF(ndim, ndim, dmat, ndim, ipiv, ierror)
     if ( ierror /= 0 ) then
         call s_print_exception('s_det_d','error in lapack subroutine dgetrf')
     endif ! back if ( ierror /= 0 ) block

! calculate determinant
     ddet = one
     do i=1,ndim
         if ( ipiv(i) == i ) then
             ddet = ddet * ( +dmat(i,i) )
         else
             ddet = ddet * ( -dmat(i,i) )
         endif ! back if ( ipiv(i) == i ) block
     enddo ! over i={1,ndim} loop

! everything is ok!
     if ( ierror == 0 ) RETURN

!-------------------------------------------------------------------------
! method B: as a backup
!-------------------------------------------------------------------------
! diagonalize amat to obtain its eigenvalues: wr and wi
     call DGEEV('N', 'N', ndim, amat, ndim, wr, wi, vl, ndim, vr, ndim, work, lwork, ierror)
     if ( ierror /= 0 ) then
         call s_print_error('s_det_d','error in lapack subroutine dgeev')
     endif ! back if ( ierror /= 0 ) block

! evaluate the final determinant
     cres = cone
     do i=1,ndim
         cres = cres * dcmplx( wr(i), wi(i) )
     enddo ! over i={1,ndim} loop
     ddet = real(cres)

! deallocate memory
     if ( allocated(ipiv) ) deallocate(ipiv)
     if ( allocated(work) ) deallocate(work)
     if ( allocated(wi  ) ) deallocate(wi  )
     if ( allocated(wr  ) ) deallocate(wr  )
     if ( allocated(vl  ) ) deallocate(vl  )
     if ( allocated(vr  ) ) deallocate(vr  )
     if ( allocated(amat) ) deallocate(amat)

     return
  end subroutine s_det_d

!!
!! @sub s_det_z
!!
!! calculate the determinant of a complex(dp) matrix
!!
  subroutine s_det_z(ndim, zmat, zdet)
     use constants, only : dp
     use constants, only : cone

     implicit none

! external arguments
! dimension of zmat matrix
     integer, intent(in)        :: ndim

! determinant of zmat matrix
     complex(dp), intent(out)   :: zdet

! object matrix, on entry, it contains the original matrix, on exit,
! it is destroyed and replaced with the L and U matrix
     complex(dp), intent(inout) :: zmat(ndim,ndim)

! local variables
! loop index
     integer :: i

! error flag
     integer :: ierror

! working arrays for lapack subroutines
     integer, allocatable :: ipiv(:)

! allocate memory
     allocate(ipiv(ndim), stat=ierror)
     if ( ierror /= 0 ) then
         call s_print_error('s_det_z','can not allocate enough memory')
     endif ! back if ( ierror /= 0 ) block

! computes the LU factorization of a general m-by-n matrix, need lapack
! package, zgetrf subroutine
     call ZGETRF(ndim, ndim, zmat, ndim, ipiv, ierror)
     if ( ierror /= 0 ) then
         call s_print_error('s_det_z','error in lapack subroutine zgetrf')
     endif ! back if ( ierror /= 0 ) block

! calculate determinant
     zdet = cone
     do i=1,ndim
         if ( ipiv(i) == i ) then
             zdet = zdet * ( +zmat(i,i) )
         else
             zdet = zdet * ( -zmat(i,i) )
         endif ! back if ( ipiv(i) == i ) block
     enddo ! over i={1,ndim} loop

! deallocate memory
     if ( allocated(ipiv) ) deallocate(ipiv)

     return
  end subroutine s_det_z

!!========================================================================
!!>>> matrix manipulation: calculate matrix's inversion                <<<
!!========================================================================

!!
!! @sub s_inv_d
!!
!! invert real(dp) matrix using lapack subroutines
!!
  subroutine s_inv_d(ndim, dmat)
     use constants, only : dp

     implicit none

! external arguments
! dimension of dmat matrix
     integer, intent(in)     :: ndim

! object matrix, on entry, it contains the original matrix, on exit,
! it is destroyed and replaced with the inversed matrix
     real(dp), intent(inout) :: dmat(ndim,ndim)

! local variables
! error flag
     integer  :: ierror

! working arrays for lapack subroutines
     integer, allocatable  :: ipiv(:)
     real(dp), allocatable :: work(:)

! allocate memory
     allocate(ipiv(ndim), stat=ierror)
     allocate(work(ndim), stat=ierror)
     if ( ierror /= 0 ) then
         call s_print_error('s_inv_d','can not allocate enough memory')
     endif ! back if ( ierror /= 0 ) block

! computes the LU factorization of a general m-by-n matrix, need lapack
! package, dgetrf subroutine
     call DGETRF(ndim, ndim, dmat, ndim, ipiv, ierror)
     if ( ierror /= 0 ) then
         call s_print_error('s_inv_d','error in lapack subroutine dgetrf')
     endif ! back if ( ierror /= 0 ) block

! computes the inverse of an LU-factored general matrix, need lapack
! package, dgetri subroutine
     call DGETRI(ndim, dmat, ndim, ipiv, work, ndim, ierror)
     if ( ierror /= 0 ) then
         call s_print_error('s_inv_d','error in lapack subroutine dgetri')
     endif ! back if ( ierror /= 0 ) block

! deallocate memory
     if ( allocated(ipiv) ) deallocate(ipiv)
     if ( allocated(work) ) deallocate(work)

     return
  end subroutine s_inv_d

!!
!! @sub s_inv_z
!!
!! invert complex(dp) matrix using lapack subroutines
!!
  subroutine s_inv_z(ndim, zmat)
     use constants, only : dp

     implicit none

! external arguments
! dimension of zmat matrix
     integer, intent(in)        :: ndim

! object matrix, on entry, it contains the original matrix, on exit,
! it is destroyed and replaced with the inversed matrix
     complex(dp), intent(inout) :: zmat(ndim,ndim)

! local variables
! error flag
     integer     :: ierror

! working arrays for lapack subroutines
     integer, allocatable     :: ipiv(:)
     complex(dp), allocatable :: work(:)

! allocate memory
     allocate(ipiv(ndim), stat=ierror)
     allocate(work(ndim), stat=ierror)
     if ( ierror /= 0 ) then
         call s_print_error('s_inv_z','can not allocate enough memory')
     endif ! back if ( ierror /= 0 ) block

! computes the LU factorization of a general m-by-n matrix, need lapack
! package, zgetrf subroutine
     call ZGETRF(ndim, ndim, zmat, ndim, ipiv, ierror)
     if ( ierror /= 0 ) then
         call s_print_error('s_inv_z','error in lapack subroutine zgetrf')
     endif ! back if ( ierror /= 0 ) block

! computes the inverse of an LU-factored general matrix, need lapack
! package, zgetri subroutine
     call ZGETRI(ndim, zmat, ndim, ipiv, work, ndim, ierror)
     if ( ierror /= 0 ) then
         call s_print_error('s_inv_z','error in lapack subroutine zgetri')
     endif ! back if ( ierror /= 0 ) block

! deallocate memory
     if ( allocated(ipiv) ) deallocate(ipiv)
     if ( allocated(work) ) deallocate(work)

     return
  end subroutine s_inv_z

!!========================================================================
!!>>> matrix manipulation: solve eigenvalues and eigenvectors problem  <<<
!!========================================================================

!!
!! @sub s_eig_dg
!!
!! diagonalize a general real(dp) matrix and return its eigenvalues and
!! eigenvectors
!!
  subroutine s_eig_dg(ldim, ndim, amat, eval, evec)
     use constants, only : dp
     use constants, only : zero

     implicit none

! external arguments
! leading dimension of matrix amat
     integer, intent(in)   :: ldim

! the order of the matrix amat
     integer, intent(in)   :: ndim

! original general real(dp) matrix to compute eigenvals and eigenvectors
     real(dp), intent(in)  :: amat(ldim,ndim)

! if info = 0, the eigenvalues in ascending order
     real(dp), intent(out) :: eval(ndim)

! if info = 0, orthonormal eigenvectors of the matrix
     real(dp), intent(out) :: evec(ldim,ndim)

! local variables
! status flag
     integer :: istat

! return information from subroutine dgeev
     integer :: info

! the length of the array work, lwork >= max(1,4*ndim)
     integer :: lwork

! workspace array
     real(dp), allocatable :: work(:)

! auxiliary real(dp) matrix: real and imaginary parts of eigenvalues
     real(dp), allocatable :: wr(:)
     real(dp), allocatable :: wi(:)

! auxiliary real(dp) matrix: left and right eigenvectors
     real(dp), allocatable :: vr(:,:)
     real(dp), allocatable :: vl(:,:)

! initialize lwork
     lwork = 4 * ndim

! allocate memory
     allocate(work(lwork),   stat=istat)
     allocate(wr(ndim),      stat=istat)
     allocate(wi(ndim),      stat=istat)
     allocate(vr(ndim,ndim), stat=istat)
     allocate(vl(ndim,ndim), stat=istat)
     if ( istat /= 0 ) then
         call s_print_error('s_eig_dg','can not allocate enough memory')
     endif ! back if ( istat /= 0 ) block

! initialize output arrays
     eval = zero
     evec = amat

! call the computational subroutine: dgeev
     call DGEEV('N', 'V', ndim, evec, ldim, wr, wi, vl, ndim, vr, ndim, work, lwork, info)

! check the status
     if ( info /= 0 ) then
         call s_print_error('s_eig_dg','error in lapack subroutine dgeev')
     endif ! back if ( info /= 0 ) block

! copy eigenvalues and eigenvectors
     eval(1:ndim) = wr(1:ndim)
     evec(1:ndim,1:ndim) = vr(1:ndim,1:ndim)

! dealloate memory for workspace array
     if ( allocated(work) ) deallocate(work)
     if ( allocated(wr  ) ) deallocate(wr  )
     if ( allocated(wi  ) ) deallocate(wi  )
     if ( allocated(vr  ) ) deallocate(vr  )
     if ( allocated(vl  ) ) deallocate(vl  )

     return
  end subroutine s_eig_dg

!!
!! @sub s_eig_zg
!!
!! diagonalize a general complex(dp) matrix and return its eigenvalues
!! and eigenvectors
!!
  subroutine s_eig_zg(ldim, ndim, zmat, zeig, zvec)
     use constants, only : dp
     use constants, only : czero

     implicit none

! external arguments
! leading dimension of matrix amat
     integer, intent(in)      :: ldim

! the order of the matrix amat
     integer, intent(in)      :: ndim

! original general complex(dp) matrix to compute eigenvals and eigenvectors
     complex(dp), intent(in)  :: zmat(ldim,ndim)

! if info = 0, the eigenvalues in ascending order
     complex(dp), intent(out) :: zeig(ndim)

! if info = 0, orthonormal eigenvectors of the matrix
     complex(dp), intent(out) :: zvec(ldim,ndim)

! local variables
! status flag
     integer :: istat

! return information from subroutine zgeev
     integer :: info

! the length of the array work, lwork >= max(1,2*ndim)
     integer :: lwork

! workspace array
     complex(dp), allocatable :: work(:)

! auxiliary real(dp) matrix
     complex(dp), allocatable :: rwork(:)

! auxiliary complex(dp) matrix: left and right eigenvectors
     complex(dp), allocatable :: vr(:,:)
     complex(dp), allocatable :: vl(:,:)

! initialize lwork
     lwork = 2 * ndim

! allocate memory
     allocate(work(lwork),   stat=istat)
     allocate(rwork(lwork),  stat=istat)
     allocate(vr(ndim,ndim), stat=istat)
     allocate(vl(ndim,ndim), stat=istat)
     if ( istat /= 0 ) then
         call s_print_error('s_eig_zg','can not allocate enough memory')
     endif ! back if ( istat /= 0 ) block

! initialize output arrays
     zeig = czero
     zvec = zmat

! call the computational subroutine: zgeev
     call ZGEEV('N', 'V', ndim, zvec, ldim, zeig, vl, ndim, vr, ndim, work, lwork, rwork, info)

! check the status
     if ( info /= 0 ) then
         call s_print_error('s_eig_zg','error in lapack subroutine zgeev')
     endif ! back if ( info /= 0 ) block

! copy eigenvectors
     zvec = vr

! dealloate memory for workspace array
     if ( allocated(work ) )  deallocate(work )
     if ( allocated(rwork) )  deallocate(rwork)
     if ( allocated(vr   ) )  deallocate(vr   )
     if ( allocated(vl   ) )  deallocate(vl   )

     return
  end subroutine s_eig_zg

!!
!! @sub s_eigvals_dg
!!
!! diagonalize a general real(dp) matrix and return its eigenvalues only
!!
  subroutine s_eigvals_dg(ldim, ndim, amat, eval)
     use constants, only : dp
     use constants, only : zero

     implicit none

! external arguments
! leading dimension of matrix amat
     integer, intent(in)   :: ldim

! the order of the matrix amat
     integer, intent(in)   :: ndim

! original general real(dp) matrix to compute eigenvals
     real(dp), intent(in)  :: amat(ldim,ndim)

! if info = 0, the eigenvalues in ascending order
     real(dp), intent(out) :: eval(ndim)

! local variables
! status flag
     integer :: istat

! return information from subroutine dgeev
     integer :: info

! the length of the array work, lwork >= max(1,4*ndim)
     integer :: lwork

! workspace array, used to store amat
     real(dp), allocatable :: evec(:,:)

! workspace array
     real(dp), allocatable :: work(:)

! auxiliary real(dp) matrix: real and imaginary parts of eigenvalues
     real(dp), allocatable :: wr(:)
     real(dp), allocatable :: wi(:)

! auxiliary real(dp) matrix: left and right eigenvectors
     real(dp), allocatable :: vr(:,:)
     real(dp), allocatable :: vl(:,:)

! initialize lwork
     lwork = 4 * ndim

! allocate memory
     allocate(evec(ldim,ndim), stat=istat)
     allocate(work(lwork),     stat=istat)
     allocate(wr(ndim),        stat=istat)
     allocate(wi(ndim),        stat=istat)
     allocate(vr(ndim,ndim),   stat=istat)
     allocate(vl(ndim,ndim),   stat=istat)
     if ( istat /= 0 ) then
         call s_print_error('s_eigvals_dg','can not allocate enough memory')
     endif ! back if ( istat /= 0 ) block

! initialize output arrays
     eval = zero
     evec = amat

! call the computational subroutine: dgeev
     call DGEEV('N', 'N', ndim, evec, ldim, wr, wi, vl, ndim, vr, ndim, work, lwork, info)

! check the status
     if ( info /= 0 ) then
         call s_print_error('s_eigvals_dg','error in lapack subroutine dgeev')
     endif ! back if ( info /= 0 ) block

! copy eigenvalues
     eval(1:ndim) = wr(1:ndim)

! dealloate memory for workspace array
     if ( allocated(evec) ) deallocate(evec)
     if ( allocated(work) ) deallocate(work)
     if ( allocated(wr  ) ) deallocate(wr  )
     if ( allocated(wi  ) ) deallocate(wi  )
     if ( allocated(vr  ) ) deallocate(vr  )
     if ( allocated(vl  ) ) deallocate(vl  )

     return
  end subroutine s_eigvals_dg

!!
!! @sub s_eigvals_zg
!!
!! diagonalize a general complex(dp) matrix and return its eigenvalues only
!!
  subroutine s_eigvals_zg(ldim, ndim, zmat, zeig)
     use constants, only : dp
     use constants, only : czero

     implicit none

! external arguments
! leading dimension of matrix amat
     integer, intent(in)      :: ldim

! the order of the matrix amat
     integer, intent(in)      :: ndim

! original general complex(dp) matrix to compute eigenvals
     complex(dp), intent(in)  :: zmat(ldim,ndim)

! if info = 0, the eigenvalues in ascending order
     complex(dp), intent(out) :: zeig(ndim)

! local variables
! status flag
     integer :: istat

! return information from subroutine zgeev
     integer :: info

! the length of the array work, lwork >= max(1,2*ndim)
     integer :: lwork

! workspace array, used to store amat
     complex(dp), allocatable :: zvec(:,:)

! workspace array
     complex(dp), allocatable :: work(:)

! auxiliary real(dp) matrix
     complex(dp), allocatable :: rwork(:)

! auxiliary complex(dp) matrix: left and right eigenvectors
     complex(dp), allocatable :: vr(:,:)
     complex(dp), allocatable :: vl(:,:)

! initialize lwork
     lwork = 2 * ndim

! allocate memory
     allocate(zvec(ldim,ndim), stat=istat)
     allocate(work(lwork),     stat=istat)
     allocate(rwork(lwork),    stat=istat)
     allocate(vr(ndim,ndim),   stat=istat)
     allocate(vl(ndim,ndim),   stat=istat)
     if ( istat /= 0 ) then
         call s_print_error('s_eigvals_zg','can not allocate enough memory')
     endif ! back if ( istat /= 0 ) block

! initialize output arrays
     zeig = czero
     zvec = zmat

! call the computational subroutine: zgeev
     call ZGEEV('N', 'N', ndim, zvec, ldim, zeig, vl, ndim, vr, ndim, work, lwork, rwork, info)

! check the status
     if ( info /= 0 ) then
         call s_print_error('s_eigvals_zg','error in lapack subroutine zgeev')
     endif ! back if ( info /= 0 ) block

! dealloate memory for workspace array
     if ( allocated(zvec ) )  deallocate(zvec )
     if ( allocated(work ) )  deallocate(work )
     if ( allocated(rwork) )  deallocate(rwork)
     if ( allocated(vr   ) )  deallocate(vr   )
     if ( allocated(vl   ) )  deallocate(vl   )

     return
  end subroutine s_eigvals_zg

!!
!! @sub s_eig_sy
!!
!! computes all eigenvalues and eigenvectors of real symmetric matrix
!!
  subroutine s_eig_sy(ldim, ndim, amat, eval, evec)
     use constants, only : dp
     use constants, only : zero

     implicit none

! external arguments
! leading dimension of matrix amat
     integer, intent(in)   :: ldim

! the order of the matrix amat
     integer, intent(in)   :: ndim

! original real symmetric matrix to compute eigenvals and eigenvectors
     real(dp), intent(in)  :: amat(ldim,ndim)

! if info = 0, the eigenvalues in ascending order
     real(dp), intent(out) :: eval(ndim)

! if info = 0, orthonormal eigenvectors of the matrix
     real(dp), intent(out) :: evec(ldim,ndim)

! local variables
! status flag
     integer :: istat

! return information from subroutine dysev
     integer :: info

! the length of the array work, lwork >= max(1,3*ndim-1)
     integer :: lwork

! workspace array
     real(dp), allocatable :: work(:)

! initialize lwork
     lwork = 3 * ndim - 1

! allocate memory
     allocate(work(lwork), stat=istat)
     if ( istat /= 0 ) then
         call s_print_error('s_eig_sy','can not allocate enough memory')
     endif ! back if ( istat /= 0 ) block

! initialize output arrays
     eval = zero
     evec = amat

! call the computational subroutine: dsyev
     call DSYEV('V', 'U', ndim, evec, ldim, eval, work, lwork, info)

! check the status
     if ( info /= 0 ) then
         call s_print_error('s_eig_sy','error in lapack subroutine dsyev')
     endif ! back if ( info /= 0 ) block

! dealloate memory for workspace array
     if ( allocated(work) ) deallocate(work)

     return
  end subroutine s_eig_sy

!!
!! @sub s_eig_he
!!
!! computes all eigenvalues and eigenvectors of complex Hermitian matrix
!!
  subroutine s_eig_he(ldim, ndim, amat, eval, evec)
     use constants, only : dp
     use constants, only : zero

     implicit none

! external arguments
! leading dimension of matrix amat
     integer, intent(in)      :: ldim

! the order of the matrix amat
     integer, intent(in)      :: ndim

! original complex Hermitian matrix to compute eigenvals and eigenvectors
     complex(dp), intent(in)  :: amat(ldim,ndim)

! if info = 0, the eigenvalues in ascending order
     real(dp), intent(out)    :: eval(ndim)

! if info = 0, orthonormal eigenvectors of the matrix
     complex(dp), intent(out) :: evec(ldim,ndim)

! local variables
! status flag
     integer :: istat

! return information from subroutine zheev
     integer :: info

! the length of the array work and rwork
! lwork >= max(1,2*ndim-1), lrwork >= max(1,3*ndim-2)
     integer :: lwork
     integer :: lrwork

! workspace array
     real(dp), allocatable    :: rwork(:)
     complex(dp), allocatable :: work(:)

! initialize lwork (lrwork)
     lwork = 2 * ndim - 1
     lrwork = 3 * ndim - 2

! allocate memory
     allocate(work(lwork),   stat=istat)
     allocate(rwork(lrwork), stat=istat)
     if ( istat /= 0 ) then
         call s_print_error('s_eig_he','can not allocate enough memory')
     endif ! back if ( istat /= 0 ) block

! initialize output arrays
     eval = zero
     evec = amat

! call the computational subroutine: zheev
     call ZHEEV('V', 'U', ndim, evec, ldim, eval, work, lwork, rwork, info)

! check the status
     if ( info /= 0 ) then
         call s_print_error('s_eig_he','error in lapack subroutine zheev')
     endif ! back if ( info /= 0 ) block

! dealloate memory for workspace array
     if ( allocated(work ) ) deallocate(work )
     if ( allocated(rwork) ) deallocate(rwork)

     return
  end subroutine s_eig_he

!!
!! @sub s_eigvals_sy
!!
!! computes all eigenvalues of real symmetric matrix
!!
  subroutine s_eigvals_sy(ldim, ndim, amat, eval)
     use constants, only : dp
     use constants, only : zero

     implicit none

! external arguments
! leading dimension of matrix amat
     integer, intent(in)   :: ldim

! the order of the matrix amat
     integer, intent(in)   :: ndim

! original real symmetric matrix to compute eigenvals
     real(dp), intent(in)  :: amat(ldim,ndim)

! if info = 0, the eigenvalues in ascending order
     real(dp), intent(out) :: eval(ndim)

! local variables
! status flag
     integer :: istat

! return information from subroutine dysev
     integer :: info

! the length of the array work, lwork >= max(1,3*ndim-1)
     integer :: lwork

! workspace array
     real(dp), allocatable :: work(:)

! workspace array, used to store amat
     real(dp), allocatable :: evec(:,:)

! initialize lwork
     lwork = 3 * ndim - 1

! allocate memory
     allocate(work(lwork),     stat=istat)
     allocate(evec(ldim,ndim), stat=istat)
     if ( istat /= 0 ) then
         call s_print_error('s_eigvals_sy','can not allocate enough memory')
     endif ! back if ( istat /= 0 ) block

! initialize output arrays
     eval = zero
     evec = amat

! call the computational subroutine: dsyev
     call DSYEV('N', 'U', ndim, evec, ldim, eval, work, lwork, info)

! check the status
     if ( info /= 0 ) then
         call s_print_error('s_eigvals_sy','error in lapack subroutine dsyev')
     endif ! back if ( info /= 0 ) block

! dealloate memory for workspace array
     if ( allocated(work) ) deallocate(work)
     if ( allocated(evec) ) deallocate(evec)

     return
  end subroutine s_eigvals_sy

!!
!! @sub s_eigvals_he
!!
!! computes all eigenvalues of complex Hermitian matrix
!!
  subroutine s_eigvals_he(ldim, ndim, amat, eval)
     use constants, only : dp
     use constants, only : zero

     implicit none

! external arguments
! leading dimension of matrix amat
     integer, intent(in)     :: ldim

! the order of the matrix amat
     integer, intent(in)     :: ndim

! original complex Hermitian matrix to compute eigenvals and eigenvectors
     complex(dp), intent(in) :: amat(ldim,ndim)

! if info = 0, the eigenvalues in ascending order
     real(dp), intent(out)   :: eval(ndim)

! local variables
! status flag
     integer :: istat

! return information from subroutine zheev
     integer :: info

! the length of the array work and rwork
! lwork >= max(1,2*ndim-1), lrwork >= max(1,3*ndim-2)
     integer :: lwork
     integer :: lrwork

! workspace array
     real(dp), allocatable    :: rwork(:)
     complex(dp), allocatable :: work(:)

! workspace array, used to store amat
     complex(dp), allocatable :: evec(:,:)

! initialize lwork (lrwork)
     lwork = 2 * ndim - 1
     lrwork = 3 * ndim - 2

! allocate memory
     allocate(work(lwork),     stat=istat)
     allocate(rwork(lrwork),   stat=istat)
     allocate(evec(ldim,ndim), stat=istat)
     if ( istat /= 0 ) then
         call s_print_error('s_eigvals_he','can not allocate enough memory')
     endif ! back if ( istat /= 0 ) block

! initialize output arrays
     eval = zero
     evec = amat

! call the computational subroutine: zheev
     call ZHEEV('N', 'U', ndim, evec, ldim, eval, work, lwork, rwork, info)

! check the status
     if ( info /= 0 ) then
         call s_print_error('s_eigvals_he','error in lapack subroutine zheev')
     endif ! back if ( info /= 0 ) block

! dealloate memory for workspace array
     if ( allocated(work ) ) deallocate(work )
     if ( allocated(rwork) ) deallocate(rwork)
     if ( allocated(evec ) ) deallocate(evec )

     return
  end subroutine s_eigvals_he

!!========================================================================
!!>>> matrix manipulation: solve linear equations                      <<<
!!========================================================================

!!
!! @sub s_solve_dg
!!
!! solve linear system AX = B, real(dp) general version
!!
  subroutine s_solve_dg(n, nrhs, A, B)
     use constants, only : dp

     implicit none

! external arguments
! the number of linear equations
     integer, intent(in)     :: n

! the number of right-hand sides
     integer, intent(in)     :: nrhs

! on entry, it is a n-by-n coefficient matrix A; on exit, it is overwritten
! by the factors L and U from the factorization of A = PLU.
     real(dp), intent(inout) :: A(n,n)

! on entry, it is a n-by-nrhs matrix of right hand side matrix B; on exit,
! it is overwritten by the solution matrix X.
     real(dp), intent(inout) :: B(n,nrhs)

! local variables
! status flag
     integer :: istat

! return information from subroutine dgesv
     integer :: info

! workspace array, its dimension is at least max(1,n)
     integer, allocatable :: ipiv(:)

! allocate memory
     allocate(ipiv(n), stat=istat)
     if ( istat /= 0 ) then
         call s_print_error('s_solve_dg','can not allocate enough memory')
     endif ! back if ( istat /= 0 ) block

! call the computational subroutine: dgesv
     call DGESV(n, nrhs, A, n, ipiv, B, n, info)

! check the status
     if ( info /= 0 ) then
         call s_print_error('s_solve_dg','error in lapack subroutine dgesv')
     endif ! back if ( info /= 0 ) block

! deallocate memory
     if ( allocated(ipiv) ) deallocate(ipiv)

     return
  end subroutine s_solve_dg

!!
!! @sub s_solve_zg
!!
!! solve linear system AX = B, complex(dp) general version
!!
  subroutine s_solve_zg(n, nrhs, A, B)
     use constants, only : dp

     implicit none

! external arguments
! the number of linear equations
     integer, intent(in)        :: n

! the number of right-hand sides
     integer, intent(in)        :: nrhs

! on entry, it is a n-by-n coefficient matrix A; on exit, it is overwritten
! by the factors L and U from the factorization of A = PLU.
     complex(dp), intent(inout) :: A(n,n)

! on entry, it is a n-by-nrhs matrix of right hand side matrix B; on exit,
! it is overwritten by the solution matrix X.
     complex(dp), intent(inout) :: B(n,nrhs)

! local variables
! status flag
     integer :: istat

! return information from subroutine zgesv
     integer :: info

! workspace array, its dimension is at least max(1,n)
     integer, allocatable :: ipiv(:)

! allocate memory
     allocate(ipiv(n), stat=istat)
     if ( istat /= 0 ) then
         call s_print_error('s_solve_zg','can not allocate enough memory')
     endif ! back if ( istat /= 0 ) block

! call the computational subroutine: zgesv
     call ZGESV(n, nrhs, A, n, ipiv, B, n, info)

! check the status
     if ( info /= 0 ) then
         call s_print_error('s_solve_zg','error in lapack subroutine zgesv')
     endif ! back if ( info /= 0 ) block

! deallocate memory
     if ( allocated(ipiv) ) deallocate(ipiv)

     return
  end subroutine s_solve_zg

!!
!! @sub s_solve_sy
!!
!! solve linear system AX = B, real(dp) symmetric version
!!
  subroutine s_solve_sy(n, nrhs, A, B)
     use constants, only : dp

     implicit none

! external arguments
! the number of linear equations
     integer, intent(in)     :: n

! the number of right-hand sides
     integer, intent(in)     :: nrhs

! on entry, it is a n-by-n coefficient matrix A; on exit, it is overwritten
! by the factors L and U from the factorization of A = PLU.
     real(dp), intent(inout) :: A(n,n)

! on entry, it is a n-by-nrhs matrix of right hand side matrix B; on exit,
! it is overwritten by the solution matrix X.
     real(dp), intent(inout) :: B(n,nrhs)

! local variables
! status flag
     integer :: istat

! return information from subroutine dsysv
     integer :: info

! workspace array, its dimension is at least max(1,n)
     integer, allocatable  :: ipiv(:)

! workspace array, its dimension is at least max(1, lwork) and lwork >= 1
     real(dp), allocatable :: work(:)

! allocate memory
     allocate(ipiv(n), stat=istat)
     allocate(work(n), stat=istat)
     if ( istat /= 0 ) then
         call s_print_error('s_solve_sy','can not allocate enough memory')
     endif ! back if ( istat /= 0 ) block

! call the computational subroutine: dsysv
     call DSYSV('U', n, nrhs, A, n, ipiv, B, n, work, n, info)

! check the status
     if ( info /= 0 ) then
         call s_print_error('s_solve_sy','error in lapack subroutine dsysv')
     endif ! back if ( info /= 0 ) block

! deallocate memory
     if ( allocated(ipiv) ) deallocate(ipiv)
     if ( allocated(work) ) deallocate(work)

     return
  end subroutine s_solve_sy

!!
!! @sub s_solve_he
!!
!! solve linear system AX = B, complex(dp) Hermitian version
!!
  subroutine s_solve_he(n, nrhs, A, B)
     use constants, only : dp

     implicit none

! external arguments
! the number of linear equations
     integer, intent(in)        :: n

! the number of right-hand sides
     integer, intent(in)        :: nrhs

! on entry, it is a n-by-n coefficient matrix A; on exit, it is overwritten
! by the factors L and U from the factorization of A = PLU.
     complex(dp), intent(inout) :: A(n,n)

! on entry, it is a n-by-nrhs matrix of right hand side matrix B; on exit,
! it is overwritten by the solution matrix X.
     complex(dp), intent(inout) :: B(n,nrhs)

! local variables
! status flag
     integer :: istat

! return information from subroutine zhesv
     integer :: info

! workspace array, its dimension is at least max(1,n)
     integer, allocatable     :: ipiv(:)

! workspace array, its dimension is at least max(1, lwork) and lwork >= 1
     complex(dp), allocatable :: work(:)

! allocate memory
     allocate(ipiv(n), stat=istat)
     allocate(work(n), stat=istat)
     if ( istat /= 0 ) then
         call s_print_error('s_solve_he','can not allocate enough memory')
     endif ! back if ( istat /= 0 ) block

! call the computational subroutine: zhesv
     call ZHESV('U', n, nrhs, A, n, ipiv, B, n, work, n, info)

! check the status
     if ( info /= 0 ) then
         call s_print_error('s_solve_he','error in lapack subroutine zhesv')
     endif ! back if ( info /= 0 ) block

! deallocate memory
     if ( allocated(ipiv) ) deallocate(ipiv)
     if ( allocated(work) ) deallocate(work)

     return
  end subroutine s_solve_he

!!========================================================================
!!>>> matrix manipulation: singular values decomposition               <<<
!!========================================================================

!!
!! @sub s_svd_dg
!!
!! perform the singular values decomposition for a general real(dp) m-by-n
!! matrix A, where A = U * SIGMA * transpose(V), return its left vectors,
!! right vectors, and singular values
!!
  subroutine s_svd_dg(m, n, min_mn, amat, umat, svec, vmat)
     use constants, only : dp

     implicit none

! external arguments
! number of rows of A matrix
     integer, intent(in)     :: m

! number of columns of A matrix
     integer, intent(in)     :: n

! minimal value of m and n
     integer, intent(in)     :: min_mn

! A matrix
     real(dp), intent(inout) :: amat(m,n)

! left vectors of svd, U
     real(dp), intent(out)   :: umat(m,min_mn)

! singular values of svd, SIGMA
     real(dp), intent(out)   :: svec(min_mn)

! right vectors of svd, transpose(V)
     real(dp), intent(out)   :: vmat(min_mn,n)

! local variables
! status flag
     integer :: istat

! return information from dgesvd
     integer :: info

! length of work array, lwork >= max(1, 3 * min_mn + max(m,n), 5 * min_mn)
     integer :: lwork

! workspace array
     real(dp), allocatable :: work(:)

! initialize lwrok
     lwork = max(1, 3 * min_mn + max(m,n), 5 * min_mn)

! allocate memory
     allocate(work(lwork), stat=istat)
     if ( istat /= 0 ) then
         call s_print_error('s_svd_dg','can not allocate enough memory')
     endif ! back if ( istat /= 0 ) block

! call the computational subroutine: dgesvd
     call DGESVD('S', 'S', m, n, amat, m, svec, umat, m, vmat, min_mn, work, lwork, info)

! check the status
     if ( info /= 0 ) then
         call s_print_error('s_svd_dg','error in lapack subroutine dgesvd')
     endif ! back if ( info /= 0 ) block

! deallocate the memory for workspace array
     if ( allocated(work) ) deallocate(work)

     return
  end subroutine s_svd_dg

!!
!! @sub s_svd_zg
!!
!! perform the singular values decomposition for a general complex(dp)
!! m-by-n matrix A, where A = U * SIGMA * conjugate-transpose(V), return
!! its left vectors, right vectors, and singular values
!!
  subroutine s_svd_zg(m, n, min_mn, amat, umat, svec, vmat)
     use constants, only : dp

     implicit none

! external arguments
! number of rows of A matrix
     integer, intent(in)        :: m

! number of columns of A matrix
     integer, intent(in)        :: n

! minimal value of m and n
     integer, intent(in)        :: min_mn

! A matrix
     complex(dp), intent(inout) :: amat(m,n)

! left vectors of svd, U
     complex(dp), intent(out)   :: umat(m,min_mn)

! singular values of svd, SIGMA
     real(dp), intent(out)      :: svec(min_mn)

! right vectors of svd, conjugate-transpose(V)
     complex(dp), intent(out)   :: vmat(min_mn,n)

! local variables
! status flag
     integer :: istat

! return information from zgesvd
     integer :: info

! length of work array, lwork >= max(1, 2 * min_mn + max(m,n))
     integer :: lwork

! workspace arrays
     complex(dp), allocatable :: work(:)
     real(dp), allocatable :: rwork(:)

! initialize lwrok
     lwork = max(1, 2 * min_mn + max(m,n))

! allocate memory
     allocate(work(lwork),     stat=istat)
     allocate(rwork(5*min_mn), stat=istat)
     if ( istat /= 0 ) then
         call s_print_error('s_svd_zg','can not allocate enough memory')
     endif ! back if ( istat /= 0 ) block

! call the computational subroutine: zgesvd
     call ZGESVD('S', 'S', m, n, amat, m, svec, umat, m, vmat, min_mn, work, lwork, rwork, info)

! check the status
     if ( info /= 0 ) then
         call s_print_error('s_svd_zg','error in lapack subroutine zgesvd')
     endif ! back if ( info /= 0 ) block

! deallocate the memory for workspace array
     if ( allocated(work ) ) deallocate(work )
     if ( allocated(rwork) ) deallocate(rwork)

     return
  end subroutine s_svd_zg
