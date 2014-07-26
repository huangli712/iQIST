
!>>> invert real(dp) matrix using lapack subroutines
  subroutine ctqmc_dmat_inv(ndim, dmat)
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
         call ctqmc_print_error('ctqmc_dmat_inv','error in lapack subroutine dgetrf')
     endif

! computes the inverse of an LU-factored general matrix, need lapack
! package, dgetri subroutine
     call dgetri(ndim, dmat, ndim, ipiv, work, ndim, ierror)
     if ( ierror /= 0 ) then
         call ctqmc_print_error('ctqmc_dmat_inv','error in lapack subroutine dgetri')
     endif

     return
  end subroutine ctqmc_dmat_inv

!>>> invert complex(dp) matrix using lapack subroutines
  subroutine ctqmc_zmat_inv(ndim, zmat)
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
         call ctqmc_print_error('ctqmc_zmat_inv','error in lapack subroutine zgetrf')
     endif

! computes the inverse of an LU-factored general matrix, need lapack
! package, zgetri subroutine
     call zgetri(ndim, zmat, ndim, ipiv, work, ndim, ierror)
     if ( ierror /= 0 ) then
         call ctqmc_print_error('ctqmc_zmat_inv','error in lapack subroutine zgetri')
     endif

     return
  end subroutine ctqmc_zmat_inv

!>>> calculate the determinant of a real(dp) matrix
  subroutine ctqmc_dmat_det(ndim, dmat, ddet)
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
     call dgetrf(ndim, ndim, dmat, ndim, ipiv, ierror)
     if ( ierror /= 0 ) then
         call ctqmc_print_exception('ctqmc_dmat_det','error in lapack subroutine dgetrf')
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
     call dgeev('N', 'N', ndim, amat, ndim, wr, wi, vl, ndim, vr, ndim, work, lwork, ierror)
     if ( ierror /= 0 ) then
         call ctqmc_print_error('ctqmc_dmat_det','error in lapack subroutine dgeev')
     endif

! evaluate the final determinant
     cres = cone
     do i=1,ndim
         cres = cres * dcmplx( wr(i), wi(i) )
     enddo ! over i={1,ndim} loop
     ddet = cres

     return
  end subroutine ctqmc_dmat_det

!>>> calculate the determinant of a complex(dp) matrix
  subroutine ctqmc_zmat_det(ndim, zmat, zdet)
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
     call zgetrf(ndim, ndim, zmat, ndim, ipiv, ierror)
     if ( ierror /= 0 ) then
         call ctqmc_print_error('ctqmc_zmat_det','error in lapack subroutine zgetrf')
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
  end subroutine ctqmc_zmat_det
