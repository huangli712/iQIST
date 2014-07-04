!>>> compute the determinant of orginal matrix
  subroutine dmat_deter( ndim, dmat, deter )

     implicit none

! external arguments
! dimension of the input square matrix
     integer, intent(in) :: ndim

! orginal matrix to compute determinant, on exit, it's replaced by L and U
     real(8), intent(inout) :: dmat(ndim, ndim)

! determinant of the orginal matrix dmat
     real(8), intent(out) :: deter

! local variables
! the pivot index, row i of the matrix was interchanged with row ipiv(i)
     integer :: ipiv(ndim)

! return information from subroutine dgetrf
     integer :: info

! determinant of original matrix. determinant = det(1) * ten**det(2)
     real(8) :: det(2)

! auxiliary constant, ten = 10.d0.
     real(8) :: ten

! loop index
     integer :: i

! compute an LU factorization using partial pivoting with row interchanges
! A = P * L * U, where P is a permutation matrix, L is lower triangular 
! with unit diagonal diagonal element, U is upper triangular
     call dgetrf( ndim, ndim, dmat, ndim, ipiv, info )
     if ( info /= 0 ) then
         print*, 'info=', info
         stop 'severe error happened in subroutine dgetrf'
     endif ! back if ( info /= 0 ) blcok

! computer determinant by just multiply the diagonal elements
! special trick to avoid numerical precision loss of determinant value
     ten    = 10.0d0
     det(1) = 1.00d0
     det(2) = 0.00d0

     do i=1,ndim
         if( ipiv(i) /= i ) det(1) = - det(1)
         det(1) = dmat(i,i) * det(1)
         if ( det(1) .eq. 0.0d0 ) exit
         do while ( abs(det(1)) < 1.0d0 )
             det(1) = ten * det(1)
             det(2) = det(2) - 1.0d0
         enddo ! over while ( abs(det(1) < 1.0d0 ) loop
         do while ( abs(det(1)) > 10.d0 )
             det(1) = det(1) / ten
             det(2) = det(2) + 1.0d0
         enddo ! over while ( abs(det(1)) > 10.d0 ) loop
     enddo ! over i={1,ndim} loop

     deter = det(1) * ten**det(2)

     return
  end subroutine dmat_deter
  subroutine dmat_dgeev( ldim, ndim, amat, eval, evec )
     implicit none

! external variables
! leading dimension of matrix amat
     integer, intent(in) :: ldim

! the order of the matrix amat
     integer, intent(in) :: ndim

! original real symmetric matrix to compute eigenval and eigenvector
     real(8), intent(in) :: amat(ldim, ndim)

! if info = 0, the eigenvalues in ascending order.
     real(8), intent(out) :: eval(ndim)

! if info = 0, orthonormal eigenvectors of the matrix A
     real(8), intent(out) :: evec(ldim, ndim)

! local variables
! status flag
     integer :: istat

! return information from subroutine dysev
     integer :: info

! the length of the array work, lwork >= max(1,4*ndim)
     integer :: lwork

! workspace array
     real(8), allocatable :: work(:)

! auxiliary real(dp) matrix
     real(8), allocatable :: wr(:)
     real(8), allocatable :: wi(:)

     real(8), allocatable :: vr(:, :)
     real(8), allocatable :: vl(:, :)

! initialize lwork and allocate memory fo array work
     lwork = 4*ndim
     allocate(work(lwork), stat=istat)
     allocate(wr(ndim), stat=istat)
     allocate(wi(ndim), stat=istat)
     allocate(vr(ndim, ndim), stat=istat)
     allocate(vl(ndim, ndim), stat=istat)
     if ( istat /= 0 ) then
         stop "allocate memory error in dmat_dsyev"
     endif ! back if ( istat /= 0 ) block

! initialize output arrays
     eval = 0.d0
     evec = amat

     call DGEEV('N', 'V', ndim, evec, ldim, wr, wi, &
                vl, ndim, vr, ndim, work, lwork, info)
     if (info /= 0) stop "Failure in subroutine dmat_dgeev"
     eval(1:ndim) = wr(1:ndim)
     evec(1:ndim, 1:ndim) = vr(1:ndim, 1:ndim)

! dealloate memory for workspace array
     if (allocated(work)) deallocate(work)
     if (allocated(wr  )) deallocate(wr  )
     if (allocated(wi  )) deallocate(wi  )
     if (allocated(vr  )) deallocate(vr  )
     if (allocated(vl  )) deallocate(vl  )

     return
  end subroutine dmat_dgeev
!>>> perform real matrix-matrix multiply operation
  subroutine dmat_dgemm( ndim1, ndim2, ndim3, amat, bmat, cmat )
     implicit none

! dimension of the input square matrix 'amat and bmat'
     integer, intent(in)  :: ndim1
     integer, intent(in)  :: ndim2
     integer, intent(in)  :: ndim3

! input square matrix 'amat and bmat'
     real(8), intent(in)  :: amat(ndim1, ndim2)
     real(8), intent(in)  :: bmat(ndim2, ndim3)

! output square matrix, cmat = amat * bmat
     real(8), intent(out) :: cmat(ndim1, ndim3)

! local variables
     real(8) :: alpha
     real(8) :: betta

     alpha = 1.0d0; betta = 0.0d0; cmat = 0.0d0
     call dgemm('N', 'N', ndim1, ndim3, ndim2, alpha, amat, ndim1, &
                                bmat, ndim2, betta, cmat, ndim1 )

     return
  end subroutine dmat_dgemm

!>>> perform real matrix-matrix multiply operation
  subroutine dmat_dgemm0( ndim, amat, bmat, cmat )
     implicit none

! dimension of the input square matrix 'amat and bmat'
     integer, intent(in)  :: ndim

! input square matrix 'amat and bmat'
     real(8), intent(in)  :: amat(ndim, ndim)
     real(8), intent(in)  :: bmat(ndim, ndim)

! output square matrix, cmat = amat * bmat
     real(8), intent(out) :: cmat(ndim, ndim)

! local variables
     real(8) :: alpha
     real(8) :: betta

     alpha = 1.0d0; betta = 0.0d0; cmat = 0.0d0
     call dgemm('N', 'N', ndim, ndim, ndim, alpha, amat, ndim, &
                                bmat, ndim, betta, cmat, ndim )

     return
  end subroutine dmat_dgemm0

!>>> perform real matrix-matrix multiply operation
  subroutine dmat_dgemm1( ndim, amat, bmat, cmat )
     implicit none

! dimension of the input square matrix 'amat and bmat'
     integer, intent(in)  :: ndim

! input square matrix 'amat and bmat'
     real(8), intent(in)  :: amat(ndim, ndim)
     real(8), intent(in)  :: bmat(ndim, ndim)

! output square matrix, cmat = amat * bmat
     real(8), intent(out) :: cmat(ndim, ndim)

! local variables
     real(8) :: alpha
     real(8) :: betta

     alpha = 1.0d0; betta = 0.0d0; cmat = 0.0d0
     call dgemm('T', 'N', ndim, ndim, ndim, alpha, amat, ndim, &
                                bmat, ndim, betta, cmat, ndim )

     return
  end subroutine dmat_dgemm1

!>>> computes all eigenvalues and eigenvectors of real symmetric matrix A.
  subroutine dmat_dsyev( ldim, ndim, amat, eval, evec )
     implicit none

! external variables
! leading dimension of matrix amat
     integer, intent(in) :: ldim

! the order of the matrix amat
     integer, intent(in) :: ndim

! original real symmetric matrix to compute eigenval and eigenvector
     real(8), intent(in) :: amat(ldim, ndim)

! if info = 0, the eigenvalues in ascending order.
     real(8), intent(out) :: eval(ndim)

! if info = 0, orthonormal eigenvectors of the matrix A
     real(8), intent(out) :: evec(ldim, ndim)

! local variables
! status flag
     integer :: istat

! return information from subroutine dysev
     integer :: info

! the length of the array work, lwork >= max(1,3*ndim-1)
     integer :: lwork

! workspace array
     real(8), allocatable :: work(:)

! initialize lwork and allocate memory fo array work
     lwork = 3*ndim-1
     allocate(work(lwork), stat=istat)
     if ( istat /= 0 ) then
         stop "allocate memory error in dmat_dsyev"
     endif ! back if ( istat /= 0 ) block

! initialize output arrays
     eval = 0.d0
     evec = amat

     call DSYEV('V', 'U', ndim, evec, ldim, eval, work, lwork, info)
     if (info /= 0) stop "Failure in subroutine cal_dsyev"

! dealloate memory for workspace array
     if (allocated(work)) deallocate(work)

     return
  end subroutine dmat_dsyev
!>>> dump a one-dimensional real(dp) matrix to user defined device
  subroutine dmat_dump1(myout, ndims, dmat)

     implicit none

! external arguments
! output file unit
     integer, intent(in) :: myout

! leading dimension of the matrix
     integer, intent(in) :: ndims

! target matrix to be export
     real(8), intent(in) :: dmat(ndims)

! local variables
! loop index over dimensions
     integer :: idim1

! dump dmat to user defined device
     do idim1=1,ndims
         if (abs(dmat(idim1)) .lt. 1.0D-9) cycle
         write(myout, '(I8, F17.10)') idim1, dmat(idim1)
     enddo ! over idim1={1,ndims} loop

     return
  end subroutine dmat_dump1

!>>> dump a two-dimensional real(dp) matrix to user defined device
  subroutine dmat_dump2(myout, ldims, ndims, dmat)

     implicit none

! external arguments
! output file unit
     integer, intent(in) :: myout

! leading dimension of the matrix
     integer, intent(in) :: ldims

! second dimension of the matrix
     integer, intent(in) :: ndims

! target matrix to be export
     real(8), intent(in) :: dmat(ldims, ndims)

! local variables
! loop index over dimensions
     integer :: idim1
     integer :: idim2

! dump dmat to user defined device
     do idim2=1,ndims
         do idim1=1,ldims
             if (abs(dmat(idim1, idim2)) .lt. 1.0D-9) cycle
             write(myout, '(2I8, F17.10)') idim1, idim2, dmat(idim1, idim2)
         enddo ! over idim1={1,ldims} loop
     enddo ! over idim2={1,ndims} loop

     return
  end subroutine dmat_dump2

!>>> dump a three-dimensional real(dp) matrix to user defined device
  subroutine dmat_dump3(myout, ldims, mdims, ndims, dmat)

     implicit none

! external arguments
! output file unit
     integer, intent(in) :: myout

! leading dimension of the matrix
     integer, intent(in) :: ldims

! second dimension of the matrix
     integer, intent(in) :: mdims

! third dimension of the matrix
     integer, intent(in) :: ndims

! target matrix to be export
     real(8), intent(in) :: dmat(ldims, mdims, ndims)

! local variables
! loop index over dimensions
     integer :: idim1
     integer :: idim2
     integer :: idim3

! dump dmat to user defined device
     do idim3=1,ndims
         do idim2=1,mdims
             do idim1=1,ldims
                 if (abs(dmat(idim1, idim2, idim3)) .lt. 1.0D-9) cycle
                 write(myout, '(3I8, F17.10)') idim1, idim2, idim3, dmat(idim1, idim2, idim3)
             enddo ! over idim1={1,ldims} loop
         enddo ! over idim2={1,mdims} loop
     enddo ! over idim3={1,ndims} loop

     return
  end subroutine dmat_dump3

!>>> dump a four-dimensional real(dp) matrix to user defined device
  subroutine dmat_dump4(myout, ldims, mdims, ndims, kdims, dmat)

     implicit none

! external arguments
! output file unit
     integer, intent(in) :: myout

! leading dimension of the matrix
     integer, intent(in) :: ldims

! second dimension of the matrix
     integer, intent(in) :: mdims

! third dimension of the matrix
     integer, intent(in) :: ndims

! fourth dimension of the matrix
     integer, intent(in) :: kdims

! target matrix to be export
     real(8), intent(in) :: dmat(ldims, mdims, ndims, kdims)

! local variables
! loop index over dimensions
     integer :: idim1
     integer :: idim2
     integer :: idim3
     integer :: idim4

! dump dmat to user defined device
     do idim4=1,kdims
         do idim3=1,ndims
             do idim2=1,mdims
                 do idim1=1,ldims
                     if (abs(dmat(idim1, idim2, idim3, idim4)) .lt. 1.0D-9) cycle
                     write(myout, '(4I8, F17.10)') idim1, idim2, idim3, idim4, dmat(idim1, idim2, idim3, idim4)
                 enddo ! over idim1={1,ldims} loop
             enddo ! over icol={1,mdims} loop
         enddo ! over irow={1,ndims} loop
     enddo ! over idim4={1,kdims} loop

     return
  end subroutine dmat_dump4

!>>> compute the inverse of a real matrix
  subroutine dmat_inver(ndim, dmat)

     implicit none

! external arguments
! dimension of the square matrix A
     integer, intent(in) :: ndim

! on entry, the ndim-by-ndim matrix will be inversed
! on exit , the inverse of the orginal matrix dmat
     real(8), intent(inout) :: dmat(ndim, ndim)

! local variables
! return flag, (info = 0) => successful exit
     integer :: info

! the pivot indices, row i of the matrix was interchanged with row ipiv(i)
     integer :: ipiv(ndim)

! real array used as workspace
     real(8) :: work(ndim)

! compute an LU factorization using partial pivoting with row interchanges
! A = P * L * U, where P is a permutation matrix, L is lower triangular 
! with unit diagonal diagonal element, U is upper triangular
     call dgetrf( ndim, ndim, dmat, ndim, ipiv, info )
     if ( info /= 0 ) then
         stop 'severe error happened in subroutine dgetrf'
     endif ! back if ( info /= 0 ) block

! compute the inverse by using LU fractorization computed by zgetrf
     call dgetri( ndim, dmat, ndim, ipiv, work, ndim, info )
     if ( info /= 0 ) then
         stop 'severe error happened in subroutine dgetri'
     endif ! back if ( info /= 0 ) block

     return
  end subroutine dmat_inver
!>>> dump a one-dimensional integer matrix to user defined device
  subroutine imat_dump1(myout, ndims, imat)

     implicit none

! external arguments
! output file unit
     integer, intent(in) :: myout

! leading dimension of the matrix
     integer, intent(in) :: ndims

! target matrix to be export
     integer, intent(in) :: imat(ndims)

! local variables
! loop index over dimensions
     integer :: idim1

! dump imat to output device
     do idim1=1,ndims
         if (imat(idim1) .eq. 0) cycle
         write(myout, '(2I8)') idim1, imat(idim1)
     enddo ! over idim1={1,ndims} loop

     return
  end subroutine imat_dump1

!>>> dump a two-dimensional integer matrix to user defined device
  subroutine imat_dump2(myout, ldims, ndims, imat)

     implicit none

! external arguments
! output file unit
     integer, intent(in) :: myout

! leading dimension of the matrix
     integer, intent(in) :: ldims

! second dimension of the matrix
     integer, intent(in) :: ndims

! target matrix to be export
     integer, intent(in) :: imat(ldims, ndims)

! local variables
! loop index over dimensions
     integer :: idim1
     integer :: idim2

! dump imat to user defined device
     do idim2=1,ndims
         do idim1=1,ndims
             if (imat(idim1, idim2) .eq. 0) cycle
             write(myout, '(3I8)') idim1, idim2, imat(idim1, idim2)
         enddo ! over idim1={1,ndims} loop
     enddo ! over idim2={1,ndims} loop

     return
  end subroutine imat_dump2

!>>> dump a three-dimensional integer matrix to user defined device
  subroutine imat_dump3(myout, ldims, mdims, ndims, imat)

     implicit none

! external arguments
! output file unit
     integer, intent(in) :: myout

! leading dimension of the matrix
     integer, intent(in) :: ldims

! second dimension of the matrix
     integer, intent(in) :: mdims

! third dimension of the matrix
     integer, intent(in) :: ndims

! target matrix to be export
     integer, intent(in) :: imat(ldims, mdims, ndims)

! local variables
! loop index over dimensions
     integer :: idim1
     integer :: idim2
     integer :: idim3

! dump imat to user defined device
     do idim3=1,ndims
         do idim2=1,ndims
             do idim1=1,ndims
                 if (imat(idim1, idim2, idim3) .eq. 0) cycle
                 write(myout, '(4I8)') idim1, idim2, idim3, imat(idim1, idim2, idim3)
             enddo ! over idim1={1,ndims} loop
         enddo ! over icol={1,ndims} loop
     enddo ! over irow={1,ndims} loop

     return
  end subroutine imat_dump3

!>>> dump a fout-dimensional integer matrix to user defined device
  subroutine imat_dump4(myout, ldims, mdims, ndims, kdims, imat)

     implicit none

! external arguments
! output file unit
     integer, intent(in) :: myout

! leading dimension of the matrix
     integer, intent(in) :: ldims

! second dimension of the matrix
     integer, intent(in) :: mdims

! third dimension of the matrix
     integer, intent(in) :: ndims

! fourth dimension of the matrix
     integer, intent(in) :: kdims

! target matrix to be export
     integer, intent(in) :: imat(ldims, mdims, ndims, kdims)

! local variables
! loop index over dimensions
     integer :: idim1
     integer :: idim2
     integer :: idim3
     integer :: idim4

! dump imat to user defined device
     do idim4=1,kdims
         do idim3=1,ndims
             do idim2=1,mdims
                 do idim1=1,ldims
                     if (imat(idim1, idim2, idim3, idim4) .eq. 0) cycle
                     write(myout, '(5I8)') idim1, idim2, idim3, idim4, imat(idim1, idim2, idim3, idim4)
                 enddo ! over idim1={1,ldims} loop
             enddo ! over idim2={1,mdims} loop
         enddo ! over idim3={1,ndims} loop
     enddo ! over idim4={1,kdims} loop

     return
  end subroutine imat_dump4

!>>> compute the determinant of orginal matrix
  subroutine zmat_deter( ndim, zmat, deter )

     implicit none

! external arguments
! dimension of the input square matrix
     integer, intent(in) :: ndim

! orginal matrix to compute determinant, on exit, it's replaced by L and U
     complex(8), intent(inout) :: zmat(ndim, ndim)

! determinant of the orginal matrix dmat
     complex(8), intent(out) :: deter

! local variables
! the pivot index, row i of the matrix was interchanged with row ipiv(i)
     integer :: ipiv(ndim)

! return information from subroutine dgetrf
     integer :: info

! determinant of original matrix. determinant = det(1) * ten**det(2)
     complex(8) :: det(2)

! auxiliary constant, ten = 10.d0.
     complex(8) :: ten

! loop index
     integer :: i

! compute an LU factorization using partial pivoting with row interchanges
! A = P * L * U, where P is a permutation matrix, L is lower triangular 
! with unit diagonal diagonal element, U is upper triangular
     call zgetrf( ndim, ndim, zmat, ndim, ipiv, info )
     if ( info .ne. 0 ) then
         print*, 'info=', info
         stop 'severe error happened in subroutine zmat_deter'
     endif ! back if ( info /= 0 ) blcok

! computer determinant by just multiply the diagonal elements
! special trick to avoid numerical precision loss of determinant value
     ten    = dcmplx(10.0d0, 0.0d0)
     det(1) = dcmplx(1.00d0, 0.0d0)
     det(2) = dcmplx(0.00d0, 0.0d0)

     do i=1,ndim
         if( ipiv(i) /= i ) det(1) = -det(1)
         det(1) = zmat(i,i)*det(1)
         if ( det(1) .eq. 0.0d0 ) exit
         do while ( cdabs(det(1)) < 1.0d0 )
             det(1) = ten * det(1)
             det(2) = det(2) - dcmplx(1.0d0, 0.0d0)
         enddo ! over while ( cdabs(det(1)) < 1.0d0 ) loop
         do while ( cdabs(det(1)) > 10.d0 )
             det(1) = det(1) / ten
             det(2) = det(2) + dcmplx(1.0d0, 0.0d0)
         enddo ! over while ( cdabs(det(1)) > 10.d0 ) loop
     enddo ! over i={1,ndim} loop

     deter = det(1) * ten**det(2)

     return
  end subroutine zmat_deter
!>>> dump a one-dimensional complex matrix to user defined device
  subroutine zmat_dump1(myout, ndims, zmat)

     implicit none

! external arguments
! output file unit
     integer, intent(in) :: myout

! leading dimension of the matrix
     integer, intent(in) :: ndims

! target matrix to be export
     complex(8), intent(in) :: zmat(ndims)

! local variables
! loop index over dimensions
     integer :: idim1

! dump zmat to user defined device
     do idim1=1,ndims
         if (abs(zmat(idim1)) .lt. 1.0D-9) cycle
         write(myout, '(I8, 2F17.10)') idim1, zmat(idim1)
     enddo ! over idim1={1,ndims} loop

     return
  end subroutine zmat_dump1

!>>> dump a two-dimensional complex matrix to user defined device
  subroutine zmat_dump2(myout, ldims, ndims, zmat)

     implicit none

! external arguments
! output file unit
     integer, intent(in) :: myout

! leading dimension of the matrix
     integer, intent(in) :: ldims

! second dimension of the matrix
     integer, intent(in) :: ndims

! target matrix to be export
     complex(8), intent(in) :: zmat(ldims, ndims)

! local variables
! loop index over dimensions
     integer :: idim1
     integer :: idim2

! dump zmat to user defined device
     do idim2=1,ndims
         do idim1=1,ldims
             if (abs(zmat(idim1, idim2)) .lt. 1.0D-9) cycle
             write(myout, '(2I8, 2F17.10)') idim1, idim2, zmat(idim1, idim2)
         enddo ! over idim1={1,ldims} loop
     enddo ! over idim2={1,ndims} loop

     return
  end subroutine zmat_dump2

!>>> dump a three-dimensional complex matrix to user defined device
  subroutine zmat_dump3(myout, ldims, mdims, ndims, zmat)

     implicit none

! external arguments
! output file unit
     integer, intent(in) :: myout

! leading dimension of the matrix
     integer, intent(in) :: ldims

! second dimension of the matrix
     integer, intent(in) :: mdims

! third dimension of the matrix
     integer, intent(in) :: ndims

! target matrix to be export
     complex(8), intent(in) :: zmat(ldims, mdims, ndims)

! local variables
! loop index over dimensions
     integer :: idim1
     integer :: idim2
     integer :: idim3

! dump zmat to user defined device
     do idim3=1,ndims
         do idim2=1,mdims
             do idim1=1,ldims
                 if (abs(zmat(idim1, idim2, idim3)) .lt. 1.0D-9) cycle
                 write(myout, '(3I8, 2F17.10)') idim1, idim2, idim3, zmat(idim1, idim2, idim3)
             enddo ! over idim1={1,ldims} loop
         enddo ! over icol={1,mdims} loop
     enddo ! over irow={1,ndims} loop

     return
  end subroutine zmat_dump3

!>>> dump a four-dimensional complex matrix to user defined device
  subroutine zmat_dump4(myout, ldims, mdims, ndims, kdims, zmat)

     implicit none

! external arguments
! output file unit
     integer, intent(in) :: myout

! leading dimension of the matrix
     integer, intent(in) :: ldims

! second dimension of the matrix
     integer, intent(in) :: mdims

! third dimension of the matrix
     integer, intent(in) :: ndims

! fourth dimension of the matrix
     integer, intent(in) :: kdims

! target matrix to be export
     complex(8), intent(in) :: zmat(ldims, mdims, ndims, kdims)

! local variables
! loop index over dimensions
     integer :: idim1
     integer :: idim2
     integer :: idim3
     integer :: idim4

! dump zmat to user defined device
     do idim4=1,kdims
         do idim3=1,ndims
             do idim2=1,mdims
                 do idim1=1,ldims
                     if (abs(zmat(idim1, idim2, idim3, idim4)) .lt. 1.0D-9) cycle
                     write(myout, '(4I8, 2F17.10)') idim1, idim2, idim3, idim4, zmat(idim1, idim2, idim3, idim4)
                 enddo ! over idim1={1,ldims} loop
             enddo ! over idim2={1,mdims} loop
         enddo ! over idim3={1,ndims} loop
     enddo ! over idim4={1,kdims} loop

     return
  end subroutine zmat_dump4

!>>> compute the inverse of a complex matrix
  subroutine zmat_inver(ndim, zmat)

     implicit none

! external arguments
! dimension of the square matrix A
     integer, intent(in) :: ndim

! on entry, the ndim-by-ndim matrix will be inversed
! on exit , the inverse of the orginal matrix zmat
     complex(8), intent(inout) :: zmat(ndim, ndim)

! local variables
! return flag, (info = 0) => successful exit
     integer :: info

! the pivot indices, row i of the matrix was interchanged with row ipiv(i)
     integer :: ipiv(ndim)

! complex array used as workspace
     complex(8) :: work(ndim)

! compute an LU factorization using partial pivoting with row interchanges
! A = P * L * U, where P is a permutation matrix, L is lower triangular 
! with unit diagonal diagonal element, U is upper triangular
     call zgetrf( ndim, ndim, zmat, ndim, ipiv, info )
     if ( info /= 0 ) then
         stop 'severe error happened in subroutine zgetrf'
     endif ! back if ( info /= 0 ) block

! compute the inverse by using LU fractorization computed by zgetrf
     call zgetri( ndim, zmat, ndim, ipiv, work, ndim, info )
     if ( info /= 0 ) then
         stop 'severe error happened in subroutine zgetri'
     endif ! back if ( info /= 0 ) block

     return
  end subroutine zmat_inver

  subroutine zmat_zgeev( ldim, ndim, zmat, zeig, zvec )
     implicit none

! external variables
! leading dimension of matrix amat
     integer, intent(in) :: ldim

! the order of the matrix amat
     integer, intent(in) :: ndim

! original real symmetric matrix to compute eigenval and eigenvector
     complex(8), intent(in) :: zmat(ldim, ndim)

! if info = 0, the eigenvalues in ascending order.
     complex(8), intent(out) :: zeig(ndim)

! if info = 0, orthonormal eigenvectors of the matrix A
     complex(8), intent(out) :: zvec(ldim, ndim)

! local variables
! status flag
     integer :: istat

! return information from subroutine dysev
     integer :: info

! the length of the array work, lwork >= max(1,4*ndim)
     integer :: lwork

! workspace array
     complex(8), allocatable :: work(:)

! auxiliary real(dp) matrix
     complex(8), allocatable :: rwork(:)

     complex(8), allocatable :: vr(:, :)
     complex(8), allocatable :: vl(:, :)

! initialize lwork and allocate memory fo array work
     lwork = 2*ndim
     allocate( work(lwork), stat=istat)
     allocate(rwork(lwork), stat=istat)
     allocate(vr(ndim, ndim), stat=istat)
     allocate(vl(ndim, ndim), stat=istat)
     if ( istat /= 0 ) then
         stop "allocate memory error in dmat_dsyev"
     endif ! back if ( istat /= 0 ) block

! initialize output arrays
     zeig = 0.d0
     zvec = zmat

     call ZGEEV('N', 'V', ndim, zvec, ldim, zeig, &
                vl, ndim, vr, ndim, work, lwork, rwork, info)
     if (info /= 0) stop "Failure in subroutine zmat_zgeev"
     zvec = vr

! dealloate memory for workspace array
     if (allocated(vr  )) deallocate(vr  )
     if (allocated(vl  )) deallocate(vl  )
     if (allocated(work)) deallocate(work)
     if (allocated(rwork)) deallocate(rwork)

     return
  end subroutine zmat_zgeev
!>>> zgemm performs one of the matrix-matrix operations <<<!
! C = alpha*op( A )*op( B ) + beta*C, where op( X ) is one of
! 'N' => op(A) = A; 'T' => op(A) = A'; 'C' => op(A) = conjg(A')
  subroutine zmat_zgemm0( ndim, amat, bmat, cmat )

     implicit none

! dimension of the input square matrix 'amat and bmat'
     integer, intent(in)  :: ndim

! input square matrix 'amat and bmat'
     complex(8), intent(in)  :: amat(ndim, ndim)
     complex(8), intent(in)  :: bmat(ndim, ndim)

! output square matrix, cmat = amat * bmat
     complex(8), intent(out) :: cmat(ndim, ndim)

! local variables
     complex(8) :: alpha
     complex(8) :: betta

     cmat  = dcmplx(0.0d0, 0.0d0)
     alpha = dcmplx(1.0d0, 0.0d0)
     betta = dcmplx(0.0d0, 0.0d0)

     call zgemm('N', 'N', ndim, ndim, ndim, alpha, amat, ndim,& 
                                bmat, ndim, betta, cmat, ndim)

     return
  end subroutine zmat_zgemm0

!>>> zgemm performs one of the matrix-matrix operations <<<!
! C = alpha*op( A )*op( B ) + beta*C, where op( X ) is one of
! 'N' => op(A) = A; 'T' => op(A) = A'; 'C' => op(A) = conjg(A')
  subroutine zmat_zgemm1( ndim, amat, bmat, cmat )

     implicit none

! dimension of the input square matrix 'amat and bmat'
     integer, intent(in)  :: ndim

! input square matrix 'amat and bmat'
     complex(8), intent(in)  :: amat(ndim, ndim)
     complex(8), intent(in)  :: bmat(ndim, ndim)

! output square matrix, cmat = amat * bmat
     complex(8), intent(out) :: cmat(ndim, ndim)

! local variables
     complex(8) :: alpha
     complex(8) :: betta

     cmat  = dcmplx(0.0d0, 0.0d0)
     alpha = dcmplx(1.0d0, 0.0d0)
     betta = dcmplx(0.0d0, 0.0d0)

     call zgemm('T', 'N', ndim, ndim, ndim, alpha, amat, ndim,& 
                                bmat, ndim, betta, cmat, ndim)

     return
  end subroutine zmat_zgemm1

!>>> zgemm performs one of the matrix-matrix operations <<<!
! C = alpha*op( A )*op( B ) + beta*C, where op( X ) is one of
! 'N' => op(A) = A; 'T' => op(A) = A'; 'C' => op(A) = conjg(A')
  subroutine zmat_zgemm2( ndim, amat, bmat, cmat )

     implicit none

! dimension of the input square matrix 'amat and bmat'
     integer, intent(in)  :: ndim

! input square matrix 'amat and bmat'
     complex(8), intent(in)  :: amat(ndim, ndim)
     complex(8), intent(in)  :: bmat(ndim, ndim)

! output square matrix, cmat = amat * bmat
     complex(8), intent(out) :: cmat(ndim, ndim)

! local variables
     complex(8) :: alpha
     complex(8) :: betta

     cmat  = dcmplx(0.0d0, 0.0d0)
     alpha = dcmplx(1.0d0, 0.0d0)
     betta = dcmplx(0.0d0, 0.0d0)

     call zgemm('C', 'N', ndim, ndim, ndim, alpha, amat, ndim,& 
                                bmat, ndim, betta, cmat, ndim)

     return
  end subroutine zmat_zgemm2
! zheev computes all eigenvalues and optionally, eigenvectors of a
! complex Hermitian matrix A.
  subroutine zmat_zheev( ldim, ndim, amat, eval, evec )
     implicit none

! external variables
! leading dimension of matrix amat
     integer, intent(in) :: ldim

! the order of the matrix amat
     integer, intent(in) :: ndim

! original real symmetric matrix to compute eigenval and eigenvector
     complex(8), intent(in) :: amat(ldim, ndim)

! if info = 0, the eigenvalues in ascending order.
     real(8), intent(out) :: eval(ndim)

! if info = 0, orthonormal eigenvectors of the matrix A
     complex(8), intent(out) :: evec(ldim, ndim)

! local variables
! status flag
     integer :: istat

! return information from subroutine dysev
     integer :: info

! the length of the array work, lwork >= max(1,2*ndim-1)
     integer :: lwork
     integer :: lrwork

! workspace array
     complex(8), allocatable :: work(:)
     real(8)   , allocatable :: rwork(:)

! initialize lwork (lrwork) and allocate memory for array work (rwork)
     lwork = 2*ndim-1
     lrwork = 3*ndim-2
     
     allocate(work(lwork), stat=istat)
     allocate(rwork(lrwork), stat=istat)
     if ( istat /= 0 ) then
         stop "allocate memory error in dmat_dsyev"
     endif ! back if ( istat /= 0 ) block

! initialize output arrays
     eval = 0.d0
     evec = amat

     call ZHEEV('V', 'L', ndim, evec, ldim, eval, work, lwork, rwork, info)
     if (info /= 0) stop "Failure in subroutine zmat_zheev"

! dealloate memory for workspace array
     if (allocated(work )) deallocate(work)
     if (allocated(rwork)) deallocate(rwork)

     return
  end subroutine zmat_zheev

! zheevd computes all eigenvalues and optionally, eigenvectoes of a 
! complex Hermitian matrix A. If eigenvectors are desired, it uses a
! devide and conquer algorithm.
  subroutine zmat_zheevd(ldim, ndim, amat, eval, evec)

! external variables
! leading dimension of matrix amat
     integer, intent(in) :: ldim

! the order of the matrix amat
     integer, intent(in) :: ndim

! original real symmetric matrix to compute eigenval and eigenvector
     complex(8), intent(in) :: amat(ldim, ndim)

! if info = 0, the eigenvalues in ascending order.
     real(8), intent(out) :: eval(ndim)

! if info = 0, orthonormal eigenvectors of the matrix A
     complex(8), intent(out) :: evec(ldim, ndim)

! local variables
! status flag
     integer :: istat

! return information from subroutine dysev
     integer :: info

! the length of the array work, lwork >= max(1,2*ndim-1)
     integer :: lwork
     integer :: liwork
     integer :: lrwork

! workspace array
     integer   , allocatable :: iwork(:)
     real(8)   , allocatable :: rwork(:)
     complex(8), allocatable ::  work(:)

! initialize lwork (lrwork) and allocate memory for array work (rwork)
     lwork = 2*ndim + ndim*ndim
     liwork = 3 + 5*ndim
     lrwork = 1 + 5*ndim + 2*ndim*ndim
     
     allocate(work(lwork), stat=istat)
     allocate(rwork(lrwork), stat=istat)
     allocate(iwork(liwork), stat=istat)
     if ( istat /= 0 ) then
         stop "allocate memory error in dmat_dsyev"
     endif ! back if ( istat /= 0 ) block

! initialize output arrays
     eval = 0.d0
     evec = amat

     call zheevd('V', 'L', ndim, evec, ldim, eval, work, lwork, &
                 rwork, lrwork, iwork, liwork, info)
     if (info /= 0) stop "Failure in subroutine zmat_zheev"

! dealloate memory for workspace array
     if (allocated(work )) deallocate(work)
     if (allocated(rwork)) deallocate(rwork)
     if (allocated(iwork)) deallocate(iwork)

     return
  end subroutine zmat_zheevd
!>>> computes all eigenvalues and eigenvectors of complex generalized 
! Hermitian-definite eigenproblem of the form: A*x=(lambda)*B*x
  subroutine zmat_zhegv( ldim, ndim, amat, bmat, eval, evec )
     implicit none

! external variables
! leading dimension of matrix amat and bmat
     integer, intent(in) :: ldim

! the order of the matrix amat and bmat
     integer, intent(in) :: ndim

! original Hermitian matrix A
     complex(8), intent(in) :: amat(ldim, ndim)

! original Hermitian positive definite matrix B
     complex(8), intent(in) :: bmat(ldim, ndim)

! if info = 0, the eigenvalues in ascending order.
     real(8), intent(out) :: eval(ndim)

! if info = 0, orthonormal eigenvectors of the matrix A
     complex(8), intent(out) :: evec(ldim, ndim)

! local variables
! status flag
     integer :: istat

! return information from subroutine dysev
     integer :: info

! the length of the array work, lwork >= max(1,2*ndim-1)
     integer :: lwork
     integer :: lrwork

! real(dp) workspace array
     real(8), allocatable :: rwork(:)

! complex(dp) workspace array
     complex(8), allocatable :: work(:)
     complex(8), allocatable :: zmat(:, :)

! initialize lwork (lrwork) and allocate memory for array work (rwork)
     lwork = 2*ndim-1
     lrwork = 3*ndim-2
     
! allocate memory for work arrays
     allocate(work(lwork), stat=istat)
     allocate(rwork(lrwork), stat=istat)
     allocate(zmat(ldim, ndim), stat=istat)
     if ( istat /= 0 ) then
         stop "allocate memory error in dmat_dsyev"
     endif ! back if ( istat /= 0 ) block

! initialize output arrays
     eval = 0.d0
     evec = amat
     zmat = bmat

     call ZHEGV(1,'V', 'L', ndim, evec, ldim, zmat, ldim, eval, work, lwork, rwork, info)
     if (info /= 0) stop "Failure in subroutine zmat_zhegv"

! dealloate memory for workspace array
     if (allocated(work )) deallocate(work)
     if (allocated(zmat )) deallocate(zmat)
     if (allocated(rwork)) deallocate(rwork)

     return
  end subroutine zmat_zhegv

