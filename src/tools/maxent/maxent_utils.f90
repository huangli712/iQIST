!---------------------------------------------------------------
! project : maxent
! program : maxent_dsyev
!         : maxent_dgesvd
!         : maxent_dgemm
!         : maxent_dgemv
!         : maxent_make_hist
! source  : maxent_utils.f90
! type    : subroutine
! author  : yilin wang (email: qhwyl2006@126.com)
! history : 05/31/2013 by yilin wang
! purpose : define some utility subroutines for maxent 
! input   :
! output  :
! status  : unstable
! comment :
!---------------------------------------------------------------

!===============================================================
! subroutine: maxent_dsyev 
! purpose   : diagonalize a symmetry real matrix, 
!           : call lapack subroutine dsyev  
!===============================================================
  subroutine maxent_dsyev(ndim, matrix, eigval, info)
      use constants
 
      implicit none

! external variables
! the dimension of the matrix
      integer, intent(in) :: ndim 

! the symmetry matrix to be diagonalized
      real(dp), intent(inout) :: matrix(ndim, ndim)

! the eigen values of matrix
      real(dp), intent(out) :: eigval(ndim) 

! the error infomation from the dsyev call
      integer, intent(out) :: info 

! local variables
      integer :: lwork 

      real(dp) :: work(5*ndim)

      lwork = 5 * ndim

! call lapack subroutine dysev
      call dsyev( 'V', 'U', ndim, matrix, ndim, eigval, work, lwork, info )

      return
  end subroutine maxent_dsyev

!===============================================================
! subroutine: maxent_dgemm 
! purpose   : multiply two real matrix, 
!           : call lapack subroutine dgemm  
!===============================================================
  subroutine maxent_dgemm(ndim1, ndim2, mat1, ndim3, mat2, mat3)
      use constants
 
      implicit none

! external variables
! the dimension of the matrix
      integer, intent(in) :: ndim1
      integer, intent(in) :: ndim2
      integer, intent(in) :: ndim3

! the matrixs to be multiplied
      real(dp), intent(in) :: mat1(ndim1, ndim2)
      real(dp), intent(in) :: mat2(ndim2, ndim3)
      real(dp), intent(inout) :: mat3(ndim1, ndim3)

! local variables

! call lapack subroutine dgemm
      call dgemm('N','N',ndim1,ndim3,ndim2,one,mat1,ndim1,mat2,ndim2,zero,mat3,ndim1)

      return
  end subroutine maxent_dgemm

!===============================================================
! subroutine: maxent_dgemv
! purpose   : a matrix multiplies a vector 
!           : call lapack subroutine dgemv  
!===============================================================
  subroutine maxent_dgemv(ndim1, ndim2, mat, vec1, vec2)
      use constants
 
      implicit none

! external variables
! the dimension of the matrix
      integer, intent(in) :: ndim1
      integer, intent(in) :: ndim2

! the matrix and the vector
      real(dp), intent(in) :: mat(ndim1, ndim2)
      real(dp), intent(in) :: vec1(ndim2)
      real(dp), intent(out) :: vec2(ndim1)

! local variables

! call lapack subroutine dgemv
      call dgemv('N',ndim1,ndim2,one,mat,ndim1,vec1,1,zero,vec2,1)

      return
  end subroutine maxent_dgemv

!===============================================================
! subtoutine: maxent_dgesvd 
! purpose   : make the singular values decomposition
!===============================================================
  subroutine maxent_dgesvd(m, n, min_mn, amat, umat, sigvec, vmatt, info)
      use constants
 
      implicit none

! external variables
! rows of A matrix
      integer, intent(in) :: m

! columns of A matrix
      integer, intent(in) :: n

! the min value of m and n
      integer, intent(in) :: min_mn

! A matrix
      real(dp), intent(inout) :: amat(m,n)

! the left vectors of svd
      real(dp), intent(out) :: umat(m,min_mn)

! the singular values of svd
      real(dp), intent(out) :: sigvec(min_mn)

! the right vectors of svd
      real(dp), intent(out) :: vmatt(min_mn,n)

! error information from dgesvd
      integer, intent(out) :: info

! local variables
! the dimension of work array
      integer :: lwork 

! work array
      real(dp), allocatable :: work(:)

      lwork = max(3 * min_mn + max(m,n), 5 * min_mn)
      allocate(work(lwork))

      call dgesvd( 'S', 'S', m, n, amat, m, sigvec, umat, m, vmatt, min_mn, work, lwork, info )

      return
  end subroutine maxent_dgesvd

!===============================================================
! subroutine maxent_make_hist is used to make a histogram for
! data bins
!===============================================================
  subroutine maxent_make_hist( grnbin, grn_ave, hmesh, hist )
      use constants
      use control
 
      implicit none

! external variables
      real(dp), intent(in)  :: grnbin(ntime, nbins)
      real(dp), intent(in)  :: grn_ave(ntime)
      real(dp), intent(out) :: hmesh(slice, ntime)
      integer, intent(out)  :: hist(slice, ntime)

! local variables
      real(dp) :: grn(ntime,nbins)
      real(dp) :: width
      real(dp) :: maxv
      real(dp) :: minv

      integer :: itime
      integer :: ibin
      integer :: islice
      integer :: indx

! subtract the average values
      do itime=1, ntime
          grn(itime,:) = grnbin(itime,:) - grn_ave(itime)
      enddo

      hist = 0
      do itime=1, ntime
! determine the max and min value, and width
          maxv = maxval(grn(itime,:))
          minv = minval(grn(itime,:))
          width = ( maxv - minv ) / real(slice)

! build mesh
          do islice=1, slice
              hmesh(islice,itime) = minv + width * ( islice - 1 ) + width / two
          enddo

! build histogram
! if width == 0, there is no histogram, we just set hist(slice/2,itime) = minv
          if ( width > eps12 ) then
              do ibin=1, nbins
                  indx = ceiling( (grn(itime, ibin) - minv) / width )
                  if ( indx == 0 ) then
                      indx = 1
                  endif
                  if ( indx == (slice + 1) ) then
                      indx = slice 
                  endif

                  hist(indx,itime) = hist(indx, itime) + 1
              enddo 
          else 
              hist(slice/2,itime) = nbins
          endif
      enddo

      return
  end subroutine maxent_make_hist

!===============================================================
! subroutine: maxent_time_builder
! purpose   : build time
! author    : lihuang
!===============================================================
  subroutine maxent_time_builder(date_time_string)
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
  end subroutine maxent_time_builder

