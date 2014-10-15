!!!---------------------------------------------------------------
!!! project : maxent
!!! program : maxent_dsyev
!!!           maxent_dgesvd 
!!!           maxent_make_hist
!!! source  : maxent_utils.f90
!!! type    : subroutine
!!! author  : yilin wang (email: qhwyl2006@126.com)
!!! history : 05/31/2013 by yilin wang
!!!         : 10/14/2014 by yilin wang
!!! purpose : define some utility subroutines for maxent 
!!! status  : unstable
!!! comment :
!!!---------------------------------------------------------------

!!>>> maxent_dsyev: diagonalize a symmetry real matrix, 
!!>>> call lapack subroutine dsyev  
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

!!>>> maxent_dgesvd: make the singular values decomposition
  subroutine maxent_dgesvd(m, n, min_mn, amat, umat, sigvec, vmatt, info)
     use constants, only : dp
 
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

!!>>> maxent_make_hist: make a histogram for data bins
  subroutine maxent_make_hist( grnbin, grn_ave, hmesh, hist )
     use constants, only : dp, two, epss
     use control, only : ntime, nbins, slice
 
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
         if ( width > epss ) then
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
