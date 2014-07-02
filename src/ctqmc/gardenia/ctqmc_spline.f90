!-------------------------------------------------------------------------
! project : gardenia
! program : ctqmc_make_htau
!           ctqmc_make_hsed
!           ctqmc_make_spline
!           ctqmc_make_splint
! source  : ctqmc_spline.f90
! type    : subroutines
! author  : li huang (email:huangli712@gmail.com)
! history : 05/05/2008 by li huang
!           02/18/2009 by li huang
!           09/23/2009 by li huang
!           09/26/2009 by li huang
!           10/03/2009 by li huang
!           11/10/2009 by li huang
!           12/18/2009 by li huang
! purpose : to provide cubic spline subroutines and wrapper functions to
!           interpolate the hybridization function in imaginary-time axis
! input   :
! output  :
! status  : unstable
! comment :
!-------------------------------------------------------------------------

!>>> evaluate the matrix elements for mmat matrix using cubic spline interpolation
  function ctqmc_make_htau(flvr, dtau) result(val)
     use constants
     use control
     use context

     implicit none

! external arguments
! current flavor channel
     integer, intent(in)  :: flvr

! delta imaginary time
     real(dp), intent(in) :: dtau

! external functions
! internal interpolation engine
     real(dp), external :: ctqmc_make_splint

! local variables
! return value
     real(dp) :: val

     val = ctqmc_make_splint(ntime, tmesh, htau(:, flvr, flvr), hsed(:, flvr, flvr), dtau)

     return
  end function ctqmc_make_htau

!>>> calculate the second order derivates of hybridization function on
! imaginary time space
  subroutine ctqmc_make_hsed(tmesh, htau, hsed)
     use constants
     use control

     implicit none

! external arguments
! imaginary time axis
     real(dp), intent(in)  :: tmesh(ntime)

! hybridization function on imaginary time axis
     real(dp), intent(in)  :: htau(ntime,norbs,norbs)

! second order derivates of hybridization function
     real(dp), intent(out) :: hsed(ntime,norbs,norbs)

! local variables
! loop index
     integer  :: i
     integer  :: j

! first derivate at start point
     real(dp) :: startu

! first derivate at end   point
     real(dp) :: startd

! \delta \tau
     real(dp) :: deltau

! second-order derivates
     real(dp) :: d2y(ntime)

! calculate deltau
     deltau = beta / real(ntime - 1)

! initialize hsed
     hsed = zero

! calculate it
     do j=1,norbs
         do i=1,norbs

! calculate first-order derivate of \Delta(0): startu
             startu = (-25.0_dp*htau(1,       i, j) +                    &
                        48.0_dp*htau(2,       i, j) -                    &
                        36.0_dp*htau(3,       i, j) +                    &
                        16.0_dp*htau(4,       i, j) -                    &
                         3.0_dp*htau(5,       i, j)) / 12.0_dp / deltau

! calculate first-order derivate of \Delta(\beta): startd
             startd = ( 25.0_dp*htau(ntime-0, i, j) -                    &
                        48.0_dp*htau(ntime-1, i, j) +                    &
                        36.0_dp*htau(ntime-2, i, j) -                    &
                        16.0_dp*htau(ntime-3, i, j) +                    &
                         3.0_dp*htau(ntime-4, i, j)) / 12.0_dp / deltau

! reinitialize d2y to zero
             d2y = zero

! call the service layer
             call ctqmc_make_spline(ntime, tmesh, htau(:,i,j), startu, startd, d2y)

! copy the results to hsed
             hsed(:,i,j) = d2y

         enddo ! over i={1,norbs} loop
     enddo ! over j={1,norbs} loop

     return
  end subroutine ctqmc_make_hsed

!>>> ctqmc_make_spline: evaluate the 2-order derivates of yval
  subroutine ctqmc_make_spline(ydim, xval, yval, startu, startd, d2y)
     use constants

     implicit none

! external arguments
! dimension of xval and yval
     integer, intent(in)   :: ydim

! first-derivate at point 1
     real(dp), intent(in)  :: startu

! first-derivate at point ydim
     real(dp), intent(in)  :: startd

! old knots
     real(dp), intent(in)  :: xval(ydim)

! old function values to be interpolated
     real(dp), intent(in)  :: yval(ydim)

! 2-order derivates
     real(dp), intent(out) :: d2y(ydim)

! local variables
! loop index
     integer  :: i
     integer  :: k

! dummy variables
     real(dp) :: p
     real(dp) :: qn
     real(dp) :: un
     real(dp) :: sig

! dummy arrays
     real(dp) :: u(ydim)

! deal with left boundary
     if ( startu > .99E30 ) then
         d2y(1) = zero
         u(1) = zero
     else
         p = xval(2) - xval(1)
         d2y(1) = -half
         u(1) = ( 3.0_dp / p ) * ( ( yval(2) - yval(1) ) / p - startu )
     endif

     do i=2,ydim-1
         sig    = ( xval(i) - xval(i-1) ) / ( xval(i+1) - xval(i-1) )
         p      = sig * d2y(i- 1) + two
         d2y(i) = ( sig - one ) / p
         u(i)   = ( 6.0_dp * ( ( yval(i+1) - yval(i) ) / &
                    ( xval(i+1) - xval(i) ) - ( yval(i) - yval(i-1) ) / &
                    ( xval(i) - xval(i-1) ) ) / &
                    ( xval(i+1) - xval(i-1) ) - sig * u(i-1) ) / p
     enddo ! over i={2,ydim-1} loop

! deal with right boundary
     if ( startd > .99E30 ) then
         qn = zero
         un = zero
     else
         p = xval(ydim) - xval(ydim-1)
         qn = half
         un = ( 3.0_dp / p ) * ( startd - ( yval(ydim) - yval(ydim-1) ) / p )
     endif

     d2y(ydim) = ( un - qn * u(ydim-1) ) / ( qn * d2y(ydim-1) + one )

     do k=ydim-1,1,-1
         d2y(k) = d2y(k) * d2y(k+1) + u(k)
     enddo ! over k={ydim-1,1} loop

     return
  end subroutine ctqmc_make_spline

!>>> ctqmc_make_splint: evaluate the spline value at x point
  function ctqmc_make_splint(xdim, xval, yval, d2y, x) result(val)
     use constants

     implicit none

! external arguments
! dimension of xval and yval
     integer, intent(in)  :: xdim

! new mesh point
     real(dp), intent(in) :: x

! old mesh
     real(dp), intent(in) :: xval(xdim)

! old function value
     real(dp), intent(in) :: yval(xdim)

! 2-order derviates of old function
     real(dp), intent(in) :: d2y(xdim)

! local variables
! lower boundary
     integer  :: khi

! higher boundary
     integer  :: klo

! distance between two successive mesh points
     real(dp) :: h

! dummy variables
     real(dp) :: a, b

! return value
     real(dp) :: val

! calculate the interval of spline zone
     h = xval(2) - xval(1)

! special trick is adopted to determine klo and kho
     klo = floor(x/h) + 1
     khi = klo + 1

! note: we do not need to check khi here, since x can not reach right boundary
! and left boundary either all
!<     if ( khi > xdim ) then
!<         klo = xdim - 1
!<         khi = xdim
!<     endif

! calculate splined parameters a and b
     a = ( xval(khi) - x ) / h
     b = ( x - xval(klo) ) / h

! spline it, obtain the fitted function value at x point
     val = a * yval(klo) + b * yval(khi) + &
               ( ( a*a*a - a ) * d2y(klo) + ( b*b*b - b ) * d2y(khi) ) * &
               ( h*h ) / 6.0_dp

     return
  end function ctqmc_make_splint
