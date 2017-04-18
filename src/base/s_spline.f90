!!!-----------------------------------------------------------------------
!!! project : CSSL (Common Service Subroutines Library)
!!! program : s_spl_deriv1
!!!           s_spl_deriv2
!!!           s_spl_funct
!!! source  : s_spline.f90
!!! type    : subroutines
!!! author  : li huang (email:lihuang.dmft@gmail.com)
!!! history : 07/10/2014 by li huang (created)
!!!           04/19/2017 by li huang (last modified)
!!! purpose : these subroutines are used to do cubic spline interpolation.
!!! status  : unstable
!!! comment :
!!!-----------------------------------------------------------------------

!!
!!
!! Introduction
!! ============
!!
!! 1. calculate 1-order derivates for a given function
!! ---------------------------------------------------
!!
!! subroutine s_spl_deriv1(...)
!!
!! 2. calculate 2-order derivates for a given function
!! ---------------------------------------------------
!!
!! subroutine s_spl_deriv2(...)
!!
!! 3. evaluate function value at a given point
!! -------------------------------------------
!!
!! function   s_spl_funct(...)
!!
!!

!!
!! @sub s_spl_deriv1
!!
!! evaluate the 1-order derivates of yval
!!
  subroutine s_spl_deriv1(ydim, xval, yval, d1y)
     use constants, only : dp

     implicit none

! external arguments
! dimension of xval and yval
     integer, intent(in)   :: ydim

! old knots
     real(dp), intent(in)  :: xval(ydim)

! old function values to be interpolated
     real(dp), intent(in)  :: yval(ydim)

! 1-order derivates
     real(dp), intent(out) :: d1y(ydim)

! local variables
! loop index
     integer :: i

! compute derivation using two-point formula
     do i=2,ydim-1
         d1y(i) = ( yval(i+1) - yval(i-1) ) / ( xval(i+1) - xval(i-1) )
     enddo ! over i={2,ydim-1} loop

! compute first and last derivation using linear extrapolation
     d1y(1) = d1y(2) + ( d1y(3) - d1y(2) ) / &
                     ( xval(3) - xval(2) ) * &
                     ( xval(1) - xval(2) )
     d1y(ydim) = d1y(ydim-1) + ( d1y(ydim-1) - d1y(ydim-2) ) / &
                             ( xval(ydim-1) - xval(ydim-2) ) * &
                             ( xval(ydim-0) - xval(ydim-1) )

     return
  end subroutine s_spl_deriv1

!!
!! @sub s_spl_deriv2
!!
!! evaluate the 2-order derivates of yval
!!
  subroutine s_spl_deriv2(ydim, xval, yval, startu, startd, d2y)
     use constants, only : dp, zero, one, two, half

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
     endif ! back if ( startu > .99E30 ) block

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
     endif ! back if ( startd > .99E30 ) block

     d2y(ydim) = ( un - qn * u(ydim-1) ) / ( qn * d2y(ydim-1) + one )

     do k=ydim-1,1,-1
         d2y(k) = d2y(k) * d2y(k+1) + u(k)
     enddo ! over k={ydim-1,1} loop

     return
  end subroutine s_spl_deriv2

!!
!! @fun s_spl_funct
!!
!! evaluate the spline value at x point
!!
  function s_spl_funct(xdim, xval, yval, d2y, x) result(val)
     use constants, only : dp

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

! note: we do not need to check khi here, since x can not reach right
! boundary and left boundary either all
!<     if ( khi > xdim ) then
!<         klo = xdim - 1
!<         khi = xdim
!<     endif ! back if ( khi > xdim ) block

! calculate splined parameters a and b
     a = ( xval(khi) - x ) / h
     b = ( x - xval(klo) ) / h

! spline it, obtain the fitted function value at x point
     val = a * yval(klo) + b * yval(khi) + &
               ( ( a*a*a - a ) * d2y(klo) + ( b*b*b - b ) * d2y(khi) ) * &
               ( h*h ) / 6.0_dp

     return
  end function s_spl_funct
