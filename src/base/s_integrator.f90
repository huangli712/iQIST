!!!-----------------------------------------------------------------------
!!! project : CSSL (Common Service Subroutines Library)
!!! program : s_int_trapezoid
!!!           s_int_simpson
!!! source  : s_integrator.f90
!!! type    : functions
!!! author  : li huang (email:lihuang.dmft@gmail.com)
!!! history : 09/20/2014 by li huang (created)
!!!           05/12/2017 by li huang (last modified)
!!! purpose : these functions are used to do numerical integration with
!!!           composite trapezoid or composite simpson algorithms.
!!! status  : unstable
!!! comment :
!!!-----------------------------------------------------------------------

!!
!!
!! Introduction
!! ============
!!
!! To use s_int_trapezoid() or s_int_simpson(), you have to define the
!! integrand at first. For example:
!!
!!  function f(x)
!!     use constants, only : dp
!!
!!     implicit none
!!
!!     real(dp) :: x
!!     real(dp) :: f
!!
!!     f = x * x
!!  end function f
!!
!! Next, you have to determine the lower bound a and upper bound b, and
!! the number of points n. Noted that now both the s_int_trapezoid() and
!! s_int_simpson() functions only support the 1-D numerical integration.
!!
!! 1. use s_int_trapezoid()
!! ------------------------
!!
!! procedure( real(dp) ) :: s_int_trapezoid
!! procedure( real(dp) ) :: f
!! real(dp) :: val
!!
!! val = s_int_trapezoid(f, a, b, n)
!!
!! 2. use s_int_simpson()
!! ----------------------
!!
!! procedure( real(dp) ) :: s_int_simpson
!! procedure( real(dp) ) :: f
!! real(dp) :: val
!!
!! val = s_int_simpson(f, a, b, n)
!!
!!

!!
!! @fun s_int_trapezoid
!!
!! numerical integration with trapezoid algorithm
!!
  function s_int_trapezoid(f, a, b, n) result(val)
     use constants, only : dp, zero, two

     implicit none

! external arguments
! number of data points
     integer, intent(in)  :: n

! boundries for numerical integration
     real(dp), intent(in) :: a
     real(dp), intent(in) :: b

! external function, it means the integrand
     real(dp), external   :: f

! local variables
! loop index
     integer  :: i

! return value
     real(dp) :: val

! step for integration
     real(dp) :: h

! sum for trapezoid rule
     real(dp) :: trapSum

! evaluate the step
     h = ( b - a ) / dble(n)

! calculate trapezoid sum
     trapSum = zero
     do i=1,n-1
         trapSum = trapSum + f(a+dble(i)*h)
     enddo ! over i={1,n-1} loop

! calculate the final value
     val = ( h / two ) * ( f(a) + f(b) + two*trapSum )

     return
  end function s_int_trapezoid

!!
!! @fun s_int_simpson
!!
!! numerical integration with simpson algorithm
!!
  function s_int_simpson(f, a, b, n) result(val)
     use constants, only : dp, zero

     implicit none

! external arguments
! number of data points
     integer, intent(in)  :: n

! boundries for numerical integration
     real(dp), intent(in) :: a
     real(dp), intent(in) :: b

! external function, it means the integrand
     real(dp), external   :: f

! local variables
! loop index
     integer  :: i

! return value
     real(dp) :: val

! step for integration
     real(dp) :: h

! sum for trapezoid rule
     real(dp) :: oddSum
     real(dp) :: evenSum

! evaluate the step
     h = ( b - a ) / dble(n)

! calculate simpson sum
     evenSum = zero
     oddSum = zero
     do i=1,n-1
         if ( mod(i,2) == 0 ) then
             evenSum = evenSum + f(a+dble(i)*h)
         else
             oddSum = oddSum + f(a+dble(i)*h)
         endif ! back if ( mod(i,2) == 0 ) block
     enddo ! over i={1,n-1} loop

! calculate the final value
     val = ( h / 3.0_dp ) * ( f(a) + f(b) + 2.0_dp * evenSum + 4.0_dp * oddsum )

     return
  end function s_int_simpson
