!!!-----------------------------------------------------------------------
!!! project : CSSL (Common Service Subroutines Library)
!!! program : s_int_trapezoid
!!!           s_int_simpson
!!! source  : s_integrator.f90
!!! type    : functions
!!! author  : li huang (email:huangli712@gmail.com)
!!! history : 09/20/2014 by li huang
!!! purpose : these subroutines are used to do numerical integration with
!!!           composite trapezoid or composite simpson algorithms.
!!! status  : unstable
!!! comment :
!!!-----------------------------------------------------------------------

!!>>> s_int_trapezoid: numerical integration with trapezoid algorithm
  function s_int_trapezoid(f, a, b, n) result(val)
     use constants, only : dp, zero, two

     implicit none

! external arguments
! number of data points
     integer, intent(in)  :: n

! boundries for numerical integration
     real(dp), intent(in) :: a
     real(dp), intent(in) :: b

! external function, it means the integrator
     real(dp) :: f

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
     h = (b-a) / dble(n)

! calculate trapezoid sum
     trapSum = zero
     do i=1,n-1
         trapSum = trapSum + f(a+dble(i)*h)
     enddo ! over i={1,n-1} loop

! calculate the final value
     val = ( h / two ) * ( f(a) + f(b) + two*trapSum )

     return 
  end function s_int_trapezoid

!!>>> s_int_simpson: numerical integration with simpson algorithm
  function s_int_simpson(f, a, b, n) result(val)
     use constants

     implicit none

! external argumenst
     double precision, intent(in) :: a,b
     integer, intent(in) :: n
     double precision :: h, oddSum, evenSum
     double precision :: f
     integer :: i

     real(dp) :: val
     h = (b-a) / dble(n)

     do i=1, n-1
         if ( mod(i,2) == 0 ) then
             evenSum = evenSum + f(a+dble(i)*h)
         else
             oddSum = oddSum + f(a+dble(i)*h)
         endif
     enddo

     val = (h/3.0d0) * ( f(a) + f(b) + 2.0d0*evenSum + 4.0d0*oddsum )

     return
  end function s_int_simpson
