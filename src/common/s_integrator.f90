!!!-----------------------------------------------------------------------
!!! project : CSSL (Common Service Subroutines Library)
!!! program : s_int_trapezoid
!!!           s_int_simpson
!!! source  : s_integrator.f90
!!! type    : subroutines
!!! author  : li huang (email:huangli712@gmail.com)
!!! history : 09/20/2014 by li huang
!!! purpose : these subroutines are used to do numerical integration with
!!!           composite trapezoid or composite simpson algorithms.
!!! status  : unstable
!!! comment :
!!!-----------------------------------------------------------------------

  function s_int_trapezoid(f, a, b, n) result(val)
     use constants

     implicit none

     double precision, intent(in) :: a, b
     integer, intent(in) :: n
     double precision :: f
     double precision :: h, trapSum
     integer :: i

     real(dp) :: val

     h = (b-a) / dble(n)

     trapSum = 0.0d0

     do i=1, n-1
         trapSum = trapSum + f(a+dble(i)*h)
     enddo

     val = (h/2.0d0) * ( f(a) + f(b) + 2.0d0*trapSum )

     return 
  end function s_int_trapezoid

  function s_int_simpson(f, a, b, n) result(val)
     use constants
     implicit none

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
