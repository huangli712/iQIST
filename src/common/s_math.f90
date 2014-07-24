!!!-----------------------------------------------------------------------
!!! project : CSSL (Common Service Subroutines Library)
!!! program : s_linspace_d
!!!           s_logspace_d
!!!           s_linspace_z
!!! source  : s_math.f90
!!! type    : subroutines
!!! author  : li huang (email:huangli712@gmail.com)
!!! history : 07/10/2014 by li huang
!!! history : 07/24/2014 by li huang
!!! purpose : these subroutines are used to
!!! status  : unstable
!!! comment :
!!!-----------------------------------------------------------------------

!!>>> s_linspace_d: create a linear mesh x in interval [xmin, xmax], real(dp) version
  subroutine s_linspace_d(xmin, xmax, n, x)
     use constants, only : dp

     implicit none

! external arguments
! left boundary
     real(dp), intent(in)  :: xmin

! right boundary
     real(dp), intent(in)  :: xmax

! size of array x
     integer,  intent(in)  :: n

! output array, containing the linear mesh
     real(dp), intent(out) :: x(n)

! local variables
! loop index
     integer :: i

     do i=1,n
         x(i) = ( xmax - xmin ) * real(i - 1, dp) / real(n - 1, dp) + xmin
     enddo ! over i={1,n} loop

     return
  end subroutine s_linspace_d

!!>>> s_logspace_d: create a log mesh x in interval [xmin, xmax], real(dp) version
  subroutine s_logspace_d(xmin, xmax, n, x)
     use constants, only : dp

     implicit none

! external arguments
! left boundary
     real(dp), intent(in)  :: xmin

! right boundary
     real(dp), intent(in)  :: xmax

! size of array x
     integer,  intent(in)  :: n

! output array, containing the linear mesh
     real(dp), intent(out) :: x(n)

! we can use the s_linspace_d() subroutine
     call s_linspace_d(log10(xmin), log10(xmax), x)
     x = 10.0_dp**x

     return
  end subroutine s_logspace_d
