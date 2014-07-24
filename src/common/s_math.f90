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

  subroutine s_linspace_d(xmin, xmax, x)
     use const, only : dp

     implicit none

! external arguments
! left boundary
     real(dp), intent(in)  :: xmin

! right boundary
     real(dp), intent(in)  :: xmax

! output array, containing the linear mesh
     real(dp), intent(out) :: x(:)

! local variables
! loop index
     integer :: i

! size of array x
     integer :: n

! get size of array x
     n = size(x)

     do i=1,n
         x(i) = ( xmax - xmin ) * real(i-1,dp) / real(n-1,dp) + xmin
     enddo ! over i={1,n} loop

     return
  end subroutine s_linspace_d

  subroutine s_logspace_d(xmin,xmax,x)
    implicit none
    real(dp),intent(in) :: xmin,xmax
    real(dp),intent(out) :: x(:)
    if (size(x) == 1 .and. xmin /= xmax) then
       write(0,'("ERROR: Cannot call logspace with n=1 and xmin /= xmax")')
       stop
    end if
    call linspace(log10(xmin),log10(xmax),x)
    x = 10._dp**x
  end subroutine s_logspace_d
