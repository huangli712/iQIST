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

  subroutine s_linspace_d(xmin,xmax,x)
    implicit none
    real(dp),intent(in) :: xmin,xmax
    real(dp),intent(out) :: x(:)
    integer :: i,n
    n = size(x)
    if (n == 1) then
       if(xmin /= xmax) then
          write(0,'("ERROR: Cannot call linspace with n=1 and xmin /= xmax")')
          stop
       else
          x = xmin
       end if
    else
       do i=1,n
          x(i) = (xmax-xmin) * real(i-1,dp) / real(n-1,dp) + xmin
       end do
    end if
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
