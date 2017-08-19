!!!-----------------------------------------------------------------------
!!! project : CSSL (Common Service Subroutines Library)
!!! program : s_linspace_i
!!!           s_linspace_d
!!!           s_linspace_z
!!!           s_cumsum_i
!!!           s_cumsum_d
!!!           s_cumsum_z
!!!           s_cumprod_i
!!!           s_cumprod_d
!!!           s_cumprod_z
!!!           s_swap_i
!!!           s_swap_d
!!!           s_swap_z
!!!           s_mix_i
!!!           s_mix_d
!!!           s_mix_z
!!!           s_vecadd_i
!!!           s_vecadd_d
!!!           s_vecadd_z
!!! source  : s_vector.f90
!!! type    : subroutines
!!! author  : li huang (email:lihuang.dmft@gmail.com)
!!! history : 07/10/2014 by li huang (created)
!!!           01/10/2018 by li huang (last modified)
!!! purpose : these subroutines are designed for vectors or arrays. they
!!!           can be used to manipulate grid and mesh.
!!! status  : unstable
!!! comment :
!!!-----------------------------------------------------------------------

!!
!!
!! Introduction
!! ============
!!
!! 1. mesh generation
!! ------------------
!!
!! subroutine s_linspace_i(...)
!! subroutine s_linspace_d(...)
!! subroutine s_linspace_z(...)
!!
!! 2. sum of vector
!! ----------------
!!
!! subroutine s_cumsum_i(...)
!! subroutine s_cumsum_d(...)
!! subroutine s_cumsum_z(...)
!!
!! 3. product of vector
!! --------------------
!!
!! subroutine s_cumprod_i(...)
!! subroutine s_cumprod_d(...)
!! subroutine s_cumprod_z(...)
!!
!! 4. swap two vectors
!! -------------------
!!
!! subroutine s_swap_i(...)
!! subroutine s_swap_d(...)
!! subroutine s_swap_z(...)
!!
!! 5. linear mixing for vectors
!! ----------------------------
!!
!! subroutine s_mix_i(...)
!! subroutine s_mix_d(...)
!! subroutine s_mix_z(...)
!!
!! 6. convert diagonal elements of matrix to vector
!! ------------------------------------------------
!!
!! subroutine s_vecadd_i(...)
!! subroutine s_vecadd_d(...)
!! subroutine s_vecadd_z(...)
!!
!!

!!========================================================================
!!>>> mesh generation                                                  <<<
!!========================================================================

!!
!! @sub s_linspace_i
!!
!! create a linear mesh x in interval [xmin, xmax], integer version
!!
  subroutine s_linspace_i(xmin, xmax, n, x)
     use constants, only : dp

     implicit none

! external arguments
! left boundary
     integer, intent(in)  :: xmin

! right boundary
     integer, intent(in)  :: xmax

! size of array x
     integer, intent(in)  :: n

! output array, containing the linear mesh
     integer, intent(out) :: x(n)

! local variables
! loop index
     integer :: i

     do i=1,n
         x(i) = ( xmax - xmin ) * real(i - 1, dp) / real(n - 1, dp) + xmin
     enddo ! over i={1,n} loop

     return
  end subroutine s_linspace_i

!!
!! @sub s_linspace_d
!!
!! create a linear mesh x in interval [xmin, xmax], real(dp) version
!!
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

!!
!! @sub s_linspace_z
!!
!! create a linear mesh x in interval [xmin, xmax], complex(dp) version
!!
  subroutine s_linspace_z(xmin, xmax, n, x)
     use constants, only : dp

     implicit none

! external arguments
! left boundary
     complex(dp), intent(in)  :: xmin

! right boundary
     complex(dp), intent(in)  :: xmax

! size of array x
     integer,  intent(in)     :: n

! output array, containing the linear mesh
     complex(dp), intent(out) :: x(n)

! local variables
! loop index
     integer :: i

     do i=1,n
         x(i) = ( xmax - xmin ) * real(i - 1, dp) / real(n - 1, dp) + xmin
     enddo ! over i={1,n} loop

     return
  end subroutine s_linspace_z

!!========================================================================
!!>>> sum operations                                                   <<<
!!========================================================================

!!
!! @sub s_cumsum_i
!!
!! return the cumsum of an integer array
!!
  subroutine s_cumsum_i(n, v, vsum)
     implicit none

! external arguments
! size of array v
     integer, intent(in)  :: n

! input integer array
     integer, intent(in)  :: v(n)

! cumsum of array v
     integer, intent(out) :: vsum(n)

! local variables
! loop index
     integer :: i

     vsum(1) = v(1)
     do i=2,n
         vsum(i) = vsum(i-1) + v(i)
     enddo ! over i={2,n} loop

     return
  end subroutine s_cumsum_i

!!
!! @sub s_cumsum_d
!!
!! return the cumsum of a real(dp) array
!!
  subroutine s_cumsum_d(n, v, vsum)
     use constants, only : dp

     implicit none

! external arguments
! size of array v
     integer, intent(in)   :: n

! input real(dp) array
     real(dp), intent(in)  :: v(n)

! cumsum of array v
     real(dp), intent(out) :: vsum(n)

! local variables
! loop index
     integer :: i

     vsum(1) = v(1)
     do i=2,n
         vsum(i) = vsum(i-1) + v(i)
     enddo ! over i={2,n} loop

     return
  end subroutine s_cumsum_d

!!
!! @sub s_cumsum_z
!!
!! return the cumsum of a complex(dp) array
!!
  subroutine s_cumsum_z(n, v, vsum)
     use constants, only : dp

     implicit none

! external arguments
! size of array v
     integer, intent(in)      :: n

! input complex(dp) array
     complex(dp), intent(in)  :: v(n)

! cumsum of array v
     complex(dp), intent(out) :: vsum(n)

! local variables
! loop index
     integer :: i

     vsum(1) = v(1)
     do i=2,n
         vsum(i) = vsum(i-1) + v(i)
     enddo ! over i={2,n} loop

     return
  end subroutine s_cumsum_z

!!========================================================================
!!>>> prod operations                                                  <<<
!!========================================================================

!!
!! @sub s_cumprod_i
!!
!! return the cumproduct of an integer array
!!
  subroutine s_cumprod_i(n, v, vprod)
     implicit none

! external arguments
! size of array v
     integer, intent(in)  :: n

! input integer array
     integer, intent(in)  :: v(n)

! cumproduct of array v
     integer, intent(out) :: vprod(n)

! local variables
! loop index
     integer :: i

     vprod(1) = v(1)
     do i=2,n
         vprod(i) = vprod(i-1) * v(i)
     enddo ! over i={2,n} loop

     return
  end subroutine s_cumprod_i

!!
!! @sub s_cumprod_d
!!
!! return the cumproduct of a real(dp) array
!!
  subroutine s_cumprod_d(n, v, vprod)
     use constants, only : dp

     implicit none

! external arguments
! size of array v
     integer, intent(in)   :: n

! input real(dp) array
     real(dp), intent(in)  :: v(n)

! cumproduct of array v
     real(dp), intent(out) :: vprod(n)

! local variables
! loop index
     integer :: i

     vprod(1) = v(1)
     do i=2,n
         vprod(i) = vprod(i-1) * v(i)
     enddo ! over i={2,n} loop

     return
  end subroutine s_cumprod_d

!!
!! @sub s_cumprod_z
!!
!! return the cumproduct of a complex(dp) array
!!
  subroutine s_cumprod_z(n, v, vprod)
     use constants, only : dp

     implicit none

! external arguments
! size of array v
     integer, intent(in)      :: n

! input complex(dp) array
     complex(dp), intent(in)  :: v(n)

! cumproduct of array v
     complex(dp), intent(out) :: vprod(n)

! local variables
! loop index
     integer :: i

     vprod(1) = v(1)
     do i=2,n
         vprod(i) = vprod(i-1) * v(i)
     enddo ! over i={2,n} loop

     return
  end subroutine s_cumprod_z

!!========================================================================
!!>>> swap operations                                                  <<<
!!========================================================================

!!
!! @sub s_swap_i
!!
!! exchange two integer vectors
!!
  subroutine s_swap_i(n, ix, iy)
     implicit none

! external arguments
! dimension of integer vector
     integer, intent(in)    :: n

! integer vector X
     integer, intent(inout) :: ix(n)

! integer vector Y
     integer, intent(inout) :: iy(n)

! local variables
! dummy integer vector
     integer :: it(n)

     it = ix
     ix = iy
     iy = it

     return
  end subroutine s_swap_i

!!
!! @sub s_swap_d
!!
!! exchange two real(dp) vectors
!!
  subroutine s_swap_d(n, dx, dy)
     use constants, only : dp

     implicit none

! external arguments
! dimension of real(dp) vector
     integer, intent(in)     :: n

! real(dp) vector X
     real(dp), intent(inout) :: dx(n)

! real(dp) vector Y
     real(dp), intent(inout) :: dy(n)

! local variables
! dummy real(dp) vector
     real(dp) :: dt(n)

     dt = dx
     dx = dy
     dy = dt

     return
  end subroutine s_swap_d

!!
!! @sub s_swap_z
!!
!! exchange two complex(dp) vectors
!!
  subroutine s_swap_z(n, zx, zy)
     use constants, only : dp

     implicit none

! external arguments
! dimension of complex(dp) vector
     integer, intent(in)        :: n

! complex(dp) vector X
     complex(dp), intent(inout) :: zx(n)

! complex(dp) vector Y
     complex(dp), intent(inout) :: zy(n)

! local variables
! dummy complex(dp) vector
     complex(dp) :: zt(n)

     zt = zx
     zx = zy
     zy = zt

     return
  end subroutine s_swap_z

!!========================================================================
!!>>> mix operations                                                   <<<
!!========================================================================

!!
!! @sub s_mix_i
!!
!! linear mixing for two integer vectors
!!
  subroutine s_mix_i(n, ix, iy, alpha)
     use constants, only : dp
     use constants, only : one

     implicit none

! external arguments
! dimension of integer vector
     integer, intent(in)    :: n

! mixing parameter
     real(dp), intent(in)   :: alpha

! integer vector X
     integer, intent(in)    :: ix(n)

! integer vector Y
     integer, intent(inout) :: iy(n)

     iy = int( real(ix) * (one - alpha) + real(iy) * alpha )

     return
  end subroutine s_mix_i

!!
!! @sub s_mix_d
!!
!! linear mixing for two real(dp) vectors
!!
  subroutine s_mix_d(n, dx, dy, alpha)
     use constants, only : dp
     use constants, only : one

     implicit none

! external arguments
! dimension of real(dp) vector
     integer, intent(in)     :: n

! mixing parameter
     real(dp), intent(in)    :: alpha

! real(dp) vector X
     real(dp), intent(in)    :: dx(n)

! real(dp) vector Y
     real(dp), intent(inout) :: dy(n)

     dy = dx * (one - alpha) + dy * alpha

     return
  end subroutine s_mix_d

!!
!! @sub s_mix_z
!!
!! linear mixing for two complex(dp) vectors
!!
  subroutine s_mix_z(n, zx, zy, alpha)
     use constants, only : dp
     use constants, only : one

     implicit none

! external arguments
! dimension of complex(dp) vector
     integer, intent(in)        :: n

! mixing parameter
     real(dp), intent(in)       :: alpha

! complex(dp) vector X
     complex(dp), intent(in)    :: zx(n)

! complex(dp) vector Y
     complex(dp), intent(inout) :: zy(n)

     zy = zx * (one - alpha) + zy * alpha

     return
  end subroutine s_mix_z

!!========================================================================
!!>>> vector add operations                                            <<<
!!========================================================================

  subroutine s_vecadd_z()
  end subroutine s_vecadd_z
