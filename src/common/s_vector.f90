!!!-----------------------------------------------------------------------
!!! project : CSSL (Common Service Subroutines Library)
!!! program : s_linspace_d
!!!           s_logspace_d
!!!           s_linspace_z
!!!           s_sum_i
!!!           s_sum_d
!!1           s_sum_z
!!!           s_cumsum_i
!!!           s_cumsum_d
!!!           s_cumsum_z
!!!           s_prod_i
!!!           s_prod_d
!!!           s_prod_z
!!!           s_cumprod_i
!!!           s_cumprod_d
!!!           s_cumprod_z
!!!           s_swap_i
!!!           s_swap_d
!!!           s_swap_z
!!!           s_legendre
!!!           s_chebyshev
!!! source  : s_math.f90
!!! type    : subroutines
!!! author  : li huang (email:huangli712@gmail.com)
!!! history : 07/10/2014 by li huang
!!!           07/24/2014 by li huang
!!!           08/01/2014 by li huang
!!! purpose : these subroutines are designed for vectors or arrays. They
!!!           can be used to manipulate grid and mesh, to generate the
!!!           Legendre polynomial and Chebyshev polynomial, etc.
!!! status  : unstable
!!! comment :
!!!-----------------------------------------------------------------------

!!========================================================================
!!>>> mesh generation                                                  <<<
!!========================================================================

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
     call s_linspace_d(log10(xmin), log10(xmax), n, x)
     x = 10.0_dp**x

     return
  end subroutine s_logspace_d

!!>>> s_linspace_z: create a linear mesh x in interval [xmin, xmax], complex(dp) version
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

!!>>> s_sum_i: return the sum of an integer array
  subroutine s_sum_i(n, v, vsum)
     implicit none

! external arguments
! size of array v
     integer, intent(in)  :: n

! sum of array v
     integer, intent(out) :: vsum

! input integer array
     integer, intent(in)  :: v(n)

     vsum = sum(v)

     return
  end subroutine s_sum_i

!!>>> s_sum_d: return the sum of a real(dp) array
  subroutine s_sum_d(n, v, vsum)
     use constants, only : dp

     implicit none

! external arguments
! size of array v
     integer, intent(in)   :: n

! sum of array v
     real(dp), intent(out) :: vsum

! input real(dp) array
     real(dp), intent(in)  :: v(n)

     vsum = sum(v)

     return
  end subroutine s_sum_d

!!>>> s_sum_z: return the sum of a complex(dp) array
  subroutine s_sum_z(n, v, vsum)
     use constants, only : dp

     implicit none

! external arguments
! size of array v
     integer, intent(in)      :: n

! sum of array v
     complex(dp), intent(out) :: vsum

! input complex(dp) array
     complex(dp), intent(in)  :: v(n)

     vsum = sum(v)

     return
  end subroutine s_sum_z

!!>>> s_cumsum_i: return the cumsum of an integer array
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

!!>>> s_cumsum_d: return the cumsum of a real(dp) array
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

!!>>> s_cumsum_z: return the cumsum of a complex(dp) array
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

!!>>> s_prod_i: return the product of an integer array
  subroutine s_prod_i(n, v, vprod)
     implicit none

! external arguments
! size of array v
     integer, intent(in)  :: n

! product of array v
     integer, intent(out) :: vprod

! input integer array
     integer, intent(in)  :: v(n)

     vprod = product(v)

     return
  end subroutine s_prod_i

!!>>> s_prod_d: return the product of a real(dp) array
  subroutine s_prod_d(n, v, vprod)
     use constants, only : dp

     implicit none

! external arguments
! size of array v
     integer, intent(in)   :: n

! product of array v
     real(dp), intent(out) :: vprod

! input real(dp) array
     real(dp), intent(in)  :: v(n)

     vprod = product(v)

     return
  end subroutine s_prod_d

!!>>> s_prod_z: return the product of a complex(dp) array
  subroutine s_prod_z(n, v, vprod)
     use constants, only : dp

     implicit none

! external arguments
! size of array v
     integer, intent(in)      :: n

! product of array v
     complex(dp), intent(out) :: vprod

! input complex(dp) array
     complex(dp), intent(in)  :: v(n)

     vprod = product(v)

     return
  end subroutine s_prod_z

!!>>> s_cumprod_i: return the cumproduct of an integer array
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

!!>>> s_cumprod_d: return the cumproduct of a real(dp) array
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

     return
  end subroutine s_cumprod_d

!!>>> s_cumprod_z: return the cumproduct of a complex(dp) array
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

     return
  end subroutine s_cumprod_z

  program test
     complex(8) :: a(4), b(4)
     a = 2
     call s_cumprod_i(4, a, b)
     print *, a
     print *, b
     
  end program test

!!========================================================================
!!>>> swap operations                                                  <<<
!!========================================================================

!!>>> s_swap_i: exchange two integer vectors
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

!!>>> s_swap_d: exchange two real(dp) vectors
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

!!>>> s_swap_z: exchange two complex(dp) vectors
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
!!>>> Legendre and Chebyshev polynomials                               <<<
!!========================================================================

!!>>> s_legendre: build legendre polynomial in [-1,1]
  subroutine s_legendre(lemax, legrd, pmesh, ppleg)
     use constants, only : dp, one

     implicit none

! external arguments
! maximum order for legendre polynomial
     integer, intent(in)   :: lemax

! number of mesh points for legendre polynomial
     integer, intent(in)   :: legrd

! mesh for legendre polynomial in [-1,1]
     real(dp), intent(in)  :: pmesh(legrd)

! legendre polynomial defined on [-1,1]
     real(dp), intent(out) :: ppleg(legrd,lemax)

! local variables
! loop index
     integer :: i
     integer :: j
     integer :: k

! check lemax
     if ( lemax <= 2 ) then
         call s_print_error('s_legendre','lemax must be larger than 2')
     endif ! back if ( lemax <= 2 ) block

! the legendre polynomials obey the three term recurrence relation known
! as Bonnet’s recursion formula:
!     $P_0(x) = 1$
!     P_1(x) = x$
!     $(n+1) P_{n+1}(x) = (2n+1) x P_n(x) - n P_{n-1}(x)$
     do i=1,legrd
         ppleg(i,1) = one
         ppleg(i,2) = pmesh(i)
         do j=3,lemax
             k = j - 1
             ppleg(i,j) = ( real(2*k-1) * pmesh(i) * ppleg(i,j-1) - real(k-1) * ppleg(i,j-2) ) / real(k)
         enddo ! over j={3,lemax} loop
     enddo ! over i={1,legrd} loop

     return
  end subroutine s_legendre

!!>>> s_chebyshev: build chebyshev polynomial in [-1,1]
!!>>> note: it is second kind chebyshev polynomial
  subroutine s_chebyshev(chmax, chgrd, qmesh, qqche)
     use constants, only : dp, one, two

     implicit none

! external arguments
! maximum order for chebyshev polynomial
     integer, intent(in)   :: chmax

! number of mesh points for chebyshev polynomial
     integer, intent(in)   :: chgrd

! mesh for chebyshev polynomial in [-1,1]
     real(dp), intent(in)  :: qmesh(chgrd)

! chebyshev polynomial defined on [-1,1]
     real(dp), intent(out) :: qqche(chgrd, chmax)

! local variables
! loop index
     integer :: i
     integer :: j

! check chmax
     if ( chmax <= 2 ) then
         call s_print_error('s_chebyshev','chmax must be larger than 2')
     endif ! back if ( chmax <= 2 ) block

! the chebyshev polynomials of the second kind can be defined by the
! following recurrence relation
!     $U_0(x) = 1$
!     $U_1(x) = 2x$
!     $U_{n+1}(x) = 2xU_n(x) - U_{n-1}(x)$
     do i=1,chgrd
         qqche(i,1) = one
         qqche(i,2) = two * qmesh(i)
         do j=3,chmax
             qqche(i,j) = two * qmesh(i) * qqche(i,j-1) - qqche(i,j-2)
         enddo ! over j={3,chmax} loop
     enddo ! over i={1,chgrd} loop

     return
  end subroutine s_chebyshev
