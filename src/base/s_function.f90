!!!-----------------------------------------------------------------------
!!! project : CSSL (Common Service Subroutines Library)
!!! program : s_legendre
!!!           s_chebyshev
!!!           s_sbessel
!!!           s_bezier
!!!           s_f_kernel
!!!           s_b_kernel
!!! source  : s_function.f90
!!! type    : subroutines & functions
!!! author  : li huang (email:lihuang.dmft@gmail.com)
!!! history : 07/10/2014 by li huang (created)
!!!           05/22/2017 by li huang (last modified)
!!! purpose : these subroutines are used to generate some auxiliary
!!!           functions, such as the Legendre orthogonal polynomial and
!!!           Chebyshev orthogonal polynomial, Bessel function, etc.
!!! status  : unstable
!!! comment :
!!!-----------------------------------------------------------------------

!!
!!
!! Introduction
!! ============
!!
!! 1. orthogonal polynomial
!! ------------------------
!!
!! subroutine s_legendre(...)
!! subroutine s_chebyshev(...)
!!
!! 2. spheric Bessel function
!! --------------------------
!!
!! subroutine s_sbessel(...)
!!
!! 3. bernstein polynomial
!! -----------------------
!!
!! subroutine s_bezier(...)
!!
!!

!!========================================================================
!!>>> Legendre and Chebyshev polynomials                               <<<
!!========================================================================

!!
!! @sub s_legendre
!!
!! build legendre polynomial in [-1,1]
!!
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
!     $P_1(x) = x$
!     $(n+1) P_{n+1}(x) = (2n+1) P_n(x) - n P_{n-1}(x)$
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

!!
!! @sub s_chebyshev
!!
!! build the second kind chebyshev polynomial in [-1,1]
!!
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

!!========================================================================
!!>>> spherical Bessel functions                                       <<<
!!========================================================================

!!
!! @sub s_sbessel
!!
!! computes the spherical Bessel functions of the first kind, j_l(x), for
!! argument x and l=0,1,\ldots,l_{max}.
!!
  subroutine s_sbessel(lmax, x, jl)
     use constants, only : dp, zero, one, two, eps8

     implicit none

! external arguments
! maximum order of spherical Bessel function
     integer, intent(in)   :: lmax

! real argument
     real(dp), intent(in)  :: x

! array of returned values
     real(dp), intent(out) :: jl(0:lmax)

! local parameters
! staring value for l above lmax (suitable for lmax < 50)
     integer, parameter  :: lst  = 25

! rescale limit
     real(dp), parameter :: rsc  = 1.0D100
     real(dp), parameter :: rsci = one / rsc

! local variables
! loop index
     integer  :: l

! real(dp) dummy variables
     real(dp) :: xi, jt
     real(dp) :: j0, j1
     real(dp) :: t1, t2

! important note: the recursion relation
!     j_{l+1}(x)=\frac{2l+1}{x}j_l(x)-j_{l-1}(x)
! is used either downwards for x < l or upwards for x >= l. for x << 1,
! the asymtotic form is used:
!     j_l(x) \approx \frac{x^l}{(2l+1)!!}
! this procedure is numerically stable and accurate to near this machine
! precision for l <= 50

! check the range of input variables
     if ( lmax < 0 .or. lmax > 50 ) then
         call s_print_error('s_sbessel','lmax is out of range')
     endif ! back if ( lmax < 0 .or. lmax > 50 ) block

     if ( x < zero .or. x > 1.0E5 ) then
         call s_print_error('s_sbessel','x is out of range')
     endif ! back if ( x < zero .or. x > 1.0E5 ) block

! treat x << 1
     if ( x < eps8 ) then
         jl(0) = one
         t1 = one; t2 = one
         do l=1,lmax
             t1 = t1 / (two * l + one)
             t2 = t2 * x
             jl(l) = t2 * t1
         enddo ! over l={1,lmax} loop
         RETURN
     endif ! back if ( x < eps8 ) block

     xi = one / x

! for x < lmax recurse down
     if ( x < lmax ) then
         if ( lmax == 0 ) then
             jl(0) = sin(x) / x; RETURN
         endif ! back if ( lmax == 0 ) block

! start from truly random numbers
         j0 = 0.6370354636449841609d0 * rsci
         j1 = 0.3532702964695481204d0 * rsci
         do l=lmax+lst,lmax+1,-1
             jt = j0 * (two * l + one) * xi - j1
             j1 = j0
             j0 = jt
! check for overflow
             if ( abs(j0) > rsc ) then
! rescale
                 jt = jt * rsci
                 j1 = j1 * rsci
                 j0 = j0 * rsci
             endif ! back if ( abs(j0) > rsc ) block
         enddo ! over l={lmax+lst,lmax+1} loop

         do l=lmax,0,-1
             jt = j0 * (two * l + one) * xi - j1
             j1 = j0
             j0 = jt
! check for overflow
             if ( abs(j0) > rsc ) then
! rescale
                 jt = jt * rsci
                 j1 = j1 * rsci
                 j0 = j0 * rsci
                 jl(l+1:lmax) = jl(l+1:lmax) * rsci
             endif ! back if ( abs(j0) > rsc ) block
             jl(l) = j1
         enddo ! over l={lmax,0} loop
! rescaling constant
         t1 = one / ( ( jl(0) - x * jl(1) ) * cos(x) + x * jl(0) * sin(x) )
         jl = t1 * jl
     else
! for large x recurse up
         jl(0) = sin(x) * xi
         if ( lmax == 0 ) RETURN
         jl(1) = ( jl(0) - cos(x) ) * xi
         if ( lmax == 1 ) RETURN
         j0 = jl(0)
         j1 = jl(1)
         do l=2,lmax
             jt = (two * l - one ) * j1 * xi - j0
             j0 = j1
             j1 = jt
             jl(l) = j1
         enddo ! over l={2,lmax} loop
     endif ! back if ( x < lmax ) block

     return
  end subroutine s_sbessel

!!========================================================================
!!>>> Bernstein polynomials                                            <<<
!!========================================================================

!!
!! @sub s_bezier
!!
!! to evaluates the bernstein polynomials at a point x
!!
  subroutine s_bezier(n, x, bern)
     use constants, only : dp, one

     implicit none

! external arguments
! the degree of the bernstein polynomials to be used. for any N, there
! is a set of N+1 bernstein polynomials, each of degree N, which form a
! basis for polynomials on [0,1].
     integer, intent(in)  :: n

! the evaluation point.
     real(dp), intent(in) :: x

! the values of the N+1 bernstein polynomials at X
     real(dp), intent(inout) :: bern(0:n)

! local variables
! loop index
     integer :: i
     integer :: j

! the bernstein polynomials are assumed to be based on [0,1].
! the formula is:
!
!    B(N,I)(X) = [N!/(I!*(N-I)!)] * (1-X)**(N-I) * X**I
!
! first values:
!
!    B(0,0)(X) = 1
!    B(1,0)(X) =      1-X
!    B(1,1)(X) =                X
!    B(2,0)(X) =     (1-X)**2
!    B(2,1)(X) = 2 * (1-X)    * X
!    B(2,2)(X) =                X**2
!    B(3,0)(X) =     (1-X)**3
!    B(3,1)(X) = 3 * (1-X)**2 * X
!    B(3,2)(X) = 3 * (1-X)    * X**2
!    B(3,3)(X) =                X**3
!    B(4,0)(X) =     (1-X)**4
!    B(4,1)(X) = 4 * (1-X)**3 * X
!    B(4,2)(X) = 6 * (1-X)**2 * X**2
!    B(4,3)(X) = 4 * (1-X)    * X**3
!    B(4,4)(X) =                X**4
!
! special values:
!
!    B(N,I)(X) has a unique maximum value at X = I/N.
!    B(N,I)(X) has an I-fold zero at 0 and and N-I fold zero at 1.
!    B(N,I)(1/2) = C(N,K) / 2**N
!    for a fixed X and N, the polynomials add up to 1:
!    sum ( 0 <= I <= N ) B(N,I)(X) = 1
!
     if ( n == 0 ) then
         bern(0) = one

     else if ( 0 < n ) then
         bern(0) = one - x
         bern(1) = x
         do i=2,n
             bern(i) = x * bern(i-1)
             do j=i-1,1,-1
                 bern(j) = x * bern(j-1) + ( one - x ) * bern(j)
             enddo ! over j={i-1,1} loop
             bern(0) = ( one - x ) * bern(0)
         enddo ! over i={2,n} loop

     endif ! back if ( n == 0 ) block

     return
  end subroutine s_bezier






  function s_f_kernel(tau, omega, beta) result(val)
     use constants, only : dp, one, two

     implicit none

! external arguments
! imaginary time point
     real(dp), intent(in) :: tau

! frequency point
     real(dp), intent(in) :: omega

! inverse temperature
     real(dp), intent(in) :: beta

! local variables
! return value
     real(dp) :: val

! dimensionless variables
     real(dp) :: x, y

     x = beta * omega / two
     y = two * tau / beta - one

     if ( x > 100.0_dp ) then
         val = exp( -x * ( y + one ) )
     else if ( x < -100.0_dp ) then
         val = exp(  x * ( one - y ) )
     else
         val = exp( -x * y ) / ( two * cosh(x) )
     endif ! back if ( x > 100.0_dp ) block

     return
  end function s_f_kernel

  function s_b_kernel(tau, omega, beta) result(val)
     use constants, only : dp, one, two

! external arguments
! imaginary time point
     real(dp), intent(in) :: tau

! frequency point
     real(dp), intent(in) :: omega

! inverse temperature
     real(dp), intent(in) :: beta

! local variables
! return value
     real(dp) :: val

! dimensionless variables
     real(dp) :: x, y

     x = beta * omega / two
     y = two * tau / beta - one

     if ( x > 100.0_dp ) then
         val = omega * exp( -x * ( y + one ) )
     else if ( x < -100.0_dp ) then
         val = -omega * exp(  x * ( one - y ) )
     else
         val = omega * exp( -x * y ) / ( two * sinh(x) )
     endif ! back if ( x > 100.0_dp ) block

     return
  end function s_b_kernel
