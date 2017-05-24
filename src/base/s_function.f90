!!!-----------------------------------------------------------------------
!!! project : CSSL (Common Service Subroutines Library)
!!! program : s_leg_basis
!!!           s_che_basis
!!!           s_svd_basis
!!!           s_sbessel
!!!           s_bezier
!!!           s_safe_exp
!!!           s_f_kernel
!!!           s_b_kernel
!!! source  : s_function.f90
!!! type    : subroutines & functions
!!! author  : li huang (email:lihuang.dmft@gmail.com)
!!! history : 07/10/2014 by li huang (created)
!!!           05/24/2017 by li huang (last modified)
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
!! 1. orthogonal polynomial basis
!! ------------------------------
!!
!! subroutine s_leg_basis(...)
!! subroutine s_che_basis(...)
!! subroutine s_svd_basis(...)
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
!! 4. some helper functions for s_svd_basis
!! ----------------------------------------
!!
!! function s_safe_exp(...)
!! function s_f_kernel(...)
!! function s_b_kernel(...)
!!

!!========================================================================
!!>>> orthogonal polynomial basis                                      <<<
!!========================================================================

!!
!! @sub s_leg_basis
!!
!! build legendre orthogonal polynomial in [-1,1] interval
!!
  subroutine s_leg_basis(lemax, legrd, lmesh, rep_l)
     use constants, only : dp, one

     implicit none

! external arguments
! maximum order for legendre orthogonal polynomial
     integer, intent(in)   :: lemax

! number of mesh points for legendre orthogonal polynomial
     integer, intent(in)   :: legrd

! mesh for legendre orthogonal polynomial in [-1,1]
     real(dp), intent(in)  :: lmesh(legrd)

! legendre orthogonal polynomial defined on [-1,1]
     real(dp), intent(out) :: rep_l(legrd,lemax)

! local variables
! loop index
     integer :: i
     integer :: j
     integer :: k

! check lemax
     if ( lemax <= 2 ) then
         call s_print_error('s_leg_basis','lemax must be larger than 2')
     endif ! back if ( lemax <= 2 ) block

! the legendre orthogonal polynomials obey the three term recurrence
! relation known as Bonnet’s recursion formula:
!     $P_0(x) = 1$
!     $P_1(x) = x$
!     $(n+1) P_{n+1}(x) = (2n+1) P_n(x) - n P_{n-1}(x)$
     do i=1,legrd
         rep_l(i,1) = one
         rep_l(i,2) = lmesh(i)
         do j=3,lemax
             k = j - 1
             rep_l(i,j) = ( real(2*k-1) * lmesh(i) * rep_l(i,j-1) - real(k-1) * rep_l(i,j-2) ) / real(k)
         enddo ! over j={3,lemax} loop
     enddo ! over i={1,legrd} loop

     return
  end subroutine s_leg_basis

!!
!! @sub s_che_basis
!!
!! build the second kind chebyshev orthogonal polynomial in [-1,1] interval
!!
  subroutine s_che_basis(chmax, chgrd, cmesh, rep_c)
     use constants, only : dp, one, two

     implicit none

! external arguments
! maximum order for chebyshev orthogonal polynomial
     integer, intent(in)   :: chmax

! number of mesh points for chebyshev orthogonal polynomial
     integer, intent(in)   :: chgrd

! mesh for chebyshev orthogonal polynomial in [-1,1]
     real(dp), intent(in)  :: cmesh(chgrd)

! chebyshev orthogonal polynomial defined on [-1,1]
     real(dp), intent(out) :: rep_c(chgrd, chmax)

! local variables
! loop index
     integer :: i
     integer :: j

! check chmax
     if ( chmax <= 2 ) then
         call s_print_error('s_che_basis','chmax must be larger than 2')
     endif ! back if ( chmax <= 2 ) block

! the chebyshev orthogonal polynomials of the second kind can be defined
! by the following recurrence relation
!     $U_0(x) = 1$
!     $U_1(x) = 2x$
!     $U_{n+1}(x) = 2xU_n(x) - U_{n-1}(x)$
     do i=1,chgrd
         rep_c(i,1) = one
         rep_c(i,2) = two * cmesh(i)
         do j=3,chmax
             rep_c(i,j) = two * cmesh(i) * rep_c(i,j-1) - rep_c(i,j-2)
         enddo ! over j={3,chmax} loop
     enddo ! over i={1,chgrd} loop

     return
  end subroutine s_che_basis

!!
!! @sub s_svd_basis
!!
!! build the svd orthogonal polynomial in [-1,1] interval
!!
  subroutine s_svd_basis(svmax, svgrd, smesh, rep_s, beta, stat)
     use constants, only : dp, zero, one

     implicit none

! external arguments
! using fermionic or bosonic kernel function
     character (len=1), intent(in) :: stat

! maximum order for svd orthogonal polynomial
     integer, intent(in)   :: svmax

! number of mesh points for svd orthogonal polynomial
     integer, intent(in)   :: svgrd

! inversion of temperature
     real(dp), intent(in)  :: beta

! mesh for svd orthogonal polynomial in [-1,1]
     real(dp), intent(in)  :: smesh(svgrd)

! svd orthogonal polynomial defined on [-1,1]
     real(dp), intent(out) :: rep_s(svgrd, svmax)

! external arguments
! used to calculate the fermionic kernel function
     procedure ( real(dp) ) :: s_f_kernel

! used to calculate the bosonic kernel function
     procedure ( real(dp) ) :: s_b_kernel

! local parameters
! number of mesh points for real axis
     integer, parameter  :: wsize = 513

! left boundary for real axis mesh, \omega_{min}
     real(dp), parameter :: w_min = -10.0_dp

! right boundary for real axis mesh, \omega_{max}
     real(dp), parameter :: w_max = +10.0_dp

! local variables
! loop index
     integer :: i
     integer :: j

! status flag
     integer :: istat

! real axis mesh
     real(dp), allocatable :: fmesh(:)

! fermionic or bosonic kernel function
     real(dp), allocatable :: fker(:,:)

! U, \Sigma, and V matrices for singular values decomposition
     real(dp), allocatable :: umat(:,:)
     real(dp), allocatable :: svec(:)
     real(dp), allocatable :: vmat(:,:)

     allocate(fmesh(wsize))
     allocate(fker(svgrd,wsize))
     allocate(umat(svgrd,wsize))
     allocate(svec(wsize))
     allocate(vmat(wsize,wsize))

! build real frequency mesh
     call s_linspace_d(rmin, rmax, wsize, fmesh)
     print *, 'hh'

! build the fermionic kernel
     do i=1,wsize
         do j=1,svgrd
             fker(j,i) = s_b_kernel(smesh(j), fmesh(i), beta)
         enddo ! over j={1,svgrd} loop
     enddo ! over i={1,wsize} loop

     print *, 'hh'

     call s_svd_dg(svgrd, wsize, wsize, fker, umat, svec, vmat)

     do i=1,wsize
         if ( umat(svgrd,i) < zero ) umat(:,i) = -one * umat(:,i)
     enddo

     !do i=1,svgrd
     !    write(*,'(i,3e16.8)') i, umat(i,1), umat(i,2), umat(i,3)
     !enddo

     do i=1,wsize
         print *, i, dot_product(umat(:,i), umat(:,i))
     enddo

     rep_s = umat(:,1:svmax)
     return
  end subroutine s_svd_basis

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





!!
!! @fun s_f_kernel
!!
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

     procedure( real(dp) ) :: s_safe_exp

! local variables
! return value
     real(dp) :: val

! dimensionless variables
     real(dp) :: x, y

     x = two * tau / beta - one
     y = beta * omega / two

     if ( y > 200.0_dp ) then
         val = s_safe_exp( -y * ( x + one ) )
     else if ( y < -200.0_dp ) then
         val = s_safe_exp(  y * ( one - x ) )
     else
         val = s_safe_exp( -x * y ) / ( two * cosh(y) )
     endif ! back if ( y > 200.0_dp ) block

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

     procedure( real(dp) ) :: s_safe_exp

! local variables
! return value
     real(dp) :: val

! dimensionless variables
     real(dp) :: x, y

     x = two * tau / beta - one
     y = beta * omega / two

     if ( abs(y) < 1E-10 ) then
         val = s_safe_exp( -x * y )
     else if ( y > 200.0_dp ) then
         val = two * y * s_safe_exp( -y * ( x + one ) )
     else if ( y < -200.0_dp ) then
         val = -two * y * s_safe_exp( y * ( one - x ) )
     else
         val = y * s_safe_exp( -x * y ) / sinh(y)
     endif ! back if ( abs(y) < 1E-10 ) block

     return
  end function s_b_kernel

  function s_safe_exp(x) result(val)
     use constants, only : dp, zero

     implicit none

     real(dp), intent(in) :: x
     real(dp) :: val

     if ( x < -60.0_dp ) then
         val = zero
     else
         val = exp(x)
     endif

     return
  end function s_safe_exp

  program test
     use constants, only : dp, zero

     integer, parameter :: svmax = 40
     integer, parameter :: svgrd = 10001
     real(dp), parameter :: beta = 10.0_dp
     real(dp) :: smesh(svgrd)
     real(dp) :: rep_s(svgrd, svmax)

! build time mesh
     call s_linspace_d(zero, beta, svgrd, smesh)
     print *, smesh
     call s_svd_basis(svmax, svgrd, beta, smesh, rep_s)
  end program test
