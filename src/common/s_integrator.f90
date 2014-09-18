! Caleb Wherry (Alone)
! Created: Oct. 24th, 2010
! Last Modified: Nov. 3, 2010
! Module containg different fortran functions for integrating functions
!
! A) 1D function to evaluate
! B) 2D function to evaluate
! c) 3D function to evaluate
!
! 1)  Composite Midpoint Rule 1D
! 2)  Composite Trapezoid Rule 1D
! 3)  Romburg's Method with Composite Trapezoid Rule 1D
! 4)  Composite Simpson's Method 1D
! 5)  Composite Wherry's Method 1D
! 6)  Gaussian Quadrature 1D w/ n=5
! 7)  Gaussian Quadrature 2D w/ n=5
! 8)  Gaussian Quadrature 3D w/ n=5
! 9)  Romburg's Method with Gaussian Quadrature 1D w/ n=5
! 10) Composite Trapezoid 2D
! 11) Romburg's Method with Composite Trapezoid 2D

module integrationMethods
  implicit none

  ! Value of PI:
  double precision, parameter ::  PI = 3.141592653589793238462643d0

  ! Roots and coeff for gaussian quadratures
  double precision, dimension(0:4) :: &
    r = (/  0.9061798459386640, &
            0.5384693101056831, &
            0.0000000000000000, &
           -0.5384693101056831, &
           -0.9061798459386640  &
        /)

  double precision, dimension(0:4) :: &
    coeff = (/  0.2369268850561891, &
                0.4786286704993665, &
                0.5688888888888889, &
                0.4786286704993665, &
                0.2369268850561891  &
            /)

  ! Interface for handling functions of different dimensions
  interface f
    module procedure f1 ! 1D
    module procedure f2 ! 2D
    module procedure f3 ! 3D
  end interface f

  contains

    ! A) 1D function
    double precision function f1(x)

      double precision, intent(in) :: x

      !f1 = (1.0d0 / sqrt(2*PI) ) * exp(-(1.0d0/2.0d0)*(x)**2)
      !f1 = -10.0d0 / (x * sqrt(x))
      f1 = (1.0d0/sqrt(2.0d0*PI)) * exp(-x**2/2.0d0)

      return

    end function f1

    ! B) 2D function
    double precision function f2(x,y)

      double precision, intent(in) :: x,y

      !f2 = sqrt(x**5 * (cos(y) + sin(y)) ) / (cos(y) + sin(y))
      !f2 = x**2 + y**2 
      f2 = ((x-0.5d0)**2 + (y-0.5d0)**2) -  

      return

    end function f2

    ! C) 3D function
    double precision function f3(x,y,z)

      double precision, intent(in) :: x,y,z

      f3 = ( (2.0d0+x**3)*sin(y+4.0d0)*(z+1.0d0)**2 ) / &
           sqrt( (x-4.0d0)**2 + (y-5.0d0)**2 + (z-6.0d0)**2 )

      return

    end function f3

    ! 1) Composite Midpoint Rule
    double precision function compMidpoint1d(a,b,n)

      double precision, intent(in) :: a, b
      integer, intent(in) :: n
      double precision :: h, evenSum
      integer :: i

      h = (b-a) / dble(n)

      do i=0, n-1
        if ( mod(i,2) == 0 ) then
          evenSum = evenSum + f(a+dble(i)*h)
        end if
      end do

      compMidpoint1d = 2.0d0 * h * evenSum

      return

    end function compMidpoint1d

    ! 2) Composite Trapezoid Rule:
    double precision function compTrapezoid1d(a,b,n)
 
      double precision, intent(in) :: a, b
      integer, intent(in) :: n
      double precision :: h, trapSum
      integer :: i

      h = (b-a) / dble(n)

      trapSum = 0.0d0

      do i=1, n-1
        trapSum = trapSum + f(a+dble(i)*h)
      end do

      compTrapezoid1d = (h/2.0d0) * ( f(a) + f(b) + 2.0d0*trapSum )

      return 

    end function compTrapezoid1d

    ! 3) Romberg Method with Composite Trapezoid Rule
    function rombergCompTrap1d(a,b,m,n)

      double precision, intent(in) :: a,b
      integer, intent(in) :: n
      double precision, dimension(0:n,0:n) :: R, rombergCompTrap1d
      integer :: i,k,m

      R = 0.0d0
  
      ! Loop to compute intial Trapezoid approximations:
      do i=0, n
        R(i,0) = compTrapezoid1d(a,b,m)
        m = 2*m
      end do
  
      ! Loop to compute Romberg approximations from Intial Trapezoids approxs:
      do k=1, n
        do i=k, n
          R(i,k) = R(i,k-1) + ( (R(i,k-1) - R(i-1,k-1)) / (4**k - 1) )
        end do
      end do

      rombergCompTrap1d = R

      return

    end function rombergCompTrap1d

    ! 4) Composite Simpson's Method:
    double precision function compSimpsons1d(a,b,n)

      double precision, intent(in) :: a,b
      integer, intent(in) :: n
      double precision :: h, oddSum, evenSum
      integer :: i

      h = (b-a) / dble(n)

      do i=1, n-1
        if ( mod(i,2) == 0 ) then
          evenSum = evenSum + f(a+dble(i)*h)
        else
          oddSum = oddSum + f(a+dble(i)*h)
        end if
      end do

      compSimpsons1d = (h/3.0d0) * ( f(a) + f(b) + 2.0d0*evenSum + 4.0d0*oddsum )

      return

    end function compSimpsons1d

    ! 5) Composite Wherry's Method
    double precision function compWherrys1d(a,b,n)

      double precision, intent(in) :: a,b
      integer, intent(in) :: n
      double precision :: h, sum1, sum2, sum3
      integer :: i

      h = (b-a) / dble(n)

      do i=1, n-1
        if ( mod(i,3) == 1 ) then
          sum1 = sum1 + f(a+dble(i)*h)
        else if ( mod(i,3) == 2) then
          sum2 = sum2 + f(a+dble(i)*h)
        else 
          sum3 = sum3 + f(a+dble(i)*h)
        end if
      end do

      compWherrys1d = (h/4.0d0) * ( f(a) + f(b) + 2.0d0*sum1 + 4.0d0*sum2 + 6.0d0*sum3 )

      return

    end function compWherrys1d

    ! 6) Gaussian Qudrature 1D w/ n=5
    double precision function gaussQuad1d5(a,b)

      double precision, intent(in) :: a, b
      double precision :: sum1
      integer :: i

      sum1 = 0.0d0
      do i=0,4
        sum1 = sum1 + coeff(i)*f( ((b-a)*r(i) + b + a) / 2.0d0 )
      end do

      gaussQuad1d5 = ( (b - a) / 2.0d0 ) * sum1

      return

    end function gaussQuad1d5

    ! 7) Gaussian Quadrature 2D w/ n=5
    double precision function gaussQuad2d5(a,b,c,d)

      double precision, intent(in) :: a, b, c, d
      double precision :: sum1
      integer :: i,j

      sum1 = 0.0d0
      do i=0,4
        do j=0,4
          sum1 = sum1 + coeff(i) * coeff(j) * &
                 f( ((b-a)*r(i) + b + a) / 2.0d0 , &
                    ((d-c)*r(j) + d + c) / 2.0d0 )
        end do
      end do

      gaussQuad2d5 = ((b - a) / 2.0d0) * ((d - c) / 2.0d0) * sum1
 
      return
 
    end function gaussQuad2d5

    ! 8) Gaussian Quadrature 3D w/ n=5
    double precision function gaussQuad3d5(a,b,c,d,e,fIn)

      double precision, intent(in) :: a, b, c, d, e, fIn
      double precision :: sum1
      integer :: i,j,k

      sum1 = 0.0d0
      do i=0,4
        do j=0,4
          do k=0,4
            sum1 = sum1 + coeff(i) * coeff(j) * coeff(k) * &
                   f( ((b-a)*r(i) + b + a) / 2.0d0 , &
                      ((d-c)*r(j) + d + c) / 2.0d0 , &
                      ((fIn-e)*r(k) + fIn + e) / 2.0d0 )
          end do
        end do
      end do

      gaussQuad3d5 = ((b - a) / 2.0d0) * ((d - c) / 2.0d0) * ((fIn - e) / 2.0d0) * sum1
 
      return
 
    end function gaussQuad3d5

    ! 9) Romberg Method with Gaussian Quadrature 1D w/ n=5
    function rombergGaussQuad1d5(a,b,n,m)

      integer, intent(in) :: n
      double precision, dimension(0:n,0:n) :: R, rombergGaussQuad1d5
      integer :: i,k,m
      double precision :: a,b,h,aNew, bNew

      ! Initialize romberg array
      R = 0.0d0

      ! Loop to compute intial Gaussian Quadrature 1D approximations:
      do i=0, n
        h = (b-a) / dble(m)

        do k=1, m
          aNew = a + (k-1)*h
          bNew = a + k*h
          R(i,0) = R(i,0) + gaussQuad1d5(aNew,bNew)
        end do
        
        m = m*2
      end do

      ! Loop to compute Romberg approximations from Intial Gaussian approxs:
      do k=1, n
        do i=k, n
          R(i,k) = R(i,k-1) + ( (R(i,k-1) - R(i-1,k-1)) / (4**k - 1) )
        end do
      end do

      rombergGaussQuad1d5 = R

      return

    end function rombergGaussQuad1d5

    ! 10 - Composite Trapezoid 2D
    double precision function compTrapezoid2d(a,b,c,d,n,m)

      double precision, intent(in) :: a, b, c, d
      integer, intent(in) :: n, m
      double precision :: h, k, sum1a,sum1b,sum2
      integer :: i,j

      h = (b-a) / dble(n)
      k = (d-c) / dble(m)

      sum1a = 0.0d0
      sum1b = 0.0d0
      sum2 = 0.0d0

      do i=1, n-1 
        sum1a = sum1a + f(a+dble(i)*h,c) + f(a+dble(i)*h,d)
      end do

      do j=1, m-1
        sum1b = sum1b + f(a,c+dble(j)*k) + f(a,c+dble(j)*k)
      end do      

      do i=1, n-1
        do j=1, m-1
          sum2 = sum2 + f(a+dble(i)*h,c+dble(j)*k)
        end do
      end do

      compTrapezoid2d = (h*k/4.0d0) * ( f(a,c) + f(b,c) + f(a,d) + f(b,d) + &
                                        2.0d0*(sum1a+sum1b) + 4.0d0*sum2)

      return 

    end function compTrapezoid2d

    ! 11 - Romberg Method with Composite Trapezoid 2D
    function rombergCompTrap2d(a,b,c,d,n,m,q)

      integer, intent(in) :: q
      double precision, dimension(0:q-1,0:q-1) :: R, rombergCompTrap2d
      integer :: i,j,n,m
      double precision , intent(in) :: a,b,c,d

      ! Initialize romberg array
      R = 0.0d0

      ! Loop to compute intial Comp Trap 2D  approximations:
      do i=0, q-1
        R(i,0) = compTrapezoid2d(a,b,c,d,n,m)
        n = n*2
        m = m*2
      end do
    
      ! Loop to compute Romberg approximations from intial CompTrap 2d approxs:
      do j=1, q-1
        do i=j, q-1
          R(i,j) = R(i,j-1) + ( (R(i,j-1) - R(i-1,j-1)) / (4**j - 1) )
        end do
      end do

      rombergCompTrap2d = R

    end function rombergCompTrap2d

end module integrationMethods
