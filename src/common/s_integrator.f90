    ! 2) Composite Trapezoid Rule:
    double precision function compTrapezoid1d(f, a, b, n)
 
      double precision, intent(in) :: a, b
      integer, intent(in) :: n
      double precision :: f
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

    ! 4) Composite Simpson's Method:
    double precision function compSimpsons1d(f, a, b, n)

      double precision, intent(in) :: a,b
      integer, intent(in) :: n
      double precision :: h, oddSum, evenSum
      double precision :: f
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
