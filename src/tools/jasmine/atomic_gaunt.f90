!-------------------------------------------------------------------------
! project : jasmine
! program : atomic_gaunt_5band
!         : atomic_gaunt_7band
! source  : atomic_gaunt.f90
! type    : subroutines
! author  : yilin wang (email: qhwyl2006@126.com)
! history : 07/09/2014 by yilin wang
! purpose : make gaunt coefficients
! input   :
! output  :
! status  : unstable
! comment :
!-------------------------------------------------------------------------

!>>> build gaunt coefficients for 5 band case
subroutine atomic_gaunt_5band(gaunt)
    use constants, only: dp, zero, one
    
    ! external variables
    real(dp), intent(out) :: gaunt(-2:2, -2:2, 0:4)

    gaunt = zero

    gaunt(-2, -2, 0) = one
    gaunt(-1, -1, 0) = one
    gaunt(0,   0, 0) = one
    gaunt(1,   1, 0) = one
    gaunt(2,   2, 0) = one

    gaunt(-2, -2, 2) = -sqrt(4.0/49.0) 
    gaunt(-2, -1, 2) =  sqrt(6.0/49.0);   gaunt(-1, -2, 2) = gaunt(-2, -1, 2) * (-1)**(-2+1) 
    gaunt(-2,  0, 2) = -sqrt(4.0/49.0);   gaunt(0,  -2, 2) = gaunt(-2,  0, 2) * (-1)**(-2-0)
    gaunt(-1, -1, 2) =  sqrt(1.0/49.0)
    gaunt(-1,  0, 2) =  sqrt(1.0/49.0);   gaunt(0,  -1, 2) = gaunt(-1,  0, 2) * (-1)**(-1-0)
    gaunt(-1,  1, 2) = -sqrt(6.0/49.0);   gaunt(1,  -1, 2) = gaunt(-1,  1, 2) * (-1)**(-1-1)
    gaunt(0,   0, 2) =  sqrt(4.0/49.0)
    gaunt(1,  -1, 2) = -sqrt(6.0/49.0);   gaunt(-1,  1, 2) = gaunt(1,  -1, 2) * (-1)**(1+1)
    gaunt(1,   0, 2) =  sqrt(1.0/49.0);   gaunt(0,   1, 2) = gaunt(1,   0, 2) * (-1)**(1-0)
    gaunt(1,   1, 2) =  sqrt(1.0/49.0)
    gaunt(2,   0, 2) = -sqrt(4.0/49.0);   gaunt(0,   2, 2) = gaunt(2,   0, 2) * (-1)**(2-0)
    gaunt(2,   1, 2) =  sqrt(6.0/49.0);   gaunt(1,   2, 2) = gaunt(2,   1, 2) * (-1)**(2-1)
    gaunt(2,   2, 2) = -sqrt(4.0/49.0)
   
    gaunt(-2, -2, 4) =  sqrt(1.0/441.0)
    gaunt(-2, -1, 4) = -sqrt(5.0/441.0);  gaunt(-1, -2, 4) = gaunt(-2, -1, 4) * (-1)**(-2+1)
    gaunt(-2,  0, 4) =  sqrt(15.0/441.0); gaunt(0,  -2, 4) = gaunt(-2,  0, 4) * (-1)**(-2-0)
    gaunt(-2,  1, 4) = -sqrt(35.0/441.0); gaunt(1,  -2, 4) = gaunt(-2,  1, 4) * (-1)**(-2-1)
    gaunt(-2,  2, 4) =  sqrt(70.0/441.0); gaunt(2,  -2, 4) = gaunt(-2,  2, 4) * (-1)**(-2-2)
    gaunt(-1, -1, 4) = -sqrt(16.0/441.0)
    gaunt(-1,  0, 4) =  sqrt(30.0/441.0); gaunt(0,  -1, 4) = gaunt(-1,  0, 4) * (-1)**(-1-0)
    gaunt(-1,  1, 4) = -sqrt(40.0/441.0); gaunt(1,  -1, 4) = gaunt(-1,  1, 4) * (-1)**(-1-1)
    gaunt( 0,  0, 4) =  sqrt(36.0/441.0)
    gaunt( 1,  0, 4) =  sqrt(30.0/441.0); gaunt(0,   1, 4) = gaunt(1,   0, 4) * (-1)**(1-0)
    gaunt( 1,  1, 4) = -sqrt(16.0/441.0)
    gaunt( 2, -1, 4) = -sqrt(35.0/441.0); gaunt(-1,  2, 4) = gaunt(2,  -1, 4) * (-1)**(2+1)
    gaunt( 2,  0, 4) =  sqrt(15.0/441.0); gaunt( 0,  2, 4) = gaunt(2,   0, 4) * (-1)**(2-0)
    gaunt( 2,  1, 4) = -sqrt(5.0/441.0);  gaunt( 1,  2, 4) = gaunt(2,   1, 4) * (-1)**(2-1)
    gaunt( 2,  2, 4) =  sqrt(1.0/441.0)
   
    return
end subroutine atomic_gaunt_5band

!>>> build gaunt coefficients for 7 band case
subroutine atomic_gaunt_7band(gaunt)
    use constants, only: dp, zero, one
 
    implicit none

    ! external variables
    real(dp), intent(out) :: gaunt(-3:3, -3:3, 0:6)
    gaunt = zero

    gaunt(-3, -3, 0) = one
    gaunt(-2, -2, 0) = one
    gaunt(-1, -1, 0) = one
    gaunt( 0,  0, 0) = one
    gaunt( 1,  1, 0) = one
    gaunt( 2,  2, 0) = one
    gaunt( 3,  3, 0) = one

    gaunt(-3, -3, 2) = -sqrt(25.0/225.0)
    gaunt(-3, -2, 2) =  sqrt(25.0/225.0);  gaunt(-2, -3, 2) = gaunt(-3, -2, 2) * (-1.0)**(-3+2)
    gaunt(-3, -1, 2) = -sqrt(10.0/225.0);  gaunt(-1, -3, 2) = gaunt(-3, -1, 2) * (-1.0)**(-3+1)
    gaunt(-2, -1, 2) =  sqrt(15.0/225.0);  gaunt(-1, -2, 2) = gaunt(-2, -1, 2) * (-1.0)**(-2+1)
    gaunt(-2,  0, 2) = -sqrt(20.0/225.0);  gaunt( 0, -2, 2) = gaunt(-2,  0, 2) * (-1.0)**(-2-0)
    gaunt(-1, -1, 2) =  sqrt( 9.0/225.0)
    gaunt(-1,  0, 2) =  sqrt( 2.0/225.0);  gaunt( 0, -1, 2) = gaunt(-1,  0, 2) * (-1.0)**(-1-0)
    gaunt(-1,  1, 2) = -sqrt(24.0/225.0);  gaunt( 1, -1, 2) = gaunt(-1,  1, 2) * (-1.0)**(-1-1)
    gaunt( 0,  0, 2) =  sqrt(16.0/225.0)
    gaunt( 1,  0, 2) =  sqrt( 2.0/225.0);  gaunt( 0,  1, 2) = gaunt( 1,  0, 2) * (-1.0)**( 1-0)
    gaunt( 1,  1, 2) =  sqrt( 9.0/225.0)
    gaunt( 2,  0, 2) = -sqrt(20.0/225.0);  gaunt( 0,  2, 2) = gaunt( 2,  0, 2) * (-1.0)**( 2-0)
    gaunt( 2,  1, 2) =  sqrt(15.0/225.0);  gaunt( 1,  2, 2) = gaunt( 2,  1, 2) * (-1.0)**( 2-1)
    gaunt( 3,  1, 2) = -sqrt(10.0/225.0);  gaunt( 1,  3, 2) = gaunt( 3,  1, 2) * (-1.0)**( 3-1)
    gaunt( 3,  2, 2) =  sqrt(25.0/225.0);  gaunt( 2,  3, 2) = gaunt( 3,  2, 2) * (-1.0)**( 3-2)
    gaunt( 3,  3, 2) = -sqrt(25.0/225.0)

    gaunt(-3, -3, 4) =  sqrt(9.0/1089.0)
    gaunt(-3, -2, 4) = -sqrt(30.0/1089.0); gaunt(-2, -3, 4) = gaunt(-3, -2, 4) * (-1.0)**(-3+2)
    gaunt(-3, -1, 4) =  sqrt(54.0/1089.0); gaunt(-1, -3, 4) = gaunt(-3, -1, 4) * (-1.0)**(-3+1)
    gaunt(-3,  0, 4) = -sqrt(63.0/1089.0); gaunt( 0, -3, 4) = gaunt(-3,  0, 4) * (-1.0)**(-3-0)
    gaunt(-3,  1, 4) =  sqrt(42.0/1089.0); gaunt( 1, -3, 4) = gaunt(-3,  1, 4) * (-1.0)**(-3-1)
    gaunt(-2, -2, 4) = -sqrt(49.0/1089.0)
    gaunt(-2, -1, 4) =  sqrt(32.0/1089.0); gaunt(-1, -2, 4) = gaunt(-2, -1, 4) * (-1.0)**(-2+1)
    gaunt(-2,  0, 4) = -sqrt( 3.0/1089.0); gaunt( 0, -2, 4) = gaunt(-2,  0, 4) * (-1.0)**(-2-0)
    gaunt(-2,  1, 4) = -sqrt(14.0/1089.0); gaunt( 1, -2, 4) = gaunt(-2,  1, 4) * (-1.0)**(-2-1)
    gaunt(-2,  2, 4) =  sqrt(70.0/1089.0); gaunt( 2, -2, 4) = gaunt(-2,  2, 4) * (-1.0)**(-2-2)
    gaunt(-1, -1, 4) =  sqrt( 1.0/1089.0)
    gaunt(-1,  0, 4) =  sqrt(15.0/1089.0); gaunt( 0, -1, 4) = gaunt(-1,  0, 4) * (-1.0)**(-1-0)
    gaunt(-1,  1, 4) = -sqrt(40.0/1089.0); gaunt( 1, -1, 4) = gaunt(-1,  1, 4) * (-1.0)**(-1-1)
    gaunt( 0,  0, 4) =  sqrt(36.0/1089.0)
    gaunt( 1,  0, 4) =  sqrt(15.0/1089.0); gaunt( 0,  1, 4) = gaunt( 1,  0, 4) * (-1.0)**( 1-0)
    gaunt( 1,  1, 4) =  sqrt( 1.0/1089.0)
    gaunt( 2, -1, 4) = -sqrt(14.0/1089.0); gaunt(-1,  2, 4) = gaunt( 2, -1, 4) * (-1.0)**( 2+1)
    gaunt( 2,  0, 4) = -sqrt( 3.0/1089.0); gaunt( 0,  2, 4) = gaunt( 2,  0, 4) * (-1.0)**( 2-0)
    gaunt( 2,  1, 4) =  sqrt(32.0/1089.0); gaunt( 1,  2, 4) = gaunt( 2,  1, 4) * (-1.0)**( 2-1)
    gaunt( 2,  2, 4) = -sqrt(49.0/1089.0)
    gaunt( 3, -1, 4) =  sqrt(42.0/1089.0); gaunt(-1,  3, 4) = gaunt( 3, -1, 4) * (-1.0)**( 3+1)
    gaunt( 3,  0, 4) = -sqrt(63.0/1089.0); gaunt( 0,  3, 4) = gaunt( 3,  0, 4) * (-1.0)**( 3-0)
    gaunt( 3,  1, 4) =  sqrt(54.0/1089.0); gaunt( 1,  3, 4) = gaunt( 3,  1, 4) * (-1.0)**( 3-1)
    gaunt( 3,  2, 4) = -sqrt(30.0/1089.0); gaunt( 2,  3, 4) = gaunt( 3,  2, 4) * (-1.0)**( 3-2)
    gaunt( 3,  3, 4) =  sqrt( 9.0/1089.0)

    gaunt(-3, -3, 6) = -sqrt(   25.0/184041)
    gaunt(-3, -2, 6) =  sqrt(  175.0/184041); gaunt(-2, -3, 6) = gaunt(-3, -2, 6) * (-1.0)**(-3+2)
    gaunt(-3, -1, 6) = -sqrt(  700.0/184041); gaunt(-1, -3, 6) = gaunt(-3, -1, 6) * (-1.0)**(-3+1)
    gaunt(-3,  0, 6) =  sqrt( 2100.0/184041); gaunt( 0, -3, 6) = gaunt(-3,  0, 6) * (-1.0)**(-3-0)
    gaunt(-3,  1, 6) = -sqrt( 5250.0/184041); gaunt( 1, -3, 6) = gaunt(-3,  1, 6) * (-1.0)**(-3-1)
    gaunt(-3,  2, 6) =  sqrt(11550.0/184041); gaunt( 2, -3, 6) = gaunt(-3,  2, 6) * (-1.0)**(-3-2)
    gaunt(-3,  3, 6) = -sqrt(23100.0/184041); gaunt( 3, -3, 6) = gaunt(-3,  3, 6) * (-1.0)**(-3-3)
    gaunt(-2, -2, 6) = -sqrt(  900.0/184041)
    gaunt(-2, -1, 6) = -sqrt( 2625.0/184041); gaunt(-1, -2, 6) = gaunt(-2, -1, 6) * (-1.0)**(-2+1)
    gaunt(-2,  0, 6) =  sqrt( 5600.0/184041); gaunt( 0, -2, 6) = gaunt(-2,  0, 6) * (-1.0)**(-2-0)
    gaunt(-2,  1, 6) = -sqrt( 9450.0/184041); gaunt( 1, -2, 6) = gaunt(-2,  1, 6) * (-1.0)**(-2-1)
    gaunt(-2,  2, 6) =  sqrt(12600.0/184041); gaunt( 2, -2, 6) = gaunt(-2,  2, 6) * (-1.0)**(-2-2)
    gaunt(-1, -1, 6) = -sqrt( 5625.0/184041)
    gaunt(-1,  0, 6) =  sqrt( 8750.0/184041); gaunt( 0, -1, 6) = gaunt(-1,  0, 6) * (-1.0)**(-1-0)
    gaunt(-1,  1, 6) = -sqrt(10500.0/184041); gaunt( 1, -1, 6) = gaunt(-1,  1, 6) * (-1.0)**(-1-1)
    gaunt( 0,  0, 6) =  sqrt(10000.0/184041)
    gaunt( 1,  0, 6) =  sqrt( 8750.0/184041); gaunt( 0,  1, 6) = gaunt( 1,  0, 6) * (-1.0)**( 1-0)
    gaunt( 1,  1, 6) = -sqrt( 5625.0/184041)
    gaunt( 2, -1, 6) = -sqrt( 9450.0/184041); gaunt(-1,  2, 6) = gaunt( 2, -1, 6) * (-1.0)**( 2+1)
    gaunt( 2,  0, 6) =  sqrt( 5600.0/184041); gaunt( 0,  2, 6) = gaunt( 2,  0, 6) * (-1.0)**( 2-0)
    gaunt( 2,  1, 6) = -sqrt( 2625.0/184041); gaunt( 1,  2, 6) = gaunt( 2,  1, 6) * (-1.0)**( 2-1)
    gaunt( 2,  2, 6) =  sqrt(  900.0/184041)
    gaunt( 3, -2, 6) =  sqrt(11550.0/184041); gaunt(-2,  3, 6) = gaunt( 3, -2, 6) * (-1.0)**( 3+2)
    gaunt( 3, -1, 6) = -sqrt( 5250.0/184041); gaunt(-1,  3, 6) = gaunt( 3, -1, 6) * (-1.0)**( 3+1)
    gaunt( 3,  0, 6) =  sqrt( 2100.0/184041); gaunt( 0,  3, 6) = gaunt( 3,  0, 6) * (-1.0)**( 3-0)
    gaunt( 3,  1, 6) = -sqrt(  700.0/184041); gaunt( 1,  3, 6) = gaunt( 3,  1, 6) * (-1.0)**( 3-1)
    gaunt( 3,  2, 6) =  sqrt(  175.0/184041); gaunt( 2,  3, 6) = gaunt( 3,  2, 6) * (-1.0)**( 3-2)
    gaunt( 3,  3, 6) = -sqrt(   25.0/184041)

    return
end subroutine atomic_gaunt_7band
