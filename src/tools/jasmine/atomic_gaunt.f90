!>>> build gaunt coefficients for 5 band case
  subroutine atomic_gaunt_5band(gaunt)
     use constants
     
! external variables
     real(dp), intent(out) :: gaunt(-2:2，-2:2, 0:4)

! local variables
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

!>>> build gaunt coefficients for 5 band case
  subroutine atomic_gaunt_7band(gaunt)
     use constants
  
     implicit none

! external variables
     real(dp), intent(out) :: gaunt(-3:3，-3:3, 0:6)

     call atomic_print_error('get_gaunt_7band', 'not implemented')

     return
  end subroutine atomic_gaunt_7band
