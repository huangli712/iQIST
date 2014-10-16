!!!-------------------------------------------------------------------------
!!! project : jasmine
!!! program : atomic_make_cdagger
!!!           atomic_make_c
!!!           atomic_make_gsz
!!!           atomic_make_gjz
!!!           atomic_make_gaunt5
!!!           atomic_make_gaunt7
!!!           atomic_make_cumatK
!!!           atomic_make_cumatS
!!!           atomic_make_socmat3
!!!           atomic_make_socmat5
!!!           atomic_make_socmat7
!!!           atomic_make_tmat_c2r
!!!           atomic_make_tmat_r2c
!!!           atomic_make_tmat_c2j
!!!           atomic_tran_cumat
!!!           atomic_tran_repr_cmpl
!!!           atomic_tran_repr_real
!!! source  : atomic_gaunt.f90
!!! type    : subroutines
!!! author  : yilin wang (email: qhwyl2006@126.com)
!!! history : 07/09/2014 by yilin wang
!!!           08/22/2014 by yilin wang
!!! purpose : make gaunt coefficients
!!! purpose : make transformation from one representation to another 
!!!           representation
!!! status  : unstable
!!! comment :
!!!-------------------------------------------------------------------------

!!>>> atomic_make_cdagger: create one electron on ipos 
!!>>> of |jold> to deduce |jnew>
  subroutine atomic_make_cdagger(ipos, jold, jnew, isgn)
     implicit none
  
! external argument
! position number (serial number of orbit)
     integer, intent(in) :: ipos
  
! old Fock state and new Fock state
     integer, intent(in ):: jold
     integer, intent(out):: jnew
  
! sgn due to anti-commute relation between fernions
     integer, intent(out):: isgn
  
! local variables
! loop index over orbit
     integer :: iorb
  
     if (btest(jold, ipos-1) .eqv. .true.) then
         call s_print_error("atomic_construct", "severe error happened")
     endif
  
     isgn = 0
     do iorb=1,ipos-1
        if (btest(jold, iorb-1)) isgn = isgn + 1
     enddo
     isgn = mod(isgn, 2)
  
     isgn = (-1)**isgn
     jnew = jold + 2**(ipos-1)
  
     return
  end subroutine atomic_make_cdagger

!!>>> atomic_make_c: destroy one electron on ipos 
!!>>> of |jold> to deduce |jnew>
  subroutine atomic_make_c(ipos, jold, jnew, isgn)
      implicit none
  
! external argument
! position number (serial number of orbit)
      integer, intent(in)  :: ipos
  
! old Fock state and new Fock state
      integer, intent(in ) :: jold
      integer, intent(out) :: jnew
  
! sgn due to anti-commute relation between fernions
      integer, intent(out) :: isgn
  
! local variables
! loop index
      integer :: iorb
  
      if (btest(jold, ipos-1) .eqv. .false.) then
          call s_print_error("atomic_eliminate", "severe error happened")
      endif 
  
      isgn = 0
      do iorb=1,ipos-1
          if (btest(jold, iorb-1)) isgn = isgn + 1
      enddo
      isgn = mod(isgn, 2)
  
      isgn = (-1)**isgn
      jnew = jold - 2**(ipos-1)
  
      return
  end subroutine atomic_make_c

!!>>> atomic_make_gsz: make sz for each orbital
  subroutine atomic_make_gsz(good_sz)
     use control, only : norbs
  
     implicit none
  
! external variables
     integer, intent(out) :: good_sz(norbs)

! local variables
     integer :: i
  
     do i=1, norbs
         if (mod(i,2) /= 0 ) then
             good_sz(i) = 1
         else
             good_sz(i) = -1
         endif
     enddo
  
     return
  end subroutine atomic_make_gsz
  
!>>> atomic_make_gjz: make jz for each orbital
  subroutine atomic_make_gjz(good_jz)
     use control, only : nband, norbs
  
     implicit none
  
! external variables
     integer, intent(out) :: good_jz(norbs)
  
     if (nband == 3) then
! j=1/2
         good_jz(1) = -1
         good_jz(2) =  1
! j=3/2
         good_jz(3) = -3
         good_jz(4) = -1
         good_jz(5) =  1
         good_jz(6) =  3
     elseif (nband == 5) then
! j=3/2
         good_jz(1) = -3
         good_jz(2) = -1
         good_jz(3) =  1
         good_jz(4) =  3
! j=5/2
         good_jz(5) = -5
         good_jz(6) = -3
         good_jz(7) = -1
         good_jz(8) =  1
         good_jz(9) =  3
         good_jz(10)=  5
     elseif (nband == 7) then
! j=5/2
         good_jz(1) = -5
         good_jz(2) = -3
         good_jz(3) = -1
         good_jz(4) =  1
         good_jz(5) =  3
         good_jz(6) =  5
! j=7/2
         good_jz(7) = -7
         good_jz(8) = -5
         good_jz(9) = -3
         good_jz(10)= -1
         good_jz(11)=  1
         good_jz(12)=  3
         good_jz(13)=  5
         good_jz(14)=  7
     else
         call s_print_error('atomic_make_good_jz', &
            'not implemented for this norbs value !')
     endif
  
     return
  end subroutine atomic_make_gjz

!!>>> atomic_make_gaunt5: build gaunt coefficients for 5 band case
  subroutine atomic_make_gaunt5(gaunt)
     use constants, only : dp, zero, one
     
! external variables
     real(dp), intent(out) :: gaunt(-2:2, -2:2, 0:4)
  
     gaunt = zero
  
     gaunt(-2, -2, 0) = one
     gaunt(-1, -1, 0) = one
     gaunt(0,   0, 0) = one
     gaunt(1,   1, 0) = one
     gaunt(2,   2, 0) = one
  
     gaunt(-2, -2, 2) = -sqrt(4.0_dp/49.0_dp) 
     gaunt(-2, -1, 2) =  sqrt(6.0_dp/49.0_dp);   gaunt(-1, -2, 2) = gaunt(-2, -1, 2) * (-1)**(-2+1) 
     gaunt(-2,  0, 2) = -sqrt(4.0_dp/49.0_dp);   gaunt(0,  -2, 2) = gaunt(-2,  0, 2) * (-1)**(-2-0)
     gaunt(-1, -1, 2) =  sqrt(1.0_dp/49.0_dp)
     gaunt(-1,  0, 2) =  sqrt(1.0_dp/49.0_dp);   gaunt(0,  -1, 2) = gaunt(-1,  0, 2) * (-1)**(-1-0)
     gaunt(-1,  1, 2) = -sqrt(6.0_dp/49.0_dp);   gaunt(1,  -1, 2) = gaunt(-1,  1, 2) * (-1)**(-1-1)
     gaunt(0,   0, 2) =  sqrt(4.0_dp/49.0_dp)
     gaunt(1,  -1, 2) = -sqrt(6.0_dp/49.0_dp);   gaunt(-1,  1, 2) = gaunt(1,  -1, 2) * (-1)**(1+1)
     gaunt(1,   0, 2) =  sqrt(1.0_dp/49.0_dp);   gaunt(0,   1, 2) = gaunt(1,   0, 2) * (-1)**(1-0)
     gaunt(1,   1, 2) =  sqrt(1.0_dp/49.0_dp)
     gaunt(2,   0, 2) = -sqrt(4.0_dp/49.0_dp);   gaunt(0,   2, 2) = gaunt(2,   0, 2) * (-1)**(2-0)
     gaunt(2,   1, 2) =  sqrt(6.0_dp/49.0_dp);   gaunt(1,   2, 2) = gaunt(2,   1, 2) * (-1)**(2-1)
     gaunt(2,   2, 2) = -sqrt(4.0_dp/49.0_dp)
     
     gaunt(-2, -2, 4) =  sqrt( 1.0_dp/441.0_dp)
     gaunt(-2, -1, 4) = -sqrt( 5.0_dp/441.0_dp); gaunt(-1, -2, 4) = gaunt(-2, -1, 4) * (-1)**(-2+1)
     gaunt(-2,  0, 4) =  sqrt(15.0_dp/441.0_dp); gaunt(0,  -2, 4) = gaunt(-2,  0, 4) * (-1)**(-2-0)
     gaunt(-2,  1, 4) = -sqrt(35.0_dp/441.0_dp); gaunt(1,  -2, 4) = gaunt(-2,  1, 4) * (-1)**(-2-1)
     gaunt(-2,  2, 4) =  sqrt(70.0_dp/441.0_dp); gaunt(2,  -2, 4) = gaunt(-2,  2, 4) * (-1)**(-2-2)
     gaunt(-1, -1, 4) = -sqrt(16.0_dp/441.0_dp)
     gaunt(-1,  0, 4) =  sqrt(30.0_dp/441.0_dp); gaunt(0,  -1, 4) = gaunt(-1,  0, 4) * (-1)**(-1-0)
     gaunt(-1,  1, 4) = -sqrt(40.0_dp/441.0_dp); gaunt(1,  -1, 4) = gaunt(-1,  1, 4) * (-1)**(-1-1)
     gaunt( 0,  0, 4) =  sqrt(36.0_dp/441.0_dp)
     gaunt( 1,  0, 4) =  sqrt(30.0_dp/441.0_dp); gaunt(0,   1, 4) = gaunt(1,   0, 4) * (-1)**(1-0)
     gaunt( 1,  1, 4) = -sqrt(16.0_dp/441.0_dp)
     gaunt( 2, -1, 4) = -sqrt(35.0_dp/441.0_dp); gaunt(-1,  2, 4) = gaunt(2,  -1, 4) * (-1)**(2+1)
     gaunt( 2,  0, 4) =  sqrt(15.0_dp/441.0_dp); gaunt( 0,  2, 4) = gaunt(2,   0, 4) * (-1)**(2-0)
     gaunt( 2,  1, 4) = -sqrt( 5.0_dp/441.0_dp); gaunt( 1,  2, 4) = gaunt(2,   1, 4) * (-1)**(2-1)
     gaunt( 2,  2, 4) =  sqrt( 1.0_dp/441.0_dp)
     
     return
  end subroutine atomic_make_gaunt5

!!>>> atomic_make_gaunt7: build gaunt coefficients for 7 band case
  subroutine atomic_make_gaunt7(gaunt)
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
  
     gaunt(-3, -3, 2) = -sqrt(25.0_dp/225.0_dp)
     gaunt(-3, -2, 2) =  sqrt(25.0_dp/225.0_dp);  gaunt(-2, -3, 2) = gaunt(-3, -2, 2) * (-1.0)**(-3+2)
     gaunt(-3, -1, 2) = -sqrt(10.0_dp/225.0_dp);  gaunt(-1, -3, 2) = gaunt(-3, -1, 2) * (-1.0)**(-3+1)
     gaunt(-2, -1, 2) =  sqrt(15.0_dp/225.0_dp);  gaunt(-1, -2, 2) = gaunt(-2, -1, 2) * (-1.0)**(-2+1)
     gaunt(-2,  0, 2) = -sqrt(20.0_dp/225.0_dp);  gaunt( 0, -2, 2) = gaunt(-2,  0, 2) * (-1.0)**(-2-0)
     gaunt(-1, -1, 2) =  sqrt( 9.0_dp/225.0_dp)
     gaunt(-1,  0, 2) =  sqrt( 2.0_dp/225.0_dp);  gaunt( 0, -1, 2) = gaunt(-1,  0, 2) * (-1.0)**(-1-0)
     gaunt(-1,  1, 2) = -sqrt(24.0_dp/225.0_dp);  gaunt( 1, -1, 2) = gaunt(-1,  1, 2) * (-1.0)**(-1-1)
     gaunt( 0,  0, 2) =  sqrt(16.0_dp/225.0_dp)
     gaunt( 1,  0, 2) =  sqrt( 2.0_dp/225.0_dp);  gaunt( 0,  1, 2) = gaunt( 1,  0, 2) * (-1.0)**( 1-0)
     gaunt( 1,  1, 2) =  sqrt( 9.0_dp/225.0_dp)
     gaunt( 2,  0, 2) = -sqrt(20.0_dp/225.0_dp);  gaunt( 0,  2, 2) = gaunt( 2,  0, 2) * (-1.0)**( 2-0)
     gaunt( 2,  1, 2) =  sqrt(15.0_dp/225.0_dp);  gaunt( 1,  2, 2) = gaunt( 2,  1, 2) * (-1.0)**( 2-1)
     gaunt( 3,  1, 2) = -sqrt(10.0_dp/225.0_dp);  gaunt( 1,  3, 2) = gaunt( 3,  1, 2) * (-1.0)**( 3-1)
     gaunt( 3,  2, 2) =  sqrt(25.0_dp/225.0_dp);  gaunt( 2,  3, 2) = gaunt( 3,  2, 2) * (-1.0)**( 3-2)
     gaunt( 3,  3, 2) = -sqrt(25.0_dp/225.0_dp)
  
     gaunt(-3, -3, 4) =  sqrt( 9.0_dp/1089.0_dp)
     gaunt(-3, -2, 4) = -sqrt(30.0_dp/1089.0_dp); gaunt(-2, -3, 4) = gaunt(-3, -2, 4) * (-1.0)**(-3+2)
     gaunt(-3, -1, 4) =  sqrt(54.0_dp/1089.0_dp); gaunt(-1, -3, 4) = gaunt(-3, -1, 4) * (-1.0)**(-3+1)
     gaunt(-3,  0, 4) = -sqrt(63.0_dp/1089.0_dp); gaunt( 0, -3, 4) = gaunt(-3,  0, 4) * (-1.0)**(-3-0)
     gaunt(-3,  1, 4) =  sqrt(42.0_dp/1089.0_dp); gaunt( 1, -3, 4) = gaunt(-3,  1, 4) * (-1.0)**(-3-1)
     gaunt(-2, -2, 4) = -sqrt(49.0_dp/1089.0_dp)
     gaunt(-2, -1, 4) =  sqrt(32.0_dp/1089.0_dp); gaunt(-1, -2, 4) = gaunt(-2, -1, 4) * (-1.0)**(-2+1)
     gaunt(-2,  0, 4) = -sqrt( 3.0_dp/1089.0_dp); gaunt( 0, -2, 4) = gaunt(-2,  0, 4) * (-1.0)**(-2-0)
     gaunt(-2,  1, 4) = -sqrt(14.0_dp/1089.0_dp); gaunt( 1, -2, 4) = gaunt(-2,  1, 4) * (-1.0)**(-2-1)
     gaunt(-2,  2, 4) =  sqrt(70.0_dp/1089.0_dp); gaunt( 2, -2, 4) = gaunt(-2,  2, 4) * (-1.0)**(-2-2)
     gaunt(-1, -1, 4) =  sqrt( 1.0_dp/1089.0_dp)
     gaunt(-1,  0, 4) =  sqrt(15.0_dp/1089.0_dp); gaunt( 0, -1, 4) = gaunt(-1,  0, 4) * (-1.0)**(-1-0)
     gaunt(-1,  1, 4) = -sqrt(40.0_dp/1089.0_dp); gaunt( 1, -1, 4) = gaunt(-1,  1, 4) * (-1.0)**(-1-1)
     gaunt( 0,  0, 4) =  sqrt(36.0_dp/1089.0_dp)
     gaunt( 1,  0, 4) =  sqrt(15.0_dp/1089.0_dp); gaunt( 0,  1, 4) = gaunt( 1,  0, 4) * (-1.0)**( 1-0)
     gaunt( 1,  1, 4) =  sqrt( 1.0_dp/1089.0_dp)
     gaunt( 2, -1, 4) = -sqrt(14.0_dp/1089.0_dp); gaunt(-1,  2, 4) = gaunt( 2, -1, 4) * (-1.0)**( 2+1)
     gaunt( 2,  0, 4) = -sqrt( 3.0_dp/1089.0_dp); gaunt( 0,  2, 4) = gaunt( 2,  0, 4) * (-1.0)**( 2-0)
     gaunt( 2,  1, 4) =  sqrt(32.0_dp/1089.0_dp); gaunt( 1,  2, 4) = gaunt( 2,  1, 4) * (-1.0)**( 2-1)
     gaunt( 2,  2, 4) = -sqrt(49.0_dp/1089.0_dp)
     gaunt( 3, -1, 4) =  sqrt(42.0_dp/1089.0_dp); gaunt(-1,  3, 4) = gaunt( 3, -1, 4) * (-1.0)**( 3+1)
     gaunt( 3,  0, 4) = -sqrt(63.0_dp/1089.0_dp); gaunt( 0,  3, 4) = gaunt( 3,  0, 4) * (-1.0)**( 3-0)
     gaunt( 3,  1, 4) =  sqrt(54.0_dp/1089.0_dp); gaunt( 1,  3, 4) = gaunt( 3,  1, 4) * (-1.0)**( 3-1)
     gaunt( 3,  2, 4) = -sqrt(30.0_dp/1089.0_dp); gaunt( 2,  3, 4) = gaunt( 3,  2, 4) * (-1.0)**( 3-2)
     gaunt( 3,  3, 4) =  sqrt( 9.0_dp/1089.0_dp)
  
     gaunt(-3, -3, 6) = -sqrt(   25.0_dp/184041_dp)
     gaunt(-3, -2, 6) =  sqrt(  175.0_dp/184041_dp); gaunt(-2, -3, 6) = gaunt(-3, -2, 6) * (-1.0)**(-3+2)
     gaunt(-3, -1, 6) = -sqrt(  700.0_dp/184041_dp); gaunt(-1, -3, 6) = gaunt(-3, -1, 6) * (-1.0)**(-3+1)
     gaunt(-3,  0, 6) =  sqrt( 2100.0_dp/184041_dp); gaunt( 0, -3, 6) = gaunt(-3,  0, 6) * (-1.0)**(-3-0)
     gaunt(-3,  1, 6) = -sqrt( 5250.0_dp/184041_dp); gaunt( 1, -3, 6) = gaunt(-3,  1, 6) * (-1.0)**(-3-1)
     gaunt(-3,  2, 6) =  sqrt(11550.0_dp/184041_dp); gaunt( 2, -3, 6) = gaunt(-3,  2, 6) * (-1.0)**(-3-2)
     gaunt(-3,  3, 6) = -sqrt(23100.0_dp/184041_dp); gaunt( 3, -3, 6) = gaunt(-3,  3, 6) * (-1.0)**(-3-3)
     gaunt(-2, -2, 6) =  sqrt(  900.0_dp/184041_dp)
     gaunt(-2, -1, 6) = -sqrt( 2625.0_dp/184041_dp); gaunt(-1, -2, 6) = gaunt(-2, -1, 6) * (-1.0)**(-2+1)
     gaunt(-2,  0, 6) =  sqrt( 5600.0_dp/184041_dp); gaunt( 0, -2, 6) = gaunt(-2,  0, 6) * (-1.0)**(-2-0)
     gaunt(-2,  1, 6) = -sqrt( 9450.0_dp/184041_dp); gaunt( 1, -2, 6) = gaunt(-2,  1, 6) * (-1.0)**(-2-1)
     gaunt(-2,  2, 6) =  sqrt(12600.0_dp/184041_dp); gaunt( 2, -2, 6) = gaunt(-2,  2, 6) * (-1.0)**(-2-2)
     gaunt(-1, -1, 6) = -sqrt( 5625.0_dp/184041_dp)
     gaunt(-1,  0, 6) =  sqrt( 8750.0_dp/184041_dp); gaunt( 0, -1, 6) = gaunt(-1,  0, 6) * (-1.0)**(-1-0)
     gaunt(-1,  1, 6) = -sqrt(10500.0_dp/184041_dp); gaunt( 1, -1, 6) = gaunt(-1,  1, 6) * (-1.0)**(-1-1)
     gaunt( 0,  0, 6) =  sqrt(10000.0_dp/184041_dp)
     gaunt( 1,  0, 6) =  sqrt( 8750.0_dp/184041_dp); gaunt( 0,  1, 6) = gaunt( 1,  0, 6) * (-1.0)**( 1-0)
     gaunt( 1,  1, 6) = -sqrt( 5625.0_dp/184041_dp)
     gaunt( 2, -1, 6) = -sqrt( 9450.0_dp/184041_dp); gaunt(-1,  2, 6) = gaunt( 2, -1, 6) * (-1.0)**( 2+1)
     gaunt( 2,  0, 6) =  sqrt( 5600.0_dp/184041_dp); gaunt( 0,  2, 6) = gaunt( 2,  0, 6) * (-1.0)**( 2-0)
     gaunt( 2,  1, 6) = -sqrt( 2625.0_dp/184041_dp); gaunt( 1,  2, 6) = gaunt( 2,  1, 6) * (-1.0)**( 2-1)
     gaunt( 2,  2, 6) =  sqrt(  900.0_dp/184041_dp)
     gaunt( 3, -2, 6) =  sqrt(11550.0_dp/184041_dp); gaunt(-2,  3, 6) = gaunt( 3, -2, 6) * (-1.0)**( 3+2)
     gaunt( 3, -1, 6) = -sqrt( 5250.0_dp/184041_dp); gaunt(-1,  3, 6) = gaunt( 3, -1, 6) * (-1.0)**( 3+1)
     gaunt( 3,  0, 6) =  sqrt( 2100.0_dp/184041_dp); gaunt( 0,  3, 6) = gaunt( 3,  0, 6) * (-1.0)**( 3-0)
     gaunt( 3,  1, 6) = -sqrt(  700.0_dp/184041_dp); gaunt( 1,  3, 6) = gaunt( 3,  1, 6) * (-1.0)**( 3-1)
     gaunt( 3,  2, 6) =  sqrt(  175.0_dp/184041_dp); gaunt( 2,  3, 6) = gaunt( 3,  2, 6) * (-1.0)**( 3-2)
     gaunt( 3,  3, 6) = -sqrt(   25.0_dp/184041_dp)

     return
  end subroutine atomic_make_gaunt7

!!>>> atomic_make_cumatK: make Coulomb interaction U according to
!!>>> Kanamori parameterized Hamiltonian
  subroutine atomic_make_cumatK()
     use constants, only : dp, czero, zero
     use control, only : norbs, Uc, Uv, Jz, Js, Jp

     use m_spmat, only : cumat
  
     implicit none
  
! local varibales
! orbital index
     integer :: alpha, betta
     integer :: delta, gamma

! band index and spin index
     integer :: aband, bband 
     integer :: dband, gband 
     integer :: aspin, bspin 
     integer :: dspin, gspin 

! dummy variables
     real(dp) :: dtmp
  
! initialize cumat to zero
     cumat = czero
  
! loop for creation operators
     alphaloop: do alpha=1,norbs-1
     bettaloop: do betta=alpha+1,norbs
  
! loop for annihilation operators
        gammaloop: do gamma=1,norbs-1
        deltaloop: do delta=gamma+1,norbs
            aband = (alpha+1)/2; aspin = mod(alpha,2)
            bband = (betta+1)/2; bspin = mod(betta,2)
            gband = (gamma+1)/2; gspin = mod(gamma,2)
            dband = (delta+1)/2; dspin = mod(delta,2)
  
            dtmp = zero
  
! intraorbital Coulomb interaction
            if ((alpha.eq.gamma) .and. (betta.eq.delta)) then
                if ((aband.eq.bband) .and. (aspin.ne.bspin)) then
                    dtmp = dtmp + Uc
                endif
            endif
  
! interorbital Coulomb interaction
            if ((alpha.eq.gamma) .and. (betta.eq.delta)) then
                if (aband .ne. bband) then
                    dtmp = dtmp + Uv
                endif
            endif
  
! Hund's exchange interaction 
            if ((alpha.eq.gamma) .and. (betta.eq.delta)) then
                if ((aband.ne.bband) .and. (aspin.eq.bspin)) then
                    dtmp = dtmp - Jz
                endif
            endif
           
! spin flip term
            if ((aband.eq.gband) .and. (bband.eq.dband)) then
                if ((aspin.ne.gspin) .and. (bspin.ne.dspin) .and. (aspin.ne.bspin)) then
                    dtmp = dtmp - Js
                endif
            endif
          
! pair hopping term
            if ((aband.eq.bband) .and. (dband.eq.gband) .and. (aband.ne.dband)) then
                if ((aspin.ne.bspin) .and. (dspin.ne.gspin) .and. (aspin.eq.gspin)) then
                    dtmp = dtmp + Jp
                endif
            endif
                 
            cumat(alpha, betta, delta, gamma) = dtmp
  
        enddo deltaloop ! over delta={gamma+1,norbs} loop
        enddo gammaloop ! over gamma={1,norbs-1} loop
     enddo bettaloop ! over betta={alpha+1,norbs} loop
     enddo alphaloop ! over alpha={1,norbs-1} loop
  
     return
  end subroutine atomic_make_cumatK

!!>>> atomic_make_cumatS: make Coulomb interation U, according to 
!!>>> Slater-Cordon parameterized Hamiltonian
  subroutine atomic_make_cumatS()
     use constants, only : dp, zero, half
     use control, only : nband, norbs, F0, F2, F4, F6
     use m_spmat, only : cumat
  
     implicit none
  
! local variables
! Slater-Cordon parameters
     real(dp), allocatable :: slater_cordon(:)

! gaunt coefficients
     real(dp), allocatable :: gaunt(:,:,:)
  
! orbital momentum quantum number
     integer :: l

! loop index
     integer :: i
     integer :: alpha, betta
     integer :: delta, gamma
     integer :: aband, aspin
     integer :: bband, bspin
     integer :: dband, dspin
     integer :: gband, gspin

! dummy variables
     real(dp) :: res
  
  
! allocate memory for slater_cordon and gaunt and then build them
     if (nband == 5) then
         l = 2
         allocate(slater_cordon(0:2*l))     
         slater_cordon = zero
         slater_cordon(0) = F0
         slater_cordon(2) = F2
         slater_cordon(4) = F4
         allocate(gaunt(-l:l, -l:l, 0:2*l))
         call atomic_make_gaunt5(gaunt) 
     elseif(nband == 7) then
         l = 3
         allocate(slater_cordon(0:2*l))     
         slater_cordon = zero
         slater_cordon(0) = F0
         slater_cordon(2) = F2
         slater_cordon(4) = F4
         slater_cordon(6) = F6
         allocate(gaunt(-l:l, -l:l, 0:2*l))
         call atomic_make_gaunt7(gaunt)
     else
         call s_print_error('atomic_make_cumat_slater', 'not implemented for this nband!')
     endif
  
! make Coulomb interaction U matrix
     do alpha=1,norbs
     do betta=1,norbs
         aband = (alpha-1)/2-l
         bband = (betta-1)/2-l
         aspin = mod(alpha, 2)
         bspin = mod(betta, 2)
  
         do delta=1,norbs
         do gamma=1,norbs
             dband = (delta-1)/2-l
             gband = (gamma-1)/2-l
             dspin = mod(delta, 2)
             gspin = mod(gamma, 2)
  
             if ((alpha .eq. betta) .or. (delta .eq. gamma)) cycle
  
             if ((aband + bband) .ne. (dband + gband)) cycle
             if ((aspin .ne. gspin) .or. (bspin .ne. dspin)) cycle
  
             res = zero
             do i=0, 2*l, 2
                 res = res + gaunt(aband, gband, i) * gaunt(dband, bband, i) * slater_cordon(i)
             enddo
             cumat(alpha, betta, delta, gamma) = res
         enddo ! over gamma={1,norbs} loop
         enddo ! over delta={1,norbs} loop
  
     enddo ! over betta={1,norbs} loop
     enddo ! over alpha={1,norbs} loop
  
     cumat = half * cumat
  
! deallocate memory
     if (allocated(slater_cordon)) deallocate(slater_cordon) 
     if (allocated(gaunt))         deallocate(gaunt)
  
     return
  end subroutine atomic_make_cumatS

!>>> atomic_make_socmat3: make spin-orbit coupling matrix for 3 bands
  subroutine atomic_make_socmat3(socmat)
     use constants, only : dp

     implicit none
  
! external variables
     complex(dp), intent(out) :: socmat(6,6)
  
! local variables
     real(dp) :: sqrt2
  
     sqrt2 = sqrt(2.0_dp)
     
! make SOC on complex orbital basis, the orbital order is:
! |-1,up>, |-1,dn>, |0,up>, |0,dn>, |1,up>, |1,dn>      
     socmat = dcmplx(0.0_dp, 0.0_dp)
  
     socmat(1,1) = -1.0_dp
     socmat(4,1) = sqrt2 
     socmat(2,2) =  1.0_dp
     socmat(6,3) = sqrt2
     socmat(1,4) = sqrt2
     socmat(5,5) = 1.0_dp
     socmat(3,6) = sqrt2
     socmat(6,6) = -1.0_dp
  
     return
  end subroutine atomic_make_socmat3

!!>>> atomic_make_socmat5: make spin-orbit coupling matrix for 5 bands
  subroutine atomic_make_socmat5(socmat)
     use constants, only : dp

     implicit none
  
! external variables
     complex(dp), intent(out) :: socmat(10,10)
  
! local variables
     real(dp) :: sqrt6
  
     sqrt6 = sqrt(6.0_dp)
! make SOC on complex orbital basis, the orbital order is:
! |-2,up>, |-2,dn>, |-1,up>, |-1,dn>, |0,up>, |0,dn>, |1,up>, |1,dn>, |2,up>, |2,dn>      
  
     socmat = dcmplx(0.0_dp, 0.0_dp)
  
     socmat(1,1) = -2.0_dp 
     socmat(4,1) =  2.0_dp
     socmat(2,2) =  2.0_dp
     socmat(3,3) = -1.0_dp
     socmat(6,3) =  sqrt6
     socmat(1,4) =  2.0_dp
     socmat(4,4) =  1.0_dp
     socmat(8,5) =  sqrt6
     socmat(3,6) =  sqrt6
     socmat(7,7) =  1.0_dp
     socmat(10,7)=  2.0_dp
     socmat(5,8) =  sqrt6  
     socmat(8,8) = -1.0_dp
     socmat(9,9) =  2.0_dp
     socmat(7,10)=  2.0_dp 
     socmat(10,10)= -2.0_dp
  
     return 
  end subroutine atomic_make_socmat5

!!>>> atomic_make_socmat7: make spin-orbit coupling matrix for 7 bands
  subroutine atomic_make_socmat7(socmat)
     use constants, only : dp

     implicit none
  
! local variables
     real(dp) :: sqrt6
     real(dp) :: sqrt10
     real(dp) :: sqrt12
  
! external variables
     complex(dp), intent(out) :: socmat(14,14)    
  
     sqrt6  = sqrt(6.0_dp)
     sqrt10 = sqrt(10.0_dp)
     sqrt12 = sqrt(12.0_dp)
  
     socmat = dcmplx(0.0_dp, 0.0_dp)
  
     socmat(1,1) = -3.0_dp
     socmat(4,1) = sqrt6
     socmat(2,2) = 3.0_dp
     socmat(3,3) = -2.0_dp
     socmat(6,3) = sqrt10 
     socmat(1,4) = sqrt6
     socmat(4,4) = 2.0_dp
     socmat(5,5) = -1.0_dp
     socmat(8,5) = sqrt12
     socmat(3,6) = sqrt10 
     socmat(6,6) = 1.0_dp
     socmat(10,7) = sqrt12
     socmat(5,8) = sqrt12
     socmat(9,9) = 1.0_dp
     socmat(12,9) = sqrt10
     socmat(7,10) = sqrt12
     socmat(10,10) = -1.0_dp
     socmat(11,11) =  2.0_dp
     socmat(14,11) =  sqrt6
     socmat(9,12) = sqrt10
     socmat(12,12) =  -2.0_dp
     socmat(13,13) =  3.0_dp
     socmat(11,14) =  sqrt6
     socmat(14,14) =  -3.0_dp
  
     return
  end subroutine atomic_make_socmat7

!!>>> atomic_make_tmat_c2r: make transformation matrix from 
!!>>> complex orbital basis to real orbital 
  subroutine atomic_make_tmat_c2r( umat_c2r )
     use constants, only : dp, czero, cone, czi
     use control, only : nband, norbs
  
     implicit none
  
! external variables
! the transformation matrix from real orbitals to complex orbitals
     complex(dp), intent(out) :: umat_c2r( norbs, norbs )
  
! sqrt(2)
     real(dp) :: sqrt2
  
     sqrt2 = sqrt(2.0_dp)
     umat_c2r = czero
  
     if ( nband == 3 ) then
! the real orbital order (t2g) is:  
! dxzup, dxzdn, dyzup, dyzdn, dxyup, dxydn
! the corresponding p orbital order is:
! pyup, pydn, pxup, pxdn, pzup, pzdn
! the complex orbital |Lz,Sz> order is
! -1up, -1dn, 0up, 0dn, 1up, 1dn 
         umat_c2r(1,1) =  czi/sqrt2
         umat_c2r(5,1) =  czi/sqrt2
         umat_c2r(2,2) =  czi/sqrt2
         umat_c2r(6,2) =  czi/sqrt2
         umat_c2r(1,3) =  cone/sqrt2
         umat_c2r(5,3) =  -cone/sqrt2
         umat_c2r(2,4) =  cone/sqrt2
         umat_c2r(6,4) =  -cone/sqrt2
         umat_c2r(3,5) =  cone
         umat_c2r(4,6) =  cone
  
     elseif ( nband == 5 ) then
! the real orbital order is: 
! dz2up, dz2dn, dxzup, dxzdn, dyzup, dyzdn, dx2-y2up, dx2-y2dn, dxyup, dxydn 
! the complex orbital |Lz,Sz> order is:
! -2up, -2dn, -1up, -1dn, 0up, 0dn, 1up, 1dn, 2up, 2dn
         umat_c2r(5,1) = cone
         umat_c2r(6,2) = cone
         umat_c2r(3,3) =  cone/sqrt2 
         umat_c2r(7,3) =  -cone/sqrt2
         umat_c2r(4,4) =  cone/sqrt2 
         umat_c2r(8,4) =  -cone/sqrt2
         umat_c2r(3,5) =    czi/sqrt2
         umat_c2r(7,5) =    czi/sqrt2
         umat_c2r(4,6) =    czi/sqrt2
         umat_c2r(8,6) =    czi/sqrt2
         umat_c2r(1,7) =  cone/sqrt2
         umat_c2r(9,7) =  cone/sqrt2
         umat_c2r(2,8) =  cone/sqrt2
         umat_c2r(10,8) =  cone/sqrt2
         umat_c2r(1,9) =    czi/sqrt2
         umat_c2r(9,9) =   -czi/sqrt2
         umat_c2r(2,10) =    czi/sqrt2
         umat_c2r(10,10) =  -czi/sqrt2
     elseif ( nband == 7 ) then
! the real orbital order is:
! fz3up, fz3dn, fxz2up, fxz2dn, fyz2up, fyz2dn, fz(x2-y2)up, fz(x2-y2)dn, fxyzup, fxyzdn,
! fx(x2-3y2)up, fx(x2-3y2)dn, fy(3x2-y2)up, fy(3x2-y2)dn
! the complex orbital order is:
! -3up, -3dn, -2up, -2dn, -1up, -1dn, 0up, 0dn, 1up, 1dn, 2up, 2dn, 3up, 3dn    
         umat_c2r( 7, 1) = cone 
         umat_c2r( 8, 2) = cone 
         umat_c2r( 5, 3) = cone/sqrt2
         umat_c2r( 9, 3) = -cone/sqrt2
         umat_c2r( 6, 4) = cone/sqrt2
         umat_c2r(10, 4) = -cone/sqrt2
         umat_c2r( 5, 5) = czi/sqrt2
         umat_c2r( 9, 5) = czi/sqrt2
         umat_c2r( 6, 6) = czi/sqrt2
         umat_c2r(10, 6) = czi/sqrt2
         umat_c2r( 3, 7) = cone/sqrt2
         umat_c2r(11, 7) = cone/sqrt2
         umat_c2r( 4, 8) = cone/sqrt2
         umat_c2r(12, 8) = cone/sqrt2
         umat_c2r( 3, 9) = czi/sqrt2
         umat_c2r(11, 9) = -czi/sqrt2
         umat_c2r( 4,10) = czi/sqrt2
         umat_c2r(12,10) = -czi/sqrt2
         umat_c2r( 1,11) = cone/sqrt2
         umat_c2r(13,11) = -cone/sqrt2
         umat_c2r( 2,12) = cone/sqrt2
         umat_c2r(14,12) = -cone/sqrt2
         umat_c2r( 1,13) = czi/sqrt2
         umat_c2r(13,13) = czi/sqrt2
         umat_c2r( 2,14) = czi/sqrt2
         umat_c2r(14,14) = czi/sqrt2
     else
         call s_print_error('atomic_make_umat_c2r', &
                    'not implemented for this nband!')
     endif
  
     return
  end subroutine atomic_make_tmat_c2r

!!>>> atomic_make_tmat_r2c: make umat from real orbital 
!!>>> basis to complex orbital basis
  subroutine atomic_make_tmat_r2c(umat_r2c)
     use constants, only : dp, czero
     use control, only : norbs
  
! external variables
     complex(dp), intent(out) :: umat_r2c(norbs, norbs)
   
! local variables
     complex(dp) :: umat_c2r(norbs, norbs)
  
     umat_c2r = czero
     call atomic_make_tmat_c2r(umat_c2r)
  
     umat_r2c = transpose(dconjg(umat_c2r))
  
     return
  end subroutine atomic_make_tmat_r2c

!!>>> atomic_make_tmat_c2j: make CG coefficients
  subroutine atomic_make_tmat_c2j( umat_c2j ) 
     use constants, only : dp, czero
     use control, only : nband, norbs
  
     implicit none
     
! the transformation matrix from complex orbitals |lz,sz> to |j2,jz>
     complex(dp), intent(out) :: umat_c2j( norbs, norbs )
     
     umat_c2j = czero
  
     if (nband == 3) then
! the |lz,sz> order is:
! |-1,up>, |-1,dn>, |0,up>, |0,dn>, |1,up>, |1,dn>
! the |j2,jz> order is:
! |1/2,-1/2>, |1/2,1/2>, |3/2,-3/2>, |3/2, -1/2>, |3/2, 1/2>, |3/2,3/2>
         umat_c2j(1,1) = -sqrt(2.0_dp/3.0_dp) 
         umat_c2j(4,1) =  sqrt(1.0_dp/3.0_dp) 
         umat_c2j(3,2) = -sqrt(1.0_dp/3.0_dp) 
         umat_c2j(6,2) =  sqrt(2.0_dp/3.0_dp) 
         umat_c2j(2,3) =  1.0_dp
         umat_c2j(1,4) =  sqrt(1.0_dp/3.0_dp) 
         umat_c2j(4,4) =  sqrt(2.0_dp/3.0_dp) 
         umat_c2j(3,5) =  sqrt(2.0_dp/3.0_dp) 
         umat_c2j(6,5) =  sqrt(1.0_dp/3.0_dp) 
         umat_c2j(5,6) =  1.0_dp
     elseif ( nband == 5 ) then
! the |lz,sz> order is:
! |-2,up>, |-2,dn>, |-1,up>, |-1,dn>, |0,up>, |0,dn>, |1,up>, |1,dn>, |2,up>, |2,dn>
! the |j2,jz> order is:
! |3/2,-3/2>, |3/2,-1/2>, |3/2,1/2>, |3/2,3/2>
! |5/2,-5/2>, |5/2,-3/2>, |5/2,-1/2>, |5/2,1/2>, |5/2,3/2>, |5/2,5/2>
         umat_c2j(1,1) = -sqrt(4.0_dp/5.0_dp) 
         umat_c2j(4,1) =  sqrt(1.0_dp/5.0_dp) 
         umat_c2j(3,2) = -sqrt(3.0_dp/5.0_dp) 
         umat_c2j(6,2) =  sqrt(2.0_dp/5.0_dp) 
         umat_c2j(5,3) = -sqrt(2.0_dp/5.0_dp) 
         umat_c2j(8,3) =  sqrt(3.0_dp/5.0_dp) 
         umat_c2j(7,4) = -sqrt(1.0_dp/5.0_dp) 
         umat_c2j(10,4)=  sqrt(4.0_dp/5.0_dp) 
         umat_c2j(2,5) = 1.0_dp 
         umat_c2j(1,6) =  sqrt(1.0_dp/5.0_dp)
         umat_c2j(4,6) =  sqrt(4.0_dp/5.0_dp)
         umat_c2j(3,7) =  sqrt(2.0_dp/5.0_dp)
         umat_c2j(6,7) =  sqrt(3.0_dp/5.0_dp)
         umat_c2j(5,8) =  sqrt(3.0_dp/5.0_dp)
         umat_c2j(8,8) =  sqrt(2.0_dp/5.0_dp)
         umat_c2j(7,9) =  sqrt(4.0_dp/5.0_dp)
         umat_c2j(10,9)=  sqrt(1.0_dp/5.0_dp)
         umat_c2j(9,10)= 1.0_dp
     elseif ( nband == 7 ) then
! the |lz,sz> order is:
! |-3,up>, |-3,dn>, |-2,up>, |-2,dn>, |-1,up>, |-1,dn>, |0,up>, 
! | 0,dn>, | 1,up>, | 1,dn>, | 2,up>, | 2,dn>, | 3,up>, |3,dn>
! the |j2,jz> order is:
! |5/2,-5/2>, |5/2,-3/2>, |5/2,-1/2>, |5/2,1/2>, |5/2,3/2>, |5/2,5/2>
! |7/2,-7/2>, |7/2,-5/2>, |7/2,-3/2>, |7/2,-1/2>, |7/2,1/2>, |7/2,3/2>, 
! |7/2,5/2>, |7/2, 7/2>
         umat_c2j(1, 1) = -sqrt(6.0_dp/7.0_dp)
         umat_c2j(4, 1) =  sqrt(1.0_dp/7.0_dp)
         umat_c2j(3, 2) = -sqrt(5.0_dp/7.0_dp)
         umat_c2j(6, 2) =  sqrt(2.0_dp/7.0_dp)
         umat_c2j(5, 3) = -sqrt(4.0_dp/7.0_dp)
         umat_c2j(8, 3) =  sqrt(3.0_dp/7.0_dp)
         umat_c2j(7, 4) = -sqrt(3.0_dp/7.0_dp)
         umat_c2j(10,4) =  sqrt(4.0_dp/7.0_dp)
         umat_c2j(9, 5) = -sqrt(2.0_dp/7.0_dp)
         umat_c2j(12,5) =  sqrt(5.0_dp/7.0_dp)
         umat_c2j(11,6) = -sqrt(1.0_dp/7.0_dp)
         umat_c2j(14,6) =  sqrt(6.0_dp/7.0_dp)
  
         umat_c2j(2, 7)  = 1.0_dp
         umat_c2j(1, 8) =  sqrt(1.0_dp/7.0_dp)
         umat_c2j(4, 8) =  sqrt(6.0_dp/7.0_dp)
         umat_c2j(3, 9) =  sqrt(2.0_dp/7.0_dp)
         umat_c2j(6, 9) =  sqrt(5.0_dp/7.0_dp)
         umat_c2j(5,10) =  sqrt(3.0_dp/7.0_dp)
         umat_c2j(8,10) =  sqrt(4.0_dp/7.0_dp)
         umat_c2j(7,11) =  sqrt(4.0_dp/7.0_dp)
         umat_c2j(10,11)=  sqrt(3.0_dp/7.0_dp)
         umat_c2j(9,12) =  sqrt(5.0_dp/7.0_dp)
         umat_c2j(12,12)=  sqrt(2.0_dp/7.0_dp)
         umat_c2j(11,13)=  sqrt(6.0_dp/7.0_dp)
         umat_c2j(14,13)=  sqrt(1.0_dp/7.0_dp)
         umat_c2j(13,14)=  1.0_dp
     else
         call s_print_error('atomic_make_umat_c2j','not implemented !')
     endif
  
     return
  end subroutine atomic_make_tmat_c2j

!!>>> atomic_tran_cumat: transform Coulomb interaction U tensor 
!!>>> from one representation to another representation
  subroutine atomic_tran_cumat(amtrx, cumat, cumat_t)
     use constants, only : dp, czero, epst
     use control, only : norbs
  
     implicit none
  
! transformation matrix from orginal basis to natural basis
     complex(dp), intent(in) :: amtrx(norbs, norbs)

! coefficents matrix for generalized interaction U in orginal basis
     complex(dp), intent(in) :: cumat(norbs, norbs, norbs, norbs)

! coefficents matrix for generalized interaction U in natural basis
     complex(dp), intent(out) :: cumat_t(norbs, norbs, norbs, norbs)
  
! local varoables
! loop index over orbits in orginal single particle basis
     integer :: alpha1, alpha2
     integer :: alpha3, alpha4
     integer :: sigma1, sigma2
     integer :: sigma3, sigma4

! auxiliary complex(dp) variables
     complex(dp) :: ctmp
  
! initialize cumat_t to be zero
     cumat_t = czero 
  
     sigma1loop: do sigma1=1,norbs
     sigma2loop: do sigma2=1,norbs
     sigma3loop: do sigma3=1,norbs
     sigma4loop: do sigma4=1,norbs
         ctmp = czero
  
         alpha1loop: do alpha1=1,norbs
         alpha2loop: do alpha2=1,norbs
         alpha3loop: do alpha3=1,norbs
         alpha4loop: do alpha4=1,norbs
             if (abs(cumat(alpha1, alpha2, alpha3, alpha4)) .lt. epst) cycle
             ctmp = ctmp + cumat(alpha1, alpha2, alpha3, alpha4)          &
                  * conjg(amtrx(alpha1, sigma1)) * amtrx(alpha3, sigma3)  &
                  * conjg(amtrx(alpha2, sigma2)) * amtrx(alpha4, sigma4)
         enddo alpha4loop ! over alpha4={1,norbs} loop
         enddo alpha3loop ! over alpha3={1,norbs} loop
         enddo alpha2loop ! over alpha2={1,norbs} loop
         enddo alpha1loop ! over alpha1={1,norbs} loop
  
         cumat_t(sigma1, sigma2, sigma3, sigma4) = ctmp
     enddo sigma4loop ! over sigma4={1,norbs} loop
     enddo sigma3loop ! over sigma3={1,norbs} loop
     enddo sigma2loop ! over sigma2={1,norbs} loop
     enddo sigma1loop ! over sigma1={1,norbs} loop
  
     return
  end subroutine atomic_tran_cumat

!!>>> atomic_tran_repr_cmpl: transformation from one representation 
!!>>> to another representation, complex version
  subroutine atomic_tran_repr_cmpl( ndim, amat, umat )
     use constants, only : dp, cone, czero
  
     implicit none
  
! external variables
     integer, intent(in) :: ndim
     complex(dp), intent(inout) :: amat(ndim, ndim)
     complex(dp), intent(in) :: umat(ndim, ndim)
  
! local variables
     complex(dp) :: tmp_mat(ndim, ndim)
     complex(dp) :: alpha
     complex(dp) :: betta
  
     alpha = cone; betta = czero
     call zgemm('N', 'N', ndim, ndim, ndim, & 
                         alpha, amat, ndim, &
                                umat, ndim, &
                      betta, tmp_mat, ndim  )
  
     alpha = cone; betta = czero
     call zgemm('C', 'N', ndim, ndim, ndim, &
                         alpha, umat, ndim, &
                             tmp_mat, ndim, &
                         betta, amat, ndim  )
  
  
     return
  end subroutine atomic_tran_repr_cmpl

!>>> atomic_tran_repr_real: transformation from one representation to 
!!>>> another representation, real version
  subroutine atomic_tran_repr_real( ndim, amat, umat )
     use constants, only : dp, zero, one
  
     implicit none
  
! external variables
     integer, intent(in) :: ndim
     real(dp), intent(inout) :: amat(ndim, ndim)
     real(dp), intent(in) :: umat(ndim, ndim)
  
! local variables
     real(dp) :: tmp_mat(ndim, ndim)
     real(dp) :: alpha
     real(dp) :: betta
  
     alpha = one; betta = zero
     call dgemm('N', 'N', ndim, ndim, ndim, &
                         alpha, amat, ndim, &
                                umat, ndim, &
                      betta, tmp_mat, ndim  )
  
     alpha = one; betta = zero 
     call dgemm('T', 'N', ndim, ndim, ndim, &
                         alpha, umat, ndim, &
                             tmp_mat, ndim, &
                         betta, amat, ndim  )
  
     return
  end subroutine atomic_tran_repr_real
