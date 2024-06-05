!!!-----------------------------------------------------------------------
!!! project : iqist @ jasmine
!!! program : atomic_make_cdagger
!!!           atomic_make_c
!!!           atomic_make_gsz
!!!           atomic_make_gjz
!!!           atomic_make_gps
!!!           atomic_make_slater3
!!!           atomic_make_slater5
!!!           atomic_make_slater7
!!!           atomic_make_gaunt3
!!!           atomic_make_gaunt5
!!!           atomic_make_gaunt7
!!!           atomic_make_hund
!!!           atomic_make_umatK
!!!           atomic_make_umatS
!!!           atomic_make_smat3
!!!           atomic_make_smat5
!!!           atomic_make_smat7
!!!           atomic_make_tmat_c2r
!!!           atomic_make_tmat_r2c
!!!           atomic_make_tmat_c2j
!!!           atomic_tran_fmat
!!!           atomic_tran_umat
!!!           atomic_tran_repr_cmpl
!!!           atomic_tran_repr_real
!!!           atomic_natural_basis1
!!!           atomic_natural_basis2
!!!           atomic_natural_basis3
!!!           atomic_natural_basis4
!!! source  : atomic_util.f90
!!! type    : subroutines
!!! author  : yilin wang (email:qhwyl2006@126.com)
!!! history : 07/09/2014 by yilin wang (created)
!!!           06/04/2024 by li huang (last modified)
!!! purpose : provide the utility subroutines for the atomic eigenvalue
!!!           problem solver, such as the Dirac algebra, calculations of
!!!           gaunt coefficients, spin-orbit coupling matrix, Coulomb
!!!           interaction tensor, etc.
!!! status  : unstable
!!! comment :
!!!-----------------------------------------------------------------------

!!========================================================================
!!>>> simulate creation and destroy operators                          <<<
!!========================================================================

!!
!! @sub atomic_make_cdagger
!!
!! simulate a creation operator. create one electron on ipos of |jold>
!! to obtain new Fock state |jnew>
!!
  subroutine atomic_make_cdagger(ipos, jold, jnew, isgn)
     implicit none

!! external arguments
     ! position number (serial number of orbital)
     integer, intent(in) :: ipos

     ! old Fock state
     integer, intent(in ):: jold

     ! new Fock state
     integer, intent(out):: jnew

     ! sign due to anti-commutation relation between fermions
     integer, intent(out):: isgn

!! local variables
     ! loop index over orbitals
     integer :: iorb

!! [body

     ! it is already occupied at ipos
     ! we can not violate the Pauli principle
     if ( btest(jold, ipos-1) .eqv. .true. ) then
         call s_print_error('atomic_make_cdagger', &
             & 'severe error happened')
     endif ! back if ( btest(jold, ipos-1) .eqv. .true. ) block

     ! evaluate the sign
     isgn = 0
     !
     do iorb=1,ipos-1
        if ( btest(jold, iorb-1) .eqv. .true. ) isgn = isgn + 1
     enddo ! over iorb={1,ipos-1} loop
     !
     isgn = mod(isgn,2)
     isgn = (-1)**isgn

     ! get the final Fock state
     jnew = jold + 2**(ipos-1)

!! body]

     return
  end subroutine atomic_make_cdagger

!!
!! @sub atomic_make_c
!!
!! simulate an annihilation operator. destroy one electron on ipos of
!! |jold> to obtain new Fock state |jnew>
!!
  subroutine atomic_make_c(ipos, jold, jnew, isgn)
     implicit none

!! external arguments
     ! position number (serial number of orbital)
     integer, intent(in)  :: ipos

     ! old Fock state
     integer, intent(in ) :: jold

     ! new Fock state
     integer, intent(out) :: jnew

     ! sign due to anti-commutation relation between fermions
     integer, intent(out) :: isgn

!! local variables
     ! loop index
     integer :: iorb

!! [body

     ! it is already unoccupied at ipos
     ! we can not violate the Pauli principle
     if ( btest(jold, ipos-1) .eqv. .false. ) then
         call s_print_error('atomic_make_c', &
             & 'severe error happened')
     endif ! back if ( btest(jold, ipos-1) .eqv. .false. ) block

     ! evaluate the sign
     isgn = 0
     !
     do iorb=1,ipos-1
         if ( btest(jold, iorb-1) .eqv. .true. ) isgn = isgn + 1
     enddo ! back iorb={1,ipos-1} loop
     !
     isgn = mod(isgn,2)
     isgn = (-1)**isgn

     ! get the final Fock state
     jnew = jold - 2**(ipos-1)

!! body]

     return
  end subroutine atomic_make_c

!!========================================================================
!!>>> calculate good quantum numbers                                   <<<
!!========================================================================

!!
!! @sub atomic_make_gsz
!!
!! calculate Sz quantum number for each orbital
!!
  subroutine atomic_make_gsz(good_sz)
     use control, only : nband, norbs

     implicit none

!! external arguments
     ! good quantum number: Sz
     integer, intent(out) :: good_sz(norbs)

!! local variables
     ! loop index over orbitals
     integer :: i

!! [body

     ! the orbital order is up up up ... dn dn dn ...
     do i=1,norbs
         if ( i <= nband ) then
             good_sz(i) = +1
         else
             good_sz(i) = -1
         endif ! back if ( i <= nband ) block
     enddo ! over i={1,norbs} loop

!! body]

     return
  end subroutine atomic_make_gsz

!!
!! @sub atomic_make_gjz
!!
!! calculate Jz quantum number for each orbital
!!
  subroutine atomic_make_gjz(good_jz)
     use control, only : nband, norbs

     implicit none

!! external arguments
     ! good quantum number: Jz
     integer, intent(out) :: good_jz(norbs)

!! [body

     select case (nband)

         case (3)
             good_jz( 1) = -1 ! j = 1/2
             good_jz( 2) =  1
             good_jz( 3) = -3 ! j = 3/2
             good_jz( 4) = -1
             good_jz( 5) =  1
             good_jz( 6) =  3

         case (5)
             good_jz( 1) = -3 ! j = 3/2
             good_jz( 2) = -1
             good_jz( 3) =  1
             good_jz( 4) =  3
             good_jz( 5) = -5 ! j = 5/2
             good_jz( 6) = -3
             good_jz( 7) = -1
             good_jz( 8) =  1
             good_jz( 9) =  3
             good_jz(10) =  5

         case (7)
             good_jz( 1) = -5 ! j = 5/2
             good_jz( 2) = -3
             good_jz( 3) = -1
             good_jz( 4) =  1
             good_jz( 5) =  3
             good_jz( 6) =  5
             good_jz( 7) = -7 ! j = 7/2
             good_jz( 8) = -5
             good_jz( 9) = -3
             good_jz(10) = -1
             good_jz(11) =  1
             good_jz(12) =  3
             good_jz(13) =  5
             good_jz(14) =  7

         case default
             call s_print_error('atomic_make_gjz', &
                 & 'not implemented for this norbs value!')

     end select

!! body]

     return
  end subroutine atomic_make_gjz

!!
!! @sub atomic_make_gps
!!
!! calculate PS quantum number for each band. note that only the pure
!! band-dependent part is calculated in this subroutine
!!
  subroutine atomic_make_gps(good_ps)
     use control, only : nband

     implicit none

!! external arguments
     ! good quantum number: PS
     integer, intent(out) :: good_ps(nband)

!! local variables
     ! loop index
     integer :: i

!! [body

     do i=1,nband
         good_ps(i) = 2**i
     enddo ! over i={1,nband} loop

!! body]

     return
  end subroutine atomic_make_gps

!!========================================================================
!!>>> determine Slater integrals                                       <<<
!!========================================================================

!!
!! @sub atomic_make_slater3
!!
!! build Slater integrals (radial integrals) F_k for 3 band case (l = 1)
!!
  subroutine atomic_make_slater3(Fk)
     use constants, only : dp
     use constants, only : zero

     use control, only : Ud, Jh

     implicit none

!! external arguments
     ! Slater integrals F_k for l = 1
     real(dp), intent(out) :: Fk(0:2)
     
!! [body

     Fk = zero
     !
     Fk(0) = Ud
     Fk(2) = Jh * 5.0_dp

!! body]

     return
  end subroutine atomic_make_slater3

!!
!! @sub atomic_make_slater5
!!
!! build Slater integrals (radial integrals) F_k for 5 band case (l = 2)
!!
  subroutine atomic_make_slater5(Fk)
     use constants, only : dp
     use constants, only : zero

     use control, only : Ud, Jh

     implicit none

!! external arguments
     ! Slater integrals F_k for l = 2
     real(dp), intent(out) :: Fk(0:4)
     
!! [body

     Fk = zero
     !
     Fk(0) = Ud
     Fk(2) = Jh * 14.0_dp / (1.0_dp + 0.625_dp)
     Fk(4) = 0.625_dp * Fk(2)

     ! the triqs code uses 0.63, insead of 0.625

!! body]

     return
  end subroutine atomic_make_slater5

!!
!! @sub atomic_make_slater7
!!
!! build Slater integrals (radial integrals) F_k for 7 band case (l = 3)
!!
  subroutine atomic_make_slater7(Fk)
     use constants, only : dp
     use constants, only : zero

     use control, only : Ud, Jh

     implicit none

!! external arguments
     ! Slater integrals F_k for l = 3
     real(dp), intent(out) :: Fk(0:6)

!! [body

     Fk = zero
     !
     Fk(0) = Ud
     Fk(2) = 286.0_dp
     Fk(2) = Fk(2) + 195.0_dp * 451.0_dp / 675.0_dp
     Fk(2) = Fk(2) + 250.0_dp * 1001.0_dp / 2025.0_dp
     Fk(2) = Jh * 6435.0_dp / Fk(2)
     Fk(4) = 451.0_dp / 675.0_dp * Fk(2)
     Fk(6) = 1001.0_dp / 2025.0_dp * Fk(2)

!! body]

     return
  end subroutine atomic_make_slater7

!!========================================================================
!!>>> determine Gaunt coefficients                                     <<<
!!========================================================================

!!
!! @sub atomic_make_gaunt3
!!
!! build c^{k}_{l}(m_1,m_2) coefficients for 3 band case (l = 1)
!!
  subroutine atomic_make_gaunt3(gaunt)
     use constants, only : dp
     use constants, only : zero

     implicit none

!! external arguments
     ! gaunt coefficients, gaunt(m_1, m_2, k)
     real(dp), intent(out) :: gaunt(-1:1,-1:1,0:2)

!! [body

     gaunt = zero
     !
     ! for k = 0
     gaunt( -1 , -1 ,  0 ) =  -1.0_dp * ( -1.0 )
     gaunt(  0 ,  0 ,  0 ) =   1.0_dp * (  1.0 )
     gaunt(  1 ,  1 ,  0 ) =  -1.0_dp * ( -1.0 )
     !
     ! for k = 2
     gaunt( -1 , -1 ,  2 ) =   sqrt(1.0_dp) / 5.0_dp * ( -1.0 )
     gaunt( -1 ,  0 ,  2 ) =  -sqrt(3.0_dp) / 5.0_dp * ( -1.0 )
     gaunt( -1 ,  1 ,  2 ) =   sqrt(6.0_dp) / 5.0_dp * ( -1.0 )
     gaunt(  0 , -1 ,  2 ) =  -sqrt(3.0_dp) / 5.0_dp * (  1.0 )
     gaunt(  0 ,  0 ,  2 ) =   sqrt(4.0_dp) / 5.0_dp * (  1.0 )
     gaunt(  0 ,  1 ,  2 ) =  -sqrt(3.0_dp) / 5.0_dp * (  1.0 )
     gaunt(  1 , -1 ,  2 ) =   sqrt(6.0_dp) / 5.0_dp * ( -1.0 )
     gaunt(  1 ,  0 ,  2 ) =  -sqrt(3.0_dp) / 5.0_dp * ( -1.0 )
     gaunt(  1 ,  1 ,  2 ) =   sqrt(1.0_dp) / 5.0_dp * ( -1.0 )

!! body]

     return
  end subroutine atomic_make_gaunt3

!!
!! @sub atomic_make_gaunt5
!!
!! build c^{k}_{l}(m_1,m_2) coefficients for 5 band case (l = 2)
!!
  subroutine atomic_make_gaunt5(gaunt)
     use constants, only : dp
     use constants, only : zero

     implicit none

!! external arguments
     ! gaunt coefficients, gaunt(m_1, m_2, k)
     real(dp), intent(out) :: gaunt(-2:2,-2:2,0:4)

!! [body

     gaunt = zero
     !
     ! for k = 0
     gaunt( -2 , -2 ,  0 ) =  1.0_dp * (  1.0 )
     gaunt( -1 , -1 ,  0 ) = -1.0_dp * ( -1.0 )
     gaunt(  0 ,  0 ,  0 ) =  1.0_dp * (  1.0 )
     gaunt(  1 ,  1 ,  0 ) = -1.0_dp * ( -1.0 )
     gaunt(  2 ,  2 ,  0 ) =  1.0_dp * (  1.0 )
     !
     ! for k = 2
     gaunt( -2 , -2 ,  2 ) = -sqrt(4.0_dp) / 7.0_dp * (  1.0 )
     gaunt( -2 , -1 ,  2 ) =  sqrt(6.0_dp) / 7.0_dp * (  1.0 )
     gaunt( -2 ,  0 ,  2 ) = -sqrt(4.0_dp) / 7.0_dp * (  1.0 )
     gaunt( -1 , -2 ,  2 ) =  sqrt(6.0_dp) / 7.0_dp * ( -1.0 )
     gaunt( -1 , -1 ,  2 ) = -sqrt(1.0_dp) / 7.0_dp * ( -1.0 )
     gaunt( -1 ,  0 ,  2 ) = -sqrt(1.0_dp) / 7.0_dp * ( -1.0 )
     gaunt( -1 ,  1 ,  2 ) =  sqrt(6.0_dp) / 7.0_dp * ( -1.0 )
     gaunt(  0 , -2 ,  2 ) = -sqrt(4.0_dp) / 7.0_dp * (  1.0 )
     gaunt(  0 , -1 ,  2 ) = -sqrt(1.0_dp) / 7.0_dp * (  1.0 )
     gaunt(  0 ,  0 ,  2 ) =  sqrt(4.0_dp) / 7.0_dp * (  1.0 )
     gaunt(  0 ,  1 ,  2 ) = -sqrt(1.0_dp) / 7.0_dp * (  1.0 )
     gaunt(  0 ,  2 ,  2 ) = -sqrt(4.0_dp) / 7.0_dp * (  1.0 )
     gaunt(  1 , -1 ,  2 ) =  sqrt(6.0_dp) / 7.0_dp * ( -1.0 )
     gaunt(  1 ,  0 ,  2 ) = -sqrt(1.0_dp) / 7.0_dp * ( -1.0 )
     gaunt(  1 ,  1 ,  2 ) = -sqrt(1.0_dp) / 7.0_dp * ( -1.0 )
     gaunt(  1 ,  2 ,  2 ) =  sqrt(6.0_dp) / 7.0_dp * ( -1.0 )
     gaunt(  2 ,  0 ,  2 ) = -sqrt(4.0_dp) / 7.0_dp * (  1.0 )
     gaunt(  2 ,  1 ,  2 ) =  sqrt(6.0_dp) / 7.0_dp * (  1.0 )
     gaunt(  2 ,  2 ,  2 ) = -sqrt(4.0_dp) / 7.0_dp * (  1.0 )
     !
     ! for k = 4
     gaunt( -2 , -2 ,  4 ) =  sqrt( 1.0_dp) / 21.0_dp * (  1.0 )
     gaunt( -2 , -1 ,  4 ) = -sqrt( 5.0_dp) / 21.0_dp * (  1.0 )
     gaunt( -2 ,  0 ,  4 ) =  sqrt(15.0_dp) / 21.0_dp * (  1.0 )
     gaunt( -2 ,  1 ,  4 ) = -sqrt(35.0_dp) / 21.0_dp * (  1.0 )
     gaunt( -2 ,  2 ,  4 ) =  sqrt(70.0_dp) / 21.0_dp * (  1.0 )
     gaunt( -1 , -2 ,  4 ) = -sqrt( 5.0_dp) / 21.0_dp * ( -1.0 )
     gaunt( -1 , -1 ,  4 ) =  sqrt(16.0_dp) / 21.0_dp * ( -1.0 )
     gaunt( -1 ,  0 ,  4 ) = -sqrt(30.0_dp) / 21.0_dp * ( -1.0 )
     gaunt( -1 ,  1 ,  4 ) =  sqrt(40.0_dp) / 21.0_dp * ( -1.0 )
     gaunt( -1 ,  2 ,  4 ) = -sqrt(35.0_dp) / 21.0_dp * ( -1.0 )
     gaunt(  0 , -2 ,  4 ) =  sqrt(15.0_dp) / 21.0_dp * (  1.0 )
     gaunt(  0 , -1 ,  4 ) = -sqrt(30.0_dp) / 21.0_dp * (  1.0 )
     gaunt(  0 ,  0 ,  4 ) =  sqrt(36.0_dp) / 21.0_dp * (  1.0 )
     gaunt(  0 ,  1 ,  4 ) = -sqrt(30.0_dp) / 21.0_dp * (  1.0 )
     gaunt(  0 ,  2 ,  4 ) =  sqrt(15.0_dp) / 21.0_dp * (  1.0 )
     gaunt(  1 , -2 ,  4 ) = -sqrt(35.0_dp) / 21.0_dp * ( -1.0 )
     gaunt(  1 , -1 ,  4 ) =  sqrt(40.0_dp) / 21.0_dp * ( -1.0 )
     gaunt(  1 ,  0 ,  4 ) = -sqrt(30.0_dp) / 21.0_dp * ( -1.0 )
     gaunt(  1 ,  1 ,  4 ) =  sqrt(16.0_dp) / 21.0_dp * ( -1.0 )
     gaunt(  1 ,  2 ,  4 ) = -sqrt( 5.0_dp) / 21.0_dp * ( -1.0 )
     gaunt(  2 , -2 ,  4 ) =  sqrt(70.0_dp) / 21.0_dp * (  1.0 )
     gaunt(  2 , -1 ,  4 ) = -sqrt(35.0_dp) / 21.0_dp * (  1.0 )
     gaunt(  2 ,  0 ,  4 ) =  sqrt(15.0_dp) / 21.0_dp * (  1.0 )
     gaunt(  2 ,  1 ,  4 ) = -sqrt( 5.0_dp) / 21.0_dp * (  1.0 )
     gaunt(  2 ,  2 ,  4 ) =  sqrt( 1.0_dp) / 21.0_dp * (  1.0 )

!! body]

     return
  end subroutine atomic_make_gaunt5

!!
!! @sub atomic_make_gaunt7
!!
!! build c^{k}_{l}(m_1,m_2) coefficients for 7 band case (l = 3)
!!
  subroutine atomic_make_gaunt7(gaunt)
     use constants, only : dp
     use constants, only : zero

     implicit none

!! external arguments
     ! gaunt coefficients, gaunt(m_1, m_2, k)
     real(dp), intent(out) :: gaunt(-3:3,-3:3,0:6)

!! [body

     gaunt = zero
     !
     ! for k = 0
     gaunt( -3 , -3 ,  0 ) = -1 * ( -1.0 )
     gaunt( -2 , -2 ,  0 ) =  1 * (  1.0 )
     gaunt( -1 , -1 ,  0 ) = -1 * ( -1.0 )
     gaunt(  0 ,  0 ,  0 ) =  1 * (  1.0 )
     gaunt(  1 ,  1 ,  0 ) = -1 * ( -1.0 )
     gaunt(  2 ,  2 ,  0 ) =  1 * (  1.0 )
     gaunt(  3 ,  3 ,  0 ) = -1 * ( -1.0 )
     !
     ! for k = 2
     gaunt( -3 , -3 ,  2 ) =  sqrt( 1.0_dp) /  3.0_dp * ( -1.0 )
     gaunt( -3 , -2 ,  2 ) = -sqrt( 1.0_dp) /  3.0_dp * ( -1.0 )
     gaunt( -3 , -1 ,  2 ) =  sqrt(10.0_dp) / 15.0_dp * ( -1.0 )
     gaunt( -2 , -3 ,  2 ) = -sqrt( 1.0_dp) /  3.0_dp * (  1.0 )
     gaunt( -2 , -1 ,  2 ) =  sqrt(15.0_dp) / 15.0_dp * (  1.0 )
     gaunt( -2 ,  0 ,  2 ) = -sqrt(20.0_dp) / 15.0_dp * (  1.0 )
     gaunt( -1 , -3 ,  2 ) =  sqrt(10.0_dp) / 15.0_dp * ( -1.0 )
     gaunt( -1 , -2 ,  2 ) =  sqrt(15.0_dp) / 15.0_dp * ( -1.0 )
     gaunt( -1 , -1 ,  2 ) = -sqrt( 1.0_dp) /  5.0_dp * ( -1.0 )
     gaunt( -1 ,  0 ,  2 ) = -sqrt( 2.0_dp) / 15.0_dp * ( -1.0 )
     gaunt( -1 ,  1 ,  2 ) =  sqrt(24.0_dp) / 15.0_dp * ( -1.0 )
     gaunt(  0 , -2 ,  2 ) = -sqrt(20.0_dp) / 15.0_dp * (  1.0 )
     gaunt(  0 , -1 ,  2 ) = -sqrt( 2.0_dp) / 15.0_dp * (  1.0 )
     gaunt(  0 ,  0 ,  2 ) =  sqrt(16.0_dp) / 15.0_dp * (  1.0 )
     gaunt(  0 ,  1 ,  2 ) = -sqrt( 2.0_dp) / 15.0_dp * (  1.0 )
     gaunt(  0 ,  2 ,  2 ) = -sqrt(20.0_dp) / 15.0_dp * (  1.0 )
     gaunt(  1 , -1 ,  2 ) =  sqrt(24.0_dp) / 15.0_dp * ( -1.0 )
     gaunt(  1 ,  0 ,  2 ) = -sqrt( 2.0_dp) / 15.0_dp * ( -1.0 )
     gaunt(  1 ,  1 ,  2 ) = -sqrt( 1.0_dp) /  5.0_dp * ( -1.0 )
     gaunt(  1 ,  2 ,  2 ) =  sqrt(15.0_dp) / 15.0_dp * ( -1.0 )
     gaunt(  1 ,  3 ,  2 ) =  sqrt(10.0_dp) / 15.0_dp * ( -1.0 )
     gaunt(  2 ,  0 ,  2 ) = -sqrt(20.0_dp) / 15.0_dp * (  1.0 )
     gaunt(  2 ,  1 ,  2 ) =  sqrt(15.0_dp) / 15.0_dp * (  1.0 )
     gaunt(  2 ,  3 ,  2 ) = -sqrt( 1.0_dp) /  3.0_dp * (  1.0 )
     gaunt(  3 ,  1 ,  2 ) =  sqrt(10.0_dp) / 15.0_dp * ( -1.0 )
     gaunt(  3 ,  2 ,  2 ) = -sqrt( 1.0_dp) /  3.0_dp * ( -1.0 )
     gaunt(  3 ,  3 ,  2 ) =  sqrt( 1.0_dp) /  3.0_dp * ( -1.0 )
     !
     ! for k = 4
     gaunt( -3 , -3 ,  4 ) = -sqrt( 1.0_dp) / 11.0_dp * ( -1.0 )
     gaunt( -3 , -2 ,  4 ) =  sqrt(30.0_dp) / 33.0_dp * ( -1.0 )
     gaunt( -3 , -1 ,  4 ) = -sqrt( 6.0_dp) / 11.0_dp * ( -1.0 )
     gaunt( -3 ,  0 ,  4 ) =  sqrt( 7.0_dp) / 11.0_dp * ( -1.0 )
     gaunt( -3 ,  1 ,  4 ) = -sqrt(42.0_dp) / 33.0_dp * ( -1.0 )
     gaunt( -2 , -3 ,  4 ) =  sqrt(30.0_dp) / 33.0_dp * (  1.0 )
     gaunt( -2 , -2 ,  4 ) = -sqrt(49.0_dp) / 33.0_dp * (  1.0 )
     gaunt( -2 , -1 ,  4 ) =  sqrt(32.0_dp) / 33.0_dp * (  1.0 )
     gaunt( -2 ,  0 ,  4 ) = -sqrt( 3.0_dp) / 33.0_dp * (  1.0 )
     gaunt( -2 ,  1 ,  4 ) = -sqrt(14.0_dp) / 33.0_dp * (  1.0 )
     gaunt( -2 ,  2 ,  4 ) =  sqrt(70.0_dp) / 33.0_dp * (  1.0 )
     gaunt( -1 , -3 ,  4 ) = -sqrt( 6.0_dp) / 11.0_dp * ( -1.0 )
     gaunt( -1 , -2 ,  4 ) =  sqrt(32.0_dp) / 33.0_dp * ( -1.0 )
     gaunt( -1 , -1 ,  4 ) = -sqrt( 1.0_dp) / 33.0_dp * ( -1.0 )
     gaunt( -1 ,  0 ,  4 ) = -sqrt(15.0_dp) / 33.0_dp * ( -1.0 )
     gaunt( -1 ,  1 ,  4 ) =  sqrt(40.0_dp) / 33.0_dp * ( -1.0 )
     gaunt( -1 ,  2 ,  4 ) = -sqrt(14.0_dp) / 33.0_dp * ( -1.0 )
     gaunt( -1 ,  3 ,  4 ) = -sqrt(42.0_dp) / 33.0_dp * ( -1.0 )
     gaunt(  0 , -3 ,  4 ) =  sqrt( 7.0_dp) / 11.0_dp * (  1.0 )
     gaunt(  0 , -2 ,  4 ) = -sqrt( 3.0_dp) / 33.0_dp * (  1.0 )
     gaunt(  0 , -1 ,  4 ) = -sqrt(15.0_dp) / 33.0_dp * (  1.0 )
     gaunt(  0 ,  0 ,  4 ) =  sqrt( 4.0_dp) / 11.0_dp * (  1.0 )
     gaunt(  0 ,  1 ,  4 ) = -sqrt(15.0_dp) / 33.0_dp * (  1.0 )
     gaunt(  0 ,  2 ,  4 ) = -sqrt( 3.0_dp) / 33.0_dp * (  1.0 )
     gaunt(  0 ,  3 ,  4 ) =  sqrt( 7.0_dp) / 11.0_dp * (  1.0 )
     gaunt(  1 , -3 ,  4 ) = -sqrt(42.0_dp) / 33.0_dp * ( -1.0 )
     gaunt(  1 , -2 ,  4 ) = -sqrt(14.0_dp) / 33.0_dp * ( -1.0 )
     gaunt(  1 , -1 ,  4 ) =  sqrt(40.0_dp) / 33.0_dp * ( -1.0 )
     gaunt(  1 ,  0 ,  4 ) = -sqrt(15.0_dp) / 33.0_dp * ( -1.0 )
     gaunt(  1 ,  1 ,  4 ) = -sqrt( 1.0_dp) / 33.0_dp * ( -1.0 )
     gaunt(  1 ,  2 ,  4 ) =  sqrt(32.0_dp) / 33.0_dp * ( -1.0 )
     gaunt(  1 ,  3 ,  4 ) = -sqrt( 6.0_dp) / 11.0_dp * ( -1.0 )
     gaunt(  2 , -2 ,  4 ) =  sqrt(70.0_dp) / 33.0_dp * (  1.0 )
     gaunt(  2 , -1 ,  4 ) = -sqrt(14.0_dp) / 33.0_dp * (  1.0 )
     gaunt(  2 ,  0 ,  4 ) = -sqrt( 3.0_dp) / 33.0_dp * (  1.0 )
     gaunt(  2 ,  1 ,  4 ) =  sqrt(32.0_dp) / 33.0_dp * (  1.0 )
     gaunt(  2 ,  2 ,  4 ) = -sqrt(49.0_dp) / 33.0_dp * (  1.0 )
     gaunt(  2 ,  3 ,  4 ) =  sqrt(30.0_dp) / 33.0_dp * (  1.0 )
     gaunt(  3 , -1 ,  4 ) = -sqrt(42.0_dp) / 33.0_dp * ( -1.0 )
     gaunt(  3 ,  0 ,  4 ) =  sqrt( 7.0_dp) / 11.0_dp * ( -1.0 )
     gaunt(  3 ,  1 ,  4 ) = -sqrt( 6.0_dp) / 11.0_dp * ( -1.0 )
     gaunt(  3 ,  2 ,  4 ) =  sqrt(30.0_dp) / 33.0_dp * ( -1.0 )
     gaunt(  3 ,  3 ,  4 ) = -sqrt( 1.0_dp) / 11.0_dp * ( -1.0 )
     !
     ! for k = 6
     gaunt( -3 , -3 ,  6 ) =  sqrt(   25.0_dp) / 429.0_dp * ( -1.0 )
     gaunt( -3 , -2 ,  6 ) = -sqrt(  175.0_dp) / 429.0_dp * ( -1.0 )
     gaunt( -3 , -1 ,  6 ) =  sqrt(  700.0_dp) / 429.0_dp * ( -1.0 )
     gaunt( -3 ,  0 ,  6 ) = -sqrt( 2100.0_dp) / 429.0_dp * ( -1.0 )
     gaunt( -3 ,  1 ,  6 ) =  sqrt( 5250.0_dp) / 429.0_dp * ( -1.0 )
     gaunt( -3 ,  2 ,  6 ) = -sqrt(11550.0_dp) / 429.0_dp * ( -1.0 )
     gaunt( -3 ,  3 ,  6 ) =  sqrt(23100.0_dp) / 429.0_dp * ( -1.0 )
     gaunt( -2 , -3 ,  6 ) = -sqrt(  175.0_dp) / 429.0_dp * (  1.0 )
     gaunt( -2 , -2 ,  6 ) =  sqrt(  100.0_dp) / 143.0_dp * (  1.0 )
     gaunt( -2 , -1 ,  6 ) = -sqrt( 2625.0_dp) / 429.0_dp * (  1.0 )
     gaunt( -2 ,  0 ,  6 ) =  sqrt( 5600.0_dp) / 429.0_dp * (  1.0 )
     gaunt( -2 ,  1 ,  6 ) = -sqrt( 1050.0_dp) / 143.0_dp * (  1.0 )
     gaunt( -2 ,  2 ,  6 ) =  sqrt( 1400.0_dp) / 143.0_dp * (  1.0 )
     gaunt( -2 ,  3 ,  6 ) = -sqrt(11550.0_dp) / 429.0_dp * (  1.0 )
     gaunt( -1 , -3 ,  6 ) =  sqrt(  700.0_dp) / 429.0_dp * ( -1.0 )
     gaunt( -1 , -2 ,  6 ) = -sqrt( 2625.0_dp) / 429.0_dp * ( -1.0 )
     gaunt( -1 , -1 ,  6 ) =  sqrt(  625.0_dp) / 143.0_dp * ( -1.0 )
     gaunt( -1 ,  0 ,  6 ) = -sqrt( 8750.0_dp) / 429.0_dp * ( -1.0 )
     gaunt( -1 ,  1 ,  6 ) =  sqrt(10500.0_dp) / 429.0_dp * ( -1.0 )
     gaunt( -1 ,  2 ,  6 ) = -sqrt( 1050.0_dp) / 143.0_dp * ( -1.0 )
     gaunt( -1 ,  3 ,  6 ) =  sqrt( 5250.0_dp) / 429.0_dp * ( -1.0 )
     gaunt(  0 , -3 ,  6 ) = -sqrt( 2100.0_dp) / 429.0_dp * (  1.0 )
     gaunt(  0 , -2 ,  6 ) =  sqrt( 5600.0_dp) / 429.0_dp * (  1.0 )
     gaunt(  0 , -1 ,  6 ) = -sqrt( 8750.0_dp) / 429.0_dp * (  1.0 )
     gaunt(  0 ,  0 ,  6 ) =  sqrt(10000.0_dp) / 429.0_dp * (  1.0 )
     gaunt(  0 ,  1 ,  6 ) = -sqrt( 8750.0_dp) / 429.0_dp * (  1.0 )
     gaunt(  0 ,  2 ,  6 ) =  sqrt( 5600.0_dp) / 429.0_dp * (  1.0 )
     gaunt(  0 ,  3 ,  6 ) = -sqrt( 2100.0_dp) / 429.0_dp * (  1.0 )
     gaunt(  1 , -3 ,  6 ) =  sqrt( 5250.0_dp) / 429.0_dp * ( -1.0 )
     gaunt(  1 , -2 ,  6 ) = -sqrt( 1050.0_dp) / 143.0_dp * ( -1.0 )
     gaunt(  1 , -1 ,  6 ) =  sqrt(10500.0_dp) / 429.0_dp * ( -1.0 )
     gaunt(  1 ,  0 ,  6 ) = -sqrt( 8750.0_dp) / 429.0_dp * ( -1.0 )
     gaunt(  1 ,  1 ,  6 ) =  sqrt(  625.0_dp) / 143.0_dp * ( -1.0 )
     gaunt(  1 ,  2 ,  6 ) = -sqrt( 2625.0_dp) / 429.0_dp * ( -1.0 )
     gaunt(  1 ,  3 ,  6 ) =  sqrt(  700.0_dp) / 429.0_dp * ( -1.0 )
     gaunt(  2 , -3 ,  6 ) = -sqrt(11550.0_dp) / 429.0_dp * (  1.0 )
     gaunt(  2 , -2 ,  6 ) =  sqrt( 1400.0_dp) / 143.0_dp * (  1.0 )
     gaunt(  2 , -1 ,  6 ) = -sqrt( 1050.0_dp) / 143.0_dp * (  1.0 )
     gaunt(  2 ,  0 ,  6 ) =  sqrt( 5600.0_dp) / 429.0_dp * (  1.0 )
     gaunt(  2 ,  1 ,  6 ) = -sqrt( 2625.0_dp) / 429.0_dp * (  1.0 )
     gaunt(  2 ,  2 ,  6 ) =  sqrt(  100.0_dp) / 143.0_dp * (  1.0 )
     gaunt(  2 ,  3 ,  6 ) = -sqrt(  175.0_dp) / 429.0_dp * (  1.0 )
     gaunt(  3 , -3 ,  6 ) =  sqrt(23100.0_dp) / 429.0_dp * ( -1.0 )
     gaunt(  3 , -2 ,  6 ) = -sqrt(11550.0_dp) / 429.0_dp * ( -1.0 )
     gaunt(  3 , -1 ,  6 ) =  sqrt( 5250.0_dp) / 429.0_dp * ( -1.0 )
     gaunt(  3 ,  0 ,  6 ) = -sqrt( 2100.0_dp) / 429.0_dp * ( -1.0 )
     gaunt(  3 ,  1 ,  6 ) =  sqrt(  700.0_dp) / 429.0_dp * ( -1.0 )
     gaunt(  3 ,  2 ,  6 ) = -sqrt(  175.0_dp) / 429.0_dp * ( -1.0 )
     gaunt(  3 ,  3 ,  6 ) =  sqrt(   25.0_dp) / 429.0_dp * ( -1.0 )

!! body]

     return
  end subroutine atomic_make_gaunt7

!!========================================================================
!!>>> determine Coulomb interaction matrix                             <<<
!!========================================================================

!!
!! @sub atomic_make_hund
!!
!! make the Hund's rule coupling matrix
!!
  subroutine atomic_make_hund(hund)
     use constants, only : dp
     use constants, only : zero

     use control, only : nband
     use control, only : Jz, Js, Jp

!! external arguments
     ! Hund's rule coupling matrix
     real(dp), intent(out) :: hund(nband,nband,3)

!! local variables
     ! loop index
     integer  :: i

!! [body

     hund = zero
     !
     hund(:,:,1) = Jz
     hund(:,:,2) = Js
     hund(:,:,3) = Jp
     !
     do i=1,nband
         hund(i,i,:) = zero
     enddo ! over i={1,nband} loop

!! body]

     return
  end subroutine atomic_make_hund

!!
!! @sub atomic_make_umatK
!!
!! make Coulomb interaction U (a rank-4 tensor) according to Kanamori
!! parameterized Hamiltonian
!!
  subroutine atomic_make_umatK()
     use constants, only : dp
     use constants, only : zero, two
     use constants, only : czero

     use control, only : nband, norbs
     use control, only : Uc

     use m_spmat, only : umat

     implicit none

!! local varibales
     ! orbital index
     integer  :: m, n
     integer  :: q, p

     ! band index and spin index
     integer  :: m_b, n_b
     integer  :: q_b, p_b
     integer  :: m_s, n_s
     integer  :: q_s, p_s

     ! dummy variables
     real(dp) :: dtmp

     ! Hund's rule matrix
     real(dp) :: hund(nband,nband,3)

!! [body

     ! initialize hund to zero
     hund = zero

     ! build Hund's rule coupling matrix
     call atomic_make_hund(hund)

     ! initialize umat to zero
     umat = czero

     ! loop for creation operators
     alpha: do m=1,norbs-1
         beta: do n=m+1,norbs

             ! loop for annihilation operators
             gamma: do p=1,norbs-1
                 delta: do q=p+1,norbs

                     !
                     ! notice:
                     !
                     ! here we just assume the orbital order is
                     !
                     ! up dn up dn up dn ...
                     !

                     ! get the band and spin indices
                     m_b = ( m + 1 ) / 2; m_s = mod(m,2)
                     n_b = ( n + 1 ) / 2; n_s = mod(n,2)
                     p_b = ( p + 1 ) / 2; p_s = mod(p,2)
                     q_b = ( q + 1 ) / 2; q_s = mod(q,2)

                     dtmp = zero

                     ! intraorbital Coulomb interaction
                     if ( ( m == p ) .and. ( n == q ) ) then
                         if ( ( m_b == n_b ) .and. ( m_s /= n_s ) ) then
                             dtmp = dtmp + Uc
                         endif ! back if block
                     endif ! back if block

                     ! interorbital Coulomb interaction
                     if ( ( m == p ) .and. ( n == q ) ) then
                         if ( m_b /= n_b ) then
                             dtmp = dtmp + (Uc - two * hund(m_b,n_b,1))
                         endif ! back if ( m_b /= n_b ) block
                     endif ! back if block

                     ! Hund's exchange interaction
                     if ( ( m == p ) .and. ( n == q ) ) then
                         if ( ( m_b /= n_b ) .and. ( m_s == n_s ) ) then
                             dtmp = dtmp - hund(m_b,n_b,1)
                         endif ! back if block
                     endif ! back if block

                     ! spin flip term
                     if ( ( m_b == p_b ) .and. ( n_b == q_b ) ) then
                         if ( ( m_s /= p_s ) .and. ( n_s /= q_s ) .and. ( m_s /= n_s ) ) then
                             dtmp = dtmp - hund(m_b,n_b,2)
                         endif ! back if block
                     endif ! back if block

                     ! pair hopping term
                     if ( ( m_b == n_b ) .and. ( q_b == p_b ) .and. ( m_b /= q_b ) ) then
                         if ( ( m_s /= n_s ) .and. ( q_s /= p_s ) .and. ( m_s == p_s ) ) then
                             dtmp = dtmp + hund(m_b,p_b,3)
                         endif ! back if block
                     endif ! back if block

                     !
                     ! notice:
                     !
                     ! now we change the orbital order to
                     !
                     ! up up up ... dn dn dn ...
                     !
                     umat( m_b + nband * (1 - m_s), &
                           n_b + nband * (1 - n_s), &
                           q_b + nband * (1 - q_s), &
                           p_b + nband * (1 - p_s) ) = dtmp

                 enddo delta ! over q={p+1,norbs} loop
             enddo gamma ! over p={1,norbs-1} loop
         enddo beta ! over n={m+1,norbs} loop
     enddo alpha ! over m={1,norbs-1} loop

!! body]

     return
  end subroutine atomic_make_umatK

!!
!! @sub atomic_make_umatS
!!
!! make Coulomb interation U rank-4 tensor, according to Slater-Cordon
!! parameterized Hamiltonian
!!
  subroutine atomic_make_umatS()
     use constants, only : dp
     use constants, only : zero, half
     use constants, only : czero

     use control, only : nband, norbs

     use m_spmat, only : umat

     implicit none

!! local variables
     ! orbital momentum quantum number
     integer  :: l

     ! loop index
     integer  :: i

     ! orbital index
     integer  :: m, n
     integer  :: q, p

     ! band index and spin index
     integer  :: m_b, m_s
     integer  :: n_b, n_s
     integer  :: q_b, q_s
     integer  :: p_b, p_s

     ! dummy variables
     real(dp) :: res

     ! c^k_l(m_1,m_2) coefficients
     real(dp), allocatable :: gaunt(:,:,:)

     ! Slater-Cordon parameters: F0, F2, F4, and F6
     real(dp), allocatable :: slater_cordon(:)

!! [body

     ! allocate memory for slater_cordon and gaunt and then build them
     select case (nband)

         case (3)
             l = 1
             !
             allocate(slater_cordon(0:2*l))
             call atomic_make_slater3(slater_cordon)
             !
             allocate(gaunt(-l:l,-l:l,0:2*l))
             call atomic_make_gaunt3(gaunt)

         case (5)
             l = 2
             !
             allocate(slater_cordon(0:2*l))
             call atomic_make_slater5(slater_cordon)
             !
             allocate(gaunt(-l:l,-l:l,0:2*l))
             call atomic_make_gaunt5(gaunt)

         case (7)
             l = 3
             !
             allocate(slater_cordon(0:2*l))
             call atomic_make_slater7(slater_cordon)
             !
             allocate(gaunt(-l:l,-l:l,0:2*l))
             call atomic_make_gaunt7(gaunt)

         case default
             call s_print_error('atomic_make_umatS', &
                 & 'not implemented for this nband!')

     end select

     ! make Coulomb interaction U matrix
     do m=1,norbs
     do n=1,norbs

         ! determine band index and spin index
         ! the orbital order is : up dn up dn up dn ...
         m_b = ( m - 1 ) / 2 - l
         n_b = ( n - 1 ) / 2 - l
         m_s = mod(m,2)
         n_s = mod(n,2)

         do p=1,norbs
         do q=1,norbs

             ! determine band index and spin index
             ! the orbital order is : up dn up dn up dn ...
             q_b = ( q - 1 ) / 2 - l
             p_b = ( p - 1 ) / 2 - l
             q_s = mod(q,2)
             p_s = mod(p,2)

             if ( ( m_b + n_b ) /= ( q_b + p_b ) ) CYCLE
             if ( ( m_s /= p_s ) .or. ( n_s /= q_s ) ) CYCLE

             res = zero
             do i=0,2*l,2
                 res = res + gaunt(m_b,p_b,i) * gaunt(q_b,n_b,i) * slater_cordon(i)
             enddo ! over i={0,2*l} loop

             ! transform the orbital order to: up up up ... dn dn dn ...
             umat( (m + 1) / 2 + nband * (1 - m_s), &
                   (n + 1) / 2 + nband * (1 - n_s), &
                   (q + 1) / 2 + nband * (1 - q_s), &
                   (p + 1) / 2 + nband * (1 - p_s) ) = res

         enddo ! over q={1,norbs} loop
         enddo ! over p={1,norbs} loop

     enddo ! over n={1,norbs} loop
     enddo ! over m={1,norbs} loop
     !
     umat = half * umat

     ! deallocate memory
     if ( allocated(gaunt) )         deallocate(gaunt)
     if ( allocated(slater_cordon) ) deallocate(slater_cordon)

!! body]

     return
  end subroutine atomic_make_umatS

!!========================================================================
!!>>> determine spin-orbit coupling matrix                             <<<
!!========================================================================

!!
!! @sub atomic_make_smat3
!!
!! make spin-orbit coupling matrix for 3-band case. it is defined in the
!! complex orbital basis: Y^{m}_{l}(\theta,\phi).
!!
  subroutine atomic_make_smat3(smat)
     use constants, only : dp
     use constants, only : one, two
     use constants, only : czero

     implicit none

!! external arguments
     ! spin-orbit coupling matrix
     complex(dp), intent(out) :: smat(6,6)

!! local parameters
     ! \sqrt{2}
     real(dp), parameter :: sqrt2 = sqrt(two)

!! [body

     ! make SOC on complex orbital basis, the orbital order is:
     !     | l,  m, up or down >
     !
     !     | 1, -1, up >,
     !     | 1,  0, up >,
     !     | 1,  1, up >,
     !     | 1, -1, dn >,
     !     | 1,  0, dn >,
     !     | 1,  1, dn >
     !
     smat = czero
     !
     smat( 1, 1 ) = -one
     smat( 1, 5 ) = sqrt2
     smat( 2, 6 ) = sqrt2
     smat( 3, 3 ) = one
     smat( 4, 4 ) = one
     smat( 5, 1 ) = sqrt2
     smat( 6, 2 ) = sqrt2
     smat( 6, 6 ) = -one

!! body]

     return
  end subroutine atomic_make_smat3

!!
!! @sub atomic_make_smat5
!!
!! make spin-orbit coupling matrix for 5-band case. it is defined in the
!! complex orbital basis: Y^{m}_{l}(\theta,\phi).
!!
  subroutine atomic_make_smat5(smat)
     use constants, only : dp
     use constants, only : one, two
     use constants, only : czero

     implicit none

!! external arguments
     ! spin-orbit coupling matrix
     complex(dp), intent(out) :: smat(10,10)

!! local parameters
     ! \sqrt{6}
     real(dp), parameter :: sqrt6 = sqrt(6.0_dp)

!! [body

     ! make SOC on complex orbital basis, the orbital order is:
     !     | l,  m, up or down >
     !
     !     | 2, -2, up >,
     !     | 2, -1, up >,
     !     | 2,  0, up >,
     !     | 2,  1, up >,
     !     | 2,  2, up >,
     !     | 2, -2, dn >,
     !     | 2, -1, dn >,
     !     | 2,  0, dn >,
     !     | 2,  1, dn >,
     !     | 2,  2, dn >
     !
     smat = czero
     !
     smat(  1,  1 ) = -two
     smat(  1,  7 ) = two
     smat(  2,  2 ) = -one
     smat(  2,  8 ) = sqrt6
     smat(  3,  9 ) = sqrt6
     smat(  4,  4 ) = one
     smat(  4, 10 ) = two
     smat(  5,  5 ) = two
     smat(  6,  6 ) = two
     smat(  7,  1 ) = two
     smat(  7,  7 ) = one
     smat(  8,  2 ) = sqrt6
     smat(  9,  3 ) = sqrt6
     smat(  9,  9 ) = -one
     smat( 10,  4 ) = two
     smat( 10, 10 ) = -two

!! body]

     return
  end subroutine atomic_make_smat5

!!
!! @sub atomic_make_smat7
!!
!! make spin-orbit coupling matrix for 7-band case, it is defined in the
!! complex orbital basis: Y^{m}_{l}(\theta,\phi).
!!
  subroutine atomic_make_smat7(smat)
     use constants, only : dp
     use constants, only : one, two
     use constants, only : czero

     implicit none

!! external arguments
     ! spin-orbit coupling matrix
     complex(dp), intent(out) :: smat(14,14)

!! local parameters
     ! three, \sqrt{6}, \sqrt{10}, and \sqrt{12}
     real(dp), parameter :: three  = 3.0_dp
     real(dp), parameter :: sqrt6  = sqrt( 6.0_dp)
     real(dp), parameter :: sqrt10 = sqrt(10.0_dp)
     real(dp), parameter :: sqrt12 = sqrt(12.0_dp)

!! [body

     ! make SOC on complex orbital basis, the orbital order is:
     !     | l,  m, up or down >
     !
     !     | 3, -3, up >,
     !     | 3, -2, up >,
     !     | 3, -1, up >,
     !     | 3,  0, up >,
     !     | 3,  1, up >,
     !     | 3,  2, up >,
     !     | 3,  3, up >,
     !     | 3, -3, dn >,
     !     | 3, -2, dn >,
     !     | 3, -1, dn >,
     !     | 3,  0, dn >,
     !     | 3,  1, dn >,
     !     | 3,  2, dn >,
     !     | 3,  3, dn >
     !
     smat = czero
     !
     smat(  1,  1 ) = -three
     smat(  1,  9 ) = sqrt6
     smat(  2,  2 ) = -two
     smat(  2, 10 ) = sqrt10
     smat(  3,  3 ) = -one
     smat(  3, 11 ) = sqrt12
     smat(  4, 12 ) = sqrt12
     smat(  5,  5 ) = one
     smat(  5, 13 ) = sqrt10
     smat(  6,  6 ) = two
     smat(  6, 14 ) = sqrt6
     smat(  7,  7 ) = three
     smat(  8,  8 ) = three
     smat(  9,  1 ) = sqrt6
     smat(  9,  9 ) = two
     smat( 10,  2 ) = sqrt10
     smat( 10, 10 ) = one
     smat( 11,  3 ) = sqrt12
     smat( 12,  4 ) = sqrt12
     smat( 12, 12 ) = -one
     smat( 13,  5 ) = sqrt10
     smat( 13, 13 ) = -two
     smat( 14,  6 ) = sqrt6
     smat( 14, 14 ) = -three

!! body]

     return
  end subroutine atomic_make_smat7

!!========================================================================
!!>>> determine transformation matrix for single particle basis        <<<
!!========================================================================

!!
!! @sub atomic_make_tmat_c2r
!!
!! make transformation matrix from complex orbital basis Y^{m}_{l} to
!! real orbital basis Y_{lm}
!!
  subroutine atomic_make_tmat_c2r(tmat_c2r)
     use constants, only : dp
     use constants, only : czero, cone, czi

     use control, only : nband, norbs

     implicit none

!! external arguments
     ! transformation matrix from complex orbitals to real orbitals
     complex(dp), intent(out) :: tmat_c2r(norbs,norbs)

!! local parameters
     ! \sqrt{2}
     real(dp), parameter :: sqrt2 = sqrt(2.0_dp)

!! [body

     tmat_c2r = czero
     !
     select case (nband)

         ! the real orbital order (t2g) is:
         !     | dxz, up >,
         !     | dxy, up >,
         !     | dyz, up >,
         !     | dxz, dn >,
         !     | dxy, dn >,
         !     | dyz, dn >
         !
         ! the corresponding p orbital order is:
         !     | py, up >,
         !     | pz, up >,
         !     | px, up >,
         !     | py, dn >,
         !     | pz, dn >,
         !     | px, dn >
         !
         ! the complex orbital order is
         !     | 1, -1, up >,
         !     | 1,  0, up >,
         !     | 1,  1, up >,
         !     | 1, -1, dn >,
         !     | 1,  0, dn >,
         !     | 1,  1, dn >
         !
         case (3)
             tmat_c2r( 1, 1 ) = czi / sqrt2
             tmat_c2r( 1, 3 ) = cone / sqrt2
             tmat_c2r( 2, 2 ) = cone
             tmat_c2r( 3, 1 ) = czi / sqrt2
             tmat_c2r( 3, 3 ) = -cone / sqrt2
             !
             tmat_c2r( 4, 4 ) = czi / sqrt2
             tmat_c2r( 4, 6 ) = cone / sqrt2
             tmat_c2r( 5, 5 ) = cone
             tmat_c2r( 6, 4 ) = czi / sqrt2
             tmat_c2r( 6, 6 ) = -cone / sqrt2

         ! the real orbital order is:
         !     | dxy, up >,
         !     | dyz, up >,
         !     | dz2, up >,
         !     | dxz, up >,
         !     | dx2-y2, up >,
         !     | dxy, dn >,
         !     | dyz, dn >,
         !     | dz2, dn >,
         !     | dxz, dn >,
         !     | dx2-y2, dn >
         !
         ! the complex orbital order is:
         !     | 2, -2, up >,
         !     | 2, -1, up >,
         !     | 2,  0, up >,
         !     | 2,  1, up >,
         !     | 2,  2, up >,
         !     | 2, -2, dn >,
         !     | 2, -1, dn >,
         !     | 2,  0, dn >,
         !     | 2,  1, dn >,
         !     | 2,  2, dn >
         !
         case (5)
             tmat_c2r(  1,  1 ) = czi / sqrt2
             tmat_c2r(  1,  5 ) = cone / sqrt2
             tmat_c2r(  2,  2 ) = czi / sqrt2
             tmat_c2r(  2,  4 ) = cone / sqrt2
             tmat_c2r(  3,  3 ) = cone
             tmat_c2r(  4,  2 ) = czi / sqrt2
             tmat_c2r(  4,  4 ) = -cone / sqrt2
             tmat_c2r(  5,  1 ) = -czi /sqrt2
             tmat_c2r(  5,  5 ) = cone / sqrt2
             !
             tmat_c2r(  6,  6 ) = czi / sqrt2
             tmat_c2r(  6, 10 ) = cone / sqrt2
             tmat_c2r(  7,  7 ) = czi / sqrt2
             tmat_c2r(  7,  9 ) = cone / sqrt2
             tmat_c2r(  8,  8 ) = cone
             tmat_c2r(  9,  7 ) = czi / sqrt2
             tmat_c2r(  9,  9 ) = -cone / sqrt2
             tmat_c2r( 10,  6 ) = -czi / sqrt2
             tmat_c2r( 10, 10 ) = cone / sqrt2

         ! the real orbital order is:
         !     | fy(3x2-y2), up >,
         !     | fxyz, up >,
         !     | fyz2, up >,
         !     | fz3, up >,
         !     | fxz2, up >,
         !     | fz(x2-y2), up >,
         !     | fx(x2-3y2), up >,
         !     | fy(3x2-y2), dn >,
         !     | fxyz, dn >,
         !     | fyz2, dn >,
         !     | fz3, dn >,
         !     | fxz2, dn >,
         !     | fz(x2-y2), dn >,
         !     | fx(x2-3y2), dn >
         !
         ! the complex orbital order is:
         !     | 3, -3, up >,
         !     | 3, -2, up >,
         !     | 3, -1, up >,
         !     | 3,  0, up >,
         !     | 3,  1, up >,
         !     | 3,  2, up >,
         !     | 3,  3, up >,
         !     | 3, -3, dn >,
         !     | 3, -2, dn >,
         !     | 3, -1, dn >,
         !     | 3,  0, dn >,
         !     | 3,  1, dn >,
         !     | 3,  2, dn >,
         !     | 3,  3, dn >
         !
         case (7)
             tmat_c2r(  1,  1 ) = czi/sqrt2
             tmat_c2r(  1,  7 ) = cone/sqrt2
             tmat_c2r(  2,  2 ) = czi/sqrt2
             tmat_c2r(  2,  6 ) = cone/sqrt2
             tmat_c2r(  3,  3 ) = czi/sqrt2
             tmat_c2r(  3,  5 ) = cone/sqrt2
             tmat_c2r(  4,  4 ) = cone
             tmat_c2r(  5,  3 ) = czi/sqrt2
             tmat_c2r(  5,  5 ) = -cone/sqrt2
             tmat_c2r(  6,  2 ) = -czi/sqrt2
             tmat_c2r(  6,  6 ) = cone/sqrt2
             tmat_c2r(  7,  1 ) = czi/sqrt2
             tmat_c2r(  7,  7 ) = -cone/sqrt2
             !
             tmat_c2r(  8,  8 ) = czi/sqrt2
             tmat_c2r(  8, 14 ) = cone/sqrt2
             tmat_c2r(  9,  9 ) = czi/sqrt2
             tmat_c2r(  9, 13 ) = cone/sqrt2
             tmat_c2r( 10, 10 ) = czi/sqrt2
             tmat_c2r( 10, 12 ) = cone/sqrt2
             tmat_c2r( 11, 11 ) = cone
             tmat_c2r( 12, 10 ) = czi/sqrt2
             tmat_c2r( 12, 12 ) = -cone/sqrt2
             tmat_c2r( 13,  9 ) = -czi/sqrt2
             tmat_c2r( 13, 13 ) = cone/sqrt2
             tmat_c2r( 14,  8 ) = czi/sqrt2
             tmat_c2r( 14, 14 ) = -cone/sqrt2

         case default
             call s_print_error('atomic_make_tmat_c2r', &
                 & 'not implemented for this nband!')

     end select

!! body]

     return
  end subroutine atomic_make_tmat_c2r

!!
!! @sub atomic_make_tmat_r2c
!!
!! make transformation matrix from real orbital basis Y_{lm} to
!! complex orbital basis Y^{m}_{l}
!!
  subroutine atomic_make_tmat_r2c(tmat_r2c)
     use constants, only : dp
     use constants, only : czero

     use control, only : norbs

!! external arguments
     ! transformation matrix from real orbitals to complex orbitals
     complex(dp), intent(out) :: tmat_r2c(norbs,norbs)

!! local variables
     ! dummy array, transpose of tmat_r2c
     complex(dp) :: tmat_c2r(norbs,norbs)

!! [body

     tmat_c2r = czero
     call atomic_make_tmat_c2r(tmat_c2r)
     tmat_r2c = transpose( dconjg(tmat_c2r) )

!! body]

     return
  end subroutine atomic_make_tmat_r2c

!!
!! @sub atomic_make_tmat_c2j
!!
!! make transformation matrix from complex orbital basis Y^{m}_{l} to
!! j^2-j_z orbital basis
!!
  subroutine atomic_make_tmat_c2j(tmat_c2j)
     use constants, only : dp
     use constants, only : one, two
     use constants, only : czero

     use control, only : nband, norbs

     implicit none

!! external arguments
     ! transformation matrix from complex orbitals to j^2-j_z basis
     complex(dp), intent(out) :: tmat_c2j(norbs,norbs)

!! local parameters
     real(dp), parameter :: three = 3.0_dp
     real(dp), parameter :: four  = 4.0_dp
     real(dp), parameter :: five  = 5.0_dp
     real(dp), parameter :: six   = 6.0_dp
     real(dp), parameter :: seven = 7.0_dp

!! [body

     tmat_c2j = czero
     !
     select case (nband)

         ! the complex orbital order is:
         !     | 1, -1, up >,
         !     | 1,  0, up >,
         !     | 1,  1, up >,
         !     | 1, -1, dn >,
         !     | 1,  0, dn >,
         !     | 1,  1, dn >
         !
         ! the j^2-j_z basis order is:
         !     | 1/2, -1/2 >,
         !     | 1/2,  1/2 >,
         !     | 3/2, -3/2 >,
         !     | 3/2, -1/2 >,
         !     | 3/2,  1/2 >,
         !     | 3/2,  3/2 >
         !
         case (3)
             tmat_c2j( 1, 1 ) = -sqrt(two/three)
             tmat_c2j( 1, 4 ) = sqrt(one/three)
             tmat_c2j( 2, 2 ) = -sqrt(one/three)
             tmat_c2j( 2, 5 ) = sqrt(two/three)
             tmat_c2j( 3, 6 ) = one
             tmat_c2j( 4, 3 ) = one
             tmat_c2j( 5, 1 ) = sqrt(one/three)
             tmat_c2j( 5, 4 ) = sqrt(two/three)
             tmat_c2j( 6, 2 ) = sqrt(two/three)
             tmat_c2j( 6, 5 ) = sqrt(one/three)

         ! the complex orbital order is:
         !     | 2, -2, up >,
         !     | 2, -1, up >,
         !     | 2,  0, up >,
         !     | 2,  1, up >,
         !     | 2,  2, up >,
         !     | 2, -2, dn >,
         !     | 2, -1, dn >,
         !     | 2,  0, dn >,
         !     | 2,  1, dn >,
         !     | 2,  2, dn >
         !
         ! the j^2-j_z basis order is:
         !     | 3/2, -3/2 >,
         !     | 3/2, -1/2 >,
         !     | 3/2,  1/2 >,
         !     | 3/2,  3/2 >,
         !     | 5/2, -5/2 >,
         !     | 5/2, -3/2 >,
         !     | 5/2, -1/2 >,
         !     | 5/2,  1/2 >,
         !     | 5/2,  3/2 >,
         !     | 5/2,  5/2 >
         !
         case (5)
             tmat_c2j(  1,  1 ) = -sqrt(four/five)
             tmat_c2j(  1,  6 ) = sqrt(one/five)
             tmat_c2j(  2,  2 ) = -sqrt(three/five)
             tmat_c2j(  2,  7 ) = sqrt(two/five)
             tmat_c2j(  3,  3 ) = -sqrt(two/five)
             tmat_c2j(  3,  8 ) = sqrt(three/five)
             tmat_c2j(  4,  4 ) = -sqrt(one/five)
             tmat_c2j(  4,  9 ) = sqrt(four/five)
             tmat_c2j(  5, 10 ) = one
             tmat_c2j(  6,  5 ) = one
             tmat_c2j(  7,  1 ) = sqrt(one/five)
             tmat_c2j(  7,  6 ) = sqrt(four/five)
             tmat_c2j(  8,  2 ) = sqrt(two/five)
             tmat_c2j(  8,  7 ) = sqrt(three/five)
             tmat_c2j(  9,  3 ) = sqrt(three/five)
             tmat_c2j(  9,  8 ) = sqrt(two/five)
             tmat_c2j( 10,  4 ) = sqrt(four/five)
             tmat_c2j( 10,  9 ) = sqrt(one/five)

         ! the complex orbital order is:
         !     | 3, -3, up >,
         !     | 3, -2, up >,
         !     | 3, -1, up >,
         !     | 3,  0, up >,
         !     | 3,  1, up >,
         !     | 3,  2, up >,
         !     | 3,  3, up >,
         !     | 3, -3, dn >,
         !     | 3, -2, dn >,
         !     | 3, -1, dn >,
         !     | 3,  0, dn >,
         !     | 3,  1, dn >,
         !     | 3,  2, dn >,
         !     | 3,  3, dn >
         !
         ! the j^2-j_z basis order is:
         !     | 5/2, -5/2 >,
         !     | 5/2, -3/2 >,
         !     | 5/2, -1/2 >,
         !     | 5/2,  1/2 >,
         !     | 5/2,  3/2 >,
         !     | 5/2,  5/2 >,
         !     | 7/2, -7/2 >,
         !     | 7/2, -5/2 >,
         !     | 7/2, -3/2 >,
         !     | 7/2, -1/2 >,
         !     | 7/2,  1/2 >,
         !     | 7/2,  3/2 >,
         !     | 7/2,  5/2 >,
         !     | 7/2,  7/2 >
         !
         case (7)
             tmat_c2j(  1,  1 ) = -sqrt(six/seven)
             tmat_c2j(  1,  8 ) = sqrt(one/seven)
             tmat_c2j(  2,  2 ) = -sqrt(five/seven)
             tmat_c2j(  2,  9 ) = sqrt(two/seven)
             tmat_c2j(  3,  3 ) = -sqrt(four/seven)
             tmat_c2j(  3, 10 ) = sqrt(three/seven)
             tmat_c2j(  4,  4 ) = -sqrt(three/seven)
             tmat_c2j(  4, 11 ) = sqrt(four/seven)
             tmat_c2j(  5,  5 ) = -sqrt(two/seven)
             tmat_c2j(  5, 12 ) = sqrt(five/seven)
             tmat_c2j(  6,  6 ) = -sqrt(one/seven)
             tmat_c2j(  6, 13 ) = sqrt(six/seven)
             tmat_c2j(  7, 14 ) = one
             tmat_c2j(  8,  7 ) = one
             tmat_c2j(  9,  1 ) = sqrt(one/seven)
             tmat_c2j(  9,  8 ) = sqrt(six/seven)
             tmat_c2j( 10,  2 ) = sqrt(two/seven)
             tmat_c2j( 10,  9 ) = sqrt(five/seven)
             tmat_c2j( 11,  3 ) = sqrt(three/seven)
             tmat_c2j( 11, 10 ) = sqrt(four/seven)
             tmat_c2j( 12,  4 ) = sqrt(four/seven)
             tmat_c2j( 12, 11 ) = sqrt(three/seven)
             tmat_c2j( 13,  5 ) = sqrt(five/seven)
             tmat_c2j( 13, 12 ) = sqrt(two/seven)
             tmat_c2j( 14,  6 ) = sqrt(six/seven)
             tmat_c2j( 14, 13 ) = sqrt(one/seven)

         case default
             call s_print_error('atomic_make_tmat_c2j', &
                 & 'not implemented for this nband!')

     end select

!! body]

     return
  end subroutine atomic_make_tmat_c2j

!!========================================================================
!!>>> perform representation transformation                            <<<
!!========================================================================

!!
!! @sub atomic_tran_fmat
!!
!! rotate annihilation or creation operator matrix (fmat) from Fock basis
!! to atomic eigenbasis
!!
  subroutine atomic_tran_fmat(ndimx, ndimy, amat, fmat, cmat)
     use constants, only : dp
     use constants, only : zero, one

     implicit none

!! external arguments
     ! x dimension of matrix
     integer, intent(in)  :: ndimx

     ! y dimension of matrix
     integer, intent(in)  :: ndimy

     ! left transformation matrix
     real(dp), intent(in) :: amat(ndimx,ndimx)

     ! right transformation matrix
     real(dp), intent(in) :: cmat(ndimy,ndimy)

     ! annihilation or creation operator matrix
     real(dp), intent(inout) :: fmat(ndimx,ndimy)

!! local variables
     ! dummy array
     real(dp), allocatable :: tmp_mat(:,:)

!! [body

     ! allocate memory
     allocate(tmp_mat(ndimx,ndimy))

     tmp_mat = zero
     call dgemm('N', 'N', ndimx, ndimy, ndimy, &
                             one, fmat, ndimx, &
                                  cmat, ndimy, &
                         zero, tmp_mat, ndimx  )

     call dgemm('T', 'N', ndimx, ndimy, ndimx, &
                             one, amat, ndimx, &
                               tmp_mat, ndimx, &
                            zero, fmat, ndimx  )

     ! deallocate memory
     deallocate(tmp_mat)

!! body]

     return
  end subroutine atomic_tran_fmat

!!
!! @sub atomic_tran_umat
!!
!! transform Coulomb interaction U tensor from one representation
!! to another representation
!!
  subroutine atomic_tran_umat(amtrx, umat, umat_t)
     use constants, only : dp
     use constants, only : czero
     use constants, only : epst

     use control, only : norbs

     implicit none

!! external arguments
     ! transformation matrix from orginal basis to natural eigenbasis
     !
     ! the original basis could be real orbial basis or complex
     ! orbital basis. it depends on how to determine the Coulomb
     ! interaction matrix
     complex(dp), intent(in)  :: amtrx(norbs,norbs)

     ! coefficents matrix for general interaction U in orginal basis
     complex(dp), intent(in)  :: umat(norbs,norbs,norbs,norbs)

     ! coefficents matrix for general interaction U in natural eigenbasis
     complex(dp), intent(out) :: umat_t(norbs,norbs,norbs,norbs)

!! local variables
     ! loop index over orbitals in orginal single particle basis
     integer :: a1, a2
     integer :: a3, a4
     integer :: s1, s2
     integer :: s3, s4

     ! auxiliary complex(dp) variables
     complex(dp) :: ctmp

!! [body

     ! initialize umat_t to be zero
     umat_t = czero
     !
     s1loop: do s1=1,norbs
     s2loop: do s2=1,norbs
     s3loop: do s3=1,norbs
     s4loop: do s4=1,norbs
         !
         ctmp = czero
         !
         a1loop: do a1=1,norbs
         a2loop: do a2=1,norbs
         a3loop: do a3=1,norbs
         a4loop: do a4=1,norbs

             if ( abs( umat(a1,a2,a3,a4) ) < epst ) CYCLE

             ctmp = ctmp + umat(a1,a2,a3,a4)                   &
                         * conjg(amtrx(a1,s1)) * amtrx(a3,s3)  &
                         * conjg(amtrx(a2,s2)) * amtrx(a4,s4)

         enddo a4loop ! over a4={1,norbs} loop
         enddo a3loop ! over a3={1,norbs} loop
         enddo a2loop ! over a2={1,norbs} loop
         enddo a1loop ! over a1={1,norbs} loop
         !
         umat_t(s1,s2,s3,s4) = ctmp
         !
     enddo s4loop ! over s4={1,norbs} loop
     enddo s3loop ! over s3={1,norbs} loop
     enddo s2loop ! over s2={1,norbs} loop
     enddo s1loop ! over s1={1,norbs} loop

!! body]

     return
  end subroutine atomic_tran_umat

!!
!! @sub atomic_tran_repr_cmpl
!!
!! transformation from one representation to another representation,
!! complex version
!!
  subroutine atomic_tran_repr_cmpl(ndim, amat, tmat)
     use constants, only : dp
     use constants, only : cone, czero

     implicit none

!! external arguments
     ! size of the matrix
     integer, intent(in) :: ndim

     ! transformation matrix
     complex(dp), intent(in) :: tmat(ndim,ndim)

     ! physical quantities
     complex(dp), intent(inout) :: amat(ndim,ndim)

!! local variables
     ! dummy matrix
     complex(dp) :: tmp_mat(ndim,ndim)

!! [body

     call zgemm('N', 'N', ndim, ndim, ndim, &
                          cone, amat, ndim, &
                                tmat, ndim, &
                      czero, tmp_mat, ndim  )

     call zgemm('C', 'N', ndim, ndim, ndim, &
                          cone, tmat, ndim, &
                             tmp_mat, ndim, &
                         czero, amat, ndim  )

!! body]

     return
  end subroutine atomic_tran_repr_cmpl

!!
!! @sub atomic_tran_repr_real
!!
!! transformation from one representation to another representation,
!! real version
!!
  subroutine atomic_tran_repr_real(ndim, amat, tmat)
     use constants, only : dp
     use constants, only : zero, one

     implicit none

!! external arguments
     ! size of the matrix
     integer, intent(in) :: ndim

     ! transformation matrix
     real(dp), intent(in) :: tmat(ndim,ndim)

     ! physical quantities
     real(dp), intent(inout) :: amat(ndim,ndim)

!! local variables
     ! dummy matrix
     real(dp) :: tmp_mat(ndim,ndim)

!! [body

     call dgemm('N', 'N', ndim, ndim, ndim, &
                           one, amat, ndim, &
                                tmat, ndim, &
                       zero, tmp_mat, ndim  )

     call dgemm('T', 'N', ndim, ndim, ndim, &
                           one, tmat, ndim, &
                             tmp_mat, ndim, &
                          zero, amat, ndim  )

!! body]

     return
  end subroutine atomic_tran_repr_real

!!========================================================================
!!>>> generate natural eigenbasis                                      <<<
!!========================================================================

!!
!! @sub atomic_natural_basis1
!!
!! make natural eigenbasis, the onsite energy of impurity (emat) and
!! the transformation matrix from original basis to natural eigenbasis
!! are determined as well.
!!
!! case 1
!!
!! no crystal field splitting or it is diagonal
!! no spin-orbit coupling
!!
  subroutine atomic_natural_basis1()
     use control, only : norbs
     use control, only : mune

     use m_spmat, only : cmat, emat, tmat

     implicit none

!! local variables
     ! loop index
     integer :: i

!! [body

     ! setup emat to crystal field splitting
     ! since smat is zero, so emat is equal to cmat
     emat = cmat

     ! add chemical potential to emat
     do i=1,norbs
         emat(i,i) = emat(i,i) + mune
     enddo ! over i={1,norbs} loop

     ! for this case, the natural eigenbasis is the real orbital basis
     ! so, the tmat is a unity matrix
     call s_identity_z(norbs, tmat)

!! body]

     return
  end subroutine atomic_natural_basis1

!!
!! @sub atomic_natural_basis2
!!
!! make natural eigenbasis, the onsite energy of impurity (emat) and
!! the transformation matrix from original basis to natural eigenbasis
!! are determined as well.
!!
!! case 2
!!
!! non-diagonal crystal field splitting
!! no spin-orbit coupling
!!
  subroutine atomic_natural_basis2()
     use constants, only : dp

     use control, only : nband, norbs
     use control, only : mune

     use m_spmat, only : cmat, emat, tmat

     implicit none

!! local variables
     ! loop index
     integer  :: i
     integer  :: j

     ! eigenvalues
     real(dp) :: eigval(nband)

     ! eigenvectors
     real(dp) :: eigvec(nband,nband)

     ! emat matrix for no spin freedom
     complex(dp) :: emat_nospin(nband,nband)

     ! tmat for no spin freedom
     complex(dp) :: tmat_nospin(nband,nband)

!! [body

     ! setup emat to crystal field splitting
     ! since smat is zero, so emat is equal to cmat
     emat = cmat

     ! get emat for no spin freedom
     do i=1,nband
         do j=1,nband
             emat_nospin(i,j) = emat(i,j)
         enddo ! over j={1,nband} loop
     enddo ! over i={1,nband}

     ! diagonalize emat_nospin to get natural eigenbasis
     call s_eig_sy(nband, nband, real(emat_nospin), eigval, eigvec)

     ! get diagonal emat for no spin freedom
     call s_diag_z(nband, dcmplx(eigval), emat_nospin)

     ! get tmat for no spin freedom
     tmat_nospin = dcmplx(eigvec)

     ! build emat and tmat with spin freedom
     do i=1,nband
         do j=1,nband
             emat(i,j) = emat_nospin(i,j)
             emat(i+nband,j+nband) = emat_nospin(i,j)
             tmat(i,j) = tmat_nospin(i,j)
             tmat(i+nband,j+nband) = tmat_nospin(i,j)
         enddo ! over j={1,nband} loop
     enddo ! over i={1,nband} loop

     ! add chemical potential to emat
     do i=1,norbs
         emat(i,i) = emat(i,i) + mune
     enddo ! over i={1,norbs} loop

!! body]

     return
  end subroutine atomic_natural_basis2

!!
!! @sub atomic_natural_basis3
!!
!! make natural eigenbasis, the onsite energy of impurity (emat) and
!! the transformation matrix from original basis to natural eigenbasis
!! are determined as well.
!!
!! case 3
!!
!! no crystal field splitting
!! spin-orbit coupling
!!
!! for this special case, the natural eigenbasis is |j^2,j_z>
!!
  subroutine atomic_natural_basis3()
     use constants, only : dp

     use control, only : norbs
     use control, only : mune

     use m_spmat, only : emat, smat, tmat

     implicit none

!! local variables
     ! loop inex
     integer :: i

     ! transformation matrix from complex orbital basis to |j^2,j_z> basis
     complex(dp) :: tmat_c2j(norbs,norbs)

!! [body

     ! setup emat to spin-orbit coupling
     ! since cmat is zero, so emat is equal to smat
     emat = smat

     ! evaluate transformation matrix tmat_c2j
     call atomic_make_tmat_c2j(tmat_c2j)

     ! the transformation matrix is from complex orbital
     ! basis to natural eigenbasis (|j^2,j_z> basis)
     tmat = tmat_c2j

     ! transform emat to natural eigenbasis
     call atomic_tran_repr_cmpl(norbs, emat, tmat)

     ! add chemical potential to emat
     do i=1,norbs
         emat(i,i) = emat(i,i) + mune
     enddo ! over i={1,norbs} loop

!! body]

     return
  end subroutine atomic_natural_basis3

!!
!! @sub atomic_natural_basis4
!!
!! make natural eigenbasis, the onsite energy of impurity (emat) and
!! the transformation matrix from original basis to natural eigenbasis
!! are determined as well.
!!
!! case 4
!!
!! crystal field splitting
!! spin-orbit coupling
!!
  subroutine atomic_natural_basis4()
     use constants, only : dp
     use constants, only : eps6

     use control, only : norbs
     use control, only : mune

     use m_spmat, only : cmat, smat, emat, tmat

     implicit none

!! local variables
     ! loop index
     integer  :: i

     ! eigenvalues
     real(dp) :: eigval(norbs)

     ! eigenvectors
     real(dp) :: eigvec(norbs,norbs)

     ! transformation matrix from real orbital basis
     ! to complex orbital basis
     complex(dp) :: tmat_r2c(norbs,norbs)

     ! transformation matrix from complex orbital basis
     ! to natural eigenbasis
     complex(dp) :: tmat_c2n(norbs,norbs)

!! [body

     ! build tmat_r2c
     call atomic_make_tmat_r2c(tmat_r2c)

     ! transfrom crystal field splitting (cmat) to complex orbital basis
     call atomic_tran_repr_cmpl(norbs, cmat, tmat_r2c)

     ! check whether cmat is real, if not, we cann't make natural eigenbasis
     if ( any( abs( aimag(cmat) ) > eps6 ) ) then
         call s_print_error('atomic_natural_basis4', &
             & 'crystal field on complex orbital basis should be real!')
     endif ! back if ( any( abs( aimag(cmat) ) > eps6 ) ) block

     ! setup emat: crystal field splitting + spin-orbit coupling
     emat = smat + cmat

     ! diagonalize real(emat)
     call s_eig_sy(norbs, norbs, real(emat), eigval, eigvec)

     ! build transformation matrix from complex orbital basis
     ! to natural eigenbasis
     tmat_c2n = eigvec
     tmat = tmat_c2n

     ! transform emat from complex orbital basis to natural eigenbasis
     call atomic_tran_repr_cmpl(norbs, emat, tmat_c2n)

     ! add chemical poential to emat
     do i=1,norbs
         emat(i,i) = emat(i,i) + mune
     enddo ! over i={1,norbs} loop

!! body]

     return
  end subroutine atomic_natural_basis4
