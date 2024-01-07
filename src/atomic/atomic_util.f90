!!!-----------------------------------------------------------------------
!!! project : iqist @ jasmine
!!! program : atomic_make_cdagger
!!!           atomic_make_c
!!!           atomic_make_gsz
!!!           atomic_make_gjz
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
!!!           atomic_2natural_case1
!!!           atomic_2natural_case2
!!!           atomic_2natural_case3
!!!           atomic_2natural_case4
!!! source  : atomic_util.f90
!!! type    : subroutines
!!! author  : yilin wang (email:qhwyl2006@126.com)
!!! history : 07/09/2014 by yilin wang (created)
!!!           01/08/2024 by li huang (last modified)
!!! purpose : provide the utility subroutines for the atomic eigenvalue
!!!           problem solver, such as the Dirac algebra, calculations of
!!!           gaunt coefficients and Coulomb interaction matrices, etc.
!!! status  : unstable
!!! comment :
!!!-----------------------------------------------------------------------

!!========================================================================
!!>>> simulate creation and destroy operators                          <<<
!!========================================================================

!!
!! @sub atomic_make_cdagger
!!
!! create one electron on ipos of |jold> to obtain new Fock state |jnew>
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

     ! sign due to anti-commute relation between fermions
     integer, intent(out):: isgn

!! local variables
     ! loop index over orbital
     integer :: iorb

!! [body

     ! it is occupied at ipos
     if ( btest(jold, ipos-1) .eqv. .true. ) then
         call s_print_error('atomic_make_cdagger','severe error happened')
     endif ! back if ( btest(jold, ipos-1) .eqv. .true. ) block

     ! evaluate the sign
     isgn = 0
     do iorb=1,ipos-1
        if ( btest(jold, iorb-1) .eqv. .true. ) isgn = isgn + 1
     enddo ! over iorb={1,ipos-1} loop
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
!! destroy one electron on ipos of |jold> to obtain new Fock state |jnew>
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

      ! sign due to anti-commute relation between fermions
      integer, intent(out) :: isgn

!! local variables
      ! loop index
      integer :: iorb

      ! it is unoccupied at ipos
      if ( btest(jold, ipos-1) .eqv. .false. ) then
          call s_print_error('atomic_make_c','severe error happened')
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
!!>>> determine good quantum numbers                                   <<<
!!========================================================================

!!
!! @sub atomic_make_gsz
!!
!! make Sz quantum number for each orbital
!!
  subroutine atomic_make_gsz(good_sz)
     use control, only : norbs

     implicit none

!! external arguments
     ! good quantum numbers: Sz
     integer, intent(out) :: good_sz(norbs)

!! local variables
     ! loop index
     integer :: i

!! [body

     ! note: the spin arrangement likes up dn up dn up dn ...
     do i=1,norbs
         if ( mod(i,2) /= 0 ) then
             good_sz(i) = +1
         else
             good_sz(i) = -1
         endif ! back if ( mod(i,2) /= 0 ) block
     enddo ! over i={1,norbs} loop

!! body]

     return
  end subroutine atomic_make_gsz

!>>> atomic_make_gjz: make Jz quantum number for each orbital
  subroutine atomic_make_gjz(good_jz)
     use control, only : nband, norbs

     implicit none

! external arguments
! good quantum numbers: Jz
     integer, intent(out) :: good_jz(norbs)

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
             call s_print_error('atomic_make_gjz','not implemented for this norbs value!')

     end select

     return
  end subroutine atomic_make_gjz

!!========================================================================
!!>>> determine gaunt coefficients                                     <<<
!!========================================================================

!!>>> atomic_make_gaunt5: build gaunt coefficients for 5 band case
  subroutine atomic_make_gaunt5(gaunt)
     use constants, only : dp, zero, one

! external arguments
! gaunt coefficients
     real(dp), intent(out) :: gaunt(-2:2,-2:2,0:4)

     gaunt = zero

     gaunt(-2, -2, 0) = one
     gaunt(-1, -1, 0) = one
     gaunt( 0,  0, 0) = one
     gaunt( 1,  1, 0) = one
     gaunt( 2,  2, 0) = one

     gaunt(-2, -2, 2) = -sqrt(4.0_dp/49.0_dp)
     gaunt(-2, -1, 2) =  sqrt(6.0_dp/49.0_dp);   gaunt(-1, -2, 2) = gaunt(-2, -1, 2) * (-1)**(-2+1)
     gaunt(-2,  0, 2) = -sqrt(4.0_dp/49.0_dp);   gaunt( 0, -2, 2) = gaunt(-2,  0, 2) * (-1)**(-2-0)
     gaunt(-1, -1, 2) =  sqrt(1.0_dp/49.0_dp)
     gaunt(-1,  0, 2) =  sqrt(1.0_dp/49.0_dp);   gaunt( 0, -1, 2) = gaunt(-1,  0, 2) * (-1)**(-1-0)
     gaunt(-1,  1, 2) = -sqrt(6.0_dp/49.0_dp);   gaunt( 1, -1, 2) = gaunt(-1,  1, 2) * (-1)**(-1-1)
     gaunt( 0,  0, 2) =  sqrt(4.0_dp/49.0_dp)
     gaunt( 1, -1, 2) = -sqrt(6.0_dp/49.0_dp);   gaunt(-1,  1, 2) = gaunt( 1, -1, 2) * (-1)**( 1+1)
     gaunt( 1,  0, 2) =  sqrt(1.0_dp/49.0_dp);   gaunt( 0,  1, 2) = gaunt( 1,  0, 2) * (-1)**( 1-0)
     gaunt( 1,  1, 2) =  sqrt(1.0_dp/49.0_dp)
     gaunt( 2,  0, 2) = -sqrt(4.0_dp/49.0_dp);   gaunt( 0,  2, 2) = gaunt( 2,  0, 2) * (-1)**( 2-0)
     gaunt( 2,  1, 2) =  sqrt(6.0_dp/49.0_dp);   gaunt( 1,  2, 2) = gaunt( 2,  1, 2) * (-1)**( 2-1)
     gaunt( 2,  2, 2) = -sqrt(4.0_dp/49.0_dp)

     gaunt(-2, -2, 4) =  sqrt( 1.0_dp/441.0_dp)
     gaunt(-2, -1, 4) = -sqrt( 5.0_dp/441.0_dp); gaunt(-1, -2, 4) = gaunt(-2, -1, 4) * (-1)**(-2+1)
     gaunt(-2,  0, 4) =  sqrt(15.0_dp/441.0_dp); gaunt( 0, -2, 4) = gaunt(-2,  0, 4) * (-1)**(-2-0)
     gaunt(-2,  1, 4) = -sqrt(35.0_dp/441.0_dp); gaunt( 1, -2, 4) = gaunt(-2,  1, 4) * (-1)**(-2-1)
     gaunt(-2,  2, 4) =  sqrt(70.0_dp/441.0_dp); gaunt( 2, -2, 4) = gaunt(-2,  2, 4) * (-1)**(-2-2)
     gaunt(-1, -1, 4) = -sqrt(16.0_dp/441.0_dp)
     gaunt(-1,  0, 4) =  sqrt(30.0_dp/441.0_dp); gaunt( 0, -1, 4) = gaunt(-1,  0, 4) * (-1)**(-1-0)
     gaunt(-1,  1, 4) = -sqrt(40.0_dp/441.0_dp); gaunt( 1, -1, 4) = gaunt(-1,  1, 4) * (-1)**(-1-1)
     gaunt( 0,  0, 4) =  sqrt(36.0_dp/441.0_dp)
     gaunt( 1,  0, 4) =  sqrt(30.0_dp/441.0_dp); gaunt( 0,  1, 4) = gaunt( 1,  0, 4) * (-1)**( 1-0)
     gaunt( 1,  1, 4) = -sqrt(16.0_dp/441.0_dp)
     gaunt( 2, -1, 4) = -sqrt(35.0_dp/441.0_dp); gaunt(-1,  2, 4) = gaunt( 2, -1, 4) * (-1)**( 2+1)
     gaunt( 2,  0, 4) =  sqrt(15.0_dp/441.0_dp); gaunt( 0,  2, 4) = gaunt( 2,  0, 4) * (-1)**( 2-0)
     gaunt( 2,  1, 4) = -sqrt( 5.0_dp/441.0_dp); gaunt( 1,  2, 4) = gaunt( 2,  1, 4) * (-1)**( 2-1)
     gaunt( 2,  2, 4) =  sqrt( 1.0_dp/441.0_dp)

     return
  end subroutine atomic_make_gaunt5

!!>>> atomic_make_gaunt7: build gaunt coefficients for 7 band case
  subroutine atomic_make_gaunt7(gaunt)
     use constants, only: dp, zero, one

     implicit none

! external arguments
! gaunt coefficients
     real(dp), intent(out) :: gaunt(-3:3,-3:3,0:6)

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

!!========================================================================
!!>>> determine Coulomb interaction matrix                             <<<
!!========================================================================

!!>>> atomic_make_hund: make the Hund's rule coupling matrix
  subroutine atomic_make_hund(hund)
     use constants, only : dp, zero

     use control, only : icu
     use control, only : nband
     use control, only : Jz, Js, Jp

! external arguments
! Hund's rule coupling matrix
     real(dp), intent(out) :: hund(nband,nband,3)

! local variables
! loop index
     integer  :: i

! dummy variables
     real(dp) :: ff2(3)
     real(dp) :: ff4(3)
     real(dp) :: jj1(3)
     real(dp) :: jj2(3)
     real(dp) :: jj3(3)
     real(dp) :: jj4(3)
     real(dp) :: jzsp(3)

! initialize hund to be zero
     hund = zero

! set jzsp
     jzsp(1) = Jz
     jzsp(2) = Js
     jzsp(3) = Jp

! for isotropic Hund's rule coupling
     if ( icu == 1 ) then
         hund(:,:,1) = Jz
         hund(:,:,2) = Js
         hund(:,:,3) = Jp
         do i=1,nband
             hund(i,i,:) = zero
         enddo ! over i={1,nband} loop
     endif ! back if ( icu == 1 ) block

! for anisotropic Hund's rule coupling
     if ( icu == 3 ) then
         if ( nband == 5 ) then
! J(dxy,dxz) = J(dxy,dyz) = J(dxz,dyz) = J(dxz,dx2) = J(dyz,dx2) = 3/49 * F^2 + 20/441 * F^4
! J(dxy,dz2) = J(dx2,dz2) = 4/49 * F^2 + 15/441 * F^4
! J(dxz,dz2) = J(dyz,dz2) = 1/49 * F^2 + 30/441 * F^4
! J(dxy,dx2) = 35/441 * F^4
! the averaged Hund's rule coupling is J_{ave} = 5/98 * (F^2 + F^4)
! and F^4 = 0.625 * F^2.
             ff2 = (98.0 * jzsp) / (1.625 * 5.0)
             ff4 = 0.625 * ff2
             jj1 = 3.0 / 49.0 * ff2 + 20.0 / 441.0 * ff4
             jj2 = 4.0 / 49.0 * ff2 + 15.0 / 441.0 * ff4
             jj3 = 1.0 / 49.0 * ff2 + 30.0 / 441.0 * ff4
             jj4 = 35.0 / 441.0 * ff4

! orbital order is: (1) dz2, (2) dxz, (3) dyz, (4) dx2, (5) dxy
             hund(2,3,:) = jj1;   hund(3,2,:) = jj1
             hund(2,5,:) = jj1;   hund(5,2,:) = jj1
             hund(3,5,:) = jj1;   hund(5,3,:) = jj1
             hund(2,4,:) = jj1;   hund(4,2,:) = jj1
             hund(3,4,:) = jj1;   hund(4,3,:) = jj1

             hund(1,4,:) = jj2;   hund(4,1,:) = jj2
             hund(1,5,:) = jj2;   hund(5,1,:) = jj2

             hund(1,2,:) = jj3;   hund(2,1,:) = jj3
             hund(1,3,:) = jj3;   hund(3,1,:) = jj3

             hund(4,5,:) = jj4;   hund(5,4,:) = jj4
         else
             call s_print_error('atomic_make_hund','not implemented!')
         endif ! back if ( nband == 5 ) block
     endif ! back if ( icu == 3 ) block

  end subroutine atomic_make_hund

!!>>> atomic_make_umatK: make Coulomb interaction U according to Kanamori
!!>>> parameterized Hamiltonian
  subroutine atomic_make_umatK()
     use constants, only : dp, zero, two, czero

     use control, only : nband, norbs
     use control, only : Uc
     use m_spmat, only : umat

     implicit none

! local varibales
! orbital index
     integer  :: alpha, betta
     integer  :: delta, gamma

! band index and spin index
     integer  :: aband, bband
     integer  :: dband, gband
     integer  :: aspin, bspin
     integer  :: dspin, gspin

! dummy variables
     real(dp) :: dtmp

! Hund's rule matrix
     real(dp) :: hund(nband,nband,3)

! initialize hund to zero
     hund = zero

! build Hund's rule coupling matrix
     call atomic_make_hund(hund)

! initialize umat to zero
     umat = czero

! loop for creation operators
     alphaloop: do alpha=1,norbs-1
         bettaloop: do betta=alpha+1,norbs

! loop for annihilation operators
             gammaloop: do gamma=1,norbs-1
                 deltaloop: do delta=gamma+1,norbs

! get the band and spin indices
                     aband = ( alpha + 1 ) / 2; aspin = mod(alpha,2)
                     bband = ( betta + 1 ) / 2; bspin = mod(betta,2)
                     gband = ( gamma + 1 ) / 2; gspin = mod(gamma,2)
                     dband = ( delta + 1 ) / 2; dspin = mod(delta,2)

                     dtmp = zero

! intraorbital Coulomb interaction
                     if ( ( alpha == gamma ) .and. ( betta == delta ) ) then
                         if ( ( aband == bband ) .and. ( aspin /= bspin ) ) then
                             dtmp = dtmp + Uc
                         endif ! back if ( ( aband == bband ) .and. ( aspin /= bspin ) ) block
                     endif ! back if ( ( alpha == gamma ) .and. ( betta == delta ) ) block

! interorbital Coulomb interaction
                     if ( ( alpha == gamma ) .and. ( betta == delta ) ) then
                         if ( aband /= bband ) then
                             dtmp = dtmp + (Uc - two * hund(aband,bband,1))
                         endif ! back if ( aband /= bband ) block
                     endif ! back if ( ( alpha == gamma ) .and. ( betta == delta ) ) block

! Hund's exchange interaction
                     if ( ( alpha == gamma ) .and. ( betta == delta ) ) then
                         if ( ( aband /= bband ) .and. ( aspin == bspin ) ) then
                             dtmp = dtmp - hund(aband,bband,1)
                         endif ! back if ( ( aband /= bband ) .and. ( aspin == bspin ) ) block
                     endif ! back if ( ( alpha == gamma ) .and. ( betta == delta ) ) block

! spin flip term
                     if ( ( aband == gband ) .and. ( bband == dband ) ) then
                         if ( ( aspin /= gspin ) .and. ( bspin /= dspin ) .and. ( aspin /= bspin ) ) then
                             dtmp = dtmp - hund(aband,bband,2)
                         endif ! back if ( ( aspin /= gspin ) .and. ( bspin /= dspin ) .and. ( aspin /= bspin ) ) block
                     endif ! back if ( ( aband == gband ) .and. ( bband == dband ) ) block

! pair hopping term
                     if ( ( aband == bband ) .and. ( dband == gband ) .and. ( aband /= dband ) ) then
                         if ( ( aspin /= bspin ) .and. ( dspin /= gspin ) .and. ( aspin == gspin ) ) then
                             dtmp = dtmp + hund(aband,gband,3)
                         endif ! back if ( ( aspin /= bspin ) .and. ( dspin /= gspin ) .and. ( aspin == gspin ) ) block
                     endif ! back if ( ( aband == bband ) .and. ( dband == gband ) .and. ( aband /= dband ) ) block

                     umat(alpha,betta,delta,gamma) = dtmp

                 enddo deltaloop ! over delta={gamma+1,norbs} loop
             enddo gammaloop ! over gamma={1,norbs-1} loop
         enddo bettaloop ! over betta={alpha+1,norbs} loop
     enddo alphaloop ! over alpha={1,norbs-1} loop

     return
  end subroutine atomic_make_umatK

!!>>> atomic_make_umatS: make Coulomb interation U, according to
!!>>> Slater-Cordon parameterized Hamiltonian
  subroutine atomic_make_umatS()
     use constants, only : dp, zero, half

     use control, only : nband, norbs
     use control, only : Ud, Jh
     use m_spmat, only : umat

     implicit none

! local variables
! orbital momentum quantum number
     integer  :: l

! loop index
     integer  :: i

! orbital index
     integer  :: alpha, betta
     integer  :: delta, gamma

! band index and spin index
     integer  :: aband, aspin
     integer  :: bband, bspin
     integer  :: dband, dspin
     integer  :: gband, gspin

! dummy variables
     real(dp) :: res

! gaunt coefficients
     real(dp), allocatable :: gaunt(:,:,:)

! Slater-Cordon parameters: F0, F2, F4, and F6
     real(dp), allocatable :: slater_cordon(:)

! allocate memory for slater_cordon and gaunt and then build them
     select case (nband)

         case (5)
             l = 2
             allocate(slater_cordon(0:2*l))
             slater_cordon = zero
             slater_cordon(0) = Ud
             slater_cordon(2) = Jh * 14.0_dp / 1.625_dp
             slater_cordon(4) = 0.625_dp * slater_cordon(2)
             allocate(gaunt(-l:l,-l:l,0:2*l))
             call atomic_make_gaunt5(gaunt)

         case (7)
             l = 3
             allocate(slater_cordon(0:2*l))
             slater_cordon = zero
             slater_cordon(0) = Ud
             slater_cordon(2) = Jh * 6435.0_dp / ( 286.0_dp + ( 195.0_dp * &
                 & 451.0_dp / 675.0_dp ) + ( 250.0_dp * 1001.0_dp / 2025.0_dp ) )
             slater_cordon(4) = 451.0_dp / 675.0_dp * slater_cordon(2)
             slater_cordon(6) = 1001.0_dp / 2025.0_dp * slater_cordon(2)
             allocate(gaunt(-l:l,-l:l,0:2*l))
             call atomic_make_gaunt7(gaunt)

         case default
             call s_print_error('atomic_make_umatS','not implemented for this nband!')

     end select

! make Coulomb interaction U matrix
     do alpha=1,norbs
         do betta=1,norbs
             aband = ( alpha - 1 ) / 2 - l
             bband = ( betta - 1 ) / 2 - l
             aspin = mod(alpha,2)
             bspin = mod(betta,2)

             do gamma=1,norbs
                 do delta=1,norbs
                     dband = ( delta - 1 ) / 2 - l
                     gband = ( gamma - 1 ) / 2 - l
                     dspin = mod(delta,2)
                     gspin = mod(gamma,2)

                     if ( ( aband + bband ) /= ( dband + gband ) ) CYCLE
                     if ( ( aspin /= gspin ) .or. ( bspin /= dspin ) ) CYCLE

                     res = zero
                     do i=0,2*l,2
                         res = res + gaunt(aband,gband,i) * gaunt(dband,bband,i) * slater_cordon(i)
                     enddo ! over i={0,2*l} loop
                     umat(alpha,betta,delta,gamma) = res
                 enddo ! over gamma={1,norbs} loop
             enddo ! over delta={1,norbs} loop

         enddo ! over betta={1,norbs} loop
     enddo ! over alpha={1,norbs} loop
     umat = half * umat

! deallocate memory
     if ( allocated(gaunt) )         deallocate(gaunt)
     if ( allocated(slater_cordon) ) deallocate(slater_cordon)

     return
  end subroutine atomic_make_umatS

!!========================================================================
!!>>> determine spin-orbital coupling matrix                           <<<
!!========================================================================

!>>> atomic_make_smat3: make spin-orbit coupling matrix for 3-band case
  subroutine atomic_make_smat3(smat)
     use constants, only : dp, czero

     implicit none

! external arguments
! SOC matrix
     complex(dp), intent(out) :: smat(6,6)

! local parameters
! \sqrt{2}
     real(dp), parameter :: sqrt2 = sqrt(2.0_dp)

! make SOC on complex orbital basis, the orbital order is:
! |-1,up>, |-1,dn>,
! | 0,up>, | 0,dn>,
! | 1,up>, | 1,dn>
     smat = czero

     smat( 1, 1) = -1.0_dp
     smat( 4, 1) =  sqrt2
     smat( 2, 2) = +1.0_dp
     smat( 6, 3) =  sqrt2
     smat( 1, 4) =  sqrt2
     smat( 5, 5) = +1.0_dp
     smat( 3, 6) =  sqrt2
     smat( 6, 6) = -1.0_dp

     return
  end subroutine atomic_make_smat3

!!>>> atomic_make_smat5: make spin-orbit coupling matrix for 5-band case
  subroutine atomic_make_smat5(smat)
     use constants, only : dp, czero

     implicit none

! external arguments
! SOC matrix
     complex(dp), intent(out) :: smat(10,10)

! local parameters
! \sqrt{6}
     real(dp), parameter :: sqrt6 = sqrt(6.0_dp)

! make SOC on complex orbital basis, the orbital order is:
! |-2,up>, |-2,dn>,
! |-1,up>, |-1,dn>,
! | 0,up>, | 0,dn>,
! | 1,up>, | 1,dn>,
! | 2,up>, | 2,dn>
     smat = czero

     smat( 1, 1) = -2.0_dp
     smat( 4, 1) = +2.0_dp
     smat( 2, 2) = +2.0_dp
     smat( 3, 3) = -1.0_dp
     smat( 6, 3) =  sqrt6
     smat( 1, 4) = +2.0_dp
     smat( 4, 4) = +1.0_dp
     smat( 8, 5) =  sqrt6
     smat( 3, 6) =  sqrt6
     smat( 7, 7) = +1.0_dp
     smat(10, 7) = +2.0_dp
     smat( 5, 8) =  sqrt6
     smat( 8, 8) = -1.0_dp
     smat( 9, 9) = +2.0_dp
     smat( 7,10) = +2.0_dp
     smat(10,10) = -2.0_dp

     return
  end subroutine atomic_make_smat5

!!>>> atomic_make_smat7: make spin-orbit coupling matrix for 7-band case
  subroutine atomic_make_smat7(smat)
     use constants, only : dp, czero

     implicit none

! external arguments
! SOC matrix
     complex(dp), intent(out) :: smat(14,14)

! local parameters
! \sqrt{6}, \sqrt{10}, and \sqrt{12}
     real(dp), parameter :: sqrt6  = sqrt( 6.0_dp)
     real(dp), parameter :: sqrt10 = sqrt(10.0_dp)
     real(dp), parameter :: sqrt12 = sqrt(12.0_dp)

! make SOC on complex orbital basis, the orbital order is:
! |-3,up>, |-3,dn>,
! |-2,up>, |-2,dn>,
! |-1,up>, |-1,dn>,
! | 0,up>, | 0,dn>,
! | 1,up>, | 1,dn>,
! | 2,up>, | 2,dn>
! | 3,up>, | 3,dn>
     smat = czero

     smat( 1, 1) = -3.0_dp
     smat( 4, 1) =  sqrt6
     smat( 2, 2) = +3.0_dp
     smat( 3, 3) = -2.0_dp
     smat( 6, 3) =  sqrt10
     smat( 1, 4) =  sqrt6
     smat( 4, 4) = +2.0_dp
     smat( 5, 5) = -1.0_dp
     smat( 8, 5) =  sqrt12
     smat( 3, 6) =  sqrt10
     smat( 6, 6) = +1.0_dp
     smat(10, 7) =  sqrt12
     smat( 5, 8) =  sqrt12
     smat( 9, 9) = +1.0_dp
     smat(12, 9) =  sqrt10
     smat( 7,10) =  sqrt12
     smat(10,10) = -1.0_dp
     smat(11,11) = +2.0_dp
     smat(14,11) =  sqrt6
     smat( 9,12) =  sqrt10
     smat(12,12) = -2.0_dp
     smat(13,13) = +3.0_dp
     smat(11,14) =  sqrt6
     smat(14,14) = -3.0_dp

     return
  end subroutine atomic_make_smat7

!!========================================================================
!!>>> determine representation transformation matrix                   <<<
!!========================================================================

!!>>> atomic_make_tmat_c2r: make transformation matrix from complex
!!>>> orbital basis to real orbital basis
  subroutine atomic_make_tmat_c2r(tmat_c2r)
     use constants, only : dp, czero, cone, czi

     use control, only : nband, norbs

     implicit none

! external arguments
! the transformation matrix from complex orbitals to real orbitals
     complex(dp), intent(out) :: tmat_c2r(norbs,norbs)

! local parameters
! \sqrt{2}
     real(dp), parameter :: sqrt2 = sqrt(2.0_dp)

     tmat_c2r = czero
     select case (nband)

! the real orbital order (t2g) is:
! dxzup, dxzdn, dyzup, dyzdn, dxyup, dxydn
!
! the corresponding p orbital order is:
! pyup, pydn, pxup, pxdn, pzup, pzdn
!
! the complex orbital |lz,sz> order is
! |-1,up>, |-1,dn>, |0,up>, |0,dn>, |1,up>, |1,dn>
         case (3)
             tmat_c2r( 1, 1) =  czi/sqrt2
             tmat_c2r( 5, 1) =  czi/sqrt2
             tmat_c2r( 2, 2) =  czi/sqrt2
             tmat_c2r( 6, 2) =  czi/sqrt2
             tmat_c2r( 1, 3) =  cone/sqrt2
             tmat_c2r( 5, 3) = -cone/sqrt2
             tmat_c2r( 2, 4) =  cone/sqrt2
             tmat_c2r( 6, 4) = -cone/sqrt2
             tmat_c2r( 3, 5) =  cone
             tmat_c2r( 4, 6) =  cone

! the real orbital order is:
! dz2up, dz2dn, dxzup, dxzdn, dyzup, dyzdn, dx2-y2up, dx2-y2dn, dxyup, dxydn
!
! the complex orbital |lz,sz> order is:
! |-2,up>, |-2,dn>, |-1,up>, |-1,dn>, |0,up>, |0,dn>, |1,up>, |1,dn>, |2,up>, |2,dn>
         case (5)
             tmat_c2r( 5, 1) =  cone
             tmat_c2r( 6, 2) =  cone
             tmat_c2r( 3, 3) =  cone/sqrt2
             tmat_c2r( 7, 3) = -cone/sqrt2
             tmat_c2r( 4, 4) =  cone/sqrt2
             tmat_c2r( 8, 4) = -cone/sqrt2
             tmat_c2r( 3, 5) =  czi/sqrt2
             tmat_c2r( 7, 5) =  czi/sqrt2
             tmat_c2r( 4, 6) =  czi/sqrt2
             tmat_c2r( 8, 6) =  czi/sqrt2
             tmat_c2r( 1, 7) =  cone/sqrt2
             tmat_c2r( 9, 7) =  cone/sqrt2
             tmat_c2r( 2, 8) =  cone/sqrt2
             tmat_c2r(10, 8) =  cone/sqrt2
             tmat_c2r( 1, 9) =  czi/sqrt2
             tmat_c2r( 9, 9) = -czi/sqrt2
             tmat_c2r( 2,10) =  czi/sqrt2
             tmat_c2r(10,10) = -czi/sqrt2

! the real orbital order is:
! fz3up, fz3dn, fxz2up, fxz2dn, fyz2up, fyz2dn, fz(x2-y2)up, fz(x2-y2)dn, fxyzup, fxyzdn,
! fx(x2-3y2)up, fx(x2-3y2)dn, fy(3x2-y2)up, fy(3x2-y2)dn
!
! the complex orbital |lz,sz> order is:
! |-3,up>, |-3,dn>, |-2,up>, |-2,dn>, |-1,up>, |-1,dn>, |0,up>,
! | 0,dn>, | 1,up>, | 1,dn>, | 2,up>, | 2,dn>, | 3,up>, |3,dn>
         case (7)
             tmat_c2r( 7, 1) =  cone
             tmat_c2r( 8, 2) =  cone
             tmat_c2r( 5, 3) =  cone/sqrt2
             tmat_c2r( 9, 3) = -cone/sqrt2
             tmat_c2r( 6, 4) =  cone/sqrt2
             tmat_c2r(10, 4) = -cone/sqrt2
             tmat_c2r( 5, 5) =  czi/sqrt2
             tmat_c2r( 9, 5) =  czi/sqrt2
             tmat_c2r( 6, 6) =  czi/sqrt2
             tmat_c2r(10, 6) =  czi/sqrt2
             tmat_c2r( 3, 7) =  cone/sqrt2
             tmat_c2r(11, 7) =  cone/sqrt2
             tmat_c2r( 4, 8) =  cone/sqrt2
             tmat_c2r(12, 8) =  cone/sqrt2
             tmat_c2r( 3, 9) =  czi/sqrt2
             tmat_c2r(11, 9) = -czi/sqrt2
             tmat_c2r( 4,10) =  czi/sqrt2
             tmat_c2r(12,10) = -czi/sqrt2
             tmat_c2r( 1,11) =  cone/sqrt2
             tmat_c2r(13,11) = -cone/sqrt2
             tmat_c2r( 2,12) =  cone/sqrt2
             tmat_c2r(14,12) = -cone/sqrt2
             tmat_c2r( 1,13) =  czi/sqrt2
             tmat_c2r(13,13) =  czi/sqrt2
             tmat_c2r( 2,14) =  czi/sqrt2
             tmat_c2r(14,14) =  czi/sqrt2

         case default
             call s_print_error('atomic_make_tmat_c2r','not implemented for this nband!')

     end select

     return
  end subroutine atomic_make_tmat_c2r

!!>>> atomic_make_tmat_r2c: make transformation matrix from real orbital
!!>>> basis to complex orbital basis
  subroutine atomic_make_tmat_r2c(tmat_r2c)
     use constants, only : dp, czero

     use control, only : norbs

! external arguments
! the transformation matrix from real orbitals to complex orbitals
     complex(dp), intent(out) :: tmat_r2c(norbs,norbs)

! local variables
! dummy array, transpose of tmat_r2c
     complex(dp) :: tmat_c2r(norbs,norbs)

     tmat_c2r = czero
     call atomic_make_tmat_c2r(tmat_c2r)
     tmat_r2c = transpose( dconjg(tmat_c2r) )

     return
  end subroutine atomic_make_tmat_r2c

!!>>> atomic_make_tmat_c2j: make transformation matrix from complex
!!>>> orbital basis (|lz,sz>) to j2-jz orbital basis (|j2,jz>), i.e.,
!!>>> the Cordon-Gaunt (CG) coefficients
  subroutine atomic_make_tmat_c2j(tmat_c2j)
     use constants, only : dp, czero

     use control, only : nband, norbs

     implicit none

! external arguments
! the transformation matrix from complex orbitals |lz,sz> to |j2,jz>
     complex(dp), intent(out) :: tmat_c2j(norbs,norbs)

     tmat_c2j = czero
     select case (nband)

! the |lz,sz> order is:
! |-1,up>, |-1,dn>, |0,up>, |0,dn>, |1,up>, |1,dn>
!
! the |j2,jz> order is:
! |1/2,-1/2>, |1/2,1/2>,
! |3/2,-3/2>, |3/2,-1/2>, |3/2,1/2>, |3/2,3/2>
         case (3)
             tmat_c2j( 1, 1) = -sqrt(2.0_dp/3.0_dp)
             tmat_c2j( 4, 1) =  sqrt(1.0_dp/3.0_dp)
             tmat_c2j( 3, 2) = -sqrt(1.0_dp/3.0_dp)
             tmat_c2j( 6, 2) =  sqrt(2.0_dp/3.0_dp)
             tmat_c2j( 2, 3) =  1.0_dp
             tmat_c2j( 1, 4) =  sqrt(1.0_dp/3.0_dp)
             tmat_c2j( 4, 4) =  sqrt(2.0_dp/3.0_dp)
             tmat_c2j( 3, 5) =  sqrt(2.0_dp/3.0_dp)
             tmat_c2j( 6, 5) =  sqrt(1.0_dp/3.0_dp)
             tmat_c2j( 5, 6) =  1.0_dp

! the |lz,sz> order is:
! |-2,up>, |-2,dn>, |-1,up>, |-1,dn>, |0,up>, |0,dn>, |1,up>, |1,dn>, |2,up>, |2,dn>
!
! the |j2,jz> order is:
! |3/2,-3/2>, |3/2,-1/2>, |3/2,1/2>, |3/2,3/2>
! |5/2,-5/2>, |5/2,-3/2>, |5/2,-1/2>, |5/2,1/2>, |5/2,3/2>, |5/2,5/2>
         case (5)
             tmat_c2j( 1, 1) = -sqrt(4.0_dp/5.0_dp)
             tmat_c2j( 4, 1) =  sqrt(1.0_dp/5.0_dp)
             tmat_c2j( 3, 2) = -sqrt(3.0_dp/5.0_dp)
             tmat_c2j( 6, 2) =  sqrt(2.0_dp/5.0_dp)
             tmat_c2j( 5, 3) = -sqrt(2.0_dp/5.0_dp)
             tmat_c2j( 8, 3) =  sqrt(3.0_dp/5.0_dp)
             tmat_c2j( 7, 4) = -sqrt(1.0_dp/5.0_dp)
             tmat_c2j(10, 4) =  sqrt(4.0_dp/5.0_dp)
             tmat_c2j( 2, 5) =  1.0_dp
             tmat_c2j( 1, 6) =  sqrt(1.0_dp/5.0_dp)
             tmat_c2j( 4, 6) =  sqrt(4.0_dp/5.0_dp)
             tmat_c2j( 3, 7) =  sqrt(2.0_dp/5.0_dp)
             tmat_c2j( 6, 7) =  sqrt(3.0_dp/5.0_dp)
             tmat_c2j( 5, 8) =  sqrt(3.0_dp/5.0_dp)
             tmat_c2j( 8, 8) =  sqrt(2.0_dp/5.0_dp)
             tmat_c2j( 7, 9) =  sqrt(4.0_dp/5.0_dp)
             tmat_c2j(10, 9) =  sqrt(1.0_dp/5.0_dp)
             tmat_c2j( 9,10) =  1.0_dp

! the |lz,sz> order is:
! |-3,up>, |-3,dn>, |-2,up>, |-2,dn>, |-1,up>, |-1,dn>, |0,up>,
! | 0,dn>, | 1,up>, | 1,dn>, | 2,up>, | 2,dn>, | 3,up>, |3,dn>
!
! the |j2,jz> order is:
! |5/2,-5/2>, |5/2,-3/2>, |5/2,-1/2>, |5/2,1/2>, |5/2,3/2>, |5/2,5/2>
! |7/2,-7/2>, |7/2,-5/2>, |7/2,-3/2>, |7/2,-1/2>, |7/2,1/2>, |7/2,3/2>, |7/2,5/2>, |7/2,7/2>
         case (7)
             tmat_c2j( 1, 1) = -sqrt(6.0_dp/7.0_dp)
             tmat_c2j( 4, 1) =  sqrt(1.0_dp/7.0_dp)
             tmat_c2j( 3, 2) = -sqrt(5.0_dp/7.0_dp)
             tmat_c2j( 6, 2) =  sqrt(2.0_dp/7.0_dp)
             tmat_c2j( 5, 3) = -sqrt(4.0_dp/7.0_dp)
             tmat_c2j( 8, 3) =  sqrt(3.0_dp/7.0_dp)
             tmat_c2j( 7, 4) = -sqrt(3.0_dp/7.0_dp)
             tmat_c2j(10, 4) =  sqrt(4.0_dp/7.0_dp)
             tmat_c2j( 9, 5) = -sqrt(2.0_dp/7.0_dp)
             tmat_c2j(12, 5) =  sqrt(5.0_dp/7.0_dp)
             tmat_c2j(11, 6) = -sqrt(1.0_dp/7.0_dp)
             tmat_c2j(14, 6) =  sqrt(6.0_dp/7.0_dp)
             tmat_c2j( 2, 7) =  1.0_dp
             tmat_c2j( 1, 8) =  sqrt(1.0_dp/7.0_dp)
             tmat_c2j( 4, 8) =  sqrt(6.0_dp/7.0_dp)
             tmat_c2j( 3, 9) =  sqrt(2.0_dp/7.0_dp)
             tmat_c2j( 6, 9) =  sqrt(5.0_dp/7.0_dp)
             tmat_c2j( 5,10) =  sqrt(3.0_dp/7.0_dp)
             tmat_c2j( 8,10) =  sqrt(4.0_dp/7.0_dp)
             tmat_c2j( 7,11) =  sqrt(4.0_dp/7.0_dp)
             tmat_c2j(10,11) =  sqrt(3.0_dp/7.0_dp)
             tmat_c2j( 9,12) =  sqrt(5.0_dp/7.0_dp)
             tmat_c2j(12,12) =  sqrt(2.0_dp/7.0_dp)
             tmat_c2j(11,13) =  sqrt(6.0_dp/7.0_dp)
             tmat_c2j(14,13) =  sqrt(1.0_dp/7.0_dp)
             tmat_c2j(13,14) =  1.0_dp

         case default
             call s_print_error('atomic_make_tmat_c2j','not implemented for this nband!')

     end select

     return
  end subroutine atomic_make_tmat_c2j

!!========================================================================
!!>>> perform representation transformation                            <<<
!!========================================================================

!!>>> atomic_tran_fmat: rotate F-matrix (fmat) from Fock basis to eigen-
!!>>> states basis
  subroutine atomic_tran_fmat(ndimx, ndimy, amat, bmat, cmat)
     use constants, only: dp, zero, one

     implicit none

! external arguments
! x dimension of matrix
     integer, intent(in)  :: ndimx

! y dimension of matrix
     integer, intent(in)  :: ndimy

! left transformation matrix
     real(dp), intent(in) :: amat(ndimx,ndimx)

! right transformation matrix
     real(dp), intent(in) :: cmat(ndimy,ndimy)

! F-matrix
     real(dp), intent(inout) :: bmat(ndimx,ndimy)

! local variables
! dummy array
     real(dp), allocatable :: tmp_mat(:,:)

! allocate memory
     allocate(tmp_mat(ndimx,ndimy))

     tmp_mat = zero
     call dgemm('N', 'N', ndimx, ndimy, ndimy, &
                             one, bmat, ndimx, &
                                  cmat, ndimy, &
                         zero, tmp_mat, ndimx  )

     call dgemm('T', 'N', ndimx, ndimy, ndimx, &
                             one, amat, ndimx, &
                               tmp_mat, ndimx, &
                            zero, bmat, ndimx  )

! deallocate memory
     deallocate(tmp_mat)

     return
  end subroutine atomic_tran_fmat

!!>>> atomic_tran_umat: transform Coulomb interaction U tensor from one
!!>>> representation to another representation
  subroutine atomic_tran_umat(amtrx, umat, umat_t)
     use constants, only : dp, czero, epst

     use control, only : norbs

     implicit none

! external arguments
! transformation matrix from orginal basis to natural basis
     complex(dp), intent(in)  :: amtrx(norbs,norbs)

! coefficents matrix for general interaction U in orginal basis
     complex(dp), intent(in)  :: umat(norbs,norbs,norbs,norbs)

! coefficents matrix for general interaction U in natural basis
     complex(dp), intent(out) :: umat_t(norbs,norbs,norbs,norbs)

! local variables
! loop index over orbitals in orginal single particle basis
     integer :: alpha1, alpha2
     integer :: alpha3, alpha4
     integer :: sigma1, sigma2
     integer :: sigma3, sigma4

! auxiliary complex(dp) variables
     complex(dp) :: ctmp

! initialize umat_t to be zero
     umat_t = czero

     sigma1loop: do sigma1=1,norbs
         sigma2loop: do sigma2=1,norbs
             sigma3loop: do sigma3=1,norbs
                 sigma4loop: do sigma4=1,norbs
                     ctmp = czero

                     alpha1loop: do alpha1=1,norbs
                         alpha2loop: do alpha2=1,norbs
                             alpha3loop: do alpha3=1,norbs
                                 alpha4loop: do alpha4=1,norbs
                                     if ( abs( umat(alpha1,alpha2,alpha3,alpha4) ) < epst ) CYCLE
                                     ctmp = ctmp + umat(alpha1,alpha2,alpha3,alpha4)                   &
                                                 * conjg(amtrx(alpha1,sigma1)) * amtrx(alpha3,sigma3)  &
                                                 * conjg(amtrx(alpha2,sigma2)) * amtrx(alpha4,sigma4)
                                 enddo alpha4loop ! over alpha4={1,norbs} loop
                             enddo alpha3loop ! over alpha3={1,norbs} loop
                         enddo alpha2loop ! over alpha2={1,norbs} loop
                     enddo alpha1loop ! over alpha1={1,norbs} loop

                     umat_t(sigma1,sigma2,sigma3,sigma4) = ctmp
                 enddo sigma4loop ! over sigma4={1,norbs} loop
             enddo sigma3loop ! over sigma3={1,norbs} loop
         enddo sigma2loop ! over sigma2={1,norbs} loop
     enddo sigma1loop ! over sigma1={1,norbs} loop

     return
  end subroutine atomic_tran_umat

!!>>> atomic_tran_repr_cmpl: transformation from one representation
!!>>> to another representation, complex version
  subroutine atomic_tran_repr_cmpl(ndim, amat, tmat)
     use constants, only : dp, cone, czero

     implicit none

! external arguments
! size of the matrix
     integer, intent(in) :: ndim

! transformation matrix
     complex(dp), intent(in) :: tmat(ndim,ndim)

! physical quantities
     complex(dp), intent(inout) :: amat(ndim,ndim)

! local variables
! dummy matrix
     complex(dp) :: tmp_mat(ndim,ndim)

     call zgemm('N', 'N', ndim, ndim, ndim, &
                          cone, amat, ndim, &
                                tmat, ndim, &
                      czero, tmp_mat, ndim  )

     call zgemm('C', 'N', ndim, ndim, ndim, &
                          cone, tmat, ndim, &
                             tmp_mat, ndim, &
                         czero, amat, ndim  )

     return
  end subroutine atomic_tran_repr_cmpl

!!>>> atomic_tran_repr_real: transformation from one representation to
!!>>> another representation, real version
  subroutine atomic_tran_repr_real(ndim, amat, tmat)
     use constants, only : dp, zero, one

     implicit none

! external arguments
! size of the matrix
     integer, intent(in) :: ndim

! transformation matrix
     real(dp), intent(in) :: tmat(ndim,ndim)

! physical quantities
     real(dp), intent(inout) :: amat(ndim,ndim)

! local variables
! dummy matrix
     real(dp) :: tmp_mat(ndim,ndim)

     call dgemm('N', 'N', ndim, ndim, ndim, &
                           one, amat, ndim, &
                                tmat, ndim, &
                       zero, tmp_mat, ndim  )

     call dgemm('T', 'N', ndim, ndim, ndim, &
                           one, tmat, ndim, &
                             tmp_mat, ndim, &
                          zero, amat, ndim  )

     return
  end subroutine atomic_tran_repr_real



!!>>> atomic_2natural_case1: make natural basis for no crystal field or
!!>>> diagonal crystal field, without spin-orbital coupling
  subroutine atomic_2natural_case1()
     use control, only : norbs
     use control, only : mune
     use m_spmat, only : cmat, emat, tmat

     implicit none

! local variables
! loop index
     integer :: i

! set emat
! since smat is zero, so emat is equal to cmat
     emat = cmat

! add chemical potential to eimpmat
     do i=1,norbs
         emat(i,i) = emat(i,i) + mune
     enddo ! over i={1,norbs} loop

! for this case, the natural basis is the real orbital basis
! so, the tmat is a unity matrix
     call s_identity_z(norbs, tmat)

     return
  end subroutine atomic_2natural_case1

!!>>> atomic_2natural_case2: make natural basis for non-diagonal
!!>>> crystal field without spin-orbital coupling
  subroutine atomic_2natural_case2()
     use constants, only : dp

     use control, only : nband, norbs
     use control, only : mune
     use m_spmat, only : cmat, emat, tmat

     implicit none

! local variables
! loop index
     integer  :: i
     integer  :: j

! eigenvalue
     real(dp) :: eigval(nband)

! eigen vector
     real(dp) :: eigvec(nband,nband)

! emat matrix for no spin freedom
     complex(dp) :: emat_nospin(nband,nband)

! tmat for no spin freedom
     complex(dp) :: tmat_nospin(nband,nband)

! set emat to crystal field
! since smat is zero, so emat is equal to cmat
     emat = cmat

! get emat for no spin freedom
     do i=1,nband
         do j=1,nband
             emat_nospin(j,i) = emat(2*j-1,2*i-1)
         enddo ! over j={1,nband} loop
     enddo ! over i={1,nband}

! diagonalize emat_nospin to get natural basis
     call s_eig_sy(nband, nband, real(emat_nospin), eigval, eigvec)

! get diagonal emat for no spin freedom
     call s_diag_z(nband, dcmplx(eigval), emat_nospin)

! get tmat for no spin freedom
     tmat_nospin = dcmplx(eigvec)

! build emat and tmat with spin freedom
     do i=1,nband
         do j=1,nband
             emat(2*j-1,2*i-1) = emat_nospin(j,i)
             emat(2*j,2*i)     = emat_nospin(j,i)
             tmat(2*j-1,2*i-1) = tmat_nospin(j,i)
             tmat(2*j,2*i)     = tmat_nospin(j,i)
         enddo ! over j={1,nband} loop
     enddo ! over i={1,nband} loop

! add chemical potential to emat
     do i=1,norbs
         emat(i,i) = emat(i,i) + mune
     enddo ! over i={1,norbs} loop

     return
  end subroutine atomic_2natural_case2

!!>>> atomic_2natural_case3: make natural basis for the case without
!!>>> crystal field and with spin-orbital coupling
!!>>> for this special case, the natural basis is |j2,jz>
  subroutine atomic_2natural_case3()
     use constants, only : dp

     use control, only : norbs
     use control, only : mune
     use m_spmat, only : emat, smat, tmat

     implicit none

! local variables
! loop inex
     integer :: i

! transformation matrix from complex orbital basis to |j2,jz> basis
     complex(dp) :: tmat_c2j(norbs,norbs)

! set emat
! since cmat is zero, so emat is equal to smat
     emat = smat

! evaluate transformation matrix tmat_c2j
     call atomic_make_tmat_c2j(tmat_c2j)

! for this case, the transformation matrix is from complex orbital basis
! to natural basis (|j2,jz> basis)
     tmat = tmat_c2j

! transform emat to natural basis
     call atomic_tran_repr_cmpl(norbs, emat, tmat)

! add chemical potential to emat
     do i=1,norbs
         emat(i,i) = emat(i,i) + mune
     enddo ! over i={1,norbs} loop

     return
  end subroutine atomic_2natural_case3

!!>>> atomic_2natural_case4: make natural basis for the case with
!!>>> crystal field and with spin-orbital coupling
  subroutine atomic_2natural_case4()
     use constants, only : dp, eps6

     use control, only : norbs
     use control, only : mune
     use m_spmat, only : cmat, smat, emat, tmat

     implicit none

! local variables
! loop index
     integer  :: i

! eigenvalue
     real(dp) :: eigval(norbs)

! eigenvector
     real(dp) :: eigvec(norbs,norbs)

! transformation matrix from real orbital basis to complex orbital basis
     complex(dp) :: tmat_r2c(norbs,norbs)

! transformation matrix from complex orbital basis to natural basis
     complex(dp) :: tmat_c2n(norbs,norbs)

! build tmat_r2c
     call atomic_make_tmat_r2c(tmat_r2c)

! transfrom crystal field (cmat) to complex orbital basis
     call atomic_tran_repr_cmpl(norbs, cmat, tmat_r2c)

! check whether cmat is real, if not, we cann't make natural basis
     if ( any( abs( aimag(cmat) ) > eps6 ) ) then
         call s_print_error('atomic_2natural_case4','crystal field on complex orbital basis should be real!')
     endif ! back if ( any( abs( aimag(cmat) ) > eps6 ) ) block

! set emat: CF + SOC
     emat = smat + cmat

! diagonalize real(emat)
     call s_eig_sy(norbs, norbs, real(emat), eigval, eigvec)

! get the transformation matrix from complex orbital basis to natural basis
     tmat_c2n = eigvec
     tmat = tmat_c2n

! transform emat from complex orbital basis to natural basis
     call atomic_tran_repr_cmpl(norbs, emat, tmat_c2n)

! add chemical poential to emat
     do i=1,norbs
         emat(i,i) = emat(i,i) + mune
     enddo ! over i={1,norbs} loop

     return
  end subroutine atomic_2natural_case4
