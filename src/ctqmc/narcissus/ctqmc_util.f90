!!!-----------------------------------------------------------------------
!!! project : narcissus
!!! program : ctqmc_four_htau
!!!           ctqmc_four_hybf
!!!           ctqmc_eval_htau
!!!           ctqmc_eval_hsed
!!!           ctqmc_eval_ktau
!!!           ctqmc_eval_ksed
!!!           ctqmc_make_uumat
!!!           ctqmc_make_state
!!!           ctqmc_make_shift
!!!           ctqmc_eval_shift
!!! source  : ctqmc_util.f90
!!! type    : functions & subroutines
!!! author  : li huang (email:lihuang.dmft@gmail.com)
!!! history : 10/01/2008 by li huang (created)
!!!           01/28/2017 by li huang (last modified)
!!! purpose : to provide utility functions and subroutines for hybridization
!!!           expansion version continuous time quantum Monte Carlo (CTQMC)
!!!           quantum impurity solver
!!! status  : unstable
!!! comment :
!!!-----------------------------------------------------------------------

!!========================================================================
!!>>> cubic spline interpolation                                       <<<
!!========================================================================

!! To provide cubic spline subroutines and wrapper functions to interpolate
!! the hybridization function in imaginary-time axis.

!!>>> ctqmc_make_htau: evaluate the matrix elements for mmat matrix using
!!>>> cubic spline interpolation
  function ctqmc_make_htau(flvr, dtau) result(val)
     use constants, only : dp

     use control, only : ntime
     use context, only : tmesh
     use context, only : htau, hsed

     implicit none

! external arguments
! current flavor channel
     integer, intent(in)  :: flvr

! delta imaginary time
     real(dp), intent(in) :: dtau

! external functions
! internal interpolation engine
     procedure( real(dp) ) :: s_spl_funct

! local variables
! return value
     real(dp) :: val

     val = s_spl_funct(ntime, tmesh, htau(:, flvr, flvr), hsed(:, flvr, flvr), dtau)

     return
  end function ctqmc_make_htau

!!>>> ctqmc_make_hsed: calculate the second order derivates of hybridization
!!>>> function on imaginary time space
  subroutine ctqmc_make_hsed(tmesh, htau, hsed)
     use constants, only : dp, zero

     use control, only : norbs
     use control, only : ntime
     use control, only : beta

     implicit none

! external arguments
! imaginary time axis
     real(dp), intent(in)  :: tmesh(ntime)

! hybridization function on imaginary time axis
     real(dp), intent(in)  :: htau(ntime,norbs,norbs)

! second order derivates of hybridization function
     real(dp), intent(out) :: hsed(ntime,norbs,norbs)

! local variables
! loop index
     integer  :: i
     integer  :: j

! first derivate at start point
     real(dp) :: startu

! first derivate at end   point
     real(dp) :: startd

! \delta \tau
     real(dp) :: deltau

! second-order derivates
     real(dp) :: d2y(ntime)

! calculate deltau
     deltau = beta / real(ntime - 1)

! initialize hsed
     hsed = zero

! calculate it
     do j=1,norbs
         do i=1,norbs

! calculate first-order derivate of \Delta(0): startu
             startu = (-25.0_dp*htau(1,       i, j) +                    &
                        48.0_dp*htau(2,       i, j) -                    &
                        36.0_dp*htau(3,       i, j) +                    &
                        16.0_dp*htau(4,       i, j) -                    &
                         3.0_dp*htau(5,       i, j)) / 12.0_dp / deltau

! calculate first-order derivate of \Delta(\beta): startd
             startd = ( 25.0_dp*htau(ntime-0, i, j) -                    &
                        48.0_dp*htau(ntime-1, i, j) +                    &
                        36.0_dp*htau(ntime-2, i, j) -                    &
                        16.0_dp*htau(ntime-3, i, j) +                    &
                         3.0_dp*htau(ntime-4, i, j)) / 12.0_dp / deltau

! reinitialize d2y to zero
             d2y = zero

! call the service layer
             call s_spl_deriv2(ntime, tmesh, htau(:,i,j), startu, startd, d2y)

! copy the results to hsed
             hsed(:,i,j) = d2y

         enddo ! over i={1,norbs} loop
     enddo ! over j={1,norbs} loop

     return
  end subroutine ctqmc_make_hsed

!! To provide cubic spline subroutines and wrapper functions to interpolate
!! the retarded interaction integrated function (i.e., screening function)
!! in imaginary-time axis.

!!>>> ctqmc_make_ktau: evaluate the intermediate elements for K(\tau)
!!>>> using cubic spline interpolation
!!>>> note: this function can be used to interpolate K'(\tau) as well
  function ctqmc_make_ktau(mode, dtau) result(val)
     use constants, only : dp

     use control, only : ntime
     use context, only : tmesh
     use context, only : ktau, ksed, ptau, psed

     implicit none

! external arguments
! order for derivates
! if mode = 1, K(\tau), ktau is considered
! if mode = 2, K'(\tau), ptau is considered
     integer, intent(in)  :: mode

! current imaginary time
     real(dp), intent(in) :: dtau

! external functions
! internal interpolation engine
     procedure( real(dp) ) :: s_spl_funct

! local variables
! return value
     real(dp) :: val

! using cubic spline interpolation for K(\tau)
     if ( mode == 1 ) then
         val = s_spl_funct(ntime, tmesh, ktau, ksed, dtau)
! using cubic spline interpolation for K'(\tau)
     else
         val = s_spl_funct(ntime, tmesh, ptau, psed, dtau)
     endif ! back if ( mode == 1 ) block

     return
  end function ctqmc_make_ktau

!!>>> ctqmc_make_ksed: calculate the second order derivates of screening
!!>>> function K(\tau) on imaginary time space
!!>>> note: this subroutine can be used to calculate the second order
!!>>> derivates of K'(\tau) as well. What you have to do is to transfer
!!>>> ptau and psed to this subroutine
  subroutine ctqmc_make_ksed(tmesh, ktau, ksed)
     use constants, only : dp, zero

     use control, only : ntime
     use control, only : beta

     implicit none

! external arguments
! imaginary time axis
     real(dp), intent(in)  :: tmesh(ntime)

! screening function on imaginary time axis
     real(dp), intent(in)  :: ktau(ntime)

! second order derivates of screening function
     real(dp), intent(out) :: ksed(ntime)

! local variables
! first derivate at start point
     real(dp) :: startu

! first derivate at end   point
     real(dp) :: startd

! \delta \tau
     real(dp) :: deltau

! calculate deltau
     deltau = beta / real(ntime - 1)

! initialize ksed
     ksed = zero

! calculate it
! calculate first-order derivate of K(0): startu
     startu = (-25.0_dp*ktau(1      ) + &
                48.0_dp*ktau(2      ) - &
                36.0_dp*ktau(3      ) + &
                16.0_dp*ktau(4      ) - &
                 3.0_dp*ktau(5      )) / 12.0_dp / deltau

! calculate first-order derivate of K(\beta): startd
     startd = ( 25.0_dp*ktau(ntime-0) - &
                48.0_dp*ktau(ntime-1) + &
                36.0_dp*ktau(ntime-2) - &
                16.0_dp*ktau(ntime-3) + &
                 3.0_dp*ktau(ntime-4)) / 12.0_dp / deltau

! call the service layer
     call s_spl_deriv2(ntime, tmesh, ktau, startu, startd, ksed)

     return
  end subroutine ctqmc_make_ksed

!!========================================================================
!!>>> fast fourier transformation                                      <<<
!!========================================================================

!! Here are forward and backward fourier transformation subroutines for
!! hybridization function. Nominally, the following subroutines are only
!! suitable for the hybridization functions, but in principle, we can also
!! apply them to the impurity green's function and bath weiss's function.

!!>>> ctqmc_four_htau: fourier htau to hybf, from imaginary time to
!!>>> matsubara frequency
  subroutine ctqmc_four_htau(htau, hybf)
     use constants, only : dp, zero, czero

     use control, only : norbs
     use control, only : mfreq
     use control, only : ntime
     use context, only : tmesh, rmesh

     implicit none

! external arguments
! hybridization function on imaginary time axis
     real(dp), intent(in) :: htau(ntime,norbs,norbs)

! hybridization function on matsubara frequency axis
     complex(dp), intent(out) :: hybf(mfreq,norbs,norbs)

! local variables
! loop index over orbitals
     integer  :: i
     integer  :: j

! dummy arrays
     real(dp) :: raux(ntime)
     complex(dp) :: caux(mfreq)

! initialize them
     raux = zero
     caux = czero

     do i=1,norbs
         do j=1,norbs

! copy the imaginary-time data to raux
             raux = htau(:,j,i)

! call the service layer
             call s_fft_forward(ntime, tmesh, raux, mfreq, rmesh, caux)

! copy the matsubara frequency data to hybf
             hybf(:,j,i) = caux

         enddo ! over j={1,norbs} loop
     enddo ! over i={1,norbs} loop

     return
  end subroutine ctqmc_four_htau

!!>>> ctqmc_four_hybf: fourier hybf to htau, from matsubara frequency to
!!>>> imaginary time
  subroutine ctqmc_four_hybf(hybf, htau)
     use constants, only : dp, zero, czero, eps6

     use control, only : norbs
     use control, only : mfreq
     use control, only : ntime
     use control, only : beta
     use context, only : tmesh, rmesh

     implicit none

! external arguments
! hybridization function on imaginary time axis
     real(dp), intent(out) :: htau(ntime,norbs,norbs)

! hybridization function on matsubara frequency axis
     complex(dp), intent(in) :: hybf(mfreq,norbs,norbs)

! local variables
! loop index over orbitals
     integer  :: i
     integer  :: j

! used to determine the bottom region of hybridiaztion function
     integer  :: start
     integer  :: last

! dummy arrays
     real(dp) :: raux(ntime)
     complex(dp) :: caux(mfreq)

! initialize them
     raux = zero
     caux = czero

     do i=1,norbs
         do j=1,norbs

! copy matsubara frequency data to caux
             caux = hybf(:,j,i)

! call the service layer
             call s_fft_backward(mfreq, rmesh, caux, ntime, tmesh, raux, beta)

! copy imaginary time data to htau
             htau(:,j,i) = raux

         enddo ! over j={1,norbs} loop
     enddo ! over i={1,norbs} loop

! checks for diagonal htau to be causal. htau should be concave. hence,
! if it becomes very small at two points, it should remain zero in all
! points between the two points. this is very important in insulators,
! because htau can overshoot to positive values multiple times and kinks
! can be trapped in the range between the two crossing points, where htau
! is causal, but should be zero.
     start = 1
     last = 1
     do i=1,norbs
         do j=1,ntime    ! search forward
             if ( htau(j,i,i) > -eps6 ) then
                 start = j; EXIT
             endif ! back if ( htau(j,i,i) > -eps6 ) block
         enddo ! over j={1,ntime} loop

         do j=ntime,1,-1 ! search backward
             if ( htau(j,i,i) > -eps6 ) then
                 last = j; EXIT
             endif ! back if ( htau(j,i,i) > -eps6 ) block
         enddo ! over j={ntime,1,-1} loop

!-------------------------------------------------------------------------
!<         if ( start > 1 .and. last > 1 ) then
!<             do j=start,last
!<                 htau(j,i,i) = -eps6
!<             enddo ! over j={start,last} loop
!<         endif ! back if ( start > 1 .and. last > 1 ) block
!-------------------------------------------------------------------------
     enddo ! over i={1,norbs} loop

! enforce hybridization function less than zero to ensure the causality
     do i=1,norbs
         do j=1,ntime
             if ( htau(j,i,i) > zero ) htau(j,i,i) = -eps6
         enddo ! over j={1,ntime} loop
     enddo ! over i={1,norbs} loop

     return
  end subroutine ctqmc_four_hybf

!!========================================================================
!!>>> Coulomb interaction matrix                                       <<<
!!========================================================================

!! note: do not support spin-flip and pair-hopping term so far.
!!
!! note: only Uc and Jz are need, the other Coulomb interaction parameters
!! are used as backup

!!>>> ctqmc_make_uumat: to build general U interaction matrix: uumat, using
!!>>> my own style
  subroutine ctqmc_make_uumat(uumat)
     use constants, only : dp, zero

     use control, only : nband, norbs
     use control, only : Uc, Jz

     implicit none

! external arguments
! Coulomb interaction matrix
     real(dp), intent(out) :: uumat(norbs,norbs)

! local variables
! loop index
     integer  :: i
     integer  :: j
     integer  :: k
     integer  :: m

! dummy u vector
     real(dp) :: ut(nband*(norbs-1))

! initialize it
     uumat = zero

! calculate it
     k = 0
     do i=1,norbs-1
         do j=i+1,norbs
             k = k + 1
             if ( i <= nband .and. j > nband ) then
                 m = j - nband
                 if ( m == i ) then
                     ut(k) = Uc
                 else
                     ut(k) = Uc - 2.0_dp * Jz
                 endif ! back if ( m == i ) block
             else
                 ut(k) = Uc - 3.0_dp * Jz
             endif ! back if ( i <= nband .and. j > nband ) block

             uumat(i,j) = ut(k)
             uumat(j,i) = ut(k)
         enddo ! over j={i+1,norbs} loop
     enddo ! over i={1,norbs-1} loop

     return
  end subroutine ctqmc_make_uumat

!!========================================================================
!!>>> atomic state converter                                           <<<
!!========================================================================

!!>>> ctqmc_make_state: convert current atomic state array into a decimal
!!>>> number (state index)
  subroutine ctqmc_make_state(norbs, pstat, state)
     implicit none

! external arguments
! index of atomic state
     integer, intent(out) :: pstat

! number of orbitals
     integer, intent(in)  :: norbs

! atomic state array
     integer, intent(in)  :: state(norbs)

! local variables
! loop index
     integer :: i

! init pstat
     pstat = 1

! evaluate pstat, for example, 0101 = 0*2^0 + 1*2^1 + 0*2^2 + 1*2^3 = 10
     do i=1,norbs
         if ( state(i) > 0 ) pstat = pstat + ishft(1, i-1)
     enddo ! over i={1,norbs} loop

     return
  end subroutine ctqmc_make_state

!!========================================================================
!!>>> auxiliary subroutines for retarded interaction                   <<<
!!========================================================================

!!>>> ctqmc_make_shift: to shift the Coulomb interaction matrix and the
!!>>> chemical potential if retarded interaction is considered
  subroutine ctqmc_make_shift(uumat, ssign)
     use constants, only : dp, two

     use control, only : norbs
     use control, only : mune

     implicit none

! external arguments
! Coulomb interaction matrix
     real(dp), intent(inout) :: uumat(norbs,norbs)

! sign for the shift, it should be 1.0_dp or -1.0_dp
     real(dp), intent(in)    :: ssign

! local variables
! loop index
     integer  :: i
     integer  :: j

! Coulomb interaction shift introduced by dynamical screening effect
     real(dp) :: shift

! evaluate the shift at first
     shift = 0.0_dp; call ctqmc_eval_shift(shift)

! multiple the shift with sign
     call s_assert( abs(ssign) == 1.0_dp ); shift = shift * ssign

! shift the Coulomb interaction matrix (skip the diagonal elements)
     do i=1,norbs-1
         do j=i+1,norbs
             uumat(i,j) = uumat(i,j) - shift
             uumat(j,i) = uumat(j,i) - shift
         enddo ! over j={i+1,norbs} loop
     enddo ! over i={1,norbs-1} loop

! shift chemical potential as a byproduct
     mune = mune - shift / two

     return
  end subroutine ctqmc_make_shift

!!>>> ctqmc_eval_shift: evaluate the shift for the Coulomb interaction and
!!>>> the chemical potential. in fact, shift = 2 K'(\tau = 0)
  subroutine ctqmc_eval_shift(shift)
     use constants, only : dp, zero, two

     use control, only : isscr
     use control, only : lc, wc
     use context, only : ptau

     implicit none

! external arguments
! the shift value for U and mune
     real(dp), intent(out) :: shift

! evaluate Coulomb interaction shift
     select case ( isscr )

         case (1)
             shift = zero               ! normal model, recover azalea

         case (2)
             shift = two * lc * lc / wc ! holstein-hubbard model

         case (3)
             shift = two * lc * lc / wc ! plasmon pole model

         case (4)
             shift = two * lc * wc      ! ohmic model

         case (99)
             shift = lc                 ! realistic materials
             call s_assert( shift /= ptau(1) )

     end select

     return
  end subroutine ctqmc_eval_shift
