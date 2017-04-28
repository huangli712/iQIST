!!!-----------------------------------------------------------------------
!!! project : narcissus
!!! program : ctqmc_four_htau
!!!           ctqmc_four_hybf <<<---
!!!           ctqmc_eval_htau
!!!           ctqmc_eval_hsed <<<---
!!!           ctqmc_eval_ktau
!!!           ctqmc_eval_ksed <<<---
!!!           ctqmc_symm_nmat
!!!           ctqmc_symm_gtau
!!!           ctqmc_symm_grnf <<<---
!!!           ctqmc_smth_sigf <<<---
!!!           ctqmc_make_umat <<<---
!!!           ctqmc_make_fock <<<---
!!!           ctqmc_make_lift
!!!           ctqmc_prep_lift <<<---
!!!           ctqmc_make_gtau
!!!           ctqmc_make_ftau <<<---
!!!           ctqmc_make_iret
!!!           ctqmc_make_pref <<<---
!!!           ctqmc_make_prod <<<---
!!!           ctqmc_make_hub1
!!!           ctqmc_make_hub2 <<<---
!!! source  : ctqmc_util.f90
!!! type    : functions & subroutines
!!! author  : li huang (email:lihuang.dmft@gmail.com)
!!! history : 10/01/2008 by li huang (created)
!!!           04/28/2017 by li huang (last modified)
!!! purpose : provide utility functions and subroutines for hybridization
!!!           expansion version continuous time quantum Monte Carlo (CTQMC)
!!!           quantum impurity solver
!!! status  : unstable
!!! comment :
!!!-----------------------------------------------------------------------

!!========================================================================
!!>>> fast fourier transformation                                      <<<
!!========================================================================

!!
!! note:
!!
!! here are forward and backward fourier transformation subroutines for
!! hybridization function. nominally, the following subroutines are only
!! suitable for the hybridization functions, but in principle, we can also
!! apply them to the impurity green's function and bath weiss's function
!!

!!
!! @sub ctqmc_four_htau
!!
!! fourier htau to hybf, from imaginary time to matsubara frequency
!!
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

! copy the imaginary time data to raux
             raux = htau(:,j,i)

! call the service layer
             call s_fft_forward(ntime, tmesh, raux, mfreq, rmesh, caux)

! copy the matsubara frequency data to hybf
             hybf(:,j,i) = caux

         enddo ! over j={1,norbs} loop
     enddo ! over i={1,norbs} loop

     return
  end subroutine ctqmc_four_htau

!!
!! @sub ctqmc_four_hybf
!!
!! fourier hybf to htau, from matsubara frequency to imaginary time
!!
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
!!>>> cubic spline interpolation 1                                     <<<
!!========================================================================

!!
!! note:
!!
!! to provide cubic spline subroutine and wrapper function to interpolate
!! the hybridization function in imaginary time axis
!!

!!
!! @fun ctqmc_eval_htau
!!
!! evaluate the matrix elements for mmat matrix using cubic spline
!! interpolation method
!!
  function ctqmc_eval_htau(flvr, dtau) result(val)
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
  end function ctqmc_eval_htau

!!
!! @sub ctqmc_eval_hsed
!!
!! calculate the second order derivates of hybridization function on
!! imaginary time space
!!
  subroutine ctqmc_eval_hsed(tmesh, htau, hsed)
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
  end subroutine ctqmc_eval_hsed

!!========================================================================
!!>>> cubic spline interpolation 2                                     <<<
!!========================================================================

!!
!! note:
!!
!! to provide cubic spline subroutine and wrapper function to interpolate
!! the retarded interaction integrated function (i.e., screening function)
!! in imaginary time axis
!!

!!
!! @fun ctqmc_eval_ktau
!!
!! evaluate the intermediate elements for K(\tau) using cubic spline
!! interpolation. this function can be used to interpolate K'(\tau) too
!!
  function ctqmc_eval_ktau(mode, dtau) result(val)
     use constants, only : dp

     use control, only : ntime
     use context, only : tmesh
     use context, only : ktau, ksed
     use context, only : ptau, psed

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
  end function ctqmc_eval_ktau

!!
!! @sub ctqmc_eval_ksed
!!
!! calculate the second order derivates of screening function K(\tau)
!! on imaginary time space. this subroutine can be used to calculate
!! the second order derivates of K'(\tau) as well. what you have to do
!! is to transfer ptau and psed to this subroutine
!!
  subroutine ctqmc_eval_ksed(tmesh, ktau, ksed)
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
  end subroutine ctqmc_eval_ksed

!!========================================================================
!!>>> symmetrize physical observables                                  <<<
!!========================================================================

!!>>> ctqmc_symm_nmat: symmetrize the nmat according to symm vector
  subroutine ctqmc_symm_nmat(symm, nmat)
     use constants, only : dp, zero, two

     use control, only : issun, isspn
     use control, only : nband, norbs

     implicit none

! external arguments
! symmetry vector
     integer, intent(in) :: symm(norbs)

! occupation number
     real(dp), intent(inout) :: nmat(norbs)

! local variables
! loop index over bands
     integer  :: ibnd
     integer  :: jbnd

! dummy variables
     real(dp) :: raux

! histogram vector
! note: it is NOT the global one
     integer  :: hist(norbs)

! build histogram
     hist = 0
     do ibnd=1,norbs
         hist(symm(ibnd)) = hist(symm(ibnd)) + 1
     enddo ! over ibnd={1,norbs} loop

! perform symmetrization for those orbitals which symm index are identity
     if ( issun == 2 ) then
         do ibnd=1,norbs
             if ( hist(ibnd) > 0 ) then         ! need to enforce symmetry
                 raux = zero

                 do jbnd=1,norbs                ! gather the data
                     if ( symm(jbnd) == ibnd ) then
                         raux = raux + nmat(jbnd)
                     endif ! back if ( symm(jbnd) == ibnd ) block
                 enddo ! over jbnd={1,norbs} loop

                 raux = raux / real(hist(ibnd)) ! calculate average value

                 do jbnd=1,norbs                ! setup it
                     if ( symm(jbnd) == ibnd ) then
                         nmat(jbnd) = raux
                     endif ! back if ( symm(jbnd) == ibnd ) block
                 enddo ! over jbnd={1,norbs} loop
             endif ! back if ( hist(ibnd) > 0 ) block
         enddo ! over ibnd={1,norbs} loop
     endif ! back if ( issun == 2 ) block

! symmetrize nmat over spin
     if ( isspn == 1 ) then
         do jbnd=1,nband
             raux = ( nmat(jbnd) + nmat(jbnd+nband) ) / two
             nmat(jbnd) = raux
             nmat(jbnd+nband) = raux
         enddo ! over jbnd={1,nband} loop
     endif ! back if ( isspn == 1 ) block

     return
  end subroutine ctqmc_symm_nmat

!!>>> ctqmc_symm_gtau: symmetrize the gtau according to symm vector
!!>>> only the diagonal elements are taken into considerations
  subroutine ctqmc_symm_gtau(symm, gtau)
     use constants, only : dp, zero, two

     use control, only : issun, isspn
     use control, only : nband, norbs
     use control, only : ntime

     implicit none

! external arguments
! symmetry vector
     integer, intent(in) :: symm(norbs)

! impurity green's function
     real(dp), intent(inout) :: gtau(ntime,norbs,norbs)

! local variables
! loop index over bands
     integer  :: ibnd
     integer  :: jbnd

! loop index over imaginary time points
     integer  :: ktau

! dummy variables
     real(dp) :: raux

! histogram vector
! note: it is NOT the global one
     integer  :: hist(norbs)

! build histogram
     hist = 0
     do ibnd=1,norbs
         hist(symm(ibnd)) = hist(symm(ibnd)) + 1
     enddo ! over ibnd={1,norbs} loop

! perform symmetrization for those orbitals which symm index are identity
     if ( issun == 2 ) then
         do ktau=1,ntime
             do ibnd=1,norbs
                 if ( hist(ibnd) > 0 ) then         ! need to enforce symmetry
                     raux = zero

                     do jbnd=1,norbs                ! gather the data
                         if ( symm(jbnd) == ibnd ) then
                             raux = raux + gtau(ktau,jbnd,jbnd)
                         endif ! back if ( symm(jbnd) == ibnd ) block
                     enddo ! over jbnd={1,norbs} loop

                     raux = raux / real(hist(ibnd)) ! calculate average value

                     do jbnd=1,norbs                ! setup it
                         if ( symm(jbnd) == ibnd ) then
                             gtau(ktau,jbnd,jbnd) = raux
                         endif ! back if ( symm(jbnd) == ibnd ) block
                     enddo ! over jbnd={1,norbs} loop
                 endif ! back if ( hist(ibnd) > 0 ) block
             enddo ! over ibnd={1,norbs} loop
         enddo ! over ktau={1,ntime} loop
     endif ! back if ( issun == 2 ) block

! symmetrize gtau over spin
     if ( isspn == 1 ) then
         do ktau=1,ntime
             do jbnd=1,nband
                 raux = ( gtau(ktau,jbnd,jbnd) + gtau(ktau,jbnd+nband,jbnd+nband) ) / two
                 gtau(ktau,jbnd,jbnd) = raux
                 gtau(ktau,jbnd+nband,jbnd+nband) = raux
             enddo ! over jbnd={1,nband} loop
         enddo ! over ktau={1,ntime} loop
     endif ! back if ( isspn == 1 ) block

     return
  end subroutine ctqmc_symm_gtau

!!>>> ctqmc_symm_grnf: symmetrize the grnf according to symm vector
!!>>> only the diagonal elements are taken into considerations
  subroutine ctqmc_symm_grnf(symm, grnf)
     use constants, only : dp, two, czero

     use control, only : issun, isspn
     use control, only : nband, norbs
     use control, only : mfreq

     implicit none

! external arguments
! symmetry vector
     integer, intent(in) :: symm(norbs)

! impurity green's function
     complex(dp), intent(inout) :: grnf(mfreq,norbs,norbs)

! local variables
! loop index over bands
     integer :: ibnd
     integer :: jbnd

! loop index over matsubara frequencies
     integer :: kfrq

! dummy variables
     complex(dp) :: caux

! histogram vector
! note: it is NOT the global one
     integer :: hist(norbs)

! build histogram
     hist = 0
     do ibnd=1,norbs
         hist(symm(ibnd)) = hist(symm(ibnd)) + 1
     enddo ! over ibnd={1,norbs} loop

! perform symmetrization for those orbitals which symm index are identity
     if ( issun == 2 ) then
         do kfrq=1,mfreq
             do ibnd=1,norbs
                 if ( hist(ibnd) > 0 ) then         ! need to enforce symmetry
                     caux = czero

                     do jbnd=1,norbs                ! gather the data
                         if ( symm(jbnd) == ibnd ) then
                             caux = caux + grnf(kfrq,jbnd,jbnd)
                         endif ! back if ( symm(jbnd) == ibnd ) block
                     enddo ! over jbnd={1,norbs} loop

                     caux = caux / real(hist(ibnd)) ! calculate average value

                     do jbnd=1,norbs                ! setup it
                         if ( symm(jbnd) == ibnd ) then
                             grnf(kfrq,jbnd,jbnd) = caux
                         endif ! back if ( symm(jbnd) == ibnd ) block
                     enddo ! over jbnd={1,norbs} loop
                 endif ! back if ( hist(ibnd) > 0 ) block
             enddo ! over ibnd={1,norbs} loop
         enddo ! over kfrq={1,mfreq} loop
     endif ! back if ( issun == 2 ) block

! symmetrize grnf over spin
     if ( isspn == 1 ) then
         do kfrq=1,mfreq
             do jbnd=1,nband
                 caux = ( grnf(kfrq,jbnd,jbnd) + grnf(kfrq,jbnd+nband,jbnd+nband) ) / two
                 grnf(kfrq,jbnd,jbnd) = caux
                 grnf(kfrq,jbnd+nband,jbnd+nband) = caux
             enddo ! over jbnd={1,nband} loop
         enddo ! over kfrq={1,mfreq} loop
     endif ! back if ( isspn == 1 ) block

     return
  end subroutine ctqmc_symm_grnf

!!========================================================================
!!>>> smooth physical observables                                      <<<
!!========================================================================

!!>>> ctqmc_smth_sigf: smooth impurity self-energy function in low
!!>>> frequency region
  subroutine ctqmc_smth_sigf(sigf)
     use constants, only : dp, czero

     use control, only : nfreq

     implicit none

! external arguments
! impurity self-energy function to be smoothen
     complex(dp), intent(inout) :: sigf(nfreq)

! local variables
! loop index
     integer  :: i
     integer  :: j
     integer  :: k

! smooth radius
     integer  :: lrad
     integer  :: srad

! imaginary part of self-energy function
     real(dp) :: ti

! real part of self-energy function
     real(dp) :: tr

! dummy variables for addition
     complex(dp) :: saux

! dummy self-energy function
     complex(dp) :: stmp(nfreq)

! determine smooth radius
     lrad = nfreq / 4  ! large radius
     srad = nfreq / 16 ! small radius

! |---------|---------|----------------------|
! 1         lrad      2*lrad                 nfreq
! deal with [1,lrad], head part
     do k=1,lrad
         stmp(k) = sigf(k)
     enddo ! over k={1,lrad} loop

! deal with [lrad+1,2*lrad], intermediate part
     do i=1,lrad
         k = lrad + i
         saux = czero
         do j=-srad,srad
             saux = saux + sigf(k+j)
         enddo ! over j={-srad,srad} loop
         stmp(k) = saux / real(2 * srad + 1)
         stmp(k) = ( (lrad - i) * sigf(k) + i * stmp(k) ) / real(lrad)
     enddo ! over i={1,lrad} loop

! deal with [nfreq-2*lrad+1,nfreq], tail part
     do k=nfreq-2*lrad+1,nfreq
         tr =  real( stmp(nfreq-2*lrad) )
         ti = aimag( stmp(nfreq-2*lrad) ) * real(nfreq - 2 * lrad) / real(k)
         stmp(k) = dcmplx(tr,ti)
     enddo ! over k={nfreq-2*lrad+1,nfreq} loop

! copy stmp to sigf
     do k=1,nfreq
         sigf(k) = stmp(k)
     enddo ! over k={1,nfreq} loop

     return
  end subroutine ctqmc_smth_sigf

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
!!>>> atomic eigenstate converter                                      <<<
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
     shift = 0.0_dp; call ctqmc_prep_shift(shift)

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

!!>>> ctqmc_prep_shift: evaluate the shift for the Coulomb interaction and
!!>>> the chemical potential. in fact, shift = 2 K'(\tau = 0)
  subroutine ctqmc_prep_shift(shift)
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
  end subroutine ctqmc_prep_shift

!!========================================================================
!!>>> postprocess physical observables                                 <<<
!!========================================================================

!!>>> ctqmc_make_gtau: build imaginary green's function using orthogonal
!!>>> polynomial representation
  subroutine ctqmc_make_gtau(tmesh, gtau, gaux)
     use constants, only : dp, zero, one, two, pi

     use control, only : isort
     use control, only : norbs
     use control, only : lemax, legrd, chmax, chgrd
     use control, only : ntime
     use control, only : beta
     use context, only : ppleg, qqche

     implicit none

! external arguments
! imaginary time mesh
     real(dp), intent(in)  :: tmesh(ntime)

! impurity green's function/orthogonal polynomial coefficients
     real(dp), intent(in)  :: gtau(ntime,norbs,norbs)

! calculated impurity green's function
     real(dp), intent(out) :: gaux(ntime,norbs,norbs)

! local parameters
! scheme of integral kernel used to damp the Gibbs oscillation
! damp = 0, Dirichlet   mode
! damp = 1, Jackson     mode, preferred
! damp = 2, Lorentz     mode
! damp = 3, Fejer       mode
! damp = 4, Wang-Zunger mode
     integer, parameter :: damp = 0

! local variables
! loop index
     integer  :: i
     integer  :: j

! loop index for legendre polynomial
     integer  :: fleg

! loop index for chebyshev polynomial
     integer  :: fche

! index for imaginary time \tau
     integer  :: curr

! interval for imaginary time slice
     real(dp) :: step

! dummy variables
     real(dp) :: raux

! initialize gaux
     gaux = zero

! select calculation method
     select case ( isort )

         case (1, 4)
             call cat_make_gtau1()

         case (2, 5)
             call cat_make_gtau2()

         case (3, 6)
             call cat_make_gtau3()

     end select

     return

  contains

!!>>> cat_make_kpm: build the integral kernel function
  subroutine cat_make_kpm(kdim, kern)
     implicit none

! external arguments
! dimension of integral kernel function
     integer, intent(in)   :: kdim

! integral kernel function
     real(dp), intent(out) :: kern(kdim)

! local variables
! loop index
     integer :: kcur

     kern = zero
     do kcur=1,kdim
         select case ( damp )

! Dirichlet mode
             case (0)
                 kern(kcur) = one

! Jackson mode
             case (1)
                 i = kcur - 1
                 curr = kdim + 1
                 raux = pi * i / curr
                 kern(kcur) = ( (curr-i) * cos(raux) + sin(raux) / tan(pi/curr) ) / curr

! Lorentz mode
             case (2)
                 kern(kcur) = sinh( one - (kcur - one) / real( kdim ) ) / sinh(one)

! Fejer mode
             case (3)
                 kern(kcur) = one - ( kcur - one ) / real( kdim )

! Wang-Zunger mode
             case (4)
                 kern(kcur) = exp( - ( (kcur - one) / real( kdim ) )**4 )

         end select
     enddo ! over kcur={1,kdim} loop

     return
  end subroutine cat_make_kpm

!!>>> cat_make_gtau1: build impurity green's function using normal
!!>>> representation
  subroutine cat_make_gtau1()
     implicit none

     raux = real(ntime) / (beta * beta)
     do i=1,norbs
         do j=1,ntime
             gaux(j,i,i) = gtau(j,i,i) * raux
         enddo ! over j={1,ntime} loop
     enddo ! over i={1,norbs} loop

     return
  end subroutine cat_make_gtau1

!!>>> cat_make_gtau2: build impurity green's function using legendre
!!>>> polynomial representation
  subroutine cat_make_gtau2()
     implicit none

! integral kernel
     real(dp) :: ker1(lemax)

! build kernel function at first
     ker1 = one; call cat_make_kpm(lemax, ker1)

! reconstruct green's function
     step = real(legrd - 1) / two
     do i=1,norbs
         do j=1,ntime
             raux = two * tmesh(j) / beta
             curr = nint(raux * step) + 1
             do fleg=1,lemax
                 raux = sqrt(two * fleg - 1) / (beta * beta) * ker1(fleg)
                 gaux(j,i,i) = gaux(j,i,i) + raux * gtau(fleg,i,i) * ppleg(curr,fleg)
             enddo ! over fleg={1,lemax} loop
         enddo ! over j={1,ntime} loop
     enddo ! over i={1,norbs} loop

     return
  end subroutine cat_make_gtau2

!!>>> cat_make_gtau3: build impurity green's function using chebyshev
!!>>> polynomial representation
  subroutine cat_make_gtau3()
     implicit none

! integral kernel
     real(dp) :: ker2(chmax)

! build kernel function at first
     ker2 = one; call cat_make_kpm(chmax, ker2)

! reconstruct green's function
     step = real(chgrd - 1) / two
     do i=1,norbs
         do j=1,ntime
             raux = two * tmesh(j) / beta
             curr = nint(raux * step) + 1
             do fche=1,chmax
                 raux = two / (beta * beta) * ker2(fche)
                 gaux(j,i,i) = gaux(j,i,i) + raux * gtau(fche,i,i) * qqche(curr,fche)
             enddo ! over fche={1,chmax} loop
         enddo ! over j={1,ntime} loop
     enddo ! over i={1,norbs} loop

     return
  end subroutine cat_make_gtau3
  end subroutine ctqmc_make_gtau

!!>>> ctqmc_make_ftau: build auxiliary correlation function using
!!>>> orthogonal polynomial representation, F(\tau)
  subroutine ctqmc_make_ftau(tmesh, ftau, faux)
     use constants, only : dp, zero, two

     use control, only : isort
     use control, only : norbs
     use control, only : lemax, legrd, chmax, chgrd
     use control, only : ntime
     use control, only : beta
     use context, only : ppleg, qqche

     implicit none

! external arguments
! imaginary time mesh
     real(dp), intent(in)  :: tmesh(ntime)

! auxiliary correlation function/orthogonal polynomial coefficients
     real(dp), intent(in)  :: ftau(ntime,norbs,norbs)

! calculated auxiliary correlation function
     real(dp), intent(out) :: faux(ntime,norbs,norbs)

! local variables
! loop index
     integer  :: i
     integer  :: j
     integer  :: k

! loop index for legendre polynomial
     integer  :: fleg

! loop index for chebyshev polynomial
     integer  :: fche

! index for imaginary time \tau
     integer  :: curr

! interval for imaginary time slice
     real(dp) :: step

! dummy variables
     real(dp) :: raux

! initialize faux
     faux = zero

! select calculation method
     select case ( isort )

         case (4)
             call cat_make_ftau1()

         case (5)
             call cat_make_ftau2()

         case (6)
             call cat_make_ftau3()

     end select

     return

  contains

!!>>> cat_make_ftau1: build auxiliary correlation function using normal
!!>>> representation
  subroutine cat_make_ftau1()
     implicit none

     raux = real(ntime) / (beta * beta)
     do i=1,norbs
         do j=1,norbs
             do k=1,ntime
                 faux(k,j,i) = ftau(k,j,i) * raux
             enddo ! over k={1,ntime} loop
         enddo ! over j={1,norbs} loop
     enddo ! over i={1,norbs} loop

     return
  end subroutine cat_make_ftau1

!!>>> cat_make_ftau2: build auxiliary correlation function using legendre
!!>>> polynomial representation
  subroutine cat_make_ftau2()
     implicit none

     step = real(legrd - 1) / two
     do i=1,norbs
         do j=1,norbs
             do k=1,ntime
                 raux = two * tmesh(k) / beta
                 curr = nint(raux * step) + 1
                 do fleg=1,lemax
                     raux = sqrt(two * fleg - 1) / (beta * beta)
                     faux(k,j,i) = faux(k,j,i) + raux * ftau(fleg,j,i) * ppleg(curr,fleg)
                 enddo ! over fleg={1,lemax} loop
             enddo ! over k={1,ntime} loop
         enddo ! over j={1,norbs} loop
     enddo ! over i={1,norbs} loop

     return
  end subroutine cat_make_ftau2

!!>>> cat_make_ftau3: build auxiliary correlation function using chebyshev
!!>>> polynomial representation
  subroutine cat_make_ftau3()
     implicit none

     step = real(chgrd - 1) / two
     do i=1,norbs
         do j=1,norbs
             do k=1,ntime
                 raux = two * tmesh(k) / beta
                 curr = nint(raux * step) + 1
                 raux = two / (beta * beta)
                 do fche=1,chmax
                     faux(k,j,i) = faux(k,j,i) + raux * ftau(fche,j,i) * qqche(curr,fche)
                 enddo ! over fche={1,chmax} loop
             enddo ! over k={1,ntime} loop
         enddo ! over j={1,norbs} loop
     enddo ! over i={1,norbs} loop

     return
  end subroutine cat_make_ftau3
  end subroutine ctqmc_make_ftau

!!========================================================================
!!>>> build prefactor for improved estimator                           <<<
!!========================================================================

!!>>> ctqmc_make_iret: to calculate the integral I(\tau_end) which is very
!!>>> important when retarded interaction is included
  subroutine ctqmc_make_iret(time, iret)
     use constants, only : dp, zero, two

     use control, only : norbs
     use context, only : index_s, index_e, time_s, time_e
     use context, only : rank

     implicit none

! external arguments
! imaginary time point, in principle, it is \tau_end
     real(dp), intent(in)  :: time

! integral value for I(\tau_end)
     real(dp), intent(out) :: iret

! local variables
! loop index for start and end points
     integer  :: it

! loop index for flavor channel
     integer  :: flvr

! length betweem two time points
     real(dp) :: dtau
     real(dp) :: daux

! calculate integral I(\tau_end), the equation is
! Eq. (39) in Phys. Rev. B 89, 235128 (2014)
     iret = zero

     do flvr=1,norbs
         do it=1,rank(flvr)
! contribution from create operators
             dtau = time_s( index_s(it, flvr), flvr ) - time
             if ( dtau >= zero ) then
                 call cat_weight_kernel(2, +dtau, daux)
                 iret = iret + daux
             else
                 call cat_weight_kernel(2, -dtau, daux)
                 iret = iret - daux
             endif ! back if ( dtau >= zero ) block

! contribution from destroy operators
             dtau = time_e( index_e(it, flvr), flvr ) - time
             if ( dtau >= zero ) then
                 call cat_weight_kernel(2, +dtau, daux)
                 iret = iret - daux
             else
                 call cat_weight_kernel(2, -dtau, daux)
                 iret = iret + daux
             endif ! back if ( dtau >= zero ) block
         enddo ! over it={1,rank(flvr)} loop
     enddo ! over flvr={1,norbs} loop

! add additional term
     call cat_weight_kernel(2, zero, daux)
     iret = -iret - two * daux

     return
  end subroutine ctqmc_make_iret

!!>>> ctqmc_make_pref: to calculate the prefactor used by the improved
!!>>> estimator for self-energy function and vertex function
  subroutine ctqmc_make_pref()
     use constants, only : dp, zero, half

     use control, only : isscr
     use control, only : norbs
     use context, only : index_e, time_e
     use context, only : rank, pref, uumat

     implicit none

! local variables
! loop index for start and end points
     integer  :: it

! loop index for flavor channel
     integer  :: flvr

! loop index for colour channel
     integer  :: clur

! occupation number at \tau_end
     real(dp) :: occu

! integral value for I(\tau_end)
     real(dp) :: iret

! if it is holstein-hubbard model, the improved estimator algorithm
! can not be used now
     call s_assert2(isscr /= 2,'sorry, isscr = 2 is not compatible with improved estimator')

     do flvr=1,norbs
         do it=1,rank(flvr)

! reset the prefactor
             pref(it,flvr) = zero

! calculate normal contribution
             do clur=1,norbs
                 call cat_occupy_status(clur, time_e( index_e(it, flvr), flvr ), occu)
                 pref(it,flvr) = pref(it,flvr) + half * ( uumat(flvr,clur) + uumat(clur,flvr) ) * occu
             enddo ! over clur={1,norbs} loop

! if retarded interaction is considered, we have to include the
! contribution from I(\tau_end)
             if ( isscr > 1 ) then
                 call ctqmc_make_iret(time_e( index_e(it, flvr), flvr ), iret)
                 pref(it,flvr) = pref(it,flvr) + iret
             endif ! back if ( isscr > 1 ) block
         enddo ! over it={1,rank(flvr)} loop
     enddo ! over flvr={1,norbs} loop

     return
  end subroutine ctqmc_make_pref

!!========================================================================
!!>>> build auxiliary two-particle related variables                   <<<
!!========================================================================

!!>>> ctqmc_make_prod: try to calculate the product of matsubara
!!>>> frequency exponents exp(i \omega_n \tau)
  subroutine ctqmc_make_prod(flvr, nfaux, mrank, caux1, caux2)
     use constants, only : dp, two, pi, czi

     use control, only : nffrq
     use control, only : beta
     use context, only : index_s, index_e, time_s, time_e
     use context, only : rank

     implicit none

! external arguments
! current flavor channel
     integer, intent(in) :: flvr

! combination of nffrq and nbfrq
     integer, intent(in) :: nfaux

! maximum number of operators in different flavor channels
     integer, intent(in) :: mrank

! matsubara frequency exponents for create operators
     complex(dp), intent(out) :: caux1(nfaux,mrank)

! matsubara frequency exponents for destroy operators
     complex(dp), intent(out) :: caux2(nfaux,mrank)

! local variables
! loop indices for start and end points
     integer  :: is
     integer  :: ie

! imaginary time for start and end points
     real(dp) :: taus
     real(dp) :: taue

! for create operators
     do is=1,rank(flvr)
         taus = time_s( index_s(is, flvr), flvr )
         caux1(:,is) = exp(-two * czi * pi * taus / beta)
         call s_cumprod_z(nfaux, caux1(:,is), caux1(:,is))
         caux1(:,is) = caux1(:,is) * exp(+(nffrq + 1) * czi * pi * taus / beta)
     enddo ! over is={1,rank(flvr)} loop

! for destroy operators
     do ie=1,rank(flvr)
         taue = time_e( index_e(ie, flvr), flvr )
         caux2(:,ie) = exp(+two * czi * pi * taue / beta)
         call s_cumprod_z(nfaux, caux2(:,ie), caux2(:,ie))
         caux2(:,ie) = caux2(:,ie) * exp(-(nffrq + 1) * czi * pi * taue / beta)
     enddo ! over ie={1,rank(flvr)} loop

     return
  end subroutine ctqmc_make_prod

!!========================================================================
!!>>> build self-energy function                                       <<<
!!========================================================================

!!>>> ctqmc_make_hub1: build atomic green's function and self-energy
!!>>> function using improved Hubbard-I approximation, and then make
!!>>> interpolation for self-energy function between low frequency QMC
!!>>> data and high frequency Hubbard-I approximation data, the full
!!>>> impurity green's function can be obtained by using dyson's equation
!!>>> finally
  subroutine ctqmc_make_hub1()
     use constants, only : dp, zero, one, two, czi, czero

     use control, only : norbs, ncfgs
     use control, only : mfreq
     use control, only : nfreq
     use control, only : mune
     use control, only : myid, master
     use context, only : rmesh
     use context, only : prob, nmat
     use context, only : eimp, uumat
     use context, only : grnf
     use context, only : hybf
     use context, only : sig2

     implicit none

! local parameters
! maximum allowable number of non-zero elements in F matrix
     integer, parameter :: nzero = 1024

! local variables
! loop index
     integer  :: i
     integer  :: j
     integer  :: k
     integer  :: m
     integer  :: n

! dummy integer variables, used to build F matrix
     integer  :: start
     integer  :: value
     integer  :: permute

! dummy real variables, used to interpolate self-energy function
     real(dp) :: ob, oe
     real(dp) :: d0, d1
     real(dp) :: shift

! dummy complex variables, used to interpolate self-energy function
     complex(dp) :: cb, ce
     complex(dp) :: sinf

! dummy atomic states: alpha, beta, gamma
     integer  :: sa(norbs)
     integer  :: sb(norbs)
     integer  :: sc(norbs)

! atomic basis sets
     integer  :: basis(ncfgs,norbs)

! F matrix, < alpha | f_{n} | beta >
     integer  :: fcounter(norbs)
     integer  :: fa(nzero,norbs)
     integer  :: fb(nzero,norbs)
     integer  :: fv(nzero,norbs)

! eigenvalues for local hmailtonian
     real(dp) :: eaux(ncfgs)

! dummy imurity green's function: G^{-1}
     complex(dp) :: gaux(norbs,norbs)

! atomic green's function and self-energy function in Hubbard-I approximation
     complex(dp) :: ghub(mfreq,norbs)
     complex(dp) :: shub(mfreq,norbs)

! evaluate the shift for the Coulomb interaction and chemical potential
! if the retarded interaction is used
     call ctqmc_prep_shift(shift)

! build atomic basis set, we do not order them according to their
! occupation numbers
     do i=1,ncfgs
         do j=1,norbs
             if ( btest(i-1,j-1) .eqv. .true. ) then
                 basis(i,j) = 1
             else
                 basis(i,j) = 0
             endif ! back if ( btest(i-1,j-1) .eqv. .true. ) block
         enddo ! over j={1,norbs} loop
     enddo ! over i={1,ncfgs} loop

! evaluate atomic eigenvalues directly
     eaux = zero
     do i=1,ncfgs
         do j=1,norbs
             eaux(i) = eaux(i) + ( eimp(j) - mune ) * basis(i,j)
         enddo ! over j={1,norbs} loop
         do j=1,norbs-1
             do k=j+1,norbs
                 if ( basis(i,j) == 1 .and. basis(i,k) == 1 ) then
                     eaux(i) = eaux(i) + uumat(j,k)
                 endif ! back if ( basis(i,j) == 1 .and. basis(i,k) == 1 ) block
             enddo ! over k={j+1,norbs} loop
         enddo ! over j={1,norbs-1} loop
     enddo ! over i={1,ncfgs} loop

! build F matrix < alpha | f_{n} | beta >
! note 1: to save the memory and accelerate the computation, we only store
! the non-zero element of F matrix
! note 2: it is crucial to check whether the number of non-zero elements
! exceed limit (nzero)
     fcounter = 0
     alpha_loop: do i=1,ncfgs
         sa = basis(i,:)
         beta_loop: do j=1,ncfgs
             sb = basis(j,:)

             orbital_loop: do m=1,norbs
                 sc = sb

                 if ( sc(m) == 1 ) then
                     permute = 1
                     do n=1,m-1
                         if ( sc(n) == 1 ) permute = -permute
                     enddo ! over n={1,m-1} loop
                     sc(m) = 0

                     value = 1
                     do n=1,norbs
                         if ( sa(n) /= sc(n) ) value = 0
                     enddo ! over n={1,norbs} loop
                     value = value * permute
                 else
                     value = 0
                 endif ! back if ( sc(m) == 1 ) block

                 if ( value /= 0 ) then
                     fcounter(m) = fcounter(m) + 1
                     if ( fcounter(m) > nzero ) then
                         call s_print_error('ctqmc_make_hub1','non-zero elements exceed limit')
                     endif ! back if ( fcounter(m) > nzero ) block
                     fa(fcounter(m),m) = i
                     fb(fcounter(m),m) = j
                     fv(fcounter(m),m) = value
                 endif ! back if ( value /= 0 ) block
             enddo orbital_loop ! over m={1,norbs} loop

         enddo beta_loop ! over j={1,ncfgs} loop
     enddo alpha_loop ! over i={1,ncfgs} loop

! calculate atomic green's function using Hubbard-I approximation
     do i=1,norbs
         do k=1,mfreq
             ghub(k,i) = czero
             do m=1,fcounter(i)
                 ob = fv(m,i) * fv(m,i) * ( prob(fa(m,i)) + prob(fb(m,i)) )
                 cb = czi * rmesh(k) + eaux(fa(m,i)) - eaux(fb(m,i))
                 ghub(k,i) = ghub(k,i) + ob / cb
             enddo ! over m={1,fcounter(i)} loop
         enddo ! over k={1,mfreq} loop
     enddo ! over i={1,norbs} loop

! calculate atomic self-energy function using dyson's equation
     do i=1,norbs
         do k=1,mfreq
             shub(k,i) = czi * rmesh(k) + mune - eimp(i) - one / ghub(k,i)
             shub(k,i) = shub(k,i) - nmat(i) * shift
         enddo ! over k={1,mfreq} loop
     enddo ! over i={1,norbs} loop

! dump the ghub and shub, only for reference, only the master node can do it
     if ( myid == master ) then
         call ctqmc_dump_hub1(rmesh, ghub, shub)
     endif ! back if ( myid == master ) block

! build self-energy function at low frequency region
!-------------------------------------------------------------------------
! filter grnf to suppress the fluctuation of its real part
!-------------------------------------------------------------------------
!<     do k=1,nfreq
!<         do i=1,norbs
!<             ob =  real( grnf(k,i,i) )
!<             oe = aimag( grnf(k,i,i) )
!<             grnf(k,i,i) = dcmplx( zero, oe )
!<         enddo ! over i={1,norbs} loop
!<     enddo ! over k={1,nfreq} loop
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
     do k=1,nfreq
         gaux = grnf(k,:,:)
         call s_inv_z(norbs, gaux)
         do i=1,norbs
             sig2(k,i,i) = czi * rmesh(k) + mune - eimp(i) - gaux(i,i) - hybf(k,i,i)
             sig2(k,i,i) = sig2(k,i,i) - shift / two
         enddo ! over i={1,norbs} loop
     enddo ! over k={1,nfreq} loop
!-------------------------------------------------------------------------
! filter sig2 to suppress the fluctuation of its imaginary part
!-------------------------------------------------------------------------
     do k=1,16
         do i=1,norbs
             call ctqmc_smth_sigf( sig2(1:nfreq,i,i) ) ! smooth it 16 times
         enddo ! over i={1,norbs} loop
     enddo ! over k={1,16} loop
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

! interpolates self-energy function between low energy QMC data and high
! energy Hubbard-I approximation
     do i=1,norbs

! determine the base point, its value is calculated by using five points
         cb = czero
         do k=nfreq-4,nfreq
             cb = cb + sig2(k,i,i)
         enddo ! over k={nfreq-4,nfreq} loop

         cb = cb / real(5)
         ob = rmesh(nfreq-2)

! step A: for the imaginary part
! determine the intermediate region [nfreq+1,start] at first
         start = 0
         do k=nfreq+1,mfreq
             start = k
             d0 = aimag( shub(k,i) - cb ) / ( rmesh(k) - ob )
             d1 = aimag( shub(k,i) - shub(k-1,i) ) / ( rmesh(k) - rmesh(k-1) )
             if ( abs( d0 - d1 ) < 0.02_dp ) EXIT
         enddo ! over k={nfreq+1,mfreq} loop

! we just constrain start \in [nfreq + 32, nfreq + 128]
         if ( start - nfreq <  32 ) start = nfreq +  32
         if ( start - nfreq > 128 ) start = nfreq + 128

         ce = shub(start,i)
         oe = rmesh(start)

! deal with the intermediate region, using linear interpolation
         do k=nfreq+1,start
             sig2(k,i,i) = dcmplx( zero, aimag(cb) + aimag( ce - cb ) * ( rmesh(k) - ob ) / ( oe - ob ) )
         enddo ! over k={nfreq+1,start} loop

! deal with the tail region, using atomic self-energy function directly
         do k=start+1,mfreq
             sig2(k,i,i) = dcmplx( zero, aimag( shub(k,i) ) )
         enddo ! over k={start+1,mfreq} loop

! step B: for the real part
         sinf = shub(mfreq,i)
         do k=nfreq+1,mfreq
             sig2(k,i,i) = sig2(k,i,i) + real(sinf) + ( ob / rmesh(k) )**2 * real( cb - sinf )
         enddo ! over k={nfreq+1,mfreq} loop

     enddo ! over i={1,norbs} loop

! calculate final impurity green's function using dyson's equation
     do k=1,mfreq
         gaux = czero
         do i=1,norbs
             gaux(i,i) = czi * rmesh(k) + mune - eimp(i) - sig2(k,i,i) - hybf(k,i,i)
             gaux(i,i) = gaux(i,i) - shift / two
         enddo ! over i={1,norbs} loop
         call s_inv_z(norbs, gaux)
         grnf(k,:,:) = gaux
     enddo ! over k={1,mfreq} loop

     return
  end subroutine ctqmc_make_hub1

!!>>> ctqmc_make_hub2: build atomic green's function and self-energy
!!>>> function using improved Hubbard-I approximation, and then make
!!>>> forward fourier transformation for impurity green's function and
!!>>> auxiliary correlation function. then the final self-energy function
!!>>> is obtained by analytical formula.
  subroutine ctqmc_make_hub2()
     use constants, only : dp, zero, one, two, pi, czi, czero

     use control, only : isort
     use control, only : norbs, ncfgs
     use control, only : lemax
     use control, only : mfreq
     use control, only : ntime
     use control, only : mune, beta
     use control, only : myid, master
     use context, only : tmesh, rmesh
     use context, only : prob, nmat
     use context, only : eimp, uumat
     use context, only : gtau, ftau, grnf, frnf
     use context, only : sig2

     implicit none

! local parameters
! maximum allowable number of non-zero elements in F matrix
     integer, parameter :: nzero = 1024

! local variables
! loop index
     integer  :: i
     integer  :: j
     integer  :: k
     integer  :: m
     integer  :: n

! dummy integer variables, used to build F matrix
     integer  :: value
     integer  :: permute

! dummy real variables, used to build atomic green's function
     real(dp) :: ob
     real(dp) :: shift

! dummy complex variables, used to build atomic green's function
     complex(dp) :: cb

! dummy atomic states: alpha, beta, gamma
     integer  :: sa(norbs)
     integer  :: sb(norbs)
     integer  :: sc(norbs)

! atomic basis sets
     integer  :: basis(ncfgs,norbs)

! F matrix, < alpha | f_{n} | beta >
     integer  :: fcounter(norbs)
     integer  :: fa(nzero,norbs)
     integer  :: fb(nzero,norbs)
     integer  :: fv(nzero,norbs)

! eigenvalues for local hmailtonian
     real(dp) :: eaux(ncfgs)

! spherical Bessel functions
     real(dp) :: jaux(mfreq,lemax)

! imaginary time green's function
     real(dp) :: gaux(ntime,norbs,norbs)

! auxiliary correlation function on imaginary time axis
     real(dp) :: faux(ntime,norbs,norbs)

! unitary transformation matrix for legendre polynomial
     complex(dp) :: taux(mfreq,lemax)

! atomic green's function and self-energy function in Hubbard-I approximation
     complex(dp) :: ghub(mfreq,norbs)
     complex(dp) :: shub(mfreq,norbs)

! evaluate the shift for the Coulomb interaction and chemical potential
! if the retarded interaction is used
     call ctqmc_prep_shift(shift)

! build atomic basis set, we do not order them according to their
! occupation numbers
     do i=1,ncfgs
         do j=1,norbs
             if ( btest(i-1,j-1) .eqv. .true. ) then
                 basis(i,j) = 1
             else
                 basis(i,j) = 0
             endif ! back if ( btest(i-1,j-1) .eqv. .true. ) block
         enddo ! over j={1,norbs} loop
     enddo ! over i={1,ncfgs} loop

! evaluate atomic eigenvalues directly
     eaux = zero
     do i=1,ncfgs
         do j=1,norbs
             eaux(i) = eaux(i) + ( eimp(j) - mune ) * basis(i,j)
         enddo ! over j={1,norbs} loop
         do j=1,norbs-1
             do k=j+1,norbs
                 if ( basis(i,j) == 1 .and. basis(i,k) == 1 ) then
                     eaux(i) = eaux(i) + uumat(j,k)
                 endif ! back if ( basis(i,j) == 1 .and. basis(i,k) == 1 ) block
             enddo ! over k={j+1,norbs} loop
         enddo ! over j={1,norbs-1} loop
     enddo ! over i={1,ncfgs} loop

! build F matrix < alpha | f_{n} | beta >
! note 1: to save the memory and accelerate the computation, we only store
! the non-zero element of F matrix
! note 2: it is crucial to check whether the number of non-zero elements
! exceed limit (nzero)
     fcounter = 0
     alpha_loop: do i=1,ncfgs
         sa = basis(i,:)
         beta_loop: do j=1,ncfgs
             sb = basis(j,:)

             orbital_loop: do m=1,norbs
                 sc = sb

                 if ( sc(m) == 1 ) then
                     permute = 1
                     do n=1,m-1
                         if ( sc(n) == 1 ) permute = -permute
                     enddo ! over n={1,m-1} loop
                     sc(m) = 0

                     value = 1
                     do n=1,norbs
                         if ( sa(n) /= sc(n) ) value = 0
                     enddo ! over n={1,norbs} loop
                     value = value * permute
                 else
                     value = 0
                 endif ! back if ( sc(m) == 1 ) block

                 if ( value /= 0 ) then
                     fcounter(m) = fcounter(m) + 1
                     if ( fcounter(m) > nzero ) then
                         call s_print_error('ctqmc_make_hub2','non-zero elements exceed limit')
                     endif ! back if ( fcounter(m) > nzero ) block
                     fa(fcounter(m),m) = i
                     fb(fcounter(m),m) = j
                     fv(fcounter(m),m) = value
                 endif ! back if ( value /= 0 ) block
             enddo orbital_loop ! over m={1,norbs} loop

         enddo beta_loop ! over j={1,ncfgs} loop
     enddo alpha_loop ! over i={1,ncfgs} loop

! calculate atomic green's function using Hubbard-I approximation
     do i=1,norbs
         do k=1,mfreq
             ghub(k,i) = czero
             do m=1,fcounter(i)
                 ob = fv(m,i) * fv(m,i) * ( prob(fa(m,i)) + prob(fb(m,i)) )
                 cb = czi * rmesh(k) + eaux(fa(m,i)) - eaux(fb(m,i))
                 ghub(k,i) = ghub(k,i) + ob / cb
             enddo ! over m={1,fcounter(i)} loop
         enddo ! over k={1,mfreq} loop
     enddo ! over i={1,norbs} loop

! calculate atomic self-energy function using dyson's equation
     do i=1,norbs
         do k=1,mfreq
             shub(k,i) = czi * rmesh(k) + mune - eimp(i) - one / ghub(k,i)
             shub(k,i) = shub(k,i) - nmat(i) * shift
         enddo ! over k={1,mfreq} loop
     enddo ! over i={1,norbs} loop

! dump the ghub and shub, only for reference, only the master node can do it
     if ( myid == master ) then
         call ctqmc_dump_hub1(rmesh, ghub, shub)
     endif ! back if ( myid == master ) block

! build final impurity green's function and then transform them into
! matsubara frequency axis
! note: only for isort == 4 .or. isort == 6 cases
     if ( isort /= 5 ) then
         call ctqmc_make_gtau(tmesh, gtau, gaux)
         call ctqmc_four_htau(gaux, grnf)
     endif ! back if ( isort /= 5 ) block

! build final auxiliary correlation function and then transform them into
! matsubara frequency axis
! note: only for isort == 4 .or. isort == 6 cases
     if ( isort /= 5 ) then
         call ctqmc_make_ftau(tmesh, ftau, faux)
         call ctqmc_four_htau(faux, frnf)
     endif ! back if ( isort /= 5 ) block

! special consideration must be taken for legendre representation, we can
! calculate grnf and frnf directly by using legendre coefficients, instead
! of performing fourier transformation
     if ( isort == 5 ) then
! build spherical Bessel functions: jaux
         jaux = zero
         do k=1,mfreq
             ob = (two * k - one) * pi / two
             call s_sbessel(lemax-1, ob, jaux(k,:))
         enddo ! over k={1,mfreq} loop

! build unitary transformation matrix: taux
         taux = czero
         do i=1,lemax
             do k=1,mfreq
                 ob = (-one)**(k - 1) * sqrt(two * i - one)
                 cb = czi**i
                 taux(k,i) = jaux(k,i) * ob * cb
             enddo ! over k={1,mfreq} loop
         enddo ! over i={1,lemax} loop

! rebuild impurity green's function on matsubara frequency (grnf) using
! orthogonal polynomial representation, G(i\omega)
! rebuild auxiliary correlation function on matsubara frequency (frnf)
! using orthogonal polynomial representation, F(i\omega)
         grnf = czero
         frnf = czero
         do i=1,norbs
             do j=1,lemax
                 do k=1,mfreq
                     grnf(k,i,i) = grnf(k,i,i) + taux(k,j) * gtau(j,i,i) / beta
                     frnf(k,i,i) = frnf(k,i,i) + taux(k,j) * ftau(j,i,i) / beta
                 enddo ! over k={1,mfreq} loop
             enddo ! over j={1,lemax} loop
         enddo ! over i={1,norbs} loop
     endif ! back if ( isort == 5 ) block

! build full self-energy function by using frnf and grnf
     do i=1,norbs
         do k=1,mfreq
             sig2(k,i,i) = frnf(k,i,i) / grnf(k,i,i)
         enddo ! over k={1,nfreq} loop
     enddo ! over i={1,norbs} loop

     return
  end subroutine ctqmc_make_hub2
