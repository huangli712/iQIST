!!!-----------------------------------------------------------------------
!!! project : gardenia
!!! program : ctqmc_four_htau
!!!           ctqmc_four_hybf
!!!           ctqmc_eval_htau
!!!           ctqmc_eval_hsed
!!!           ctqmc_make_uumat
!!!           ctqmc_make_state
!!!           ctqmc_symm_nmat
!!!           ctqmc_symm_gtau
!!!           ctqmc_symm_grnf
!!!           ctqmc_smth_sigf   <<<---
!!!           ctqmc_make_gtau
!!!           ctqmc_make_ftau   <<<---
!!!           ctqmc_make_iret
!!!           ctqmc_make_pref   <<<---
!!!           ctqmc_make_prod   <<<---
!!!           ctqmc_make_hub1
!!!           ctqmc_make_hub2   <<<---
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
!!>>> cubic spline interpolation                                       <<<
!!========================================================================

!! To provide cubic spline subroutines and wrapper functions to interpolate
!! the hybridization function in imaginary-time axis.

!!>>> ctqmc_eval_htau: evaluate the matrix elements for mmat matrix using
!!>>> cubic spline interpolation
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

!!>>> ctqmc_eval_hsed: calculate the second order derivates of hybridization
!!>>> function on imaginary time space
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

! loop index over imaginary-time points
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


