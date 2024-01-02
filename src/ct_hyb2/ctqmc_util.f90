!!!-----------------------------------------------------------------------
!!! project : manjushaka
!!! program : ctqmc_four_htau
!!!           ctqmc_four_hybf <<<---
!!!           ctqmc_eval_htau
!!!           ctqmc_eval_hsed <<<---
!!!           ctqmc_symm_nimp
!!!           ctqmc_symm_gtau
!!!           ctqmc_symm_grnf <<<---
!!!           ctqmc_make_gtau
!!!           ctqmc_make_ftau <<<---
!!!           ctqmc_make_prod <<<---
!!!           ctqmc_make_hub2 <<<---
!!! source  : ctqmc_util.f90
!!! type    : functions & subroutines
!!! author  : li huang (email:lihuang.dmft@gmail.com)
!!! history : 10/01/2008 by li huang (created)
!!!           05/18/2017 by li huang (last modified)
!!! purpose : provide utility functions and subroutines for hybridization
!!!           expansion version continuous time quantum Monte Carlo (CTQMC)
!!!           quantum impurity solver.
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
!!>>> cubic spline interpolation                                       <<<
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
  subroutine ctqmc_eval_hsed(htau, hsed)
     use constants, only : dp, zero

     use control, only : norbs
     use control, only : ntime
     use control, only : beta

     use context, only : tmesh

     implicit none

! external arguments
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
!!>>> symmetry operation                                               <<<
!!========================================================================

!!
!! @sub ctqmc_symm_nimp
!!
!! symmetrize the occupation number array, nimp, according to symm vector
!!
  subroutine ctqmc_symm_nimp(symm, nimp)
     use constants, only : dp, zero, two

     use control, only : isbnd, isspn
     use control, only : nband, norbs

     implicit none

! external arguments
! symmetry vector
     integer, intent(in) :: symm(norbs)

! occupation number
     real(dp), intent(inout) :: nimp(norbs)

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

! perform symmetrization for those orbitals with the same symmetry
     if ( isbnd == 2 ) then
         do ibnd=1,norbs
             if ( hist(ibnd) > 0 ) then         ! need to enforce symmetry
                 raux = zero

                 do jbnd=1,norbs                ! gather the data
                     if ( symm(jbnd) == ibnd ) then
                         raux = raux + nimp(jbnd)
                     endif ! back if ( symm(jbnd) == ibnd ) block
                 enddo ! over jbnd={1,norbs} loop

                 raux = raux / real(hist(ibnd)) ! calculate average value

                 do jbnd=1,norbs                ! setup it
                     if ( symm(jbnd) == ibnd ) then
                         nimp(jbnd) = raux
                     endif ! back if ( symm(jbnd) == ibnd ) block
                 enddo ! over jbnd={1,norbs} loop
             endif ! back if ( hist(ibnd) > 0 ) block
         enddo ! over ibnd={1,norbs} loop
     endif ! back if ( isbnd == 2 ) block

! symmetrize nimp over spin
     if ( isspn == 2 ) then
         do jbnd=1,nband
             raux = ( nimp(jbnd) + nimp(jbnd+nband) ) / two
             nimp(jbnd) = raux
             nimp(jbnd+nband) = raux
         enddo ! over jbnd={1,nband} loop
     endif ! back if ( isspn == 2 ) block

     return
  end subroutine ctqmc_symm_nimp

!!
!! @sub ctqmc_symm_gtau
!!
!! symmetrize the gtau according to symm vector. only the diagonal
!! elements are taken into considerations
!!
  subroutine ctqmc_symm_gtau(symm, gtau)
     use constants, only : dp, zero, two

     use control, only : isbnd, isspn
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

! perform symmetrization for those orbitals with the same symmetry
     if ( isbnd == 2 ) then
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
     endif ! back if ( isbnd == 2 ) block

! symmetrize gtau over spin
     if ( isspn == 2 ) then
         do ktau=1,ntime
             do jbnd=1,nband
                 raux = ( gtau(ktau,jbnd,jbnd) + gtau(ktau,jbnd+nband,jbnd+nband) ) / two
                 gtau(ktau,jbnd,jbnd) = raux
                 gtau(ktau,jbnd+nband,jbnd+nband) = raux
             enddo ! over jbnd={1,nband} loop
         enddo ! over ktau={1,ntime} loop
     endif ! back if ( isspn == 2 ) block

     return
  end subroutine ctqmc_symm_gtau

!!
!! @sub ctqmc_symm_grnf
!!
!! symmetrize the grnf according to symm vector. only the diagonal
!! elements are taken into considerations
!!
  subroutine ctqmc_symm_grnf(symm, grnf)
     use constants, only : dp, two, czero

     use control, only : isbnd, isspn
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

! perform symmetrization for those orbitals with the same symmetry
     if ( isbnd == 2 ) then
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
     endif ! back if ( isbnd == 2 ) block

! symmetrize grnf over spin
     if ( isspn == 2 ) then
         do kfrq=1,mfreq
             do jbnd=1,nband
                 caux = ( grnf(kfrq,jbnd,jbnd) + grnf(kfrq,jbnd+nband,jbnd+nband) ) / two
                 grnf(kfrq,jbnd,jbnd) = caux
                 grnf(kfrq,jbnd+nband,jbnd+nband) = caux
             enddo ! over jbnd={1,nband} loop
         enddo ! over kfrq={1,mfreq} loop
     endif ! back if ( isspn == 2 ) block

     return
  end subroutine ctqmc_symm_grnf

!!========================================================================
!!>>> advanced representation                                          <<<
!!========================================================================

!!
!! @sub ctqmc_make_gtau
!!
!! build imaginary time green's function using different representation
!!
  subroutine ctqmc_make_gtau(tmesh, gtau, gaux)
     use constants, only : dp, zero, two

     use control, only : isort
     use control, only : norbs
     use control, only : lemax, legrd
     use control, only : ntime
     use control, only : beta

     use context, only : rep_l

     implicit none

! external arguments
! imaginary time mesh
     real(dp), intent(in)  :: tmesh(ntime)

! impurity green's function/orthogonal polynomial coefficients
     real(dp), intent(in)  :: gtau(ntime,norbs,norbs)

! calculated impurity green's function
     real(dp), intent(out) :: gaux(ntime,norbs,norbs)

! local variables
! loop index
     integer  :: i
     integer  :: j

! loop index for legendre polynomial
     integer  :: fleg

! index for imaginary time \tau
     integer  :: curr

! interval for imaginary time slice
     real(dp) :: step

! dummy variables
     real(dp) :: raux

! initialize gaux
     gaux = zero

!-------------------------------------------------------------------------
! using normal representation
!-------------------------------------------------------------------------
     STD_BLOCK: if ( isort == 1 ) then
         raux = real(ntime) / (beta * beta)
         do i=1,norbs
             do j=1,ntime
                 gaux(j,i,i) = gtau(j,i,i) * raux
             enddo ! over j={1,ntime} loop
         enddo ! over i={1,norbs} loop
     endif STD_BLOCK ! back if ( isort == 1 ) block
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

!-------------------------------------------------------------------------
! using legendre polynomial representation
!-------------------------------------------------------------------------
     LEG_BLOCK: if ( isort == 2 ) then
         step = real(legrd - 1) / two
         do i=1,norbs
             do j=1,ntime
                 raux = two * tmesh(j) / beta
                 curr = nint(raux * step) + 1
                 do fleg=1,lemax
                     raux = sqrt(two * fleg - 1) / (beta * beta) * rep_l(curr,fleg)
                     gaux(j,i,i) = gaux(j,i,i) + raux * gtau(fleg,i,i)
                 enddo ! over fleg={1,lemax} loop
             enddo ! over j={1,ntime} loop
         enddo ! over i={1,norbs} loop
     endif LEG_BLOCK ! back if ( isort == 2 ) block
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

     return
  end subroutine ctqmc_make_gtau

!!
!! @sub ctqmc_make_ftau
!!
!! build auxiliary correlation function using different representation
!!
  subroutine ctqmc_make_ftau(tmesh, ftau, faux)
     use constants, only : dp, zero, two

     use control, only : isort
     use control, only : norbs
     use control, only : lemax, legrd
     use control, only : ntime
     use control, only : beta

     use context, only : rep_l

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

! loop index for legendre polynomial
     integer  :: fleg

! index for imaginary time \tau
     integer  :: curr

! interval for imaginary time slice
     real(dp) :: step

! dummy variables
     real(dp) :: raux

! initialize faux
     faux = zero

!-------------------------------------------------------------------------
! using normal representation
!-------------------------------------------------------------------------
     STD_BLOCK: if ( isort == 1 ) then
         raux = real(ntime) / (beta * beta)
         do i=1,norbs
             do j=1,ntime
                 faux(j,i,i) = ftau(j,i,i) * raux
             enddo ! over j={1,ntime} loop
         enddo ! over i={1,norbs} loop
     endif STD_BLOCK ! back if ( isort == 1 ) block
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

!-------------------------------------------------------------------------
! using legendre polynomial representation
!-------------------------------------------------------------------------
     LEG_BLOCK: if ( isort == 2 ) then
         step = real(legrd - 1) / two
         do i=1,norbs
             do j=1,ntime
                 raux = two * tmesh(j) / beta
                 curr = nint(raux * step) + 1
                 do fleg=1,lemax
                     raux = sqrt(two * fleg - 1) / (beta * beta) * rep_l(curr,fleg)
                     faux(j,i,i) = faux(j,i,i) + raux * ftau(fleg,i,i)
                 enddo ! over fleg={1,lemax} loop
             enddo ! over j={1,ntime} loop
         enddo ! over i={1,norbs} loop
     endif LEG_BLOCK ! back if ( isort == 2 ) block
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

     return
  end subroutine ctqmc_make_ftau

!!========================================================================
!!>>> two-particle green's function                                    <<<
!!========================================================================

!!
!! @sub ctqmc_make_prod
!!
!! calculate product of matsubara frequency exponents exp(i \omega_n \tau)
!!
  subroutine ctqmc_make_prod(flvr, nfaux, mrank, caux1, caux2)
     use constants, only : dp, two, pi, czi

     use control, only : nffrq
     use control, only : beta

     use context, only : index_s, index_e
     use context, only : time_s, time_e
     use context, only : rank

     implicit none

! external arguments
! current flavor channel
     integer, intent(in) :: flvr

! combination of nffrq and nbfrq
     integer, intent(in) :: nfaux

! maximum number of operators in different flavor channels
     integer, intent(in) :: mrank

! matsubara frequency exponents for creation operators
     complex(dp), intent(out) :: caux1(nfaux,mrank)

! matsubara frequency exponents for annihilation operators
     complex(dp), intent(out) :: caux2(nfaux,mrank)

! local variables
! loop indices for start and end points
     integer  :: is
     integer  :: ie

! imaginary time for start and end points
     real(dp) :: taus
     real(dp) :: taue

! for creation operators
     do is=1,rank(flvr)
         taus = time_s( index_s(is, flvr), flvr )
         caux1(:,is) = exp(-two * czi * pi * taus / beta)
         call s_cumprod_z(nfaux, caux1(:,is), caux1(:,is))
         caux1(:,is) = caux1(:,is) * exp(+(nffrq + 1) * czi * pi * taus / beta)
     enddo ! over is={1,rank(flvr)} loop

! for annihilation operators
     do ie=1,rank(flvr)
         taue = time_e( index_e(ie, flvr), flvr )
         caux2(:,ie) = exp(+two * czi * pi * taue / beta)
         call s_cumprod_z(nfaux, caux2(:,ie), caux2(:,ie))
         caux2(:,ie) = caux2(:,ie) * exp(-(nffrq + 1) * czi * pi * taue / beta)
     enddo ! over ie={1,rank(flvr)} loop

     return
  end subroutine ctqmc_make_prod

!!========================================================================
!!>>> self-energy function                                             <<<
!!========================================================================

!!
!! @sub ctqmc_make_hub2
!!
!! first of all, build impurity green's function and auxiliary correlation
!! function via fast fourier transformation (if isort == 1) or analytical
!! formula (if isort == 2). and then, the self-energy function is obtained
!! by using the improved estimator trick
!!
  subroutine ctqmc_make_hub2()
     use constants, only : dp, zero, one, two, pi, czi, czero

     use control, only : isort
     use control, only : norbs
     use control, only : lemax
     use control, only : mfreq
     use control, only : nfreq
     use control, only : ntime
     use control, only : beta

     use context, only : tmesh
     use context, only : gtau, ftau
     use context, only : grnf, frnf
     use context, only : sig2

     implicit none

! local variables
! loop index
     integer  :: i
     integer  :: j
     integer  :: k

! dummy real variables
     real(dp) :: ob

! spherical Bessel functions
     real(dp) :: jaux(mfreq,lemax)

! imaginary time green's function
     real(dp) :: gaux(ntime,norbs,norbs)

! imaginary time auxiliary correlation function
     real(dp) :: faux(ntime,norbs,norbs)

! unitary transformation matrix for legendre orthogonal polynomial
     complex(dp) :: taux(mfreq,lemax)

! used to backup the sampled impurity green's function
     complex(dp) :: gtmp(nfreq,norbs,norbs)

! task 1: backup the sampled impurity green's function
!-------------------------------------------------------------------------
     gtmp = grnf(1:nfreq,:,:)

! task 2: build impurity green's function and auxiliary correlation function
!-------------------------------------------------------------------------
! using fast fourier transformation
     STD_BLOCK: if ( isort == 1 ) then

         call ctqmc_make_gtau(tmesh, gtau, gaux)
         call ctqmc_four_htau(gaux, grnf)
         call ctqmc_make_ftau(tmesh, ftau, faux)
         call ctqmc_four_htau(faux, frnf)

     endif STD_BLOCK ! back if ( isort == 1 ) block

! task 3: build impurity green's function and auxiliary correlation function
!-------------------------------------------------------------------------
! special consideration must be taken for legendre representation, we can
! calculate grnf and frnf directly by using legendre coefficients, instead
! of performing fourier transformation
     LEG_BLOCK: if ( isort == 2 ) then

! 3.1 build spherical Bessel functions: jaux
         jaux = zero
         do k=1,mfreq
             ob = (two * k - one) * pi / two
             call s_sph_jl(lemax-1, ob, jaux(k,:))
         enddo ! over k={1,mfreq} loop

! 3.2 build unitary transformation matrix: taux
         taux = czero
         do i=1,lemax
             do k=1,mfreq
                 ob = (-one)**(k - 1) * sqrt(two * i - one)
                 taux(k,i) = jaux(k,i) * ob * ( czi**i )
             enddo ! over k={1,mfreq} loop
         enddo ! over i={1,lemax} loop
         taux = taux / beta

! 3.3 rebuild impurity green's function on matsubara frequency
!     using orthogonal polynomial representation, G(i\omega)
!
! 3.4 rebuild auxiliary correlation function on matsubara frequency
!     using orthogonal polynomial representation, F(i\omega)
         grnf = czero
         frnf = czero
         do i=1,norbs
             do j=1,lemax
                 do k=1,mfreq
                     grnf(k,i,i) = grnf(k,i,i) + taux(k,j) * gtau(j,i,i)
                     frnf(k,i,i) = frnf(k,i,i) + taux(k,j) * ftau(j,i,i)
                 enddo ! over k={1,mfreq} loop
             enddo ! over j={1,lemax} loop
         enddo ! over i={1,norbs} loop

     endif LEG_BLOCK ! back if ( isort == 2 ) block

! task 4: build final self-energy function by using improved estimator
!-------------------------------------------------------------------------
     do i=1,norbs
         do k=1,mfreq
             sig2(k,i,i) = frnf(k,i,i) / grnf(k,i,i)
         enddo ! over k={1,nfreq} loop
     enddo ! over i={1,norbs} loop

! task 5: restore the sampled impurity green's function
!-------------------------------------------------------------------------
     grnf(1:nfreq,:,:) = gtmp(1:nfreq,:,:)

     return
  end subroutine ctqmc_make_hub2
