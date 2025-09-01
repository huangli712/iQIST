!!!-----------------------------------------------------------------------
!!! project : iqist @ narcissus
!!! program : ctqmc_four_htau
!!!           ctqmc_four_hybf
!!!           ctqmc_eval_htau
!!!           ctqmc_eval_hsed
!!!           ctqmc_eval_ktau
!!!           ctqmc_eval_ksed
!!!           ctqmc_symm_gtau
!!!           ctqmc_symm_grnf
!!!           ctqmc_symm_nimp
!!!           ctqmc_tran_gtau
!!!           ctqmc_tran_grnf
!!!           ctqmc_tran_twop
!!!           ctqmc_make_fock
!!!           ctqmc_make_umat
!!!           ctqmc_make_lift
!!!           ctqmc_make_iret
!!!           ctqmc_make_pref
!!!           ctqmc_make_fexp
!!!           ctqmc_make_bexp
!!!           ctqmc_make_hub2
!!! source  : ctqmc_util.f90
!!! type    : functions & subroutines
!!! author  : li huang (email:huangli@caep.cn)
!!! history : 10/01/2008 by li huang (created)
!!!           05/24/2025 by li huang (last modified)
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
     use constants, only : dp
     use constants, only : zero
     use constants, only : czero

     use control, only : norbs
     use control, only : mfreq
     use control, only : ntime

     use context, only : tmesh, rmesh

     implicit none

!! external arguments
     ! hybridization function on imaginary time axis
     real(dp), intent(in) :: htau(ntime,norbs,norbs)

     ! hybridization function on matsubara frequency axis
     complex(dp), intent(out) :: hybf(mfreq,norbs,norbs)

!! local variables
     ! loop index over orbitals
     integer  :: i
     integer  :: j

     ! dummy arrays
     real(dp) :: raux(ntime)
     complex(dp) :: caux(mfreq)

!! [body

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

!! body]

     return
  end subroutine ctqmc_four_htau

!!
!! @sub ctqmc_four_hybf
!!
!! fourier hybf to htau, from matsubara frequency to imaginary time
!!
  subroutine ctqmc_four_hybf(hybf, htau)
     use constants, only : dp
     use constants, only : zero
     use constants, only : czero
     use constants, only : eps6

     use control, only : norbs
     use control, only : mfreq
     use control, only : ntime
     use control, only : beta

     use context, only : tmesh, rmesh

     implicit none

!! external arguments
     ! hybridization function on imaginary time axis
     real(dp), intent(out) :: htau(ntime,norbs,norbs)

     ! hybridization function on matsubara frequency axis
     complex(dp), intent(in) :: hybf(mfreq,norbs,norbs)

!! local variables
     ! loop index over orbitals
     integer  :: i
     integer  :: j

     ! used to determine the bottom region of hybridiaztion function
     integer  :: start
     integer  :: last

     ! dummy arrays
     real(dp) :: raux(ntime)
     complex(dp) :: caux(mfreq)

!! [body

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

     ! checks for diagonal htau to be causal. htau should be concave.
     ! hence, if it becomes very small at two points, it should remain
     ! zero in all points between the two points. this is very important
     ! in insulators, because htau can overshoot to positive values
     ! multiple times and kinks can be trapped in the range between the
     ! two crossing points, where htau is causal, but should be zero.
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

!<!-----------------------------------------------------------------------
!<       if ( start > 1 .and. last > 1 ) then
!<           do j=start,last
!<               htau(j,i,i) = -eps6
!<           enddo ! over j={start,last} loop
!<       endif ! back if ( start > 1 .and. last > 1 ) block
!<!-----------------------------------------------------------------------
     enddo ! over i={1,norbs} loop

     ! enforce hybridization function less than zero to ensure the causality
     do i=1,norbs
         do j=1,ntime
             if ( htau(j,i,i) > zero ) htau(j,i,i) = -eps6
         enddo ! over j={1,ntime} loop
     enddo ! over i={1,norbs} loop

!! body]

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

!! external arguments
     ! current flavor channel
     integer, intent(in)  :: flvr

     ! delta imaginary time
     real(dp), intent(in) :: dtau

!! external functions
     ! internal interpolation engine
     procedure( real(dp) ) :: s_spl_funct

!! local variables
     ! return value
     real(dp) :: val

!! [body

     val = s_spl_funct(ntime, tmesh, htau(:, flvr, flvr), hsed(:, flvr, flvr), dtau)

!! body]

     return
  end function ctqmc_eval_htau

!!
!! @sub ctqmc_eval_hsed
!!
!! calculate the second order derivates of hybridization function on
!! imaginary time space
!!
  subroutine ctqmc_eval_hsed(htau, hsed)
     use constants, only : dp
     use constants, only : zero

     use control, only : norbs
     use control, only : ntime
     use control, only : beta

     use context, only : tmesh

     implicit none

!! external arguments
     ! hybridization function on imaginary time axis
     real(dp), intent(in)  :: htau(ntime,norbs,norbs)

     ! second order derivates of hybridization function
     real(dp), intent(out) :: hsed(ntime,norbs,norbs)

!! local variables
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

!! [body

     ! calculate deltau
     deltau = beta / real(ntime - 1)

     ! initialize hsed
     hsed = zero

     ! calculate it
     do j=1,norbs
         do i=1,norbs

             ! calculate first-order derivate of \Delta(0): startu
             startu = (-25.0_dp*htau(1,       i, j) + &
                        48.0_dp*htau(2,       i, j) - &
                        36.0_dp*htau(3,       i, j) + &
                        16.0_dp*htau(4,       i, j) - &
                         3.0_dp*htau(5,       i, j)) / 12.0_dp / deltau

             ! calculate first-order derivate of \Delta(\beta): startd
             startd = ( 25.0_dp*htau(ntime-0, i, j) - &
                        48.0_dp*htau(ntime-1, i, j) + &
                        36.0_dp*htau(ntime-2, i, j) - &
                        16.0_dp*htau(ntime-3, i, j) + &
                         3.0_dp*htau(ntime-4, i, j)) / 12.0_dp / deltau

             ! reinitialize d2y to zero
             d2y = zero

             ! call the service layer
             call s_spl_deriv2(ntime, tmesh, htau(:,i,j), startu, startd, d2y)

             ! copy the results to hsed
             hsed(:,i,j) = d2y

         enddo ! over i={1,norbs} loop
     enddo ! over j={1,norbs} loop

!! body]

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

!! external arguments
     ! order for derivates
     ! if mode = 1, K(\tau), ktau is considered
     ! if mode = 2, K'(\tau), ptau is considered
     integer, intent(in)  :: mode

     ! current imaginary time
     real(dp), intent(in) :: dtau

!! external functions
     ! internal interpolation engine
     procedure( real(dp) ) :: s_spl_funct

!! local variables
     ! return value
     real(dp) :: val

!! [body

     ! using cubic spline interpolation for K(\tau)
     if ( mode == 1 ) then
         val = s_spl_funct(ntime, tmesh, ktau, ksed, dtau)
     !
     ! using cubic spline interpolation for K'(\tau)
     else
         val = s_spl_funct(ntime, tmesh, ptau, psed, dtau)
     !
     endif ! back if ( mode == 1 ) block

!! body]

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
  subroutine ctqmc_eval_ksed(ktau, ksed)
     use constants, only : dp
     use constants, only : zero

     use control, only : ntime
     use control, only : beta

     use context, only : tmesh

     implicit none

!! external arguments
     ! screening function on imaginary time axis
     real(dp), intent(in)  :: ktau(ntime)

     ! second order derivates of screening function
     real(dp), intent(out) :: ksed(ntime)

!! local variables
     ! first derivate at start point
     real(dp) :: startu

     ! first derivate at end   point
     real(dp) :: startd

     ! \delta \tau
     real(dp) :: deltau

!! [body

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

!! body]

     return
  end subroutine ctqmc_eval_ksed

!!========================================================================
!!>>> symmetry operation                                               <<<
!!========================================================================

!!
!! @sub ctqmc_symm_gtau
!!
!! symmetrize the gtau according to symm vector. only the diagonal terms
!! are taken into considerations
!!
  subroutine ctqmc_symm_gtau(symm, gtau)
     use constants, only : dp
     use constants, only : zero, two

     use control, only : isbnd, isspn
     use control, only : nband, norbs
     use control, only : ntime

     implicit none

!! external arguments
     ! symmetry vector
     integer, intent(in) :: symm(norbs)

     ! impurity green's function
     real(dp), intent(inout) :: gtau(ntime,norbs,norbs)

!! local variables
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

!! [body

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

!! body]

     return
  end subroutine ctqmc_symm_gtau

!!
!! @sub ctqmc_symm_grnf
!!
!! symmetrize the grnf according to symm vector. only the diagonal terms
!! are taken into considerations
!!
  subroutine ctqmc_symm_grnf(symm, grnf)
     use constants, only : dp
     use constants, only : two, czero

     use control, only : isbnd, isspn
     use control, only : nband, norbs
     use control, only : mfreq

     implicit none

!! external arguments
     ! symmetry vector
     integer, intent(in) :: symm(norbs)

     ! impurity green's function
     complex(dp), intent(inout) :: grnf(mfreq,norbs,norbs)

!! local variables
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

!! [body

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

!! body]

     return
  end subroutine ctqmc_symm_grnf

!!
!! @sub ctqmc_symm_nimp
!!
!! symmetrize the occupation number array, nimp, according to symm vector
!!
  subroutine ctqmc_symm_nimp(symm, nimp)
     use constants, only : dp
     use constants, only : zero, two

     use control, only : isbnd, isspn
     use control, only : nband, norbs

     implicit none

!! external arguments
     ! symmetry vector
     integer, intent(in) :: symm(norbs)

     ! occupation number
     real(dp), intent(inout) :: nimp(norbs)

!! local variables
     ! loop index over bands
     integer  :: ibnd
     integer  :: jbnd

     ! dummy variables
     real(dp) :: raux

     ! histogram vector
     ! note: it is NOT the global one
     integer  :: hist(norbs)

!! [body

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

!! body]

     return
  end subroutine ctqmc_symm_nimp

!!========================================================================
!!>>> advanced representation                                          <<<
!!========================================================================

!!
!! @sub ctqmc_tran_gtau
!!
!! build imaginary time green's function using different representation
!!
  subroutine ctqmc_tran_gtau(gaux, gtau)
     use constants, only : dp
     use constants, only : zero, one, two

     use control, only : isort
     use control, only : norbs
     use control, only : lemax, legrd
     use control, only : svmax, svgrd
     use control, only : ntime
     use control, only : beta

     use context, only : tmesh
     use context, only : rep_l, rep_s

     implicit none

!! external arguments
     ! impurity green's function/orthogonal polynomial coefficients
     real(dp), intent(in)  :: gaux(ntime,norbs,norbs)

     ! calculated impurity green's function
     real(dp), intent(out) :: gtau(ntime,norbs,norbs)

!! local variables
     ! loop index
     integer  :: i
     integer  :: j

     ! loop index for orthogonal polynomial
     integer  :: fleg
     integer  :: fsvd

     ! index for imaginary time \tau
     integer  :: curr

     ! interval for imaginary time slice
     real(dp) :: step

     ! dummy variables
     real(dp) :: raux

!! [body

     ! initialize gtau
     gtau = zero

     !--------------------------------------------------------------------
     ! using normal representation
     !--------------------------------------------------------------------
     STD_BLOCK: if ( isort == 1 ) then
         raux = real(ntime) / (beta * beta)
         do i=1,norbs
             do j=1,ntime
                 gtau(j,i,i) = gaux(j,i,i) * raux
             enddo ! over j={1,ntime} loop
         enddo ! over i={1,norbs} loop
     endif STD_BLOCK ! back if ( isort == 1 ) block
     !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

     !--------------------------------------------------------------------
     ! using legendre orthogonal polynomial representation
     !--------------------------------------------------------------------
     ! see Eq. (1) and Eq. (C19) in Phys. Rev. B 84, 075145 (2011)
     LEG_BLOCK: if ( isort == 2 ) then
         step = real(legrd - 1) / two
         do i=1,norbs
             do j=1,ntime
                 raux = two * tmesh(j) / beta ! map tmesh to [0,2]
                 curr = nint(raux * step) + 1
                 do fleg=1,lemax
                     raux = sqrt(two * fleg - 1) / (beta * beta) * rep_l(curr,fleg)
                     gtau(j,i,i) = gtau(j,i,i) + raux * gaux(fleg,i,i)
                 enddo ! over fleg={1,lemax} loop
             enddo ! over j={1,ntime} loop
         enddo ! over i={1,norbs} loop
     endif LEG_BLOCK ! back if ( isort == 2 ) block
     !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

     !--------------------------------------------------------------------
     ! using svd orthogonal polynomial representation
     !--------------------------------------------------------------------
     ! see Eq. (14) in Phys. Rev. B 96, 035147 (2017)
     SVD_BLOCK: if ( isort == 3 ) then
         step = real(svgrd - 1) / two
         do i=1,norbs
             do j=1,ntime
                 raux = two * tmesh(j) / beta - one ! map tmesh to [-1,1]
                 call s_svd_point(raux, step, curr)
                 do fsvd=1,svmax
                     raux = two / (beta * beta) * rep_s(curr,fsvd)
                     gtau(j,i,i) = gtau(j,i,i) + raux * gaux(fsvd,i,i)
                 enddo ! over fsvd={1,svmax} loop
             enddo ! over j={1,ntime} loop
         enddo ! over i={1,norbs} loop
     endif SVD_BLOCK ! back if ( isort == 3 ) block
     !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

!! body]

     return
  end subroutine ctqmc_tran_gtau

!!
!! @sub ctqmc_tran_grnf
!!
!! build matsubara green's function using different representation
!!
  subroutine ctqmc_tran_grnf(gaux, grnf)
     use constants, only : dp
     use constants, only : zero, one, two
     use constants, only : czero, czi

     use mmpi, only : mp_allreduce
     use mmpi, only : mp_barrier

     use control, only : isort
     use control, only : norbs
     use control, only : lemax
     use control, only : svmax, svgrd
     use control, only : mfreq
     use control, only : ntime
     use control, only : beta
     use control, only : myid, nprocs

     use context, only : tmesh, rmesh
     use context, only : rep_s

     implicit none

!! external arguments
     ! orthogonal polynomial coefficients for impurity green's function
     real(dp), intent(in) :: gaux(ntime,norbs,norbs)

     ! calculated impurity green's function
     complex(dp), intent(out) :: grnf(mfreq,norbs,norbs)

!! local variables
     ! loop index
     integer  :: i
     integer  :: j
     integer  :: k

     ! index for imaginary time \tau
     integer  :: curr

     ! status flag
     integer  :: istat

     ! dummy real(dp) variable
     real(dp) :: raux

     ! step for the linear frequency mesh
     real(dp) :: step

     ! j_n(x), for legendre orthogonal polynomial representation
     real(dp), allocatable :: pfun(:,:)

     ! u_l(x(\tau)), for svd orthogonal polynomial representation
     real(dp), allocatable :: ufun(:,:)

     ! calculated impurity green's function, imaginary time axis
     real(dp), allocatable :: gtau(:,:,:)

     ! unitary transformation matrix for orthogonal polynomials
     complex(dp), allocatable :: tleg(:,:)
     complex(dp), allocatable :: tsvd(:,:)
     complex(dp), allocatable :: tmpi(:,:)

!! [body

     ! allocate memory
     allocate(pfun(mfreq,lemax), stat=istat)
     allocate(ufun(ntime,svmax), stat=istat)

     allocate(gtau(ntime,norbs,norbs), stat=istat)

     allocate(tleg(mfreq,lemax), stat=istat)
     allocate(tsvd(mfreq,svmax), stat=istat)
     allocate(tmpi(mfreq,svmax), stat=istat)

     if ( istat /= 0 ) then
         call s_print_error('ctqmc_tran_grnf','can not allocate enough memory')
     endif ! back if ( istat /= 0 ) block

     !--------------------------------------------------------------------
     ! using normal representation
     !--------------------------------------------------------------------
     STD_BLOCK: if ( isort == 1 ) then
         call ctqmc_tran_gtau(gaux, gtau)
         call ctqmc_four_htau(gtau, grnf)
     endif STD_BLOCK ! back if ( isort == 1 ) block
     !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

     !--------------------------------------------------------------------
     ! using legendre orthogonal polynomial representation
     !--------------------------------------------------------------------
     ! see Eq. (5) in Phys. Rev. B 84, 075145 (2011)
     LEG_BLOCK: if ( isort == 2 ) then

         ! calculate spherical Bessel functions at first
         pfun = zero
         do k=1,mfreq
             call s_sph_jl(lemax-1, rmesh(k) * beta / two, pfun(k,:))
         enddo ! over k={1,mfreq} loop

         ! build unitary transformation matrix: tleg
         tleg = czero
         do i=1,lemax
             raux = sqrt(two * i - one)
             do k=1,mfreq
                 tleg(k,i) = pfun(k,i) * (-one)**(k - 1) * raux * czi**i
             enddo ! over k={1,mfreq} loop
         enddo ! over i={1,lemax} loop

         ! normalize tleg
         ! note: the beta is from Eq. (C19) in Phys. Rev. B 84, 075145 (2011)
         tleg = tleg / beta

         ! build impurity green's function on matsubara frequency using
         ! orthogonal polynomial representation: grnf
         grnf = czero
         do i=1,norbs
             do j=1,lemax
                 do k=1,mfreq
                     grnf(k,i,i) = grnf(k,i,i) + tleg(k,j) * gaux(j,i,i)
                 enddo ! over k={1,mfreq} loop
             enddo ! over j={1,lemax} loop
         enddo ! over i={1,norbs} loop

     endif LEG_BLOCK ! back if ( isort == 2 ) block
     !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

     !--------------------------------------------------------------------
     ! using svd orthogonal polynomial representation
     !--------------------------------------------------------------------
     SVD_BLOCK: if ( isort == 3 ) then

         ! copy rep_s to ufun, prepare u_l(x(\tau))
         step = real(svgrd - 1) / two
         do i=1,ntime
             raux = two * tmesh(i) / beta - one
             call s_svd_point(raux, step, curr)
             ufun(i,:) = rep_s(curr,:)
         enddo ! over i={1,ntime} loop

         ! build unitary transformation matrix: tsvd
         ! actually, we do the fourier transformation
         tmpi = czero
         do i=1+myid,svmax,nprocs
             call s_fft_forward(ntime, tmesh, ufun(:,i), mfreq, rmesh, tmpi(:,i))
         enddo ! over i={1+myid,svmax} loop

! build tsvd, collect data from children processes
# if defined (MPI)

         ! collect data
         call mp_allreduce(tmpi, tsvd)

         ! block until all processes have reached here
         call mp_barrier()

# else  /* MPI */

         tsvd = tmpi

# endif /* MPI */

         ! normalize tsvd
         ! note: the first beta is from Eq. (C19), while the second beta
         ! is from Eq. (E1) in Phys. Rev. B 84, 075145 (2011)
         tsvd = tsvd * (two / beta / beta)

         ! build impurity green's function on matsubara frequency using
         ! orthogonal polynomial representation: grnf
         grnf = czero
         do i=1,norbs
             do j=1,svmax
                 do k=1,mfreq
                     grnf(k,i,i) = grnf(k,i,i) + tsvd(k,j) * gaux(j,i,i)
                 enddo ! over k={1,mfreq} loop
             enddo ! over j={1,svmax} loop
         enddo ! over i={1,norbs} loop

     endif SVD_BLOCK ! back if ( isort == 3 ) block
     !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

!! body]

     ! deallocate memory
     deallocate(pfun)
     deallocate(ufun)
     deallocate(gtau)
     deallocate(tleg)
     deallocate(tsvd)
     deallocate(tmpi)

     return
  end subroutine ctqmc_tran_grnf

!!
!! @sub ctqmc_tran_twop
!!
!! build two-particle green's function using different representation
!!
  subroutine ctqmc_tran_twop(gaux, grnf)
     use constants, only : dp
     use constants, only : one, two
     use constants, only : czero

     use mmpi, only : mp_allreduce
     use mmpi, only : mp_barrier

     use control, only : isort
     use control, only : norbs
     use control, only : lemax, legrd
     use control, only : svmax, svgrd
     use control, only : nffrq, nbfrq
     use control, only : ntime
     use control, only : beta
     use control, only : myid, nprocs

     use context, only : tmesh, rmesh
     use context, only : rep_l, rep_s

     implicit none

!! external arguments
     ! orthogonal polynomial coefficients for two-particle green's function
     complex(dp), intent(in)  :: gaux(nffrq,nffrq,nbfrq,norbs,norbs)

     ! calculated two-particle green's function
     complex(dp), intent(out) :: grnf(nffrq,nffrq,nbfrq,norbs,norbs)

!! local variables
     ! loop index
     integer  :: i
     integer  :: j
     integer  :: k
     integer  :: l

     ! index for imaginary time \tau
     integer  :: curr

     ! status flag
     integer  :: istat

     ! dummy real(dp) variable
     real(dp) :: raux

     ! step for the linear frequency mesh
     real(dp) :: step

     ! symmetric fermionic matsubara frequency mesh
     real(dp), allocatable :: fmesh(:)

     ! p_l(x(\tau)), for legendre orthogonal polynomial representation
     real(dp), allocatable :: pfun(:,:)

     ! u_l(x(\tau)), for svd orthogonal polynomial representation
     real(dp), allocatable :: ufun(:,:)

     ! unitary transformation matrix for orthogonal polynomials
     complex(dp), allocatable :: tleg(:,:)
     complex(dp), allocatable :: tsvd(:,:)
     complex(dp), allocatable :: tmpi(:,:)

!! [body

     ! allocate memory
     allocate(fmesh(nffrq),      stat=istat)
     allocate(pfun(ntime,lemax), stat=istat)
     allocate(ufun(ntime,svmax), stat=istat)
     allocate(tleg(nffrq,lemax), stat=istat)
     allocate(tsvd(nffrq,svmax), stat=istat)

     if ( istat /= 0 ) then
         call s_print_error('ctqmc_tran_twop','can not allocate enough memory')
     endif ! back if ( istat /= 0 ) block

     ! build symmetric fermionic matsubara frequency mesh
     do i=nffrq/2+1,nffrq
         fmesh(i) = rmesh(i-nffrq/2)  ! > 0
         fmesh(nffrq-i+1) = -fmesh(i) ! < 0
     enddo ! over i={nffrq/2+1,nffrq} loop

     !--------------------------------------------------------------------
     ! using normal representation
     !--------------------------------------------------------------------
     STD_BLOCK: if ( isort == 1 ) then
         allocate(tmpi(  1  ,  1  ), stat=istat); tmpi = czero
         grnf = gaux
     endif STD_BLOCK ! back if ( isort == 1 ) block
     !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

     !--------------------------------------------------------------------
     ! using legendre orthogonal polynomial representation
     !--------------------------------------------------------------------
     LEG_BLOCK: if ( isort == 2 ) then

         ! copy rep_l to pfun, prepare p_l(x(\tau))
         step = real(legrd - 1) / two
         do i=1,ntime
             raux = two * tmesh(i) / beta
             curr = nint( raux * step ) + 1
             pfun(i,:) = rep_l(curr,:)
         enddo ! over i={1,ntime} loop

         ! build unitary transformation matrix: tleg
         ! we do the fourier transformation directly using Eq. (E1) in
         ! Phys. Rev. B 84, 075145 (2011). the advantage is that it
         ! doesn't depend on the spherical Bessel functions any more
         allocate(tmpi(nffrq,lemax), stat=istat); tmpi = czero
         do i=1+myid,lemax,nprocs
             call s_fft_forward(ntime, tmesh, pfun(:,i), nffrq, fmesh, tmpi(:,i))
             tmpi(:,i) = tmpi(:,i) * sqrt(two * i - one)
         enddo ! over i={1+myid,lemax} loop

! build tleg, collect data from children processes
# if defined (MPI)

         ! collect data
         call mp_allreduce(tmpi, tleg)

         ! block until all processes have reached here
         call mp_barrier()

# else  /* MPI */

         tleg = tmpi

# endif /* MPI */

         ! normalize tleg
         ! note: the beta is from Eq. (E1) in Phys. Rev. B 84, 075145 (2011)
         tleg = tleg / beta

         ! build two-particle green's function on matsubara frequency
         ! using orthogonal polynomial representation: grnf
         ! see Eq. (14) in Phys. Rev. B 84, 075145 (2011)
         grnf = czero
         do i=1,nffrq             ! for v' index
             do j=1,nffrq         ! for v  index
                 do k=1,lemax     ! for l' index
                     do l=1,lemax ! for l  index
                         associate ( val => grnf(i,j,:,:,:), &
                                     gkl => gaux(k,l,:,:,:) )
                             val = val + tleg(j,l) * gkl * conjg( tleg(i,k) )
                         end associate
                     enddo ! over l={1,lemax} loop
                 enddo ! over k={1,lemax} loop
             enddo ! over j={1,nffrq} loop
         enddo ! over i={1,nffrq} loop

     endif LEG_BLOCK ! back if ( isort == 2 ) block
     !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

     !--------------------------------------------------------------------
     ! using svd orthogonal polynomial representation
     !--------------------------------------------------------------------
     SVD_BLOCK: if ( isort == 3 ) then

         ! copy rep_s to ufun, prepare u_l(x(\tau))
         step = real(svgrd - 1) / two
         do i=1,ntime
             raux = two * tmesh(i) / beta - one
             call s_svd_point(raux, step, curr)
             ufun(i,:) = rep_s(curr,:)
         enddo ! over i={1,ntime} loop

         ! build unitary transformation matrix: tsvd
         ! actually, we do the fourier transformation
         allocate(tmpi(nffrq,svmax), stat=istat); tmpi = czero
         do i=1+myid,svmax,nprocs
             call s_fft_forward(ntime, tmesh, ufun(:,i), nffrq, fmesh, tmpi(:,i))
         enddo ! over i={1+myid,svmax} loop

! build tsvd, collect data from children processes
# if defined (MPI)

         ! collect data
         call mp_allreduce(tmpi, tsvd)

         ! block until all processes have reached here
         call mp_barrier()

# else  /* MPI */

         tsvd = tmpi

# endif /* MPI */

         ! normalize tsvd
         tsvd = tsvd * (two / beta)

         ! build two-particle green's function on matsubara frequency
         ! using svd orthogonal polynomial representation: grnf
         grnf = czero
         do i=1,nffrq             ! for v' index
             do j=1,nffrq         ! for v  index
                 do k=1,svmax     ! for l' index
                     do l=1,svmax ! for l  index
                         associate ( val => grnf(nffrq-i+1,j,:,:,:), &
                                     gkl => gaux(k,l,:,:,:) )
                             val = val + tsvd(j,l) * gkl * conjg( tsvd(i,k) )
                         end associate
                     enddo ! over l={1,svmax} loop
                 enddo ! over k={1,svmax} loop
             enddo ! over j={1,nffrq} loop
         enddo ! over i={1,nffrq} loop

     endif SVD_BLOCK ! back if ( isort == 3 ) block
     !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

     ! deallocate memory
     deallocate(fmesh)
     deallocate(pfun)
     deallocate(ufun)
     deallocate(tleg)
     deallocate(tsvd)
     deallocate(tmpi)

!! body]

     return
  end subroutine ctqmc_tran_twop

!!========================================================================
!!>>> atomic eigenstates                                               <<<
!!========================================================================

!!
!! @sub ctqmc_make_fock
!!
!! convert current atomic eigenstate into a decimal number (state index)
!!
  subroutine ctqmc_make_fock(norbs, pstat, state)
     implicit none

!! external arguments
     ! index of atomic state
     integer, intent(out) :: pstat

     ! number of orbitals
     integer, intent(in)  :: norbs

     ! atomic state array
     integer, intent(in)  :: state(norbs)

!! local variables
     ! loop index
     integer :: i

!! [body

     ! init pstat
     pstat = 1

     ! evaluate pstat, for example, 0101 = 0*2^0 + 1*2^1 + 0*2^2 + 1*2^3 = 10
     do i=1,norbs
         if ( state(i) > 0 ) pstat = pstat + ishft(1, i-1)
     enddo ! over i={1,norbs} loop

!! body]

     return
  end subroutine ctqmc_make_fock

!!========================================================================
!!>>> Coulomb interaction matrix                                       <<<
!!========================================================================

!!
!! note:
!!
!! since the narcissus code is based on the segment representation, it
!! does not support the spin-flip and pair-hopping terms of course. only
!! the Uc and Jz parameters are need to build the interaction matrix
!!

!!
!! @sub ctqmc_make_umat
!!
!! build density-density two-fermions Coulomb interaction matrix: umat.
!! here the used equation is Eq. (13) in Rev. Mod. Phys. 83, 349 (2011)
!!
  subroutine ctqmc_make_umat(umat)
     use constants, only : dp
     use constants, only : zero

     use control, only : nband, norbs
     use control, only : Uc, Jz

     implicit none

!! external arguments
     ! Coulomb interaction matrix
     real(dp), intent(out) :: umat(norbs,norbs)

!! local variables
     ! loop index
     integer  :: i
     integer  :: j
     integer  :: k
     integer  :: m

     ! dummy u vector
     real(dp) :: ut(nband*(norbs-1))

!! [body

     ! initialize it
     umat = zero

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

             umat(i,j) = ut(k)
             umat(j,i) = ut(k)
         enddo ! over j={i+1,norbs} loop
     enddo ! over i={1,norbs-1} loop

!! body]

     return
  end subroutine ctqmc_make_umat

!!========================================================================
!!>>> retarded interaction                                             <<<
!!========================================================================

!!
!! @sub ctqmc_make_lift
!!
!! shift the Coulomb interaction matrix and the chemical potential if
!! retarded interaction is considered
!!
!! for plasmon pole model and ohmic model, please refer to
!!     Phys. Rev. Lett. 104, 146401 (2010)
!!
!! for general U(\omega), see
!!     Eq. (58)-(60) in J. Phys.: Condens. Matter 28, 383001 (2016)
!! or
!!     Eq. (20)-(21) in Phys. Rev. B 89, 235128 (2014)
  subroutine ctqmc_make_lift(umat, ssign)
     use constants, only : dp
     use constants, only : zero, two

     use control, only : isscr
     use control, only : norbs
     use control, only : lc, wc
     use control, only : mune

     use context, only : ptau

     implicit none

!! external arguments
     ! Coulomb interaction matrix
     real(dp), intent(inout) :: umat(norbs,norbs)

     ! sign for the shift, it should be 1.0_dp or -1.0_dp
     real(dp), intent(in)    :: ssign

!! local variables
     ! loop index
     integer  :: i
     integer  :: j

     ! Coulomb interaction shift introduced by dynamic screening effect
     real(dp) :: shift

!! [body

     ! evaluate Coulomb interaction shift
     DYNAMIC_MODEL: select case ( isscr )

         case (1) ! static interaction
             shift = zero

         case (2) ! dynamic interaction, plasmon pole model
             shift = two * lc * lc / wc

         case (3) ! dynamic interaction, ohmic model
             shift = two * lc * wc

         case (4) ! dynamic interaction, realistic materials
             shift = two * ptau(1)

     end select DYNAMIC_MODEL

     ! multiple the shift with sign
     shift = shift * ssign

     ! shift the Coulomb interaction matrix (skip the diagonal elements)
     do i=1,norbs-1
         do j=i+1,norbs
             umat(i,j) = umat(i,j) - shift
             umat(j,i) = umat(j,i) - shift
         enddo ! over j={i+1,norbs} loop
     enddo ! over i={1,norbs-1} loop

     ! shift chemical potential as a byproduct
     mune = mune - shift / two

!! body]

     return
  end subroutine ctqmc_make_lift

!!========================================================================
!!>>> improved estimator                                               <<<
!!========================================================================

!!
!! @sub ctqmc_make_iret
!!
!! calculate the integral I(\tau_end) which is very important when the
!! retarded interaction is included. here the used equation is Eq. (39)
!! in Phys. Rev. B 89, 235128 (2014)
!!
  subroutine ctqmc_make_iret(time, iret)
     use constants, only : dp
     use constants, only : zero, two

     use control, only : norbs

     use context, only : index_s, index_e
     use context, only : time_s, time_e
     use context, only : rank

     implicit none

!! external arguments
     ! imaginary time point, in principle, it is \tau_end
     real(dp), intent(in)  :: time

     ! integral value for I(\tau_end)
     real(dp), intent(out) :: iret

!! local variables
     ! loop index for start and end points
     integer  :: it

     ! loop index for flavor channel
     integer  :: flvr

     ! length betweem two time points
     real(dp) :: dtau
     real(dp) :: daux

!! [body

     ! init integral I(\tau_end)
     iret = zero

     FLVR_CYCLE: do flvr=1,norbs
         do it=1,rank(flvr)

             ! contribution from creation operators
             dtau = time_s( index_s(it, flvr), flvr ) - time
             if ( dtau >= zero ) then
                 call cat_weight_kernel(2, +dtau, daux)
                 iret = iret + daux
             else
                 call cat_weight_kernel(2, -dtau, daux)
                 iret = iret - daux
             endif ! back if ( dtau >= zero ) block

             ! contribution from annihilation operators
             dtau = time_e( index_e(it, flvr), flvr ) - time
             if ( dtau >= zero ) then
                 call cat_weight_kernel(2, +dtau, daux)
                 iret = iret - daux
             else
                 call cat_weight_kernel(2, -dtau, daux)
                 iret = iret + daux
             endif ! back if ( dtau >= zero ) block

         enddo ! over it={1,rank(flvr)} loop
     enddo FLVR_CYCLE ! over flvr={1,norbs} loop

!<   ! add additional term: -2K'(0^+)
!<   ! note: this static contribution should already be accounted for by the
!<   ! renormalization of the static U
!<   call cat_weight_kernel(2, zero, daux)
!<   iret = -iret - two * daux
     iret = -iret

!! body]

     return
  end subroutine ctqmc_make_iret

!!
!! @sub ctqmc_make_pref
!!
!! calculate the prefactor used by the improved estimator for self-energy
!! function and two-particle green's function. when retarded interaction
!! is present, an additional term is supplemented
!!
  subroutine ctqmc_make_pref()
     use constants, only : dp
     use constants, only : zero, half

     use control, only : isscr
     use control, only : norbs

     use context, only : index_e
     use context, only : time_e
     use context, only : rank, pref, umat

     implicit none

!! local variables
     ! loop index for start and end points
     integer  :: it

     ! loop index for flavor channel
     integer  :: f1
     integer  :: f2

     ! occupation number at \tau_end
     real(dp) :: occu

     ! integral value for I(\tau_end)
     real(dp) :: iret

!! [body

     FLVR_CYCLE: do f1=1,norbs
         do it=1,rank(f1)

             ! reset the prefactor
             pref(it,f1) = zero

             ! calculate contribution from static interaction
             ! see Eq. (37) in Phys. Rev. B 89, 235128 (2014)
             do f2=1,norbs
                 call cat_occupy_status(f2, time_e( index_e(it, f1), f1 ), occu)
                 pref(it,f1) = pref(it,f1) + half * ( umat(f1,f2) + umat(f2,f1) ) * occu
             enddo ! over f2={1,norbs} loop

             ! calculate contribution from retarded (dynamic) interaction
             ! see Eq. (38) in Phys. Rev. B 89, 235128 (2014)
             if ( isscr > 1 ) then
                 call ctqmc_make_iret(time_e( index_e(it, f1), f1 ), iret)
                 pref(it,f1) = pref(it,f1) + iret
             endif ! back if ( isscr > 1 ) block

         enddo ! over it={1,rank(f1)} loop
     enddo FLVR_CYCLE ! over f1={1,norbs} loop

!! body]

     return
  end subroutine ctqmc_make_pref

!!========================================================================
!!>>> two-particle green's function                                    <<<
!!========================================================================

!!
!! @sub ctqmc_make_fexp
!!
!! calculate product of matsubara frequency exponents exp(i \omega_n \tau)
!!
!! note:
!!
!!     here, we provide two versions of ctqmc_make_fexp subroutines. the
!!     difference lies in how to evaluate exp(i \omega_n \tau). one just
!!     copies data from exp_s and exp_e, which is fast. but nfreq must be
!!     larger than nfaux (= nffrq + nbfrq - 1). another one just tries to
!!     calculate the quantity directly, which is a bit slow, but safe. we
!!     generally prefer to use the first one. but if you want to use the
!!     second one, you have to comment out the relevant codes (see below)
!!     and recompile them
!!
!! version 1
!!
  subroutine ctqmc_make_fexp(flvr, nfaux, mrank, caux1, caux2)
     use constants, only : dp

     use control, only : nfreq
     use control, only : nffrq

     use context, only : index_s, index_e
     use context, only : exp_s, exp_e
     use context, only : rank

     implicit none

!! external arguments
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

!! local variables
     ! loop indices for start and end points
     integer :: is
     integer :: ie

     ! loop index for frequency
     integer :: ix

     ! index for frequency
     integer :: ir

!! [body

     ! make sure nfreq is larger than nfaux, or else this subroutine will fail
     call s_assert2( nfreq > nfaux, 'in ctqmc_make_fexp' )

     ! creation operators
     !--------------------------------------------------------------------
     ! for each \tau_s, we try to calculate
     !     exp ( i \omega_n \tau_s ) where n \in [1,nfaux]
     !     \omega_n = -(v + w), v: -v ---> +v, w: -0 ---> +w
     ! so,
     !     \omega_n = +v,   when n = 1
     !     \omega_n = -v-w, when n = nfaux
     !
     do is=1,rank(flvr)
         do ix=1,nffrq/2
             ir = nffrq / 2 + 1 - ix
             caux1(ix,is) = exp_s(ir, index_s(is, flvr), flvr)
         enddo ! over ix={1,nffrq/2} loop
         do ix=nffrq/2+1,nfaux
             ir = nffrq / 2 + 1 - ix
             ir = abs(ir) + 1
             caux1(ix,is) = dconjg( exp_s(ir, index_s(is, flvr), flvr) )
         enddo ! over ix={nffrq/2+1,nfaux} loop
     enddo ! over is={1,rank(flvr)} loop

     ! annihilation operators
     !--------------------------------------------------------------------
     ! for each \tau_e, we try to calculate
     !     exp ( i \omega_n \tau_e ) where n \in [1,nfaux]
     !     \omega_n = +(v + w), v: -v ---> +v, w: -0 ---> +w
     ! so,
     !     \omega_n = -v,   when n = 1
     !     \omega_n = +v+w, when n = nfaux
     !
     do ie=1,rank(flvr)
         do ix=1,nffrq/2
             ir = -nffrq/2 + ix
             ir = abs(ir) + 1
             caux2(ix,ie) = dconjg( exp_e(ir, index_e(ie, flvr), flvr) )
         enddo ! over ix={1,nffrq/2} loop
         do ix=nffrq/2+1,nfaux
             ir = -nffrq/2 + ix
             caux2(ix,ie) = exp_e(ir, index_e(ie, flvr), flvr)
         enddo ! over ix={nffrq/2+1,nfaux} loop
     enddo ! over ie={1,rank(flvr)} loop

!! body]

     return
  end subroutine ctqmc_make_fexp

!!
!! @sub ctqmc_make_fexp
!!
!! calculate product of matsubara frequency exponents exp(i \omega_n \tau)
!!
!! note:
!!
!!     here, we provide two versions of ctqmc_make_fexp subroutines. the
!!     difference lies in how to evaluate exp(i \omega_n \tau). one just
!!     copies data from exp_s and exp_e, which is fast. but nfreq must be
!!     larger than nfaux (= nffrq + nbfrq - 1). another one just tries to
!!     calculate the quantity directly, which is a bit slow, but safe. we
!!     generally prefer to use the first one. but if you want to use the
!!     second one, you have to comment out the relevant codes (see below)
!!     and recompile them
!!
!! version 2
!!
!<  subroutine ctqmc_make_fexp(flvr, nfaux, mrank, caux1, caux2)
!<     use constants, only : dp
!<     use constants, only : pi, two
!<     use constants, only : czi
!<
!<     use control, only : nffrq
!<     use control, only : beta
!<
!<     use context, only : index_s, index_e
!<     use context, only : time_s, time_e
!<     use context, only : rank
!<
!<     implicit none
!<
!<!! external arguments
!<     ! current flavor channel
!<     integer, intent(in) :: flvr
!<
!<     ! combination of nffrq and nbfrq
!<     integer, intent(in) :: nfaux
!<
!<     ! maximum number of operators in different flavor channels
!<     integer, intent(in) :: mrank
!<
!<     ! matsubara frequency exponents for creation operators
!<     complex(dp), intent(out) :: caux1(nfaux,mrank)
!<
!<     ! matsubara frequency exponents for annihilation operators
!<     complex(dp), intent(out) :: caux2(nfaux,mrank)
!<
!<!! local variables
!<     ! loop indices for start and end points
!<     integer :: is
!<     integer :: ie
!<
!<     ! imaginary time for start and end points
!<     ! actually, they are i\pi\tau_s/\beta and i\pi\tau_e/\beta
!<     complex(dp) :: zs
!<     complex(dp) :: ze
!<
!<!! [body
!<
!<     ! creation operators
!<     !------------------------------------------------------------------
!<     ! for each \tau_s, we try to calculate
!<     !     exp ( i \omega_n \tau_s ) where n \in [1,nfaux]
!<     !     \omega_n = -(v + w), v: -v ---> +v, w: -0 ---> +w
!<     ! so,
!<     !     \omega_n = +v,   when n = 1
!<     !     \omega_n = -v-w, when n = nfaux
!<     !
!<     do is=1,rank(flvr)
!<         zs = czi * pi * time_s( index_s(is, flvr), flvr ) / beta
!<         caux1(:,is) = exp(-two * zs)
!<         call s_cumprod_z(nfaux, caux1(:,is), caux1(:,is))
!<         caux1(:,is) = caux1(:,is) * exp(+(nffrq + 1) * zs)
!<     enddo ! over is={1,rank(flvr)} loop
!<
!<     ! annihilation operators
!<     !------------------------------------------------------------------
!<     ! for each \tau_e, we try to calculate
!<     !     exp ( i \omega_n \tau_e ) where n \in [1,nfaux]
!<     !     \omega_n = +(v + w), v: -v ---> +v, w: -0 ---> +w
!<     ! so,
!<     !     \omega_n = -v,   when n = 1
!<     !     \omega_n = +v+w, when n = nfaux
!<     !
!<     do ie=1,rank(flvr)
!<         ze = czi * pi * time_e( index_e(ie, flvr), flvr ) / beta
!<         caux2(:,ie) = exp(+two * ze)
!<         call s_cumprod_z(nfaux, caux2(:,ie), caux2(:,ie))
!<         caux2(:,ie) = caux2(:,ie) * exp(-(nffrq + 1) * ze)
!<     enddo ! over ie={1,rank(flvr)} loop
!<
!<!! body]
!<
!<     return
!<  end subroutine ctqmc_make_fexp

!!
!! @sub ctqmc_make_bexp
!!
!! calculate product of matsubara frequency exponents exp(i \omega_n \tau)
!!
!! note:
!!
!!     unlike the above ctqmc_make_fexp(), here the matsubara frequency
!!     mesh is bosonic. since the number of bosonic frequency points is
!!     usually very small (that is nbfrq << nffrq), so the calculation is
!!     very efficient
!!
  subroutine ctqmc_make_bexp(flvr, nfaux, mrank, caux1, caux2)
     use constants, only : dp
     use constants, only : pi, two
     use constants, only : czi

     use control, only : beta

     use context, only : index_s, index_e
     use context, only : time_s, time_e
     use context, only : rank

     implicit none

!! external arguments
     ! current flavor channel
     integer, intent(in) :: flvr

     ! number of frequency points, usually it is equal to nbfrq
     integer, intent(in) :: nfaux

     ! maximum number of operators in different flavor channels
     integer, intent(in) :: mrank

     ! matsubara frequency exponents for creation operators
     complex(dp), intent(out) :: caux1(nfaux,mrank)

     ! matsubara frequency exponents for annihilation operators
     complex(dp), intent(out) :: caux2(nfaux,mrank)

!! local variables
     ! loop indices for start and end points
     integer :: is
     integer :: ie

     ! loop index for matsubara frequency
     integer :: iw

     ! imaginary time for start and end points
     ! actually, they are i\pi\tau_s/\beta and i\pi\tau_e/\beta
     complex(dp) :: zs
     complex(dp) :: ze

!! [body

     ! creation operators
     !--------------------------------------------------------------------
     ! for each \tau_s, we try to calculate
     !     exp ( i \omega_n \tau_s ) where n \in [1,nfaux]
     !     \omega_n = - 2 (n - 1) \pi / beta
     !
     do is=1,rank(flvr)
         zs = czi * pi * time_s( index_s(is, flvr), flvr ) / beta
         do iw=1,nfaux
             caux1(iw,is) = exp( -two * float(iw - 1) * zs )
         enddo ! over iw={1,nfaux} loop
     enddo ! over is={1,rank(flvr)} loop

     ! annihilation operators
     !--------------------------------------------------------------------
     ! for each \tau_e, we try to calculate
     !     exp ( i \omega_n \tau_e ) where n \in [1,nfaux]
     !     \omega_n = + 2 (n - 1) \pi / beta
     !
     do ie=1,rank(flvr)
         ze = czi * pi * time_e( index_e(ie, flvr), flvr ) / beta
         do iw=1,nfaux
             caux2(iw,ie) = exp( +two * float(iw - 1) * ze )
         enddo ! over iw={1,nfaux} loop
     enddo ! over ie={1,rank(flvr)} loop

!! body]

     return
  end subroutine ctqmc_make_bexp

!!========================================================================
!!>>> self-energy function                                             <<<
!!========================================================================

!!
!! @sub ctqmc_make_hub2
!!
!! first of all, build impurity green's function and auxiliary correlation
!! function via fast fourier transformation (if isort == 1) or analytical
!! formula (if isort == 2 or isort == 3). and then, self-energy function
!! is obtained by using the improved estimator trick
!!
!! see Eq. (13) in Phys. Rev. B 85, 205106 (2012)
!! or  Eq. (30) in Phys. Rev. B 89, 235128 (2014)
!!
  subroutine ctqmc_make_hub2()
     use constants, only : dp

     use control, only : norbs
     use control, only : mfreq
     use control, only : nfreq

     use context, only : gtau, ftau
     use context, only : grnf, frnf
     use context, only : sig2

     implicit none

!! local variables
     ! loop index
     integer  :: i
     integer  :: k

     ! status flag
     integer  :: istat

     ! it is used to backup the sampled impurity green's function
     complex(dp), allocatable :: gtmp(:,:,:)

!! [body

     ! allocate memory
     allocate(gtmp(nfreq,norbs,norbs), stat=istat)

     ! backup the sampled impurity green's function
     gtmp = grnf(1:nfreq,:,:)

     ! build impurity green's function and auxiliary correlation function
     call ctqmc_tran_grnf(gtau, grnf)
     call ctqmc_tran_grnf(ftau, frnf)

     ! build final self-energy function by using improved estimator
     do i=1,norbs
         do k=1,mfreq
             sig2(k,i,i) = frnf(k,i,i) / grnf(k,i,i)
         enddo ! over k={1,nfreq} loop
     enddo ! over i={1,norbs} loop

     ! restore the sampled impurity green's function
     grnf(1:nfreq,:,:) = gtmp(1:nfreq,:,:)

     ! deallocate memory
     deallocate(gtmp)

!! body]

     return
  end subroutine ctqmc_make_hub2
