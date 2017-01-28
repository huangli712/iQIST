!!!-----------------------------------------------------------------------
!!! project : gardenia
!!! program : ctqmc_four_htau
!!!           ctqmc_four_hybf
!!!           ctqmc_eval_htau
!!!           ctqmc_eval_hsed
!!!           ctqmc_make_uumat
!!!           ctqmc_make_state
!!! source  : ctqmc_util.f90
!!! type    : functions & subroutines
!!! author  : li huang (email:lihuang.dmft@gmail.com)
!!! history : 10/01/2008 by li huang (created)
!!!           01/28/2015 by li huang (last modified)
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
