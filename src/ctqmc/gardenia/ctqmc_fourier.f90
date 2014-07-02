!-------------------------------------------------------------------------
! project : gardenia
! program : ctqmc_fourier_htau
!           ctqmc_fourier_hybf
!           ctqmc_fourier_tails
!           ctqmc_fourier_forward
!           ctqmc_fourier_backward
! source  : ctqmc_fourier.f90
! type    : subroutines
! author  : li huang (email:huangli712@gmail.com)
! history : 05/05/2008 by li huang
!           01/18/2009 by li huang
!           09/27/2009 by li huang
!           10/20/2009 by li huang
!           11/01/2009 by li huang
!           12/01/2009 by li huang
!           02/28/2010 by li huang
!           03/07/2010 by li huang
! purpose : forward and backward fourier transformation subroutines for
!           hybridization function
! input   :
! output  :
! status  : unstable
! comment : nominally, the following subroutines are only suitable for the
!           hybridization functions, but in principle, we can also apply
!           them to the impurity green's function and bath weiss's function
!-------------------------------------------------------------------------

!>>> fourier htau to hybf, from imaginary time to matsubara frequency
  subroutine ctqmc_fourier_htau(htau, hybf)
     use constants
     use control

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
             call ctqmc_fourier_forward(ntime, raux, mfreq, caux)

! copy the matsubara frequency data to hybf
             hybf(:,j,i) = caux

         enddo ! over j={1,norbs} loop
     enddo ! over i={1,norbs} loop

     return
  end subroutine ctqmc_fourier_htau

!>>> fourier hybf to htau, from matsubara frequency to imaginary time
  subroutine ctqmc_fourier_hybf(hybf, htau)
     use constants
     use control

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
             call ctqmc_fourier_backward(mfreq, caux, ntime, raux)

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
                 start = j
                 EXIT
             endif
         enddo ! over j={1,ntime} loop

         do j=ntime,1,-1 ! search backward
             if ( htau(j,i,i) > -eps6 ) then
                 last = j
                 EXIT
             endif
         enddo ! over j={ntime,1,-1} loop

!-------------------------------------------------------------------------
!<         if ( start > 1 .and. last > 1 ) then
!<             do j=start,last
!<                 htau(j,i,i) = -eps6
!<             enddo ! over j={start,last} loop
!<         endif
!-------------------------------------------------------------------------
     enddo ! over i={1,norbs} loop

! enforce hybridization function less than zero to ensure the causality
     do i=1,norbs
         do j=1,ntime
             if ( htau(j,i,i) > zero ) htau(j,i,i) = -eps6
         enddo ! over j={1,ntime} loop
     enddo ! over i={1,norbs} loop

     return
  end subroutine ctqmc_fourier_hybf

!>>> calculate high frequency tails using K. Haule's trick
  subroutine ctqmc_fourier_tails(tail, rmesh, fmat)
     use constants
     use control

     implicit none

! external arguments
! high frequency tail
     real(dp), intent(out) :: tail

! matsubara frequency grid
     real(dp), intent(in) :: rmesh(mfreq)

! function on matsubara frequency space
     complex(dp), intent(in) :: fmat(mfreq)

! local variables
! loop index
     integer  :: j

! dummy variables
     real(dp) :: Sn, Sx, Sy
     real(dp) :: Sxx, Sxy

     Sn = zero
     Sx = zero
     Sy = zero

     Sxx = zero
     Sxy = zero

     do j=mfreq - nfreq, mfreq
         Sn = Sn + one
         Sx = Sx + one / rmesh(j)**2
         Sy = Sy + aimag(fmat(j)) * rmesh(j)
         Sxx = Sxx + one / rmesh(j)**4
         Sxy = Sxy + aimag(fmat(j)) * rmesh(j) / rmesh(j)**2
     enddo ! over j={mfreq - nfreq, mfreq} loop

     tail = (Sx * Sxy - Sxx * Sy) / (Sn * Sxx - Sx * Sx)

     return
  end subroutine ctqmc_fourier_tails

!>>> fourier from imaginary time space forward to matsubara frequency space
! using linear fourier algorithm
  subroutine ctqmc_fourier_forward(ntime, ftau, mfreq, fmat)
     use constants
     use control, only : beta
     use context, only : rmesh, tmesh

     implicit none

! external arguments
! number of matsubara frequency points
     integer, intent(in) :: mfreq

! number of imaginary time points
     integer, intent(in) :: ntime

! function on imaginary time axis
     real(dp), intent(in) :: ftau(ntime)

! function on matsubara frequency axis
     complex(dp), intent(out) :: fmat(mfreq)

! local variables
! loop index
     integer  :: i
     integer  :: j

! dummy variables
     real(dp) :: sre, sim

     real(dp) :: c0, c1
     real(dp) :: s0, s1
     real(dp) :: g0, g1, dg

     do i=1,mfreq
         sre = zero
         sim = zero

         do j=1,ntime-1
             c0 = cos( tmesh(j)   * rmesh(i) )
             c1 = cos( tmesh(j+1) * rmesh(i) )
             s0 = sin( tmesh(j)   * rmesh(i) )
             s1 = sin( tmesh(j+1) * rmesh(i) )
             g0 = ftau(j)
             g1 = ftau(j+1)
             dg = ( g1 - g0 ) / ( tmesh(j+1) - tmesh(j) )
             sim = sim + ( c0 * g0 - c1 * g1 + dg * (s1 - s0) / rmesh(i) ) / rmesh(i)
             sre = sre + ( s1 * g1 - s0 * g0 + dg * (c1 - c0) / rmesh(i) ) / rmesh(i)
         enddo ! over j={1,ntime-1} loop

         fmat(i) = dcmplx(sre, sim)
     enddo ! over i={1,mfreq} loop

     return
  end subroutine ctqmc_fourier_forward

!>>> fourier from matsubara frequency space backward to imaginary time space
  subroutine ctqmc_fourier_backward(mfreq, fmat, ntime, ftau)
     use constants
     use control, only : beta
     use context, only : rmesh, tmesh

     implicit none

! external arguments
! number of matsubara frequency points
     integer, intent(in) :: mfreq

! number of imaginary time points
     integer, intent(in) :: ntime

! function on imaginary time axis
     real(dp), intent(out) :: ftau(ntime)

! function on matsubara frequency axis
     complex(dp), intent(in) :: fmat(mfreq)

! local variables
! loop index
     integer  :: i
     integer  :: j

! real(dp) dummy variables
     real(dp) :: raux
     real(dp) :: tail

! calculate high frequency tails need to be subtracted
     call ctqmc_fourier_tails(tail, rmesh, fmat)

! perform infourier transformation
     do i=1,ntime
         raux = zero
         do j=1,mfreq
             raux = raux + cos( rmesh(j) * tmesh(i) ) *   real( fmat(j) )
             raux = raux + sin( rmesh(j) * tmesh(i) ) * (aimag( fmat(j) ) + tail / rmesh(j))
         enddo ! over j={1,mfreq} loop
         ftau(i) = two * raux / beta - half * tail
     enddo ! over i={1,ntime} loop

! corrections for the boundary point
     raux = real( fmat(mfreq) ) * rmesh(mfreq) / pi
     ftau(1) = ftau(1) + raux
     ftau(ntime) = ftau(ntime) - raux

! additional corrections, may be useful for lda + dmft calculations
     ftau(1) = 3.0_dp * ftau(2) - 3.0_dp * ftau(3) + ftau(4)
     ftau(ntime) = 3.0_dp * ftau(ntime-1) - 3.0_dp * ftau(ntime-2) + ftau(ntime-3)

     return
  end subroutine ctqmc_fourier_backward
