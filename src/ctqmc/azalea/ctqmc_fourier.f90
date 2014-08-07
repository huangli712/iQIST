!-------------------------------------------------------------------------
! project : azalea
! program : ctqmc_fourier_htau
!           ctqmc_fourier_hybf
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
  end subroutine ctqmc_fourier_htau

!>>> fourier hybf to htau, from matsubara frequency to imaginary time
  subroutine ctqmc_fourier_hybf(hybf, htau)
     use constants
     use control
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
