!-------------------------------------------------------------------------
! project : azalea
! program : ctqmc_make_uumat
!           ctqmc_make_state
! source  : ctqmc_util.f90
! type    : functions & subroutines
! author  : li huang (email:huangli712@gmail.com)
! history : 10/01/2008 by li huang
!           02/08/2009 by li huang
!           09/23/2009 by li huang
!           09/26/2009 by li huang
!           11/17/2009 by li huang
!           11/21/2009 by li huang
!           12/18/2009 by li huang
!           12/22/2009 by li huang
!           12/29/2009 by li huang
!           01/12/2010 by li huang
!           02/27/2010 by li huang
!           06/08/2010 by li huang
!           06/22/2010 by li huang
! purpose : to provide utility functions and subroutines for hybridization
!           expansion version continuous time quantum Monte Carlo (CTQMC)
!           quantum impurity solver
! status  : unstable
! comment :
!-------------------------------------------------------------------------

!>>> to build general U interaction matrix: uumat, using my own style
! note: do not support spin-flip and pair-hopping term so far
! note: only Uc and Jz are need, the other Coulomb interaction parameters
! are used as backup
  subroutine ctqmc_make_uumat(uumat)
     use constants
     use control

     implicit none

! external arguments
! Coulomb interaction matrix
     real(dp), intent(out) :: uumat(norbs, norbs)

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
                 endif
             else
                 ut(k) = Uc - 3.0_dp * Jz
             endif

             uumat(i,j) = ut(k)
             uumat(j,i) = ut(k)
         enddo ! over j={i+1,norbs} loop
     enddo ! over i={1,norbs-1} loop

     return
  end subroutine ctqmc_make_uumat

!>>> convert current atomic state array into a decimal number (state index)
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
!-------------------------------------------------------------------------
! project : azalea
! program : ctqmc_make_htau
!           ctqmc_make_hsed
! source  : ctqmc_spline.f90
! type    : subroutines
! author  : li huang (email:huangli712@gmail.com)
! history : 05/05/2008 by li huang
!           02/18/2009 by li huang
!           09/23/2009 by li huang
!           09/26/2009 by li huang
!           10/03/2009 by li huang
!           11/10/2009 by li huang
!           12/18/2009 by li huang
! purpose : to provide cubic spline subroutines and wrapper functions to
!           interpolate the hybridization function in imaginary-time axis
! input   :
! output  :
! status  : unstable
! comment :
!-------------------------------------------------------------------------

!>>> evaluate the matrix elements for mmat matrix using cubic spline interpolation
  function ctqmc_make_htau(flvr, dtau) result(val)
     use constants
     use control
     use context

     implicit none

! external arguments
! current flavor channel
     integer, intent(in)  :: flvr

! delta imaginary time
     real(dp), intent(in) :: dtau

! external functions
! internal interpolation engine
     procedure(real(dp)) :: s_spl_splint

! local variables
! return value
     real(dp) :: val

     val = s_spl_splint(ntime, tmesh, htau(:, flvr, flvr), hsed(:, flvr, flvr), dtau)

     return
  end function ctqmc_make_htau

!>>> calculate the second order derivates of hybridization function on
! imaginary time space
  subroutine ctqmc_make_hsed(tmesh, htau, hsed)
     use constants
     use control

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
             call s_spl_splder(ntime, tmesh, htau(:,i,j), startu, startd, d2y)

! copy the results to hsed
             hsed(:,i,j) = d2y

         enddo ! over i={1,norbs} loop
     enddo ! over j={1,norbs} loop

     return
  end subroutine ctqmc_make_hsed
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
