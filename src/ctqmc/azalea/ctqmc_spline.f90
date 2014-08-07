!-------------------------------------------------------------------------
! project : azalea
! program : ctqmc_make_htau
!           ctqmc_make_hsed
!           ctqmc_make_spline
!           ctqmc_make_splint
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
