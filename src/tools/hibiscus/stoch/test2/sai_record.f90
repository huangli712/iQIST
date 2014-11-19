!-------------------------------------------------------------------------
! project : hibiscus
! program : sai_recording
!           sai_reducing
!           sai_make_image
! source  : sai_record.f90
! type    : subroutines
! author  : li huang (email:huangli712@yahoo.com.cn)
! history : 01/09/2011 by li huang
!           01/10/2011 by li huang
!           08/10/2011 by li huang
! purpose : measure, record, and postprocess the key observables
! input   :
! output  :
! status  : unstable
! comment :
!-------------------------------------------------------------------------

!>>> update the image function over every alpha parameters
  subroutine sai_recording()
     use control, only : nalph

     implicit none

! local variables
! loop index over alpha parameters
     integer :: i

     do i=1,nalph
         call sai_make_image(i)
     enddo ! over i={1,nalph} loop

     return
  end subroutine sai_recording

!>>> reduce image function from all children nodes
  subroutine sai_reducing()
     use constants
     use control
     use context

     use mmpi

     implicit none

! init immpi
     immpi = zero

! record image to immpi, collect data from all children processes
# if defined (MPI)

! collect data
     call mp_allreduce(image, immpi)

! block until all processes have reached here
     call mp_barrier()

# else  /* MPI */

     immpi = image

# endif /* MPI */

! calculate the average image function
     immpi = immpi / real(nprocs)

     return
  end subroutine sai_reducing

!>>> calculate image function using current configurations
  subroutine sai_make_image(ia)
     use constants
     use control
     use context

     implicit none

! external arguments
! index of alpha parameters
     integer, intent(in) :: ia

! local variables
! loop index
     integer  :: i
     integer  :: j

! integer dummy variables
     integer  :: si, curr

! real(dp) dummy variables
     real(dp) :: sr, step

! dummy image function
     real(dp) :: image_t(-nwmax:nwmax)

! init dummy image function
     image_t = zero

! evaluate the image function using \delta function
     if ( ltype == 1 ) then
         do i=1,ngamm
             sr = rgamm(ia,i)
             si = igamm(ia,i)
             do j=-nwmax,nwmax
                 image_t(j) = image_t(j) + sr * delta(j,si)
             enddo ! over j={-nwmax,nwmax} loop
         enddo ! over i={1,ngamm} loop
! evaluate the image function using legendre polynomial function
     else
         step = real(legrd - 1) / two
         do i=1,ngamm
             sr = rgamm(ia,i)
             si = igamm(ia,i)
             curr = nint( two * xgrid(si) * step ) + 1
! note: since only odd term of legendre polynomial coefficients are nonzero
             do j=1,lemax,2
                 image_t(j) = image_t(j) + sr * sqrt(two * j - one) * ppleg(curr,j)
             enddo ! over j={1,lemax} loop
         enddo ! over i={1,ngamm} loop
     endif ! back if ( ltype == 1 ) block

! add up image_t to image
     image(ia,:) = image(ia,:) + image_t

     return
  end subroutine sai_make_image
