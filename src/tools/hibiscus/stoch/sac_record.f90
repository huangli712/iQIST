!!!-----------------------------------------------------------------------
!!! project : hibiscus/stoch
!!! program : sac_recording
!!!           sac_reducing
!!! source  : sac_record.f90
!!! type    : subroutines
!!! author  : li huang (email:lihuang.dmft@gmail.com)
!!! history : 01/09/2011 by li huang
!!!           08/10/2011 by li huang
!!!           11/18/2014 by li huang
!!! purpose : measure, record, and postprocess the key observables
!!! status  : unstable
!!! comment :
!!!-----------------------------------------------------------------------

!!>>> sac_recording: update the image function over every alpha parameters
  subroutine sac_recording()
     use constants, only : dp, zero, one, two

     use control, only : nwmax, ngamm, nalph
     use control, only : ltype, lemax, legrd
     use context, only : igamm, rgamm
     use context, only : ppleg, delta, image, xgrid

     implicit none

! local variables
! loop index
     integer  :: i
     integer  :: j

! loop index of alpha parameters
     integer  :: ia

! status flag
     integer  :: istat

! integer dummy variables
     integer  :: si, curr

! real(dp) dummy variables
     real(dp) :: sr, step

! dummy image function
     real(dp), allocatable :: image_t(:)

! allocate memory
     allocate(image_t(-nwmax:nwmax), stat=istat)

     ALPH_LOOP: do ia=1,nalph

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
         image(:,ia) = image(:,ia) + image_t

     enddo ALPH_LOOP ! over ia={1,nalph} loop

! deallocate memory
     deallocate(image_t)

     return
  end subroutine sac_recording

!!>>> sac_reducing: reduce image function from all children nodes
  subroutine sac_reducing()
     use constants, only : zero
     use mmpi, only : mp_allreduce, mp_barrier

     use control, only : nprocs
     use context, only : image, immpi

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
  end subroutine sac_reducing
