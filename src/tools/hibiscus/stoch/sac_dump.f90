!!!-----------------------------------------------------------------------
!!! project : hibiscus/stoch
!!! program : sac_dump_image
!!!           sac_dump_aprob
!!! source  : sac_dump.f90
!!! type    : subroutines
!!! author  : li huang (email:huangli712@gmail.com)
!!! history : 01/08/2011 by li huang
!!!           12/14/2011 by li huang
!!!           11/18/2014 by li huang
!!! purpose : dump key observables produced by the stochastic analytic
!!!           continuation code
!!! status  : unstable
!!! comment :
!!!-----------------------------------------------------------------------

!>>> write out image function in real frequency space
! note: sai.image.dat file contain the alpha-resolved image function, and
! sai.imsum.dat file contain the averaged image function
  subroutine sai_dump_image(step)
     use constants
     use control
     use context

     implicit none

! external arguments
! current monte carlo sweep number
     real(dp), intent(in) :: step

! local variables
! loop index
     integer  :: i
     integer  :: j

! dummy image function in legendre polynomial representation
     real(dp) :: image_l(nalph,lemax)

! dummy image function in normal representation
     real(dp) :: image_t(nalph,-nwmax:nwmax)

! postprocess the image, only for legendre polynomial representation
     if ( ltype /= 1 ) then
! convert image function from legendre polynomial representation to 
! normal representation
         image_l = immpi(:,1:lemax) / step
         call sai_warp_image(image_l, image_t)

! dump the legendre polynomial coefficient
         open(mytmp, file='sai.ppleg.dat', form='formatted', status='unknown')
         do i=1,nalph
             do j=1,lemax
                 write(mytmp,*) j, image_l(i,j)
             enddo ! over j={1,lemax} loop
             write(mytmp,*)
             write(mytmp,*)
         enddo ! over i={1,nalph} loop
         close(mytmp)
     endif ! back if ( ltype /=1 ) block
 
! postprocess the image, copy it to image_t
     if ( ltype == 1 ) then
         do i=1,nalph
             do j=-nwmax,nwmax
                 image_t(i,j) = immpi(i,j) * model(j) / ( pi * step )
             enddo ! over j={-nwmax,nwmax} loop
         enddo ! over i={1,nalph} loop
     else
         do i=1,nalph
             do j=-nwmax,nwmax
                 image_t(i,j) = image_t(i,j) * model(j)
             enddo ! over j={-nwmax,nwmax} loop
         enddo ! over i={1,nalph} loop
     endif ! back if ( ltype == 1 ) block

! open data file: sai.image.dat
     open(mytmp, file='sai.image.dat', form='formatted', status='unknown')

! write it
     do i=1,nalph
         write(mytmp,'(a,i4,f16.8)') '# alpha:', i, alpha(i)
         do j=-nwmax,nwmax
             write(mytmp,'(2f16.8)') wmesh(j), image_t(i,j)
         enddo ! over j={-nwmax,nwmax} loop
         write(mytmp,*) ! write empty lines
         write(mytmp,*)
     enddo ! over i={1,nalph} loop

! close data file
     close(mytmp)

! open data file: sai.imsum.dat
     open(mytmp, file='sai.imsum.dat', form='formatted', status='unknown')

! write it
     do j=-nwmax,nwmax
         write(mytmp,'(2f16.8)') wmesh(j), sum( image_t(:,j) ) / real(nalph)
     enddo ! over j={-nwmax,nwmax} loop

! close data file
     close(mytmp)

     return
  end subroutine sai_dump_image

!>>> write out the accept/reject probability for different alpha parameters
  subroutine sai_dump_aprob(step)
     use constants
     use control
     use context

     implicit none

! external arguments
! current monte carlo sweep number
     real(dp), intent(in) :: step

! local variables
! loop index
     integer :: i

! open data file: sai.move.dat
     open(mytmp, file='sai.move.dat', form='formatted', access='append')

! write it
     write(mytmp,'(i8)', advance='no') int(step) / ndump
     do i=1,nalph
         write(mytmp,'(f12.6)', advance='no') move_accept(i) / move_tcount(i)
     enddo ! over i={1,nalph} loop
     write(mytmp,*)

! close data file
     close(mytmp)

! open data file: sai.swap.dat
     open(mytmp, file='sai.swap.dat', form='formatted', access='append')

! write it
     write(mytmp,'(i8)', advance='no') int(step) / ndump
     do i=1,nalph
         if ( nalph == 1 ) then
             write(mytmp,'(f12.6)', advance='no') zero
         else
             write(mytmp,'(f12.6)', advance='no') swap_accept(i) / swap_tcount(i)
         endif ! back if ( nalph == 1 ) block
     enddo ! over i={1,nalph} loop
     write(mytmp,*)

! close data file
     close(mytmp)

     return
  end subroutine sai_dump_aprob
