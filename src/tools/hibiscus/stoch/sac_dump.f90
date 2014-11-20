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

!!>>> sac_dump_image: write out image function in real frequency space
  subroutine sac_dump_image(step)
     use constants, only : dp, pi, mytmp

     use control, only : nwmax, nalph, ltype, lemax
     use context, only : immpi, model, wmesh, alpha

     implicit none

! external arguments
! current monte carlo sweep number
     real(dp), intent(in) :: step

! local variables
! loop index
     integer  :: i
     integer  :: j

! dummy image function in legendre polynomial representation
     real(dp) :: image_l(lemax,nalph)

! dummy image function in normal representation
     real(dp) :: image_t(-nwmax:nwmax,nalph)

! postprocess the image, only for legendre polynomial representation
     if ( ltype /= 1 ) then
! convert image function from legendre polynomial representation to
! normal representation
         image_l = immpi(1:lemax,:) / step
         call sac_warp_image(image_l, image_t)
     endif ! back if ( ltype /=1 ) block

! postprocess the image, copy it to image_t
     if ( ltype == 1 ) then
         do i=1,nalph
             do j=-nwmax,nwmax
                 image_t(j,i) = immpi(j,i) * model(j) / ( pi * step )
             enddo ! over j={-nwmax,nwmax} loop
         enddo ! over i={1,nalph} loop
     else
         do i=1,nalph
             do j=-nwmax,nwmax
                 image_t(j,i) = image_t(j,i) * model(j)
             enddo ! over j={-nwmax,nwmax} loop
         enddo ! over i={1,nalph} loop
     endif ! back if ( ltype == 1 ) block

! note: sac.image.dat file contain the alpha-resolved image function, and
! sac.imsum.dat file contain the averaged image function

! open data file: sac.image.dat
     open(mytmp, file='sac.image.dat', form='formatted', status='unknown')

! write it
     do i=1,nalph
         write(mytmp,'(a,i4,f16.8)') '# alpha:', i, alpha(i)
         do j=-nwmax,nwmax
             write(mytmp,'(2f16.8)') wmesh(j), image_t(j,i)
         enddo ! over j={-nwmax,nwmax} loop
         write(mytmp,*) ! write empty lines
         write(mytmp,*)
     enddo ! over i={1,nalph} loop

! close data file
     close(mytmp)

! open data file: sac.imsum.dat
     open(mytmp, file='sac.imsum.dat', form='formatted', status='unknown')

! write it
     do j=-nwmax,nwmax
         write(mytmp,'(2f16.8)') wmesh(j), sum( image_t(j,:) ) / real(nalph)
     enddo ! over j={-nwmax,nwmax} loop

! close data file
     close(mytmp)

     return
  end subroutine sac_dump_image

!!>>> sac_dump_aprob: write out the accept/reject probability for
!!>>> different alpha parameters
  subroutine sac_dump_aprob(step)
     use constants, only : dp, zero, mytmp

     use control, only : nalph, ndump
     use context, only : move_accept, move_tcount
     use context, only : swap_accept, swap_tcount

     implicit none

! external arguments
! current monte carlo sweep number
     real(dp), intent(in) :: step

! local variables
! loop index
     integer :: i

! open data file: sac.move.dat
     open(mytmp, file='sac.move.dat', form='formatted', access='append')

! write it
     write(mytmp,'(i8)', advance='no') int(step) / ndump
     do i=1,nalph
         write(mytmp,'(f12.6)', advance='no') move_accept(i) / move_tcount(i)
     enddo ! over i={1,nalph} loop
     write(mytmp,*)

! close data file
     close(mytmp)

! open data file: sac.swap.dat
     open(mytmp, file='sac.swap.dat', form='formatted', access='append')

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
  end subroutine sac_dump_aprob
