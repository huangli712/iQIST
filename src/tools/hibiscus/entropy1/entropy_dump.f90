!!!-----------------------------------------------------------------------
!!! project : hibiscus/entropy1
!!! program : entropy_dump_image
!!!           entropy_dump_srule
!!! source  : entropy_dump.f90
!!! type    : subroutines
!!! author  : li huang (email:huangli712@gmail.com)
!!! history : 01/08/2011 by li huang
!!!           01/26/2011 by li huang
!!!           11/17/2014 by li huang
!!! purpose : dump key observables produced by the classic maximum entropy
!!!           method code
!!! status  : unstable
!!! comment :
!!!-----------------------------------------------------------------------

!!>>> entropy_dump_image: write out image function in real frequency space
  subroutine entropy_dump_image(wmesh, image)
     use constants
     use control

     implicit none

! external arguments
! real frequency mesh
     real(dp), intent(in) :: wmesh(-nwmax:nwmax)

! image function, i.e., spectral function
     real(dp), intent(in) :: image(-nwmax:nwmax,norbs)

! local variables
! loop index
     integer :: i
     integer :: j

! open the data file: mem.dos.dat
     open(mytmp, file='mem.dos.dat', form='formatted', status='unknown')

! write it
     do i=1,norbs
         do j=-nwmax,nwmax
             write(mytmp,'(2f16.8)') wmesh(j), image(j,i)
         enddo ! over j={-nwmax,nwmax} loop
         write(mytmp,*) ! write two blank lines
         write(mytmp,*)
     enddo ! over i={1,norbs} loop

! close the data file
     close(mytmp)

     return
  end subroutine entropy_dump_image

!!>>> entropy_dump_srule: write out final sumrules values for image function
  subroutine entropy_dump_srule(srule)
     use constants
     use control

     implicit none

! external arguments
! sumrules values
     real(dp), intent(in) :: srule(3,norbs)

! local variables
! loop index
     integer :: i

! open the data file: mem.sum.dat
     open(mytmp, file='mem.sum.dat', form='formatted', status='unknown')

! write it
     do i=1,norbs
         write(mytmp,'(i4,3f16.8)') i, srule(1,i), srule(2,i), srule(3,i)
     enddo ! over i={1,norbs} loop

! close the data file
     close(mytmp)

     return
  end subroutine entropy_dump_srule
