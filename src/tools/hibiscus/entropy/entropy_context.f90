!!!-----------------------------------------------------------------------
!!! project : hibiscus/entropy
!!! program : context    module
!!! source  : entropy_context.f90
!!! type    : module
!!! author  : li huang (email:huangli712@gmail.com)
!!! history : 01/08/2011 by li huang
!!!           01/26/2011 by li huang
!!!           11/17/2014 by li huang
!!! purpose : define the key data structure and global arrays/variables
!!!           for classic maximum entropy method code
!!! status  : unstable
!!! comment :
!!!-----------------------------------------------------------------------

!!========================================================================
!!>>> module context                                                   <<<
!!========================================================================

!!>>> containing memory management subroutines and define global variables
  module context
     use constants
     use control

     implicit none

! imaginary time mesh in [0,\beta]
     real(dp), public, save, allocatable :: tmesh(:)

! real frequency mesh
     real(dp), public, save, allocatable :: wmesh(:)

! default model D(\omega)
     real(dp), public, save, allocatable :: model(:)

! sum-rules values of image function
     real(dp), public, save, allocatable :: srule(:,:)

! normalization function
     real(dp), public, save, allocatable :: fnorm(:,:)

! imaginary time green's function from QMC quantum impurity solver: G(t)
     real(dp), public, save, allocatable :: G_qmc(:,:)

! standard deviation of G(t) from QMC data: \sigma(t)
     real(dp), public, save, allocatable :: G_dev(:,:)

! fermion kernel function
     real(dp), public, save, allocatable :: fkern(:,:)

! image function, i.e., spectral function
     real(dp), public, save, allocatable :: image(:,:)

  contains ! encapsulated functionality

!!========================================================================
!!>>> allocate memory subroutines                                      <<<
!!========================================================================

!!>>> entropy_allocate_memory: allocate module memory
  subroutine entropy_allocate_memory()
     implicit none

! status flag
     integer :: istat

! allocate memory
     allocate(tmesh(ntime),              stat=istat)
     allocate(wmesh(-nwmax:nwmax),       stat=istat)
     allocate(model(-nwmax:nwmax),       stat=istat)

     allocate(srule(3,norbs),            stat=istat)
     allocate(fnorm(-nwmax:nwmax,3),     stat=istat)

     allocate(G_qmc(ntime,norbs),        stat=istat)
     allocate(G_dev(ntime,norbs),        stat=istat)

     allocate(fkern(-nwmax:nwmax,ntime), stat=istat)
     allocate(image(-nwmax:nwmax,norbs), stat=istat)

! check the status
     if ( istat /= 0 ) then
         call s_print_error('entropy_allocate_memory','can not allocate enough memory')
     endif ! back if ( istat /= 0 ) block

! initialize them
     tmesh = zero
     wmesh = zero
     model = zero

     srule = zero
     fnorm = zero

     G_qmc = zero
     G_dev = zero

     fkern = zero
     image = zero

     return
  end subroutine entropy_allocate_memory

!!========================================================================
!!>>> deallocate memory subroutines                                    <<<
!!========================================================================

!!>>> entropy_deallocate_memory: deallocate module memory
  subroutine entropy_deallocate_memory()
     implicit none

     if ( allocated(tmesh) ) deallocate(tmesh)
     if ( allocated(wmesh) ) deallocate(wmesh)
     if ( allocated(model) ) deallocate(model)

     if ( allocated(srule) ) deallocate(srule)
     if ( allocated(fnorm) ) deallocate(fnorm)

     if ( allocated(G_qmc) ) deallocate(G_qmc)
     if ( allocated(G_dev) ) deallocate(G_dev)

     if ( allocated(fkern) ) deallocate(fkern)
     if ( allocated(image) ) deallocate(image)

     return
  end subroutine entropy_deallocate_memory

  end module context
