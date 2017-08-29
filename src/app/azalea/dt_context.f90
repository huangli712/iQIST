!!!-----------------------------------------------------------------------
!!! project : azalea
!!! program : dt_mesh module
!!!           dt_dmft module
!!!           dt_dual module
!!!           dt_latt module
!!!           dt_vert module
!!!           context module
!!! source  : dmft_context.f90
!!! type    : modules
!!! author  : li huang (email:lihuang.dmft@gmail.com)
!!! history : 09/16/2009 by li huang (created)
!!!           08/29/2017 by li huang (last modified)
!!! purpose : define the key data structure and global arrays/variables
!!!           for diagrammatic framework for dynamical mean field theory.
!!! status  : unstable
!!! comment :
!!!-----------------------------------------------------------------------

!!========================================================================
!!>>> module dt_mesh                                                   <<<
!!========================================================================

!!
!! @mod dt_mesh
!!
  module dt_mesh
     use constants, only : dp

     real(dp), public, save, allocatable :: kx(:)
     real(dp), public, save, allocatable :: ky(:)
     real(dp), public, save, allocatable :: kz(:)

!!
!! @var ek
!!
!! band dispersion
!!
     real(dp), public, save, allocatable :: ek(:)

!!
!! @var fmesh
!!
!! fermionic frequencies
!!
     real(dp), public, save, allocatable :: fmesh(:)

!!
!! @var bmesh
!!
!! bosonic frequencies
!!
     real(dp), public, save, allocatable :: bmesh(:)

  end module dt_mesh

!!========================================================================
!!>>> module dt_dmft                                                   <<<
!!========================================================================

!!
!! @mod dt_dmft
!!
  module dt_dmft
     use constants, only : dp

     implicit none

!!
!! @var dmft_g
!!
!! impurity green's function from dmft
!!
     complex(dp), public, save, allocatable :: dmft_g(:,:)

!!
!! @var dmft_s
!!
!! self-energy function from dmft
!!
     complex(dp), public, save, allocatable :: dmft_s(:,:)

!!
!! @var dmft_h
!!
!! hybridization function from dmft
!!
     complex(dp), public, save, allocatable :: dmft_h(:,:)

  end module dt_dmft

!!========================================================================
!!>>> module dt_dual                                                   <<<
!!========================================================================

!!
!! @mod dt_dual
!!
  module dt_dual
     use constants, only : dp

     implicit none

!!
!! @var dual_g
!!
!! dual green's function
!!
     complex(dp), public, save, allocatable :: dual_g(:,:,:)

!!
!! @var dual_s
!!
!! dual self-energy function
!!
     complex(dp), public, save, allocatable :: dual_s(:,:,:)

!!
!! @var dual_b
!!
!! dual bare green's function
!!
     complex(dp), public, save, allocatable :: dual_b(:,:,:)

  end module dt_dual

!!========================================================================
!!>>> module dt_latt                                                   <<<
!!========================================================================

!!
!! @mod latt
!!
  module dt_latt
     use constants, only : dp

     implicit none

!!
!! @var latt_g
!!
!! lattice green's function
!!
     complex(dp), public, save, allocatable :: latt_g(:,:,:)

!!
!! @var latt_s
!!
!! lattice self-energy function
!!
     complex(dp), public, save, allocatable :: latt_s(:,:,:)

  end module dt_latt

!!========================================================================
!!>>> module dt_vert                                                   <<<
!!========================================================================

!!
!! @mod dt_vert
!!
  module dt_vert
     use constants, only : dp

     implicit none

!!
!! @var vert_d
!!
!! density vertex
!!
     complex(dp), public, save, allocatable :: vert_d(:,:,:)

!!
!! @var vert_m
!!
!! magnetic vertex
!!
     complex(dp), public, save, allocatable :: vert_m(:,:,:)

  end module dt_vert

!!========================================================================
!!>>> module context                                                   <<<
!!========================================================================

  module context
     use constants, only : dp
     use constants, only : zero, czero

     use control, only : nkp_x, nkp_y, nkp_z, nkpts
     use control, only : nffrq, nbfrq
     use control, only : norbs

     use dt_mesh
     use dt_dmft
     use dt_dual
     use dt_latt
     use dt_vert

!!========================================================================
!!>>> declare private variables                                        <<<
!!========================================================================

! status flag
     integer, private :: istat

!!========================================================================
!!>>> declare accessibility for module routines                        <<<
!!========================================================================

! declaration of module procedures: allocate memory
     public :: dt_allocate_memory_mesh
     public :: dt_allocate_memory_dmft
     public :: dt_allocate_memory_dual
     public :: dt_allocate_memory_latt
     public :: dt_allocate_memory_vert

! declaration of module procedures: deallocate memory
     public :: dt_deallocate_memory_mesh
     public :: dt_deallocate_memory_dmft
     public :: dt_deallocate_memory_dual
     public :: dt_deallocate_memory_latt
     public :: dt_deallocate_memory_vert

  contains ! encapsulated functionality

!!========================================================================
!!>>> allocate memory subroutines                                      <<<
!!========================================================================

!!
!! @sub
!!
!! allocate memory for clur-related variables
!!
  subroutine dt_allocate_memory_mesh()
     implicit none

     allocate(kx(nkp_x), stat=istat)
     allocate(ky(nkp_y), stat=istat)
     allocate(kz(nkp_z), stat=istat)
     allocate(ek(nkpts), stat=istat)

     allocate(fmesh(nffrq), stat=istat)
     allocate(bmesh(nbfrq), stat=istat)

     if ( istat /= 0 ) then
         call s_print_error('dt_allocate_memory_mesh','can not allocate enough memory')
     endif

     kx = zero
     ky = zero
     kz = zero
     ek = zero

     fmesh = zero
     bmesh = zero

     return
  end subroutine dt_allocate_memory_mesh

  subroutine dt_allocate_memory_dmft()
     implicit none

     allocate(dmft_g(nffrq,norbs), stat=istat)
     allocate(dmft_s(nffrq,norbs), stat=istat)
     allocate(dmft_h(nffrq,norbs), stat=istat)

     if ( istat /= 0 ) then
         call s_print_error('dt_allocate_memory_dmft','can not allocate enough memory')
     endif

     dmft_g = czero
     dmft_s = czero
     dmft_h = czero

     return
  end subroutine dt_allocate_memory_dmft

  subroutine dt_allocate_memory_dual()
     implicit none

     allocate(dual_g(nkpts,nffrq,norbs), stat=istat)
     allocate(dual_s(nkpts,nffrq,norbs), stat=istat)
     allocate(dual_b(nkpts,nffrq,norbs), stat=istat)

     if ( istat /= 0 ) then
         call s_print_error('dt_allocate_memory_dual','can not allocate enough memory')
     endif

     dual_g = czero
     dual_s = czero
     dual_b = czero

     return
  end subroutine dt_allocate_memory_dual

  subroutine dt_allocate_memory_latt()
     implicit none

     allocate(latt_g(nkpts,nffrq,norbs), stat=istat)
     allocate(latt_s(nkpts,nffrq,norbs), stat=istat)

     if ( istat /= 0 ) then
         call s_print_error('dt_allocate_memory_latt','can not allocate enough memory')
     endif

     latt_g = czero
     latt_s = czero

     return
  end subroutine dt_allocate_memory_latt

  subroutine dt_allocate_memory_vert()
     implicit none

     allocate(vert_d(nffrq,nffrq,nbfrq), stat=istat)
     allocate(vert_m(nffrq,nffrq,nbfrq), stat=istat)

     if ( istat /= 0 ) then
         call s_print_error('dt_allocate_memory_vert','can not allocate enough memory')
     endif

     vert_d = czero
     vert_m = czero

     return
  end subroutine dt_allocate_memory_vert

!!========================================================================
!!>>> deallocate memory subroutines                                    <<<
!!========================================================================

  subroutine dt_deallocate_memory_mesh()
     implicit none

     if ( allocated(kx) ) deallocate(kx)
     if ( allocated(ky) ) deallocate(ky)
     if ( allocated(kz) ) deallocate(kz)
     if ( allocated(ek) ) deallocate(ek)

     if ( allocated(fmesh) ) deallocate(fmesh)
     if ( allocated(bmesh) ) deallocate(bmesh)

     return
  end subroutine dt_deallocate_memory_mesh

  subroutine dt_deallocate_memory_dmft()
     implicit none

     if ( allocated(dmft_g) ) deallocate(dmft_g)
     if ( allocated(dmft_s) ) deallocate(dmft_s)
     if ( allocated(dmft_h) ) deallocate(dmft_h)

     return
  end subroutine dt_deallocate_memory_dmft

  subroutine dt_deallocate_memory_dual()
     implicit none

     if ( allocated(dual_g) ) deallocate(dual_g)
     if ( allocated(dual_s) ) deallocate(dual_s)
     if ( allocated(dual_b) ) deallocate(dual_b)

     return
  end subroutine dt_deallocate_memory_dual

  subroutine dt_deallocate_memory_latt()
     implicit none

     if ( allocated(latt_g) ) deallocate(latt_g)
     if ( allocated(latt_s) ) deallocate(latt_s)

     return
  end subroutine dt_deallocate_memory_latt

  subroutine dt_deallocate_memory_vert()
     implicit none

     if ( allocated(vert_d) ) deallocate(vert_d)
     if ( allocated(vert_m) ) deallocate(vert_m)

     return
  end subroutine dt_deallocate_memory_vert

  end module context
