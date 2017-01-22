
!!========================================================================
!!>>> module df_dmft                                                   <<<
!!========================================================================

  module df_dmft
     use constants, only : dp

     implicit none

!! dmft variables
! dmft green's function
     complex(dp), public, save, allocatable :: dmft_g(:,:)

! dmft self-energy function
     complex(dp), public, save, allocatable :: dmft_s(:,:)

! dmft hybridization function
     complex(dp), public, save, allocatable :: dmft_h(:,:)

  end module df_dmft

!!========================================================================
!!>>> module df_dual                                                   <<<
!!========================================================================

  module df_dual
     use constants, only : dp

     implicit none

!! dual variables
! dual green's function
     complex(dp), public, save, allocatable :: dual_g(:,:,:)

! dual self-energy function
     complex(dp), public, save, allocatable :: dual_s(:,:,:)

! dual bare green's function
     complex(dp), public, save, allocatable :: dual_b(:,:,:)

  end module df_dual

!!========================================================================
!!>>> module df_latt                                                   <<<
!!========================================================================

  module df_latt
     use constants, only : dp

     implicit none

!! lattice variables
! lattice green's function
     complex(dp), public, save, allocatable :: latt_g(:,:,:)

! lattice self-energy function
     complex(dp), public, save, allocatable :: latt_s(:,:,:)

  end module df_latt

!!========================================================================
!!>>> module df_vert                                                   <<<
!!========================================================================

  module df_vert
     use constants, only : dp

     implicit none

!! vertex variables
! density vertex
     complex(dp), public, save, allocatable :: vert_d(:,:,:,:)

! magnetic vertex
     complex(dp), public, save, allocatable :: vert_m(:,:,:,:)

  end module df_vert

!!========================================================================
!!>>> module df_context                                                <<<
!!========================================================================

  module df_context
     use constants

     use df_control

     use df_dmft
     use df_dual
     use df_latt
     use df_vert

!!========================================================================
!!>>> declare global variables                                         <<<
!!========================================================================

! status flag
     integer, private :: istat

!!========================================================================
!!>>> declare accessibility for module routines                        <<<
!!========================================================================

! declaration of module procedures: allocate memory
     public :: df_allocate_memory_dmft
     public :: df_allocate_memory_dual
     public :: df_allocate_memory_latt
     public :: df_allocate_memory_vert

! declaration of module procedures: deallocate memory
     public :: df_deallocate_memory_dmft
     public :: df_deallocate_memory_dual
     public :: df_deallocate_memory_latt
     public :: df_deallocate_memory_vert

  contains ! encapsulated functionality

!!========================================================================
!!>>> allocate memory subroutines                                      <<<
!!========================================================================

  subroutine df_allocate_memory_dmft()
     implicit none

     allocate(dmft_g(nffrq,norbs), stat=istat)
     allocate(dmft_s(nffrq,norbs), stat=istat)
     allocate(dmft_h(nffrq,norbs), stat=istat)

     if ( istat /= 0 ) then
         call s_print_error('df_allocate_memory_dmft','can not allocate enough memory')
     endif

     dmft_g = czero
     dmft_s = czero
     dmft_h = czero

     return
  end subroutine df_allocate_memory_dmft

  subroutine df_allocate_memory_dual()
     implicit none

     allocate(dual_g(nkpts,nffrq,norbs), stat=istat)
     allocate(dual_s(nkpts,nffrq,norbs), stat=istat)
     allocate(dual_b(nkpts,nffrq,norbs), stat=istat)

     if ( istat /= 0 ) then
         call s_print_error('df_allocate_memory_dual','can not allocate enough memory')
     endif

     dual_g = czero
     dual_s = czero
     dual_b = czero

     return
  end subroutine df_allocate_memory_dual

  subroutine df_allocate_memory_latt()
     implicit none

     allocate(latt_g(nkpts,nffrq,norbs), stat=istat)
     allocate(latt_s(nkpts,nffrq,norbs), stat=istat)

     if ( istat /= 0 ) then
         call s_print_error('df_allocate_memory_latt','can not allocate enough memory')
     endif

     latt_g = czero
     latt_s = czero

     return
  end subroutine df_allocate_memory_latt

  subroutine df_allocate_memory_vert()
     implicit none

     allocate(vert_d(nffrq,nffrq,nbfrq,norbs), stat=istat)
     allocate(vert_m(nffrq,nffrq,nbfrq,norbs), stat=istat)

     if ( istat /= 0 ) then
         call s_print_error('df_allocate_memory_vert','can not allocate enough memory')
     endif

     vert_d = czero
     vert_m = czero

     return
  end subroutine df_allocate_memory_vert

!!========================================================================
!!>>> deallocate memory subroutines                                    <<<
!!========================================================================

  subroutine df_deallocate_memory()
     implicit none

     if ( allocated(dmft_g) ) deallocate(dmft_g)
     if ( allocated(dmft_s) ) deallocate(dmft_s)
     if ( allocated(dmft_h) ) deallocate(dmft_h)

     if ( allocated(dual_g) ) deallocate(dual_g)
     if ( allocated(dual_s) ) deallocate(dual_s)
     if ( allocated(dual_b) ) deallocate(dual_b)

     if ( allocated(latt_g) ) deallocate(latt_g)
     if ( allocated(latt_s) ) deallocate(latt_s)

     if ( allocated(vertex_d) ) deallocate(vertex_d)
     if ( allocated(vertex_m) ) deallocate(vertex_m)

     return
  end subroutine df_deallocate_memory

  end module df_context
