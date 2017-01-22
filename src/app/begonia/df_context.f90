  module df_dmft
     implicit none

!! dmft variables
! dmft green's function
     complex(dp), public, save, allocatable :: dmft_g(:,:)

! dmft self-energy function
     complex(dp), public, save, allocatable :: dmft_s(:,:)

! dmft hybridization function
     complex(dp), public, save, allocatable :: dmft_h(:,:)

  end module df_dmft

  module df_dual
     implicit none

!! dual variables
! dual green's function
     complex(dp), public, save, allocatable :: dual_g(:,:,:)

! dual self-energy function
     complex(dp), public, save, allocatable :: dual_s(:,:,:)

! dual bare green's function
     complex(dp), public, save, allocatable :: dual_b(:,:,:)

  end module df_dual

  module df_latt
     implicit none

!! lattice variables
! lattice green's function
     complex(dp), public, save, allocatable :: latt_g(:,:,:)

! lattice self-energy function
     complex(dp), public, save, allocatable :: latt_s(:,:,:)

  end module df_latt

  module df_vert
     implicit none

!! vertex variables
! density vertex
     complex(dp), public, save, allocatable :: vertex_d(:,:,:,:)

! magnetic vertex
     complex(dp), public, save, allocatable :: vertex_m(:,:,:,:)

  end module df_vert

  module df_context
     use constants

     use df_control

     integer, private :: istat











     public :: df_allocate_memory
     public :: df_deallocate_memory

  contains

  subroutine df_allocate_memory()
     implicit none

     allocate(dmft_g(nffrq,norbs), stat=istat)
     allocate(dmft_s(nffrq,norbs), stat=istat)
     allocate(dmft_h(nffrq,norbs), stat=istat)

     allocate(dual_g(nkpts,nffrq,norbs), stat=istat)
     allocate(dual_s(nkpts,nffrq,norbs), stat=istat)
     allocate(dual_b(nkpts,nffrq,norbs), stat=istat)

     allocate(latt_g(nkpts,nffrq,norbs), stat=istat)
     allocate(latt_s(nkpts,nffrq,norbs), stat=istat)

     allocate(vertex_d(nffrq,nffrq,nbfrq,norbs), stat=istat)
     allocate(vertex_m(nffrq,nffrq,nbfrq,norbs), stat=istat)

     if ( istat /= 0 ) then
         call s_print_error('df_allocate_memory','can not allocate enough memory')
     endif

     dmft_g = czero
     dmft_s = czero
     dmft_h = czero

     dual_g = czero
     dual_s = czero
     dual_b = czero

     latt_g = czero
     latt_s = czero

     vertex_d = czero
     vertex_m = czero

     return
  end subroutine df_allocate_memory

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
