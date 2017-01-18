
  module df_context
     use constants

     use df_control

     integer, private :: istat

!! dmft variables
! dmft hybridization function
     complex(dp), public, save, allocatable :: dmft_h(:,:)

! dmft green's function
     complex(dp), public, save, allocatable :: dmft_g(:,:)

! dmft self-energy function
     complex(dp), public, save, allocatable :: dmft_s(:,:)


!! dual variables
! dual green's function
     complex(dp), public, save, allocatable :: dual_g(:,:)

! dual self-energy function
     complex(dp), public, save, allocatable :: dual_s(:,:)

! dual bare green's function
     complex(dp), public, save, allocatable :: dual_b(:,:)


!! lattice variables
! lattice green's function
     complex(dp), public, save, allocatable :: latt_g(:,:)

! lattice self-energy function
     complex(dp), public, save, allocatable :: latt_s(:,:)

     public :: df_allocate_memory
     public :: df_deallocate_memory

  contains

  subroutine df_allocate_memory()
     implicit none

     allocate(dmft_g(nffrq,norbs), stat=istat)
     allocate(dmft_s(nffrq,norbs), stat=istat)
     allocate(dmft_h(nffrq,norbs), stat=istat)

     if ( istat /= 0 ) then
         call s_print_error('df_allocate_memory','can not allocate enough memory')
     endif

     dmft_g = czero
     dmft_s = czero
     dmft_h = czero

     return
  end subroutine df_allocate_memory

  subroutine df_deallocate_memory()
     implicit none

     if ( allocated(dmft_g) ) deallocate(dmft_g)
     if ( allocated(dmft_s) ) deallocate(dmft_s)
     if ( allocated(dmft_h) ) deallocate(dmft_h)

     return
  end subroutine df_deallocate_memory

  end module df_context
