! parse input files, readin data
  subroutine df_config()
     implicit none

     return
  end subroutine df_config

  subroutine df_setup_array()
     use df_context ! ALL

     implicit none

! allocate memory for df_context module
     call df_allocate_memory_dmft()
     call df_allocate_memory_dual()
     call df_allocate_memory_latt()
     call df_allocate_memory_vert()
     
     return
  end subroutine df_setup_array

  subroutine df_mesh_init()
     implicit none

     return
  end subroutine df_mesh_init

  subroutine df_dmft_init()
     use constants, only : dp, mytmp

     use df_control
     use df_context

     implicit none

     integer  :: i
     real(dp) :: r1, r2
     real(dp) :: c1, c2

     open(mytmp, file = 'df.dmft_g.in', form = 'formatted', status = 'unknown')
     do i=1,nffrq
     enddo ! over i={1,nffrq} loop
     close(mytmp)

     open(mytmp, file = 'df.dmft_h.in', form = 'formatted', status = 'unknown')
     do i=1,nffrq
     enddo ! over i={1,nffrq} loop
     close(mytmp)

     return
  end subroutine df_dmft_init

  subroutine df_dual_init()
     implicit none

     return
  end subroutine df_dual_init

  subroutine df_latt_init()
     implicit none

     return
  end subroutine df_latt_init

  subroutine df_vert_init()
     implicit none

     return
  end subroutine df_vert_init

  subroutine df_final_array()
     use df_context ! ALL

     implicit none

! deallocate memory for df_context module
     call df_deallocate_memory_dmft()
     call df_deallocate_memory_dual()
     call df_deallocate_memory_latt()
     call df_deallocate_memory_vert()

     return
  end subroutine df_final_array
