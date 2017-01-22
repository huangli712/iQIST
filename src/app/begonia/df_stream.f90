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

  subroutine df_latt_init()
     implicit none

     return
  end subroutine df_latt_init

  subroutine df_dmft_init()
     implicit none

     return
  end subroutine df_dmft_init

  subroutine df_vert_init()
     implicit none

     return
  end subroutine df_vert_init

  subroutine df_final_array()
     implicit none

     return
  end subroutine df_final_array
