
  subroutine df_diagram()
     implicit none

     return
  end subroutine df_diagram

  subroutine df_static_bubble()
     use constants

     use df_control
     use df_context

     implicit none

     print *, 'static'
     return
  end subroutine df_static_bubble

  subroutine df_bubble()
     use constants

     use df_control
     use df_context

     implicit none

     print *, 'dynamic'
     return
  end subroutine df_bubble
