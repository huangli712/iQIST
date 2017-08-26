  subroutine dt_run()
     use control, only : isdia

     implicit none

     DT_CORE: &
     select case ( isdia )

         case (1)
             call dt_df_std()

         case (2)
             call dt_df_ladder()

         case default
             call s_print_error('dt_run','this feature is not implemented')

     end select DT_CORE

     return
  end subroutine dt_run
