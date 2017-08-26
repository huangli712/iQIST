!!!-----------------------------------------------------------------------
!!! project : azalea
!!! program : dt_run
!!! source  : dt_run.f90
!!! type    : subroutines
!!! author  : li huang (email:lihuang.dmft@gmail.com)
!!! history : 01/10/2018 by li huang (created)
!!!           01/10/2018 by li huang (last modified)
!!! purpose : main entry of the program.
!!! status  : unstable
!!! comment :
!!!-----------------------------------------------------------------------

!!
!! @sub dt_run
!!
!! core computational engine, it is used to dispatch the jobs
!!
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
