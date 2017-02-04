
  subroutine df_bse_solver(bubbleM, vertexM)
     use constants

     use df_control
     use df_context

     implicit none

! external arguments
     complex(dp), intent(in) :: bubbleM(nffrq,nffrq)
     complex(dp), intent(in) :: vertexM(nffrq,nffrq)

     return
  end subroutine df_bse_solver
