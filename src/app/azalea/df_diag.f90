
  subroutine df_diagram()
     implicit none

     return
  end subroutine df_diagram

  subroutine df_static_bubble()
     use constants

     use df_control
     use df_context

     implicit none

     integer :: i
     complex(dp) :: gk(nkpts), gr(nkpts)

     print *, 'static'

     do i=1,nffrq
         print *, i, fmesh(i)
         gk = dual_g(:,i,1)
         gr = czero
         call df_fft2d(+1, nkp_x, nkp_y, gk, gr)
         gr = gr * gr
         call df_fft2d(-1, nkp_x, nkp_y, gr, gk)
         print *, gk / float(nkpts*nkpts)
         print *
     enddo

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
