
  subroutine df_static_bubble(bubble, w)
     use constants

     use df_control
     use df_context

     implicit none

! external arguments
     real(dp), intent(in) :: w
     complex(dp), intent(out) :: bubble(nkpts,nffrq,norbs)

! local variables
     integer :: i
     integer :: j

     complex(dp) :: gk(nkpts)
     complex(dp) :: gr(nkpts)

     do i=1,norbs
         do j=1,nffrq
             gk = dual_g(:,j,i)
             gr = czero
             call df_fft2d(+1, nkp_x, nkp_y, gk, gr) ! gk -> gr
             gr = gr * gr
             call df_fft2d(-1, nkp_x, nkp_y, gr, gk) ! gr -> gk
             bubble(:,j,i) = -gk
         enddo ! over j={1,nffrq} loop
     enddo ! over i={1,norbs} loop
     bubble = bubble / real(nkpts * nkpts * beta)

     return
  end subroutine df_static_bubble

  subroutine df_bubble(bubble, w)
     use constants

     use df_control
     use df_context

     implicit none

! external arguments
     real(dp), intent(in) :: w
     complex(dp), intent(out) :: bubble(nkpts,nffrq,norbs)

! local variables
     integer :: i
     integer :: j

     complex(dp) :: gk(nkpts)
     complex(dp) :: gr(nkpts)

     return
  end subroutine df_bubble
