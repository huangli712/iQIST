
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
     integer :: k

     real(dp) :: fw

     complex(dp) :: gk(nkpts)
     complex(dp) :: gr(nkpts), gr1(nkpts), gr2(nkpts)
     complex(dp) :: gs(nkpts,nffrq,norbs)

     do i=1,norbs
         do j=1,nffrq
             fw = fmesh(j) + w
             if      ( fw > fmesh(nffrq) ) then
                 gs(:,j,i) = czero
             else if ( fw < fmesh(  1  ) ) then
                 gs(:,j,i) = czero
             else
                 k = ceiling( (fw * beta / pi + nffrq + one) / two )
                 gs(:,j,i) = dual_g(:,k,i)
             endif
         enddo ! over j={1,nffrq} loop
     enddo ! over i={1,norbs} loop

     do i=1,norbs
         do j=1,nffrq
             gk = dual_g(:,j,i)
             gr1 = czero
             call df_fft2d(+1, nkp_x, nkp_y, gk, gr1) ! gk -> gr

             gk = gs(:,j,i)
             gr2 = czero
             call df_fft2d(+1, nkp_x, nkp_y, gk, gr2) ! gk -> gr

             gr = gr1 * gr2
             call df_fft2d(-1, nkp_x, nkp_y, gr, gk) ! gr -> gk
             bubble(:,j,i) = -gk
         enddo ! over j={1,nffrq} loop
     enddo ! over i={1,norbs} loop
     bubble = bubble / real(nkpts * nkpts * beta)

     return
  end subroutine df_bubble
