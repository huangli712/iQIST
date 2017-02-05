
  subroutine df_run()
     use constants

     use df_control
     use df_context

     implicit none

! local variables
     integer :: i
     integer :: j
     integer :: k
     integer :: l
     integer :: n

     real(dp) :: w
     complex(dp) :: mval
     complex(dp) :: dval
     complex(dp) :: vr(nkpts)
     complex(dp) :: gr(nkpts)
     complex(dp), allocatable :: gshift(:,:,:)
     complex(dp), allocatable :: full_v(:,:,:)
     complex(dp), allocatable :: bubble(:,:,:)
     complex(dp), allocatable :: bubbleM(:,:)
     complex(dp), allocatable :: vertexM(:,:)
     complex(dp), allocatable :: vertexD(:,:)
     complex(dp), allocatable :: gammaM(:,:)
     complex(dp), allocatable :: gammaD(:,:)
     complex(dp), allocatable :: gammaM2(:,:)
     complex(dp), allocatable :: gammaD2(:,:)

     allocate(gshift(nkpts,nffrq,norbs))
     allocate(full_v(nkpts,nffrq,norbs))
     allocate(bubble(nkpts,nffrq,norbs))
     allocate(bubbleM(nffrq,nffrq))
     allocate(vertexM(nffrq,nffrq))
     allocate(vertexD(nffrq,nffrq))
     allocate(gammaM(nffrq,nffrq))
     allocate(gammaD(nffrq,nffrq))
     allocate(gammaM2(nffrq,nffrq))
     allocate(gammaD2(nffrq,nffrq))

     DF_LOOP: do i=1,ndfit
         write(mystd,'(2X,A,I3)') 'Ladder Dual Fermion Iteration:', i

         Q_LOOP: do j=2,nbfrq-1
             w = bmesh(j)
             write(mystd,'(2X,A,F12.6)') 'Bosonic Frequency:', w
             call cat_dual_shift(dual_g, gshift, w)

             if ( abs(w - zero) < epss ) then
                 call df_static_bubble(bubble, w)
             else
                 call df_bubble(bubble, w)
             endif

             vertexM = vert_m(:,:,j)
             vertexD = vert_d(:,:,j)

             K_LOOP: do k=1,nkpts
                 print *, 'K:', k
                 call s_diag_z(nffrq, bubble(k,:,1), bubbleM)
                 call df_bse_solver(bubbleM, vertexM, gammaM)
                 call df_bse_solver(bubbleM, vertexD, gammaD)

                 call df_bse_solver_iter(1, 1.0_dp, bubbleM, vertexM, gammaM2)
                 call df_bse_solver_iter(1, 1.0_dp, bubbleM, vertexD, gammaD2)

                 W_LOOP: do l=1,nffrq
                     mval = gammaM(l,l) - 0.5*gammaM2(l,l)
                     dval = gammaD(l,l) - 0.5*gammaD2(l,l)
                     full_v(k,l,1) = 0.5 * (3.0 * mval + dval) 
                     full_v(k,l,2) = 0.5 * (3.0 * mval + dval) 
                 enddo W_LOOP
             enddo K_LOOP

             do n=1,nffrq
                 call df_fft2d(+1, nkp_x, nkp_y, full_v(:,n,1), vr)
                 call df_fft2d(-1, nkp_x, nkp_y, gshift(:,n,1), gr)
                 gr = vr * gr / real(nkpts * nkpts)
                 call df_fft2d(+1, nkp_x, nkp_y, gr, vr)
                 dual_s(:,n,1) = dual_s(:,n,1) + vr / beta
                 print *, n, fmesh(n)
                 print *, dual_s(:,n,1)
             enddo

         enddo Q_LOOP
         STOP

         write(mystd,*)
     enddo DF_LOOP

     deallocate(gshift)
     deallocate(full_v)
     deallocate(bubble)
     deallocate(bubbleM)
     deallocate(vertexM)
     deallocate(vertexD)
     deallocate(gammaM)
     deallocate(gammaD)
     deallocate(gammaM2)
     deallocate(gammaD2)

     return
  end subroutine df_run

  subroutine cat_dual_shift(dual_in, dual_out, shift)
     use constants

     use df_control
     use df_context

     implicit none

     real(dp), intent(in) :: shift
     complex(dp), intent(in) :: dual_in(nkpts,nffrq,norbs)
     complex(dp), intent(out) :: dual_out(nkpts,nffrq,norbs)

     integer :: i
     integer :: j
     integer :: k
     real(dp) :: fw

     do i=1,norbs
         do j=1,nffrq
             fw = fmesh(j) + shift
             k = floor( (fw * beta / pi + nffrq + one) / two + 0.5 )
             if ( k >= 1 .and. k <= nffrq ) then
                 dual_out(:,j,i) = dual_in(:,k,i)
             else
                 dual_out(:,j,i) = czero
             endif
         enddo ! over j={1,nffrq} loop
     enddo ! over i={1,norbs} loop

     return
  end subroutine cat_dual_shift
