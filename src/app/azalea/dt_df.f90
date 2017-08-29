!!!-----------------------------------------------------------------------
!!! project : azalea
!!! program : dt_run
!!! source  : dt_df.f90
!!! type    : subroutines
!!! author  : li huang (email:lihuang.dmft@gmail.com)
!!! history : 09/16/2009 by li huang (created)
!!!           08/30/2017 by li huang (last modified)
!!! purpose :
!!! status  : unstable
!!! comment :
!!!-----------------------------------------------------------------------

  subroutine dt_run()
     use constants, only : dp
     use constants, only : zero, epss
     use constants, only : mystd

     use control, only : norbs
     use control, only : nffrq, nbfrq
     use control, only : nkpts, nkp_x, nkp_y
     use control, only : ndfit
     use control, only : beta

     use context, only : fmesh, bmesh
     use context, only : dual_g, dual_s
     use context, only : vert_d, vert_m

     implicit none

! local variables
! loop index
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

         Q_LOOP: do j=1,nbfrq-1
             w = bmesh(j)
             write(mystd,'(2X,A,F12.6)') 'Bosonic Frequency:', w
             call cat_dual_shift(dual_g, gshift, w)

             if ( abs(w - zero) < epss ) then
                 call dt_static_bubble(bubble, w)
             else
                 call dt_bubble(bubble, w)
             endif

             vertexM = vert_m(:,:,j)
             vertexD = vert_d(:,:,j)

             K_LOOP: do k=1,nkpts
                 print *, 'K:', k
                 call s_diag_z(nffrq, bubble(k,:,1), bubbleM)
                 call dt_bse_solver(bubbleM, vertexM, gammaM)
                 call dt_bse_solver(bubbleM, vertexD, gammaD)

                 call dt_bse_solver_iter(1, 1.0_dp, bubbleM, vertexM, gammaM2)
                 call dt_bse_solver_iter(1, 1.0_dp, bubbleM, vertexD, gammaD2)

                 W_LOOP: do l=1,nffrq
                     mval = gammaM(l,l) - 0.5*gammaM2(l,l)
                     dval = gammaD(l,l) - 0.5*gammaD2(l,l)
                     full_v(k,l,1) = 0.5 * (3.0 * mval + dval) 
                     full_v(k,l,2) = 0.5 * (3.0 * mval + dval) 
                 enddo W_LOOP
             enddo K_LOOP

             do n=1,nffrq
                 call dt_fft2d(+1, nkp_x, nkp_y, full_v(:,n,1), vr)
                 call dt_fft2d(-1, nkp_x, nkp_y, gshift(:,n,1), gr)
                 gr = vr * gr / real(nkpts * nkpts)
                 call dt_fft2d(+1, nkp_x, nkp_y, gr, vr)
                 dual_s(:,n,1) = dual_s(:,n,1) + vr / beta
                 print *, n, fmesh(n)
                 print *, dual_s(:,n,1)
             enddo
         STOP

         enddo Q_LOOP

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
  end subroutine dt_run

  subroutine cat_dual_shift(dual_in, dual_out, shift)
     use constants

     use control
     use context

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

  subroutine dt_static_bubble(bubble, w)
     use constants

     use control
     use context

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
             call dt_fft2d(+1, nkp_x, nkp_y, gk, gr) ! gk -> gr
             gr = gr * gr
             call dt_fft2d(-1, nkp_x, nkp_y, gr, gk) ! gr -> gk
             bubble(:,j,i) = -gk
         enddo ! over j={1,nffrq} loop
     enddo ! over i={1,norbs} loop
     bubble = bubble / real(nkpts * nkpts * beta)

     return
  end subroutine dt_static_bubble

  subroutine dt_bubble(bubble, w)
     use constants

     use control
     use context

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
             k = floor( (fw * beta / pi + nffrq + one) / two + 0.5 )
             if ( k >= 1 .and. k <= nffrq ) then
                 gs(:,j,i) = dual_g(:,k,i)
             else
                 gs(:,j,i) = czero
             endif
         enddo ! over j={1,nffrq} loop
     enddo ! over i={1,norbs} loop

     do i=1,norbs
         do j=1,nffrq
             gk = dual_g(:,j,i)
             gr1 = czero
             call dt_fft2d(+1, nkp_x, nkp_y, gk, gr1) ! gk -> gr

             gk = gs(:,j,i)
             gr2 = czero
             call dt_fft2d(+1, nkp_x, nkp_y, gk, gr2) ! gk -> gr

             gr = gr1 * gr2
             call dt_fft2d(-1, nkp_x, nkp_y, gr, gk) ! gr -> gk
             bubble(:,j,i) = -gk
         enddo ! over j={1,nffrq} loop
     enddo ! over i={1,norbs} loop
     bubble = bubble / real(nkpts * nkpts * beta)

     return
  end subroutine dt_bubble

  subroutine dt_bse_solver(bubbleM, vertexM, gammaM)
     use constants

     use control
     use context

     implicit none

! external arguments
     complex(dp), intent(in) :: bubbleM(nffrq,nffrq)
     complex(dp), intent(in) :: vertexM(nffrq,nffrq)
     complex(dp), intent(out) :: gammaM(nffrq,nffrq)

! local variables
     complex(dp) :: Imat(nffrq,nffrq)
     complex(dp) :: v4chi(nffrq,nffrq)
     complex(dp) :: zdet

     call s_identity_z(nffrq, Imat)
     v4chi = Imat - matmul(vertexM,bubbleM)
     gammaM = v4chi
     call s_det_z(nffrq, v4chi, zdet)
     call s_inv_z(nffrq, gammaM)
     gammaM = matmul(gammaM, vertexM) 

     !print *, zdet
     return
  end subroutine dt_bse_solver

  subroutine dt_bse_solver_iter(niter, mix, bubbleM, vertexM, gammaM)
     use constants

     use control
     use context

     implicit none

! external arguments
     integer, intent(in) :: niter
     real(dp), intent(in) :: mix

     complex(dp), intent(in) :: bubbleM(nffrq,nffrq)
     complex(dp), intent(in) :: vertexM(nffrq,nffrq)
     complex(dp), intent(out) :: gammaM(nffrq,nffrq)

! local variables
     integer :: it
     complex(dp) :: V4old(nffrq,nffrq)
     complex(dp) :: V4chi(nffrq,nffrq)
     complex(dp) :: diff

     gammaM = vertexM
     V4old = vertexM
     V4chi = matmul(vertexM, bubbleM)

     do it=1,niter
         if ( it == niter ) then
             gammaM = ( vertexM + matmul(V4chi,V4old) ) * mix + (1.0 - mix) * V4old
         else
             gammaM = ( matmul(V4chi,V4old) ) * mix + (1.0 - mix) * V4old
         endif
         diff = abs(sum(gammaM - V4old)) / real(nffrq * nffrq)
         V4old = gammaM
     enddo

     return
  end subroutine dt_bse_solver_iter

! calculate resulting observables
  subroutine dt_spin_susc()
     implicit none

     return
  end subroutine dt_spin_susc

  subroutine dt_char_susc()
     implicit none

     return
  end subroutine dt_char_susc
