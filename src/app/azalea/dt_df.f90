!!!-----------------------------------------------------------------------
!!! project : azalea
!!! program : dt_df_core
!!!           dt_df_schi
!!!           dt_df_cchi
!!! source  : dt_df.f90
!!! type    : subroutines
!!! author  : li huang (email:lihuang.dmft@gmail.com)
!!! history : 09/16/2009 by li huang (created)
!!!           01/04/2018 by li huang (last modified)
!!! purpose : main subroutines for the dual fermion framework.
!!! status  : unstable
!!! comment :
!!!-----------------------------------------------------------------------

!!
!! @sub dt_df_core
!!
!! implement the dual fermion framework
!!
  subroutine dt_df_core()
     use constants, only : dp
     use constants, only : zero, one, half, epss, czero, cone
     use constants, only : mystd

     use control, only : norbs
     use control, only : nffrq, nbfrq
     use control, only : nkpts, nkp_x, nkp_y
     use control, only : ndfit, dfmix
     use control, only : beta

     use context, only : fmesh, bmesh
     use context, only : dual_g, dual_s, dual_b
     use context, only : vert_d, vert_m

     implicit none

! local variables
! loop index
     integer :: it
     integer :: k
     !integer :: m
     integer :: n
     integer :: v, w

     real(dp) :: om
     complex(dp) :: mval
     complex(dp) :: dval
     complex(dp) :: vr(nkpts)
     complex(dp) :: gr(nkpts)

     complex(dp), allocatable :: gshift(:,:,:)
     complex(dp), allocatable :: dual_g_new(:,:,:)
     complex(dp), allocatable :: full_v(:,:,:)
     complex(dp), allocatable :: bubble(:,:,:)
     complex(dp), allocatable :: bubbleM(:,:)
     complex(dp), allocatable :: vertexM(:,:)
     complex(dp), allocatable :: vertexD(:,:)
     complex(dp), allocatable :: gammaM(:,:)
     complex(dp), allocatable :: gammaD(:,:)
     complex(dp), allocatable :: gammaM2(:,:)
     complex(dp), allocatable :: gammaD2(:,:)

     allocate(dual_g_new(nffrq,norbs,nkpts))
     allocate(gshift(nffrq,norbs,nkpts))
     allocate(full_v(nffrq,norbs,nkpts))
     allocate(bubble(nffrq,norbs,nkpts))
     allocate(bubbleM(nffrq,nffrq))
     allocate(vertexM(nffrq,nffrq))
     allocate(vertexD(nffrq,nffrq))
     allocate(gammaM(nffrq,nffrq))
     allocate(gammaD(nffrq,nffrq))
     allocate(gammaM2(nffrq,nffrq))
     allocate(gammaD2(nffrq,nffrq))

     DF_LOOP: do it=1,ndfit
         write(mystd,'(2X,A,I3)') 'Ladder Dual Fermion Iteration:', it

         V_LOOP: do v=1,nbfrq
             om = bmesh(v)
             write(mystd,'(2X,A,F12.6)') 'Bosonic Frequency:', om
             call cat_fill_k(dual_g, gshift, om)
             call cat_dia_2d(dual_g, gshift, bubble)

             vertexM = vert_m(:,:,v)
             vertexD = vert_d(:,:,v)

             K_LOOP: do k=1,nkpts
                 call s_diag_z(nffrq, bubble(:,1,k), bubbleM)
                 call cat_bse_solver(bubbleM, vertexM, gammaM)
                 call cat_bse_solver(bubbleM, vertexD, gammaD)

                 call cat_bse_iterator(1, one, bubbleM, vertexM, gammaM2)
                 call cat_bse_iterator(1, one, bubbleM, vertexD, gammaD2)

                 W_LOOP: do w=1,nffrq
                     mval = gammaM(w,w) - half * gammaM2(w,w)
                     dval = gammaD(w,w) - half * gammaD2(w,w)
                     full_v(w,1,k) = half * (3.0 * mval + dval) 
                     full_v(w,2,k) = half * (3.0 * mval + dval) 
                 enddo W_LOOP
             enddo K_LOOP

             do n=1,nffrq
                 call cat_fft_2d(+1, nkp_x, nkp_y, full_v(n,1,:), vr)
                 call cat_fft_2d(-1, nkp_x, nkp_y, gshift(n,1,:), gr)
                 gr = vr * gr / real(nkpts * nkpts)
                 call cat_fft_2d(+1, nkp_x, nkp_y, gr, vr)
                 dual_s(n,1,:) = dual_s(n,1,:) + vr / beta
                 dual_s(n,2,:) = dual_s(n,2,:) + vr / beta
             enddo

         enddo V_LOOP

     do k=1,nkpts
         do v=1,norbs
             do w=1,nffrq
                 dual_g_new(w,v,k) = one / ( one / dual_b(w,v,k) - dual_s(w,v,k) ) * dfmix + dual_g(w,v,k) * ( one - dfmix )
             enddo ! over i={1,nffrq} loop
         enddo ! over j={1,norbs} loop
     enddo ! over k={1,nkpts} loop

         dual_g = dual_g_new
         dual_s = czero

         write(mystd,*)
     enddo DF_LOOP

     do k=1,nkpts
         do v=1,norbs
             do w=1,nffrq
                 dual_s(w,v,k) = one / dual_b(w,v,k) - one / dual_g(w,v,k)
             enddo ! over i={1,nffrq} loop
         enddo ! over j={1,norbs} loop
     enddo ! over k={1,nkpts} loop

     do n=1,nffrq
       print *, n, fmesh(n)
       print *, dual_s(n,1,:)
     enddo

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
  end subroutine dt_df_core

!!
!! @sub dt_df_schi
!!
!! calculate the spin susceptibility within the dual fermion framework
!!
  subroutine dt_df_schi()
     implicit none

     return
  end subroutine dt_df_schi

!!
!! @sub dt_df_cchi
!!
!! calculate the charge susceptibility within the dual fermion framework
!!
  subroutine dt_df_cchi()
     implicit none

     return
  end subroutine dt_df_cchi
