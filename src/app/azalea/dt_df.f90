!!!-----------------------------------------------------------------------
!!! project : azalea
!!! program : dt_df_core
!!!           dt_df_schi
!!!           dt_df_cchi
!!! source  : dt_df.f90
!!! type    : subroutines
!!! author  : li huang (email:lihuang.dmft@gmail.com)
!!! history : 09/16/2009 by li huang (created)
!!!           09/05/2017 by li huang (last modified)
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

     DF_LOOP: do i=1,ndfit
         write(mystd,'(2X,A,I3)') 'Ladder Dual Fermion Iteration:', i

         Q_LOOP: do j=1,nbfrq-1
             w = bmesh(j)
             write(mystd,'(2X,A,F12.6)') 'Bosonic Frequency:', w
             call cat_shift_k(dual_g, gshift, w)

             if ( abs(w - zero) < epss ) then
                 call dt_static_bubble(bubble, w)
             else
                 call dt_bubble(bubble, w)
             endif

             vertexM = vert_m(:,:,j)
             vertexD = vert_d(:,:,j)

             K_LOOP: do k=1,nkpts
                 print *, 'K:', k
                 call s_diag_z(nffrq, bubble(:,1,k), bubbleM)
                 call dt_bse_solver(bubbleM, vertexM, gammaM)
                 call dt_bse_solver(bubbleM, vertexD, gammaD)

                 call dt_bse_solver_iter(1, 1.0_dp, bubbleM, vertexM, gammaM2)
                 call dt_bse_solver_iter(1, 1.0_dp, bubbleM, vertexD, gammaD2)

                 W_LOOP: do l=1,nffrq
                     mval = gammaM(l,l) - 0.5*gammaM2(l,l)
                     dval = gammaD(l,l) - 0.5*gammaD2(l,l)
                     full_v(l,1,k) = 0.5 * (3.0 * mval + dval) 
                     full_v(l,2,k) = 0.5 * (3.0 * mval + dval) 
                 enddo W_LOOP
             enddo K_LOOP

             do n=1,nffrq
                 call cat_fft2d(+1, nkp_x, nkp_y, full_v(n,1,:), vr)
                 call cat_fft2d(-1, nkp_x, nkp_y, gshift(n,1,:), gr)
                 gr = vr * gr / real(nkpts * nkpts)
                 call cat_fft2d(+1, nkp_x, nkp_y, gr, vr)
                 dual_s(n,1,:) = dual_s(n,1,:) + vr / beta
                 print *, n, fmesh(n)
                 print *, dual_s(n,1,:)
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
