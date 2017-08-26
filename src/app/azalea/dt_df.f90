!!!-----------------------------------------------------------------------
!!! project : azalea
!!! program : dt_df_core
!!!           dt_df_dual
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
     integer :: o
     integer :: v
     integer :: w

     real(dp) :: om
     complex(dp) :: vr(nkpts)
     complex(dp) :: gr(nkpts)

     complex(dp), allocatable :: gstp(:,:,:)
     complex(dp), allocatable :: gnew(:,:,:)
     complex(dp), allocatable :: gvrt(:,:,:)
     complex(dp), allocatable :: gcnv(:,:,:)

     complex(dp), allocatable :: Bmat(:,:)
     complex(dp), allocatable :: vertexM(:,:)
     complex(dp), allocatable :: vertexD(:,:)
     complex(dp), allocatable :: gammaM(:,:)
     complex(dp), allocatable :: gammaM2(:,:)

     allocate(gstp(nffrq,norbs,nkpts))
     allocate(gnew(nffrq,norbs,nkpts))
     allocate(gvrt(nffrq,norbs,nkpts))
     allocate(gcnv(nffrq,norbs,nkpts))

     allocate(Bmat(nffrq,nffrq))
     allocate(vertexM(nffrq,nffrq))
     allocate(vertexD(nffrq,nffrq))
     allocate(gammaM(nffrq,nffrq))
     allocate(gammaM2(nffrq,nffrq))

     DF_LOOP: do it=1,ndfit

         write(mystd,'(2X,A,I3)') 'Ladder Dual Fermion Iteration:', it

         V_LOOP: do v=1,nbfrq

             om = bmesh(v)
             write(mystd,'(2X,A,F12.6)') 'Bosonic Frequency:', om

             call cat_fill_k(dual_g, gstp, om)
             call cat_dia_2d(dual_g, gstp, gcnv)
             gvrt = czero

             O_LOOP: do o=1,norbs

                 vertexM = vert_m(:,:,v)
                 vertexD = vert_d(:,:,v)

                 K_LOOP: do k=1,nkpts

                     call s_diag_z(nffrq, gcnv(:,o,k), Bmat)

                     call cat_bse_solver(Bmat, vertexM, gammaM)
                     call cat_bse_iterator(1, one, Bmat, vertexM, gammaM2)
                     gammaM =  gammaM - half * gammaM2
                     call s_vecadd_z(nffrq, gvrt(:,o,k), gammaM, half * 3.0_dp)

                     call cat_bse_solver(Bmat, vertexD, gammaM)
                     call cat_bse_iterator(1, one, Bmat, vertexD, gammaM2)
                     gammaM = gammaM - half * gammaM2
                     call s_vecadd_z(nffrq, gvrt(:,o,k), gammaM, half * 1.0_dp)

                 enddo K_LOOP

                 do w=1,nffrq
                     call cat_fft_2d(+1, nkp_x, nkp_y, gvrt(w,o,:), vr)
                     call cat_fft_2d(-1, nkp_x, nkp_y, gstp(w,o,:), gr)
                     gr = vr * gr / real(nkpts * nkpts)
                     call cat_fft_2d(+1, nkp_x, nkp_y, gr, vr)
                     dual_s(w,o,:) = dual_s(w,o,:) + vr / beta
                 enddo

             enddo O_LOOP

         enddo V_LOOP

         call dt_df_dual(+1, gnew, dual_s, dual_b)
         call s_mix_z( size(gnew), dual_g, gnew, dfmix)

         dual_g = gnew
         dual_s = czero

         write(mystd,*)
     enddo DF_LOOP

     call dt_df_dual(-1, dual_g, dual_s, dual_b)

     do w=1,nffrq
       print *, w, fmesh(w)
       print *, dual_s(w,1,:)
     enddo

     deallocate(gstp)
     deallocate(gnew)
     deallocate(gvrt)
     deallocate(gcnv)
     deallocate(Bmat)
     deallocate(vertexM)
     deallocate(vertexD)
     deallocate(gammaM)
     deallocate(gammaM2)

     return
  end subroutine dt_df_core

!!
!! @sub dt_df_dual
!!
!! try to calculate the dual green's function or self-energy function by
!! using the dyson equation
!!
  subroutine dt_df_dual(op, dual_g, dual_s, dual_b)
     use constants, only : dp
     use constants, only : one

     use control, only : norbs
     use control, only : nffrq
     use control, only : nkpts

     implicit none

! external arguments
     integer, intent(in) :: op

     complex(dp), intent(inout) :: dual_g(nffrq,norbs,nkpts)
     complex(dp), intent(inout) :: dual_s(nffrq,norbs,nkpts)
     complex(dp), intent(inout) :: dual_b(nffrq,norbs,nkpts)

     if ( op == 1 ) then
         dual_g = one / ( one / dual_b - dual_s )
     else
         dual_s = one / dual_b - one / dual_g
     endif ! back if ( op == 1 ) block

     return
  end subroutine dt_df_dual

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
