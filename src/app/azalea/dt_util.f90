!!!-----------------------------------------------------------------------
!!! project : azalea
!!! program : cat_fft1d
!!!           cat_fft2d
!!!           cat_fft3d
!!!           cat_dual_shift
!!!           dt_static_bubble
!!!           dt_bubble
!!!           dt_bse_solver
!!!           dt_bse_solver_iter
!!! source  : dt_util.f90
!!! type    : subroutines
!!! author  : li huang (email:lihuang.dmft@gmail.com)
!!! history : 10/01/2008 by li huang (created)
!!!           09/05/2017 by li huang (last modified)
!!! purpose :
!!! status  : unstable
!!! comment :
!!!-----------------------------------------------------------------------

  subroutine dt_fft1d(op, nx, fin, fout)
     use iso_c_binding
     use constants, only : dp

     implicit none

     include 'fftw3.f03'

! external arguments
     integer, intent(in) :: op
     integer, intent(in) :: nx

     complex(dp), intent(inout) :: fin(nx)
     complex(dp), intent(inout) :: fout(nx)

! local variables
     type(c_ptr) :: plan

     select case (op)

         case (+1)
             plan = fftw_plan_dft_1d(nx, fin, fout, FFTW_FORWARD, FFTW_ESTIMATE)

         case (-1)
             plan = fftw_plan_dft_1d(nx, fin, fout, FFTW_BACKWARD, FFTW_ESTIMATE)

         case default
             STOP

     end select

     call fftw_execute_dft(plan, fin, fout)
     call fftw_destroy_plan(plan)

     return
  end subroutine dt_fft1d

  subroutine dt_fft2d(op, nx, ny, fin, fout)
     use iso_c_binding
     use constants, only : dp

     implicit none

     include 'fftw3.f03'

! external arguments
     integer, intent(in) :: op
     integer, intent(in) :: nx
     integer, intent(in) :: ny

     complex(dp), intent(inout) :: fin(nx,ny)
     complex(dp), intent(inout) :: fout(nx,ny)

! local variables
     type(c_ptr) :: plan

     select case (op)

         case (+1)
             plan = fftw_plan_dft_2d(nx, ny, fin, fout, FFTW_FORWARD, FFTW_ESTIMATE)

         case (-1)
             plan = fftw_plan_dft_2d(nx, ny, fin, fout, FFTW_BACKWARD, FFTW_ESTIMATE)

         case default
             STOP

     end select

     call fftw_execute_dft(plan, fin, fout)
     call fftw_destroy_plan(plan)

     return
  end subroutine dt_fft2d

  subroutine dt_fft3d(op, nx, ny, nz, fin, fout)
     use iso_c_binding
     use constants, only : dp

     implicit none

     include 'fftw3.f03'

! external arguments
     integer, intent(in) :: op
     integer, intent(in) :: nx
     integer, intent(in) :: ny
     integer, intent(in) :: nz

     complex(dp), intent(inout) :: fin(nx,ny,nz)
     complex(dp), intent(inout) :: fout(nx,ny,nz)

! local variables
     type(c_ptr) :: plan

     select case (op)

         case (+1)
             plan = fftw_plan_dft_3d(nx, ny, nz, fin, fout, FFTW_FORWARD, FFTW_ESTIMATE)

         case (-1)
             plan = fftw_plan_dft_3d(nx, ny, nz, fin, fout, FFTW_BACKWARD, FFTW_ESTIMATE)

         case default
             STOP

     end select

     call fftw_execute_dft(plan, fin, fout)
     call fftw_destroy_plan(plan)

     return
  end subroutine dt_fft3d

  subroutine cat_dual_shift(dual_in, dual_out, shift)
     use constants, only : dp
     use constants, only : one, two, pi, czero

     use control, only : norbs
     use control, only : nffrq
     use control, only : nkpts
     use control, only : beta

     use context, only : fmesh

     implicit none

! external arguments
     real(dp), intent(in) :: shift
     complex(dp), intent(in) :: dual_in(nffrq,norbs,nkpts)
     complex(dp), intent(out) :: dual_out(nffrq,norbs,nkpts)

! local variables
     integer :: i
     integer :: j
     integer :: k
     real(dp) :: fw

     do i=1,norbs
         do j=1,nffrq
             fw = fmesh(j) + shift
             k = floor( (fw * beta / pi + nffrq + one) / two + 0.5 )
             if ( k >= 1 .and. k <= nffrq ) then
                 dual_out(j,i,:) = dual_in(k,i,:)
             else
                 dual_out(j,i,:) = czero
             endif
         enddo ! over j={1,nffrq} loop
     enddo ! over i={1,norbs} loop

     return
  end subroutine cat_dual_shift

  subroutine dt_static_bubble(bubble, w)
     use constants, only : dp
     use constants, only : czero

     use control, only : norbs
     use control, only : nffrq
     use control, only : nkpts, nkp_x, nkp_y
     use control, only : beta

     use context, only : dual_g

     implicit none

! external arguments
     real(dp), intent(in) :: w
     complex(dp), intent(out) :: bubble(nffrq,norbs,nkpts)

! local variables
     integer :: i
     integer :: j

     complex(dp) :: gk(nkpts)
     complex(dp) :: gr(nkpts)

     do i=1,norbs
         do j=1,nffrq
             gk = dual_g(j,i,:)
             gr = czero
             call dt_fft2d(+1, nkp_x, nkp_y, gk, gr) ! gk -> gr
             gr = gr * gr
             call dt_fft2d(-1, nkp_x, nkp_y, gr, gk) ! gr -> gk
             bubble(j,i,:) = -gk
         enddo ! over j={1,nffrq} loop
     enddo ! over i={1,norbs} loop
     bubble = bubble / real(nkpts * nkpts * beta)

     return
  end subroutine dt_static_bubble

  subroutine dt_bubble(bubble, w)
     use constants, only : dp
     use constants, only : one, two, pi, czero

     use control, only : norbs
     use control, only : nffrq
     use control, only : nkpts, nkp_x, nkp_y
     use control, only : beta

     use context, only : fmesh
     use context, only : dual_g

     implicit none

! external arguments
     real(dp), intent(in) :: w
     complex(dp), intent(out) :: bubble(nffrq,norbs,nkpts)

! local variables
     integer :: i
     integer :: j
     integer :: k

     real(dp) :: fw

     complex(dp) :: gk(nkpts)
     complex(dp) :: gr(nkpts), gr1(nkpts), gr2(nkpts)
     complex(dp) :: gs(nffrq,norbs,nkpts)

     do i=1,norbs
         do j=1,nffrq
             fw = fmesh(j) + w
             k = floor( (fw * beta / pi + nffrq + one) / two + 0.5 )
             if ( k >= 1 .and. k <= nffrq ) then
                 gs(j,i,:) = dual_g(k,i,:)
             else
                 gs(j,i,:) = czero
             endif
         enddo ! over j={1,nffrq} loop
     enddo ! over i={1,norbs} loop

     do i=1,norbs
         do j=1,nffrq
             gk = dual_g(j,i,:)
             gr1 = czero
             call dt_fft2d(+1, nkp_x, nkp_y, gk, gr1) ! gk -> gr

             gk = gs(j,i,:)
             gr2 = czero
             call dt_fft2d(+1, nkp_x, nkp_y, gk, gr2) ! gk -> gr

             gr = gr1 * gr2
             call dt_fft2d(-1, nkp_x, nkp_y, gr, gk) ! gr -> gk
             bubble(j,i,:) = -gk
         enddo ! over j={1,nffrq} loop
     enddo ! over i={1,norbs} loop
     bubble = bubble / real(nkpts * nkpts * beta)

     return
  end subroutine dt_bubble

  subroutine dt_bse_solver(bubbleM, vertexM, gammaM)
     use constants, only : dp

     use control, only : nffrq

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
     use constants, only : dp

     use control, only : nffrq

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
