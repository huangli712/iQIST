!!!-----------------------------------------------------------------------
!!! project : azalea
!!! program : cat_fill_l
!!!           cat_fill_k <<<---
!!!           cat_fft_1d
!!!           cat_fft_2d
!!!           cat_fft_3d <<<---
!!!           cat_dia_1d
!!!           cat_dia_2d
!!!           cat_dia_3d <<<---
!!!           cat_bse_solver
!!!           cat_bse_iterator
!!! source  : dt_util.f90
!!! type    : subroutines
!!! author  : li huang (email:lihuang.dmft@gmail.com)
!!! history : 10/01/2008 by li huang (created)
!!!           09/06/2017 by li huang (last modified)
!!! purpose :
!!! status  : unstable
!!! comment :
!!!-----------------------------------------------------------------------

!!========================================================================
!!>>> shift matsubara frequency                                        <<<
!!========================================================================

!!
!! @sub cat_fill_l
!!
!! try to fill G(\nu + \omega) by G(\nu)
!!
  subroutine cat_fill_l(gin, gout, shift)
     use constants, only : dp
     use constants, only : one, two, half, pi, czero

     use control, only : norbs
     use control, only : nffrq
     use control, only : beta

     use context, only : fmesh

     implicit none

! external arguments
! shifted frequency, \omega
     real(dp), intent(in) :: shift

! input array, G(\nu)
     complex(dp), intent(in)  :: gin(nffrq,norbs)

! filled array, G(\nu + \omega)
     complex(dp), intent(out) :: gout(nffrq,norbs)

! local variables
! loop index
     integer  :: i

! resultant index for \nu + \omega
     integer  :: k

! resultant frequency, \nu + \omega
     real(dp) :: w

     do i=1,nffrq
         w = fmesh(i) + shift
         k = floor( (w * beta / pi + nffrq + one) / two + half )
         if ( k >= 1 .and. k <= nffrq ) then
             gout(i,:) = gin(k,:)
         else
             gout(i,:) = czero
         endif ! back if ( k >= 1 .and. k <= nffrq ) block
     enddo ! over i={1,nffrq} loop

     return
  end subroutine cat_fill_l

!!
!! @sub cat_fill_k
!!
!! try to fill G(\nu + \omega, K) by G(\nu, K)
!!
  subroutine cat_fill_k(gin, gout, shift)
     use constants, only : dp
     use constants, only : one, two, half, pi, czero

     use control, only : norbs
     use control, only : nffrq
     use control, only : nkpts
     use control, only : beta

     use context, only : fmesh

     implicit none

! external arguments
! shifted frequency, \omega
     real(dp), intent(in) :: shift

! input array, G(\nu, K)
     complex(dp), intent(in)  :: gin(nffrq,norbs,nkpts)

! filled array, G(\nu + \omega, K)
     complex(dp), intent(out) :: gout(nffrq,norbs,nkpts)

! local variables
! loop index
     integer  :: i

! resultant index for \nu + \omega
     integer  :: k

! resultant frequency, \nu + \omega
     real(dp) :: w

     do i=1,nffrq
         w = fmesh(i) + shift
         k = floor( (w * beta / pi + nffrq + one) / two + half )
         if ( k >= 1 .and. k <= nffrq ) then
             gout(i,:,:) = gin(k,:,:)
         else
             gout(i,:,:) = czero
         endif ! back if ( k >= 1 .and. k <= nffrq ) block
     enddo ! over i={1,nffrq} loop

     return
  end subroutine cat_fill_k

!!========================================================================
!!>>> fast fourier transformation                                      <<<
!!========================================================================

!!
!! @sub cat_fft_1d
!!
!! conduct fast fourier transformation in 1d
!!
  subroutine cat_fft_1d(op, nx, fin, fout)
     use iso_c_binding
     use constants, only : dp

     implicit none

! import fftw header file
     include 'fftw3.f03'

! external arguments
! fft direction, forward or backward
     integer, intent(in) :: op

! size of operand
     integer, intent(in) :: nx

! operand
     complex(dp), intent(inout) :: fin(nx)
     complex(dp), intent(inout) :: fout(nx)

! local variables
! fftw descriptor handler
     type(c_ptr) :: plan

     select case (op)

         case (+1)
             plan = fftw_plan_dft_1d(nx, fin, fout, FFTW_FORWARD, FFTW_ESTIMATE)

         case (-1)
             plan = fftw_plan_dft_1d(nx, fin, fout, FFTW_BACKWARD, FFTW_ESTIMATE)

         case default
             call s_print_error('cat_fft_1d','unrecognized fft operation')

     end select

     call fftw_execute_dft(plan, fin, fout)
     call fftw_destroy_plan(plan)

     return
  end subroutine cat_fft_1d

!!
!! @sub cat_fft_2d
!!
!! conduct fast fourier transformation in 2d
!!
  subroutine cat_fft_2d(op, nx, ny, fin, fout)
     use iso_c_binding
     use constants, only : dp

     implicit none

! import fftw header file
     include 'fftw3.f03'

! external arguments
! fft direction, forward or backward
     integer, intent(in) :: op

! size of operand
     integer, intent(in) :: nx
     integer, intent(in) :: ny

! operand
     complex(dp), intent(inout) :: fin(nx,ny)
     complex(dp), intent(inout) :: fout(nx,ny)

! local variables
! fftw descriptor handler
     type(c_ptr) :: plan

     select case (op)

         case (+1)
             plan = fftw_plan_dft_2d(nx, ny, fin, fout, FFTW_FORWARD, FFTW_ESTIMATE)

         case (-1)
             plan = fftw_plan_dft_2d(nx, ny, fin, fout, FFTW_BACKWARD, FFTW_ESTIMATE)

         case default
             call s_print_error('cat_fft_2d','unrecognized fft operation')

     end select

     call fftw_execute_dft(plan, fin, fout)
     call fftw_destroy_plan(plan)

     return
  end subroutine cat_fft_2d

!!
!! @sub cat_fft_3d
!!
!! conduct fast fourier transformation in 3d
!!
  subroutine cat_fft_3d(op, nx, ny, nz, fin, fout)
     use iso_c_binding
     use constants, only : dp

     implicit none

! import fftw header file
     include 'fftw3.f03'

! external arguments
! fft direction, forward or backward
     integer, intent(in) :: op

! size of operand
     integer, intent(in) :: nx
     integer, intent(in) :: ny
     integer, intent(in) :: nz

! operand
     complex(dp), intent(inout) :: fin(nx,ny,nz)
     complex(dp), intent(inout) :: fout(nx,ny,nz)

! local variables
! fftw descriptor handler
     type(c_ptr) :: plan

     select case (op)

         case (+1)
             plan = fftw_plan_dft_3d(nx, ny, nz, fin, fout, FFTW_FORWARD, FFTW_ESTIMATE)

         case (-1)
             plan = fftw_plan_dft_3d(nx, ny, nz, fin, fout, FFTW_BACKWARD, FFTW_ESTIMATE)

         case default
             call s_print_error('cat_fft_3d','unrecognized fft operation')

     end select

     call fftw_execute_dft(plan, fin, fout)
     call fftw_destroy_plan(plan)

     return
  end subroutine cat_fft_3d

!!========================================================================
!!>>> calculate bubble diagram                                         <<<
!!========================================================================

!!
!! @cat_dia_1d
!!
!! calculate the two-particle bubble diagram, 1d version
!!
  subroutine cat_dia_1d(gin, ginp, chiq)
     use constants, only : dp
     use constants, only : czero

     use control, only : norbs
     use control, only : nffrq
     use control, only : nkpts, nkp_x
     use control, only : beta

     implicit none

! external arguments
! G(\nu, K)
     complex(dp), intent(in)  :: gin(nffrq,norbs,nkpts)

! G(\nu + \omega, K + Q)
     complex(dp), intent(in)  :: ginp(nffrq,norbs,nkpts)

! two-particle bubble diagram
     complex(dp), intent(out) :: chiq(nffrq,norbs,nkpts)

! local variables
! loop index
     integer :: i
     integer :: j

! dummy complex(dp) variables
     complex(dp) :: gk(nkpts)
     complex(dp) :: gr(nkpts)
     complex(dp) :: g1(nkpts)
     complex(dp) :: g2(nkpts)

     do i=1,norbs
         do j=1,nffrq
             gk = gin(j,i,:)
             g1 = czero
             call cat_fft_1d(+1, nkp_x, gk, g1) ! gk -> gr

             gk = ginp(j,i,:)
             g2 = czero
             call cat_fft_1d(+1, nkp_x, gk, g2) ! gk -> gr

             gr = g1 * g2
             gk = czero
             call cat_fft_1d(-1, nkp_x, gr, gk) ! gr -> gk
             chiq(j,i,:) = -gk
         enddo ! over j={1,nffrq} loop
     enddo ! over i={1,norbs} loop
     chiq = chiq / real(nkpts * nkpts * beta)

     return
  end subroutine cat_dia_1d

!!
!! @cat_dia_2d
!!
!! calculate the two-particle bubble diagram, 2d version
!!
  subroutine cat_dia_2d(gin, ginp, chiq)
     use constants, only : dp
     use constants, only : czero

     use control, only : norbs
     use control, only : nffrq
     use control, only : nkpts, nkp_x, nkp_y
     use control, only : beta

     implicit none

! external arguments
! G(\nu, K)
     complex(dp), intent(in)  :: gin(nffrq,norbs,nkpts)

! G(\nu + \omega, K + Q)
     complex(dp), intent(in)  :: ginp(nffrq,norbs,nkpts)

! two-particle bubble diagram
     complex(dp), intent(out) :: chiq(nffrq,norbs,nkpts)

! local variables
! loop index
     integer :: i
     integer :: j

! dummy complex(dp) variables
     complex(dp) :: gk(nkpts)
     complex(dp) :: gr(nkpts)
     complex(dp) :: g1(nkpts)
     complex(dp) :: g2(nkpts)

     do i=1,norbs
         do j=1,nffrq
             gk = gin(j,i,:)
             g1 = czero
             call cat_fft_2d(+1, nkp_x, nkp_y, gk, g1) ! gk -> gr

             gk = ginp(j,i,:)
             g2 = czero
             call cat_fft_2d(+1, nkp_x, nkp_y, gk, g2) ! gk -> gr

             gr = g1 * g2
             gk = czero
             call cat_fft_2d(-1, nkp_x, nkp_y, gr, gk) ! gr -> gk
             chiq(j,i,:) = -gk
         enddo ! over j={1,nffrq} loop
     enddo ! over i={1,norbs} loop
     chiq = chiq / real(nkpts * nkpts * beta)

     return
  end subroutine cat_dia_2d





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
