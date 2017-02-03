
  subroutine df_fft1d(op, nx, fin, fout)
     use iso_c_binding
     use constants

     use df_control

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
  end subroutine df_fft1d

  subroutine df_fft2d(op, nx, ny, fin, fout)
     use iso_c_binding
     use constants

     use df_control

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
  end subroutine df_fft2d

  subroutine df_fft3d(op, nx, ny, nz, fin, fout)
     use iso_c_binding
     use constants

     use df_control

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
  end subroutine df_fft3d
