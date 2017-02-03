
  subroutine df_util()
  end subroutine df_util

  subroutine df_fft_forward(gk, gr)
     use iso_c_binding

     use constants
     use df_control

     implicit none

     include 'fftw3.f03'

     complex(dp), intent(inout) :: gk(nkp_x, nkp_y)
     complex(dp), intent(out) :: gr(nkp_x, nkp_y)

     type(c_ptr) :: plan

     plan = fftw_plan_dft_2d( nkp_x, nkp_y, gk, gr, FFTW_FORWARD, FFTW_ESTIMATE)
     call fftw_execute_dft(plan, gk, gr)
     call fftw_destroy_plan(plan)

     return
  end subroutine df_fft_forward
  
  subroutine df_fft_backward(gr, gk)
     use iso_c_binding

     use constants
     use df_control

     implicit none

     include 'fftw3.f03'

     complex(dp), intent(inout) :: gr(nkp_x, nkp_y)
     complex(dp), intent(out) :: gk(nkp_x, nkp_y)

     type(c_ptr) :: plan

     plan = fftw_plan_dft_2d( nkp_x, nkp_y, gr, gk, FFTW_BACKWARD, FFTW_ESTIMATE)
     call fftw_execute_dft(plan, gr, gk)
     call fftw_destroy_plan(plan)

     return
  end subroutine df_fft_backward
