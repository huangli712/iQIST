
  subroutine df_run()
     use df_control, only : ndfit, nbsit

     implicit none

     integer :: df_it
     integer :: bs_it

     call df_dual_wssf()
     do df_it=1,ndfit
         call df_bubble()

         do bs_it=1,nbsit
             call df_full_vert()
         enddo

         call df_diagram()
         call df_dual_sigf()
         call df_dual_grnf()
     enddo
     call df_latt_grnf()
     call df_latt_sigf()

     call df_dmft_grnf()
     call df_dmft_sigf()
     call df_dmft_hybf()

     call df_spin_susc()
     call df_char_susc()

     return
  end subroutine df_run



  subroutine df_dmft_grnf()
     implicit none

     return
  end subroutine df_dmft_grnf
  
  subroutine df_dmft_sigf()
     implicit none

     return
  end subroutine df_dmft_sigf

  subroutine df_dmft_hybf()
     implicit none

     return
  end subroutine df_dmft_hybf

  subroutine df_dual_grnf()
     implicit none

     return
  end subroutine df_dual_grnf

  subroutine df_dual_sigf()
     implicit none

     return
  end subroutine df_dual_sigf

  subroutine df_dual_wssf()
     implicit none

     return
  end subroutine df_dual_wssf

  subroutine df_latt_grnf()
     implicit none

     return
  end subroutine df_latt_grnf

  subroutine df_latt_sigf()
     implicit none

     return
  end subroutine df_latt_sigf

  subroutine df_full_vert()
     implicit none

     return
  end subroutine df_full_vert
