
  subroutine df_dump_dmft_grnf(rmesh, grnf)
     use constants, only : dp, zero, mytmp

     use df_control, only : norbs
     use df_control, only : nffrq

     implicit none

! external arguments
! matsubara frequency mesh
     real(dp), intent(in)    :: rmesh(nffrq)

! impurity green's function
     complex(dp), intent(in) :: grnf(nffrq,norbs)

! local variables
! loop index
     integer :: i
     integer :: j

! open data file: df.dmft_g.dat
     open(mytmp, file='df.dmft_g.dat', form='formatted', status='unknown')

! write it
     do i=1,norbs
         do j=1,nffrq
             write(mytmp,'(i6,5f16.8)') &
                 i, rmesh(j), real(grnf(j,i)), aimag(grnf(j,i)), zero, zero
         enddo ! over j={1,nffrq} loop
         write(mytmp,*) ! write empty lines
         write(mytmp,*)
     enddo ! over i={1,norbs} loop

! close data file
     close(mytmp)

     return
  end subroutine df_dump_dmft_grnf

  subroutine df_dump_dmft_sigf()
     implicit none

     return
  end subroutine df_dump_dmft_sigf

  subroutine df_dump_dmft_hybf()
     implicit none

     return
  end subroutine df_dump_dmft_hybf

  subroutine df_dump_dual_grnf()
     implicit none

     return
  end subroutine df_dump_dual_grnf

  subroutine df_dump_dual_sigf()
     implicit none

     return
  end subroutine df_dump_dual_sigf

  subroutine df_dump_dual_wssf()
     implicit none

     return
  end subroutine df_dump_dual_wssf

  subroutine df_dump_latt_grnf()
     implicit none

     return
  end subroutine df_dump_latt_grnf

  subroutine df_dump_latt_sigf()
     implicit none

     return
  end subroutine df_dump_latt_sigf

  subroutine df_dump_spin_susc()
     implicit none

     return
  end subroutine df_dump_spin_susc

  subroutine df_dump_char_susc()
     implicit none

     return
  end subroutine df_dump_char_susc

  subroutine df_dump_dens_vert()
     implicit none

     return
  end subroutine df_dump_dens_vert

  subroutine df_dump_magn_vert()
     implicit none

     return
  end subroutine df_dump_magn_vert

  subroutine df_dump_full_vert()
     implicit none

     return
  end subroutine df_dump_full_vert
