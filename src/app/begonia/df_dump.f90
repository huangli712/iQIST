
  subroutine df_dump_dmft_grnf()
     implicit none

     use constants, only : dp, mytmp

     use control, only : norbs
     use control, only : mfreq

     implicit none

! external arguments
! matsubara frequency mesh
     real(dp), intent(in)    :: rmesh(mfreq)

! impurity green's function
     complex(dp), intent(in) :: grnf(mfreq,norbs,norbs)
     complex(dp), intent(in) :: gerr(mfreq,norbs,norbs)

! local variables
! loop index
     integer :: i
     integer :: j

! open data file: solver.grn.dat
     open(mytmp, file='solver.grn.dat', form='formatted', status='unknown')

! write it
     do i=1,norbs
         do j=1,mfreq
             write(mytmp,'(i6,5f16.8)') i, rmesh(j), &
              real(grnf(j,i,i)), aimag(grnf(j,i,i)), &
              real(gerr(j,i,i)), aimag(gerr(j,i,i))
         enddo ! over j={1,mfreq} loop
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
