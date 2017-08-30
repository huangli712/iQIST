
  subroutine dt_dump_grnf(rmesh, grnf)
     use constants, only : dp
     use constants, only : czero
     use constants, only : mytmp

     use control, only : norbs
     use control, only : nffrq

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

! open data file: dt.dmft_g.dat
     open(mytmp, file='dt.dmft_g.dat', form='formatted', status='unknown')

! write it
     do i=1,norbs
         do j=1,nffrq
             write(mytmp,'(i6,5f16.8)') i, rmesh(j), grnf(j,i), czero
         enddo ! over j={1,nffrq} loop
         write(mytmp,*) ! write empty lines
         write(mytmp,*)
     enddo ! over i={1,norbs} loop

! close data file
     close(mytmp)

     return
  end subroutine dt_dump_grnf

  subroutine dt_dump_dmft_sigf(rmesh, sigf)
     use constants, only : dp, zero, mytmp

     use control, only : norbs
     use control, only : nffrq

     implicit none

! external arguments
! matsubara frequency mesh
     real(dp), intent(in)    :: rmesh(nffrq)

! impurity self-energy function
     complex(dp), intent(in) :: sigf(nffrq,norbs)

! local variables
! loop index
     integer :: i
     integer :: j

! open data file: df.dmft_s.dat
     open(mytmp, file='df.dmft_s.dat', form='formatted', status='unknown')

! write it
     do i=1,norbs
         do j=1,nffrq
             write(mytmp,'(i6,5f16.8)') &
                 i, rmesh(j), real(sigf(j,i)), aimag(sigf(j,i)), zero, zero
         enddo ! over j={1,nffrq} loop
         write(mytmp,*) ! write empty lines
         write(mytmp,*)
     enddo ! over i={1,norbs} loop

! close data file
     close(mytmp)

     return
  end subroutine dt_dump_dmft_sigf

  subroutine dt_dump_dmft_hybf(rmesh, hybf)
     use constants, only : dp, zero, mytmp

     use control, only : norbs
     use control, only : nffrq

     implicit none

! external arguments
! matsubara frequency mesh
     real(dp), intent(in)    :: rmesh(nffrq)

! impurity hybridization function
     complex(dp), intent(in) :: hybf(nffrq,norbs)

! local variables
! loop index
     integer :: i
     integer :: j

! open data file: df.dmft_h.dat
     open(mytmp, file='df.dmft_h.dat', form='formatted', status='unknown')

! write it
     do i=1,norbs
         do j=1,nffrq
             write(mytmp,'(i6,5f16.8)') &
                 i, rmesh(j), real(hybf(j,i)), aimag(hybf(j,i)), zero, zero
         enddo ! over j={1,nffrq} loop
         write(mytmp,*) ! write empty lines
         write(mytmp,*)
     enddo ! over i={1,norbs} loop

! close data file
     close(mytmp)

     return
  end subroutine dt_dump_dmft_hybf

  subroutine dt_dump_dual_grnf()
     implicit none

     return
  end subroutine dt_dump_dual_grnf

  subroutine dt_dump_dual_sigf()
     implicit none

     return
  end subroutine dt_dump_dual_sigf

  subroutine dt_dump_dual_wssf()
     implicit none

     return
  end subroutine dt_dump_dual_wssf

  subroutine dt_dump_latt_grnf()
     implicit none

     return
  end subroutine dt_dump_latt_grnf

  subroutine dt_dump_latt_sigf()
     implicit none

     return
  end subroutine dt_dump_latt_sigf

  subroutine dt_dump_spin_susc()
     implicit none

     return
  end subroutine dt_dump_spin_susc

  subroutine dt_dump_char_susc()
     implicit none

     return
  end subroutine dt_dump_char_susc

  subroutine dt_dump_dens_vert()
     implicit none

     return
  end subroutine dt_dump_dens_vert

  subroutine dt_dump_magn_vert()
     implicit none

     return
  end subroutine dt_dump_magn_vert

  subroutine dt_dump_full_vert()
     implicit none

     return
  end subroutine dt_dump_full_vert
