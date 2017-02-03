
  subroutine df_run()
     use constants

     use df_control
     use df_context

     implicit none

! local variables
     integer :: i
     integer :: j
     integer :: k

     real(dp) :: w
     complex(dp), allocatable :: bubble(:,:,:)
     allocate(bubble(nkpts,nffrq,norbs))

     DF_LOOP: do i=1,ndfit
         write(mystd,'(2X,A,I3)') 'Ladder Dual Fermion Iteration:', i

         Q_LOOP: do j=2,nbfrq-1
             w = bmesh(j)
             write(mystd,'(2X,A,F12.6)') 'Bosonic Frequency:', w
             if ( abs(w - zero) < epss ) then
                 call df_static_bubble(bubble, w)
             else
                 call df_bubble(bubble, w)
             endif

             K_LOOP: do k=1,nkpts
             enddo K_LOOP

         enddo Q_LOOP

         write(mystd,*)
     enddo DF_LOOP

     deallocate(bubble)

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
