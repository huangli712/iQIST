
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
     complex(dp), allocatable :: bubbleM(:,:)
     complex(dp), allocatable :: vertexM(:,:)
     complex(dp), allocatable :: vertexD(:,:)

     allocate(bubble(nkpts,nffrq,norbs))
     allocate(bubbleM(nffrq,nffrq))
     allocate(vertexM(nffrq,nffrq))
     allocate(vertexD(nffrq,nffrq))

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

             vertexM = vert_m(:,:,j)
             vertexD = vert_d(:,:,j)
             K_LOOP: do k=1,nkpts
                 call s_diag_z(nffrq, bubble(k,:,1), bubbleM)
             enddo K_LOOP

         enddo Q_LOOP

         write(mystd,*)
     enddo DF_LOOP

     deallocate(bubble)
     deallocate(bubbleM)
     deallocate(vertexM)
     deallocate(vertexD)

     return
  end subroutine df_run
