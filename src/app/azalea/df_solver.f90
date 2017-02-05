
  subroutine df_bse_solver(bubbleM, vertexM, gammaM)
     use constants

     use df_control
     use df_context

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
  end subroutine df_bse_solver

  subroutine df_bse_solver_iter(niter, mix, bubbleM, vertexM, gammaM)
     use constants

     use df_control
     use df_context

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
  end subroutine df_bse_solver_iter
