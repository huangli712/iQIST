!---------------------------------------------------------------
! project : maxent
! program : maxent_makeaw
!         : maxent_makemu
! source  : maxent_makeaw.f90
! type    : subroutine
! author  : yilin wang (email: qhwyl2006@126.com)
! history : 06/10/2013 by yilin wang
! purpose : for fixed $\alpha$, find the optimal spectrum $A(\omega)$. 
! input   :
! output  :
! status  : unstable
! comment :
!---------------------------------------------------------------

!=================================================================
! subroutine: maxent_makeaw
! purpose   : find the optimal spectrum $A(\omega)$
!           : for fixed parmaeter $\alpha$                                
!=================================================================
  subroutine maxent_makeaw(ialpha)
      use constants
      use control
      use context

      implicit none

! external variables
      integer, intent(in) :: ialpha

! local variables
! loop index
      integer :: itime
      integer :: iw
      integer :: iss

! the error information from the dsyev call
      integer :: info

! counter for loop times
      integer :: counter

! current alpha
      real(dp) :: alpha

! parameter $\mu$
      real(dp) :: mune

! the u vector
      real(dp) :: uvec(ns)

! the optimal spectrum function
      real(dp) :: aw(nw)

! rkern * aw
      real(dp) :: keraw(ntime)

! gardient of L
      real(dp) :: gardl(ntime)

! g vector
      real(dp) :: gvec(ns)

! $-\alpha * u - g$ vector augvec
      real(dp) :: augvec(ns)

! matrix $K = U^{T} * diag{aw} * U$
      real(dp) :: kmat(ns,ns)

! eigvalues of K
      real(dp) :: keigval(ns)

! eigvectors of K
      real(dp) :: pmat(ns,ns)

! matrix $A = diag{\sqrt{keigval}} * P^{T} * M * P * diag{\sqrt{keigval}}$
      real(dp) :: amat(ns,ns)

! eigvalues of A
      real(dp) :: aeigval(ns)

! eigvectors of A
      real(dp) :: rmat(ns,ns)

! $Y^{-1}$
      real(dp) :: yinv(ns,ns)

! $Y^{-T}$
      real(dp) :: yinvt(ns,ns)

! $M * Y^{-T}$
      real(dp) :: myinvt(ns,ns)

! $-\alpha * Y^{-1} * u - Y^{-1} * g$
      real(dp) :: yiaugvec(ns)

! $L=-\chi^2 / 2$
      real(dp) :: sum_l

! entropy S
      real(dp) :: sum_s

! -\alpha * s
      real(dp) :: alpha_s

! Q = -\alpha * S - L
      real(dp) :: sum_q

! norm of L
      real(dp) :: norml

! norm of s
      real(dp) :: norms

! norm of -\alpha *s - l
      real(dp) :: normsl

! difference
      real(dp) :: diff

! posterior probability for present alpha
      real(dp) :: post

! temp vectors
      real(dp) :: vtemp1(ns)
      real(dp) :: vtemp2(ns)
      real(dp) :: vtemp3(nw)

! temp matrixs
      real(dp) :: mtemp1(nw, ns)
      real(dp) :: mtemp2(ns, nw)
      real(dp) :: mtemp3(ns, ns)
      real(dp) :: mtemp4(ns, ns)
      real(dp) :: mtemp5(ns, ns)

! begin to search the optimal spectrum function $A(\omega)$ by using
! the quasi-Newton method

! step 1: initialize aw and the uvec 
      aw = aw_model
      uvec = zero
      alpha = amesh(ialpha) 
! the main loop
      counter = 0   
      SEARCH_LOOP: do while(.true.)
      
          counter = counter + 1

! step 2: construct the gardient of L
          call maxent_dgemv(ntime, nw, rkern, aw, keraw)

          do itime=1, ntime
              gardl(itime) = ( keraw(itime) - rgrn(itime) ) / eigcov(itime)
          enddo           

! step 3: construct g vector
          call maxent_dgemv(ns, ntime, vmatt, gardl, vtemp1)
        
          do iss=1, ns
              gvec(iss) = sigvec(iss) * vtemp1(iss)
          enddo

! step 4: construct $-\alpha * u - g$
          augvec = -alpha * uvec - gvec

! step 5: construct $K=U^{T} * diag{aw} * U$ and diagonalize it
! step 5.1: construct K
          do iw=1, nw
              mtemp1(iw,:) =  aw(iw) * umat(iw,:)
          enddo 

          mtemp2 = transpose(umat)

          call maxent_dgemm(ns, nw, mtemp2, ns, mtemp1, kmat)

! step 5.2: diagonalize K
          mtemp3 = kmat
          call maxent_dsyev(ns, mtemp3, keigval, info)
! check error from dsyev call
          if( info /= 0 ) then
              call maxent_print_advice("Kmat")
          endif
          pmat = mtemp3

! step 6: construct $A=diag{\sqrt{keigval}} * P^{T} * M * P * diag{\sqrt{keigval}} $
! and diagonalize it
! step 6.1: construct A
          do iss=1, ns
              mtemp3(:,iss) = pmat(:,iss) * sqrt( keigval(iss) )
          enddo 
          mtemp4 = transpose(mtemp3)

          call maxent_dgemm(ns, ns, mtemp4, ns, mmat, mtemp5)
          call maxent_dgemm(ns, ns, mtemp5, ns, mtemp3, mtemp4)
        

! step 6.2: diagonalize A
          amat = mtemp4
          call maxent_dsyev(ns, mtemp4, aeigval, info)
! check error from dsyev call
          if( info /= 0 ) then
              call maxent_print_advice("Amat")
          endif
          rmat = mtemp4 

! step 7: construct $Y^{-1}$ and $Y^{-T}$
! step 7.1: $Y^{-1}$
          mtemp3 = transpose(rmat)
          do iss=1, ns
              mtemp4(iss,:) = sqrt(keigval(iss)) * pmat(:,iss)
          enddo        
          call maxent_dgemm(ns, ns, mtemp3, ns, mtemp4, mtemp5)
          yinv = mtemp5

! step 7.2: $Y^{-T}$
          mtemp3 = pmat
          do iss=1, ns
              mtemp4(iss,:) = sqrt(keigval(iss)) * rmat(iss,:)
          enddo
          call maxent_dgemm(ns, ns, mtemp3, ns, mtemp4, mtemp5)
          yinvt = mtemp5

! step 8: consturct $M * Y^{-T}$
          call maxent_dgemm(ns, ns, mmat, ns, yinvt, myinvt)

! step 9: consturct $-\alpha * Y^{-1} * u - Y^{-1} * g$ 
          call maxent_dgemv(ns, ns, yinv, augvec, yiaugvec)

! step 10: adjust the parameter $\mu$
          call maxent_makemu(yiaugvec, aeigval, alpha, mune)

! step 11: construct new spectrum 
! step 11.1: construct $Y^{-1} * du$
          do iss=1, ns
              vtemp1(iss) = ( yiaugvec(iss) ) / ( alpha + mune + aeigval(iss) )            
          enddo
! step 11.2: construct $M * Y^{-T} * ( Y^{-1} * du)$
          call maxent_dgemv(ns, ns, myinvt, vtemp1, vtemp2) 

! step 11.3: construct new uvec
          do iss=1, ns
              uvec(iss) = uvec(iss) + ( augvec(iss) - vtemp2(iss) ) / ( alpha + mune )
          enddo

! step 11.4: construct new aw
          call maxent_dgemv(nw, ns, umat, uvec, vtemp3)

          do iw=1, nw
              aw(iw) = aw_model(iw) * exp( vtemp3(iw) )
          enddo

! step 12: calculate S and L
! step 12.1: L
          sum_l = zero
          call maxent_dgemv(ntime, nw, rkern, aw, keraw)
          do itime=1, ntime
              sum_l = sum_l + half * ( keraw(itime) - rgrn(itime) ) **2 / eigcov(itime)
          enddo

! step 12.2: S
          sum_s = zero
          do iw=1, nw
              sum_s = sum_s - aw(iw) * vtemp3(iw) + aw(iw) - aw_model(iw) 
          enddo

! step 12.3: \alpha * s and Q
          alpha_s = alpha * sum_s
          sum_q = alpha_s - sum_l

! step 13: check the convergence
          call maxent_dgemv(ns, ns, yinv, uvec, vtemp1)  
          call maxent_dgemv(ns, ns, yinv, gvec, vtemp2)  

          norml = zero
          norms = zero
          do iss=1, ns
             norml = norml + vtemp1(iss) **2 
             norms = norms + vtemp2(iss) **2 
          enddo
          norml = sqrt(norml)
          norms = sqrt(norms)

          call maxent_dgemv(ns, ns, yinv, augvec, vtemp1)
          normsl = zero
          do iss=1, ns
              normsl = normsl + vtemp1(iss) ** 2
          enddo
          normsl = sqrt(normsl)

          diff = two * normsl **2 / (alpha * norms + norml) **2

          if ( diff <= eps8 ) then
               write(mystd,"(4X,a,f12.6)") "Convergence For Alpha: ",amesh(ialpha)
               EXIT SEARCH_LOOP
          endif
      enddo SEARCH_LOOP

! step 14: calculate the posterior probability for the present alpha
      post = one
      do iss=1, ns
          post = post * alpha / (alpha + aeigval(iss))
      enddo
      post = sqrt(post) / alpha

! put local variables into the global variables
      weight(ialpha) = post 
      ent(ialpha) = alpha_s
      chi2(ialpha) = sum_l
      qval(ialpha) = sum_q
  
      aw_alpha(:,ialpha) = aw

! print runtime information for reference
    
      write(mystd,"(4X,a,i5)") "Loops to Find The Optimal Aw: ", counter
      write(mystd,"(4X,a,f16.6)") "alpha * S:                 ", alpha_s
      write(mystd,"(4X,a,f16.6)") "L=chi^2 / 2:               ", sum_l
      write(mystd,"(4X,a,f16.6)") "Q=alpha * S - L:           ", sum_q
      write(mystd,"(4X,a)") "---------------------------------------------------------"
      write(mystd,*)
      write(mystd,*)

      return        
  end subroutine maxent_makeaw

!===============================================================
! subroutine: maxent_makemu 
! purpose   : adjust the parameter mune to reach the condition 
!           : $|Y^{-1}\delta u|^{2}$ <= \sum m_{i}
!===============================================================
  subroutine maxent_makemu(yiaugvec, aeigval, alpha, mune)
     use constants
     use control
     use context

     implicit none 

! external variables
     real(dp), intent(in) :: yiaugvec(ns)
     real(dp), intent(in) :: aeigval(ns)
     real(dp), intent(in) :: alpha
     real(dp), intent(out) :: mune

! local variables
     real(dp) :: sum_m 
     real(dp) :: sum_yidu 

! current mune
     real(dp) :: cur_mu

! loop inex
     integer :: i


! initialize it to one
     cur_mu = half 

! calculate sum of m
     sum_m = zero
     do i=1,nw
         sum_m = sum_m + aw_model(i) 
     enddo     

     MUNE_LOOP: do while(.true.)
! step 1: calculate sum_yidu

         sum_yidu = zero
         do i=1, ns
             sum_yidu = sum_yidu + ( yiaugvec(i) / ( alpha + cur_mu + aeigval(i) ) ) ** 2
         enddo 
!>>>         print *, sum_yidu, sum_m, cur_mu
         if ( sum_yidu <= sum_m ) then
             mune = cur_mu
             EXIT  MUNE_LOOP
         endif
       
         if ( sum_yidu > sum_m ) then
! increase mune
             !cur_mu = cur_mu + 10.0_dp    
             cur_mu = cur_mu * 2.0_dp    
         endif

     enddo MUNE_LOOP

     return
  end subroutine maxent_makemu
