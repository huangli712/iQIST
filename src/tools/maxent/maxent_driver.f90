!---------------------------------------------------------------
! project : maxent
! program : maxent_driver
! source  : maxent_driver.f90
! type    : subroutine
! author  : yilin wang (email: qhwyl2006@126.com)
! history : 06/11/2013 by yilin wang
! purpose : for all the values of $\alpha$, find the optimal 
!         : spectrum $A(\omega)$. 
! input   :
! output  :
! status  : unstable
! comment :
!---------------------------------------------------------------

!=================================================================
! subroutine maxent_driver is used to find the optimal $A(\omega)$
! for all of the values of $\alpha$.
!================================================================= 
  subroutine maxent_driver()
      use constants
      use control
      use context

      implicit none
     
! local variables
! loop index
      integer :: ialpha
      integer :: iw

! the base of exp(Q)
      real(dp) :: qbase

! the distribution probability of alpha
      real(dp) :: palpha(nalpha)

! norm factor
      real(dp) :: norm 

! the final result of aw
      real(dp) :: aw(nw)

! begin the main loop of alpha
      write(mystd,"(2X,a)") "MAXENT >>> Begin The Loop For Each Alpha"
      write(mystd,"(2X,a)") "============================================================="

      do ialpha=1, nalpha
          write(mystd,"(4X,a,i5,4X,a,f12.6)") "No. ", ialpha, "Alpha: ", amesh(ialpha)
          call maxent_makeaw(ialpha)
      enddo

      write(mystd,"(2X,a)") "============================================================="
      write(mystd,"(2X,a)") "MAXENT >>> End The Loop For Each Alpha"
      write(mystd,*)

      write(mystd,"(2X,a)") "MAXENT >>> Calculate The Final Spectrum Function..."

! calculate the average spectrum over alpha
! step 1: extract the base of exp(Q), we choose the max value of 
! Q values 
      qbase =  maxval(qval) 
      do ialpha=1, nalpha
          qval(ialpha) = qval(ialpha) - qbase
      enddo
! step 2: calculate the distribution probability of alpha
      do ialpha=1, nalpha
          palpha(ialpha) =  weight(ialpha) * exp(qval(ialpha))
      enddo

! step 3: dump palpha for reference
      call maxent_dump_palpha( amesh, palpha )

! step 4: calculate the final spectrum function $A(\omega)$
      norm = zero
      do ialpha=1, nalpha-1
          norm = norm + half * ( palpha(ialpha) + palpha(ialpha+1) ) * ( amesh(1) - amesh(2) ) 
      enddo

      aw = zero
      do iw=1, nw
          do ialpha=1, nalpha-1
              aw(iw) = aw(iw) + aw_alpha(iw,ialpha) * half * ( palpha(ialpha)+palpha(ialpha+1) ) &
                                                    * ( amesh(1) - amesh(2) ) 
          enddo
          aw(iw) = aw(iw) / norm
      enddo

! rescale aw for dump
      aw = aw / step

! step 5: dump the final result of aw
      call maxent_dump_aw( fmesh, aw )

      return
  end subroutine maxent_driver
