!---------------------------------------------------------------
! project : maxent
! program : maxent_dump_hist
! source  : maxent_dump.f90
! type    : subroutine
! author  : yilin wang (email: qhwyl2006@126.com)
! history : 06/11/2013 by yilin wang
! purpose : 
! input   :
! output  :
! status  : unstable
! comment :
!---------------------------------------------------------------

!===============================================================
! subroutine maxent_dump_hist is used to dump the histogram of
! databins for reference
!===============================================================
  subroutine maxent_dump_hist( hmesh, hist )
      use constants
      use control

      implicit none
! external variables
      real(dp), intent(in) :: hmesh(slice, ntime)
      integer, intent(in) :: hist(slice, ntime)

! local variables
! loop index
      integer :: i
      integer :: j

! open file maxent.hist.dat
      open(mytmp, file="maxent.hist.dat", form="formatted", status="unknown")
      do i=1, ntime
          do j=1, slice
              write(mytmp,"(i5,f16.8,i5)") i, hmesh(j,i), hist(j,i) 
          enddo
          write(mytmp,*)
          write(mytmp,*)
      enddo
      close(mytmp)

      return
  end subroutine maxent_dump_hist

!===============================================================
! subroutine maxent_dump_model is used to dump the default model
! for reference
!===============================================================
  subroutine maxent_dump_model()
      use constants
      use control
      use context, only: fmesh, aw_model
    
      implicit none

! local variables
! loop index
      integer :: i

! open file maxent.model.dat
      open(mytmp, file="maxent.model.dat", form="formatted", status="unknown")
      do i=1, nw
          write(mytmp,"(2f16.8)") fmesh(i), aw_model(i) 
      enddo
      close(mytmp)

      return
  end subroutine maxent_dump_model

!===============================================================
! subroutine maxent_dump_eigcov is used to dump the eigvalues
! of covariance matrix
!===============================================================
  subroutine maxent_dump_eigcov()
      use constants
      use control
      use context, only: eigcov

      implicit none

! local variables
! loop index
      integer :: i

! open file maxent.eigcov.dat
      open(mytmp, file="maxent.eigcov.dat", form="formatted", status="unknown")
      do i=1, ntime
          write(mytmp,"(i5,E16.8)") i, eigcov(i)
      enddo
      close(mytmp)

      return
  end subroutine maxent_dump_eigcov

!===============================================================
! subroutine maxent_dump_palpha is used to dump the posterior
! probability of alpha for reference
!===============================================================
  subroutine maxent_dump_palpha( amesh, palpha )
      use constants
      use control
      
      implicit none

! external variables
      real(dp), intent(in) :: amesh(nalpha)
      real(dp), intent(in) :: palpha(nalpha)

! local variables
! loop index
      integer :: i

! open file maxent.palpha.dat 
      open(mytmp, file="maxent.palpha.dat", form="formatted", status="unknown")
      do i=1, nalpha
          write(mytmp,"(f16.6,5X,E20.8E8)") amesh(i), palpha(i)
      enddo
      close(mytmp)
 
      return
  end subroutine maxent_dump_palpha

!===============================================================
! subroutine maxent_dump_aw is used to dump the final result of
! aw
!===============================================================
  subroutine maxent_dump_aw( fmesh, aw )
      use constants
      use control
    
      implicit none

! external variables
      real(dp), intent(in) :: fmesh(nw)
      real(dp), intent(in) :: aw(nw)

! local variables
! loop index
      integer :: i

! open file maxent.aw.dat
      open(mytmp, file="maxent.aw.dat", form="formatted", status="unknown")
      do i=1, nw
          write(mytmp,"(2f16.8)") fmesh(i), aw(i)
      enddo
      close(mytmp)

      return
  end subroutine maxent_dump_aw

!===============================================================
! subroutine: maxent_dump_svd 
! purpose   : dump the svd sigma vector of kern^{T}
!===============================================================
  subroutine maxent_dump_svd( n, vec )
      use constants

      implicit none

! external variables
      integer, intent(in) :: n
      real(dp), intent(in) :: vec(n)

! local variables
      integer :: i

      open(mytmp, file="maxent.svd.dat", form="formatted", status="unknown")
      do i=1, n
          write(mytmp, "(i5, E16.8)") i, vec(i)
      enddo
      close(mytmp)

      return
  end subroutine maxent_dump_svd
