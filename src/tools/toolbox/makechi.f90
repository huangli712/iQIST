!=========+=========+=========+=========+=========+=========+=========+>>>
! postprocess the spin-spin correlation function (solver.schi.dat). it   !
! is calculated by continuous time quantum Monte Carlo impurity solver.  !
! author  : li huang                                                     !
! version : v2011.08.18T                                                 !
! status  : WARNING: IN TESTING STAGE, USE IT IN YOUR RISK               !
! comment : any question, please contact with huangli712@yahoo.com.cn    !
!=========+=========+=========+=========+=========+=========+=========+>>>

  program makechi
     use constants

     implicit none

!-------------------------------------------------------------------------
! local setting parameters
!-------------------------------------------------------------------------
! number of bands
     integer  :: nband = 1

! number of time slice in [0, \beta]
     integer  :: ntime = 1024

! inversion of temperature
     real(dp) :: beta  = 10.00_dp
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

!-------------------------------------------------------------------------
! local variables
!-------------------------------------------------------------------------
! loop index
     integer  :: i
     integer  :: j

! status flag
     integer  :: istat

! magnetic susceptibility \chi_{loc}
     real(dp) :: xchi

! totally spin-spin correlation function: < Sz(0) Sz(\tau) >
     real(dp), allocatable :: schi(:)
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

! print program header
     write(mystd,'(2X,a)') 'MCHI'
     write(mystd,'(2X,a)') 'making magnetic susceptibility from spin-spin correlation function'
     write(mystd,'(2X,a)') 'version: 2011.08.18T'
     write(mystd,*) ! print blank line

! setup necessary parameters
     write(mystd,'(2X,a)')   '>>> number of bands (default = 1):'
     write(mystd,'(2X,a,$)') '>>> '
     read (mystd,'(i)') nband
     write(mystd,*)

     write(mystd,'(2X,a)')   '>>> number of time slice (default = 1024):'
     write(mystd,'(2X,a,$)') '>>> '
     read (mystd,'(i)') ntime
     write(mystd,*)

     write(mystd,'(2X,a)')   '>>> inversion of temperature (default = 10.0):'
     write(mystd,'(2X,a,$)') '>>> '
     read (mystd,  *  ) beta
     write(mystd,*)

! allocate memory
     allocate(schi(ntime), stat=istat)

! read in spin-spin correlation function
     schi = zero
     open(mytmp, file='solver.schi.dat', form='formatted', status='unknown')

! skip data block for orbital-resolved spin-spin correlation function
     do j=1,nband
         read(mytmp,*)
         do i=1,ntime
             read(mytmp,*)
         enddo ! over i={1,ntime} loop
         read(mytmp,*)
         read(mytmp,*)
     enddo ! over j={1,nband} loop

     read(mytmp,*)
     do i=1,ntime
        read(mytmp,*) xchi, schi(i)
     enddo ! over i={1,ntime} loop
     close(mytmp)

! rescale spin-spin correlation function
     schi = schi * real(nband)

! calculate local magnetic susceptibility
     xchi = zero
     do i=1,ntime-1
         xchi = xchi + ( schi(i) + schi(i+1) ) / two * ( beta / real(ntime - 1) )
     enddo ! over i={1,ntime-1} loop

! write out the calculated results
     write(mystd,'(2X,a,f12.6)') 'magnetic susceptibility  :', xchi
     write(mystd,'(2X,a,f12.6)') 'effective magnetic moment:', sqrt(xchi / beta)

! deallocate memory
     deallocate(schi)

  end program makechi
