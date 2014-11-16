!!!=========+=========+=========+=========+=========+=========+=========+!
!!! HIBISCUS/toolbox/makeups @ iQIST                                     !
!!!                                                                      !
!!! This tool is used to calculate Gaussian broadening for density of    !
!!! states, which can be used to compare with PES and XAS experiments.   !
!!! author  : Li Huang (at IOP/CAS & SPCLab/CAEP & UNIFR)                !
!!! version : v2014.10.11T                                               !
!!! status  : WARNING: IN TESTING STAGE, USE IT IN YOUR RISK             !
!!! comment : any question, please contact with huangli712@gmail.com     !
!!!=========+=========+=========+=========+=========+=========+=========+!

!!
!!
!! Introduction
!! ============
!!
!! The makeups code is often used to postprocess the spectral function
!! data to compare with the XAS and UPS experiments. Now this code is
!! interfaced with entropy1 code merely. It can read the mem.dos.dat file
!! as input data. While to interface it with the stoch code is straight-
!! forward. What you need to do is to rename sai.imsum.dat to mem.dos.dat
!! file, and then supplement the lost orbital data.
!!
!! About the smearing parameter
!! The standard deviation smearing parameter was chosen to be in the same
!! range as estimates of experimental resolution (which are around 0.1
!! for high resolution PES, and approximately 0.2 to 0.4 for XAS. A good
!! test to decide if the broadening is correct is the comparison of the
!! Fermi edge in theory and experiment.
!!
!! About the beta parameter
!! The beta parameter practically plays no role if one uses the Fermi
!! function at the experimental temperature or at the temperature of the
!! QMC calculations
!!
!! Usage
!! =====
!!
!! # ./mups or bin/mups.x
!!
!! Input
!! =====
!!
!! mem.dos.dat (necessary)
!!
!! Output
!! ======
!!
!! ups.pes.dat
!! ups.xas.dat
!!
!! Documents
!! =========
!!
!! For more details, please go to iqist/doc/manual directory.
!!
!!
  program makeups
     use constants

     implicit none

!-------------------------------------------------------------------------
! local setting parameters
!-------------------------------------------------------------------------
! number of frequency grid on half plane (total number = 2*nw + 1)
     integer  :: nw   = 400

! number of orbitals, include spin degree of freedom
     integer  :: nq   = 2

! inversion temperature
     real(dp) :: beta = 10.0_dp

! broadening parameter
     real(dp) :: gamm = 0.30_dp
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

!-------------------------------------------------------------------------
! local variables
!-------------------------------------------------------------------------
! loop index for frequency grid
     integer :: iw

! loop index for orbitals
     integer :: iq

! status flag
     integer :: istat

! logical file exist flag
     logical :: fexist

! energy grid
     real(dp), allocatable :: mesh(:)

! original data of density of states
     real(dp), allocatable :: dos1(:,:)
     real(dp), allocatable :: dos3(:,:)

! broadened data of density of states (PES)
     real(dp), allocatable :: dos2(:,:)

! broadened data of density of states (XAS)
     real(dp), allocatable :: dos4(:,:)
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

! print program header
     write(mystd,'(2X,a)') 'MUPS'
     write(mystd,'(2X,a)') 'making ultraviolet photoemission spectra'
     write(mystd,'(2X,a)') 'version: 2011.08.18T'
     write(mystd,*) ! print blank line

! setup necessary parameters
     write(mystd,'(2X,a)')   '>>> number of orbitals (default = 2):'
     write(mystd,'(2X,a,$)') '>>> '
     read (mystd,'(i)') nq
     write(mystd,*)

     write(mystd,'(2X,a)')   '>>> number of frequency points (default = 400):'
     write(mystd,'(2X,a,$)') '>>> '
     read (mystd,'(i)') nw
     write(mystd,*)

     write(mystd,'(2X,a)')   '>>> inversion of temperature (default = 10.0):'
     write(mystd,'(2X,a,$)') '>>> '
     read (mystd,  *  ) beta
     write(mystd,*)

     write(mystd,'(2X,a)')   '>>> broadening parameter (default = 0.30):'
     write(mystd,'(2X,a,$)') '>>> '
     read (mystd,  *  ) gamm
     write(mystd,*)

! allocate memory
     allocate(mesh(-nw:nw),    stat=istat)

     allocate(dos1(-nw:nw,nq), stat=istat)
     allocate(dos2(-nw:nw,nq), stat=istat)
     allocate(dos3(-nw:nw,nq), stat=istat)
     allocate(dos4(-nw:nw,nq), stat=istat)

! initialize arrays
     mesh = zero

     dos1 = zero
     dos2 = zero
     dos3 = zero
     dos4 = zero

!-------------------------------------------------------------------------

! inquire data file
     inquire(file = 'mem.dos.dat', exist = fexist)
     if ( fexist == .false. ) then
         write(mystd,'(2X,a)') 'file mem.dos.dat does not exist'
         STOP ! terminate this program
     endif

! open density of states file
     write(mystd,'(2X,a)') '>>> reading mem.dos.dat ...'
     open(mytmp, file='mem.dos.dat', form='formatted', status='unknown')

! read in data
     do iq=1,nq
         do iw=-nw,nw
             read(mytmp,*) mesh(iw), dos1(iw,iq)
         enddo ! over iw={-nw,nw} loop
         read(mytmp,*) ! skip two lines
         read(mytmp,*)
     enddo ! over iq={1,nq} loop

! close original density of states file
     close(mytmp)
     write(mystd,'(2X,a)') '>>> status: OK'
     write(mystd,*)

!-------------------------------------------------------------------------

! copy data from dos1 to dos3
     dos3 = dos1

! multiply by the Fermi-Dirac function
     do iw=-nw,nw
         dos1(iw,:) = one / ( exp( +beta * mesh(iw) ) + one ) * dos1(iw,:)
         dos3(iw,:) = one / ( exp( -beta * mesh(iw) ) + one ) * dos3(iw,:)
     enddo ! over iw={-nw,nw} loop

! broadening density of states
     do iq=1,nq
         write(mystd,'(2X,a,i2)') 'PES >>> smearing band #', iq

         call broadening(nw, gamm, mesh, dos1(:,iq), dos2(:,iq))

         write(mystd,'(2X,a)') '>>> status: OK'
         write(mystd,*)
     enddo ! over iq={1,nq} loop

     do iq=1,nq
         write(mystd,'(2X,a,i2)') 'XAS >>> smearing band #', iq

         call broadening(nw, gamm, mesh, dos3(:,iq), dos4(:,iq))

         write(mystd,'(2X,a)') '>>> status: OK'
         write(mystd,*)
     enddo ! over iq={1,nq} loop

!-------------------------------------------------------------------------

! open density of states file
     write(mystd,'(2X,a)') '>>> writing ups.pes.dat ...'
     open(mytmp, file='ups.pes.dat', form='formatted', status='unknown')

! write out broadened density of states
     do iq=1,nq
         do iw=-nw,nw
             write(mytmp,'(2f16.8)') mesh(iw), dos2(iw,iq)
         enddo ! over iw={-nw,nw} loop
         write(mytmp,*) ! print two blank lines
         write(mytmp,*)
     enddo ! over iq={1,nq} loop

! close data file
     close(mytmp)
     write(mystd,'(2X,a)') '>>> status: OK'
     write(mystd,*)

!-------------------------------------------------------------------------

! open density of states file
     write(mystd,'(2X,a)') '>>> writing ups.xas.dat ...'
     open(mytmp, file='ups.xas.dat', form='formatted', status='unknown')

! write out broadened density of states
     do iq=1,nq
         do iw=-nw,nw
             write(mytmp,'(2f16.8)') mesh(iw), dos4(iw,iq)
         enddo ! over iw={-nw,nw} loop
         write(mytmp,*) ! print two blank lines
         write(mytmp,*)
     enddo ! over iq={1,nq} loop

! close data file
     close(mytmp)
     write(mystd,'(2X,a)') '>>> status: OK'
     write(mystd,*)

!-------------------------------------------------------------------------

! deallocate memory
     deallocate(mesh)

     deallocate(dos1)
     deallocate(dos2)
     deallocate(dos3)
     deallocate(dos4)

  end program makeups

!>>> calculate broadening with the Gaussian function (S = smearing):
! we define:
!     fct(x) = \frac{1}{S\sqrt{2\pi}} e^{-\frac{x^{2}}{2 S^{2}}}
! and
!     Precision = S / dx. (dx = Delta Energy between 2 points)
! then
!     Broad. Spectrum(x) =
!         \sum_{x'=-3*Precision}^{3*Precision} fct(x')*Spectrum(x+x')
! where x' = DeltaJ in the procedure.
  subroutine broadening(nw, gamm, mesh, xinp, xout)
     use constants

     implicit none

! external arguments
! number of frequency mesh (2*nw + 1)
     integer, intent(in)   :: nw

! smearing parameter
     real(dp), intent(in)  :: gamm

! frequency grid
     real(dp), intent(in)  :: mesh(-nw:nw)

! input density of states
     real(dp), intent(in)  :: xinp(-nw:nw)

! broadened density of states
     real(dp), intent(out) :: xout(-nw:nw)

! local variables
! loop index
     integer  :: iw

! define the smearing range
     integer  :: prec
     integer  :: delj

! fac1 and fac2 are used to build Gaussian function
     real(dp) :: fac1
     real(dp) :: fac2

! initialize xout
     do iw=-nw,nw
         xout(iw) = zero
     enddo ! over iw={-nw,nw} loop

! calculate broadening range
     prec = nint( abs( gamm / ( mesh(2) - mesh(1) ) ) )

! calculate broadening density of states
     if ( prec > 0 ) then

! calculate fac1 and fac2
         fac1 = one / ( real(prec) * sqrt( two * pi ) )
         fac2 = real( two * prec * prec )

         do iw=-nw,nw
             do delj=max(-3*prec,-nw-iw),min(3*prec,nw-iw)
                 xout(iw) = xout(iw) + xinp(iw+delj) * exp( -real( delj * delj ) / fac2 )
             enddo
             xout(iw) = xout(iw) * fac1
         enddo ! over iw={-nw,nw} loop

! set xinp == xout, do nothing
     else

         do iw=-nw,nw
             xout(iw) = xinp(iw)
         enddo ! over iw={-nw,nw} loop

     endif ! back if ( prec > 0 ) block

     return
  end subroutine broadening
