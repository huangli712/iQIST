!!!=========+=========+=========+=========+=========+=========+=========+!
!!! HIBISCUS/toolbox/maketau @ iQIST                                     !
!!!                                                                      !
!!! This tool is used to postprocess the bins data of the imaginary time !
!!! green's function                                                     !
!!! author  : Li Huang (at IOP/CAS & SPCLab/CAEP & UNIFR)                !
!!! version : v2015.01.06T                                               !
!!! status  : WARNING: IN TESTING STAGE, USE IT IN YOUR RISK             !
!!! comment : any question, please contact with lihuang.dmft@gmail.com   !
!!!=========+=========+=========+=========+=========+=========+=========+!

!!
!!
!! Introduction
!! ============
!!
!! The maketau code is often used to convert file solver.green.dat.*
!! or solver.green.dat to tau.grn.dat, prepare necessary input data for
!! the hibiscus/entropy or hibiscus/stoch codes.
!!
!! About solver.green.dat.* files:
!! In order to obtain solver.green.dat.* files, you have to run the
!! ctqmc/hfqmc codes for several times, and rename the solver.green.dat
!! file to solver.green.dat.* file manually.
!!
!! About nskip control parameter:
!! If ctqmc == 2 or 4, then nskip must be 1. If ctqmc == 1 or 3, then
!! nskip could be positive integer. Be careful, nskip can not be any
!! integer. Notice that mod(ntime - 1, nskip) must be 0, or else the
!! obtained tau.grn.dat should be wrong.
!!
!! Usage
!! =====
!!
!! # ./mtau or bin/mtau.x
!!
!! Input
!! =====
!!
!! solver.green.dat or solver.green.dat.* (necessary)
!!
!! Output
!! ======
!!
!! tau.grn.dat
!!
!! Documents
!! =========
!!
!! For more details, please go to iqist/doc/manual directory.
!!
!!

  program maketau
     use constants, only : dp, zero, mystd, mytmp

     implicit none

! local control parameters
! number of orbitals
     integer  :: norbs = 2

! number of time slices, 129 or 1024, in [0, \beta]
     integer  :: ntime = 129

! number of data bins
     integer  :: nbins = 1

! file type generated by quantum impurity solver
! if ctqmc == 1, ctqmc in std mode;
! if ctqmc == 2, hfqmc in std mode;
! if ctqmc == 3, ctqmc in bin mode;
! if ctqmc == 4, hfqmc in bin mode.
     integer  :: ctqmc = 1

! number of skipped points between two successive selected points
     integer  :: nskip = 1

! inversion of temperature
     real(dp) :: beta  = 10.0_dp

! local variables
! loop index
     integer  :: i
     integer  :: j

! loop index over QMC data bins
     integer  :: ibin

! dummy integer variables
     integer  :: itmp, jtmp

! status flag
     integer  :: istat

! used to check whether the input file exists
     logical  :: exists

! real(dp) dummy variables
     real(dp) :: rtmp

! current bin index, string representation
     character(len=10) :: sbin

! time slice mesh
     real(dp), allocatable :: tau(:)

! green's function data
     real(dp), allocatable :: grn(:,:), grn_bin(:,:,:)

! error bar data
     real(dp), allocatable :: err(:,:)

! print program header
     write(mystd,'(2X,a)') 'HIBISCUS/toolbox/maketau'
     write(mystd,'(2X,a)') '>>> Making tau-dependent imaginary time green''s function'
     write(mystd,*) ! print blank line

     write(mystd,'(2X,a)') 'Version: 2015.01.06T '//'(built at '//__TIME__//" "//__DATE__//')'
     write(mystd,'(2X,a)') 'Develop: by li huang (at IOP/CAS & SPCLab/CAEP & UNIFR)'
     write(mystd,'(2X,a)') 'Support: lihuang.dmft@gmail.com'
     write(mystd,'(2X,a)') 'License: GNU General Public License version 3'
     write(mystd,*) ! print blank line

! setup necessary parameters
     write(mystd,'(2X,a)')   'Number of orbitals (default = 2):'
     write(mystd,'(2X,a,$)') '>>> '
     read (*,*) norbs
     write(mystd,*)

     write(mystd,'(2X,a)')   'Number of time slices (default = 129 or 1024):'
     write(mystd,'(2X,a,$)') '>>> '
     read (*,*) ntime
     write(mystd,*)

     write(mystd,'(2X,a)')   'Number of data bins (default = 1):'
     write(mystd,'(2X,a,$)') '>>> '
     read (*,*) nbins
     write(mystd,*)

     write(mystd,'(2X,a)')   'File type generated by quantum impurity solver (default = 1):'
     write(mystd,'(2X,a)')   '  ctqmc: 1 (std mode)'
     write(mystd,'(2X,a)')   '  hfqmc: 2 (std mode)'
     write(mystd,'(2X,a)')   '  ctqmc: 3 (bin mode)'
     write(mystd,'(2X,a)')   '  hfqmc: 4 (bin mode)'
     write(mystd,'(2X,a,$)') '>>> '
     read (*,*) ctqmc
     write(mystd,*)

     write(mystd,'(2X,a)')   'Number of skipped points between two successive selected points (default = 1):'
     write(mystd,'(2X,a,$)') '>>> '
     read (*,*) nskip
     write(mystd,*)

     write(mystd,'(2X,a)')   'Inversion of temperature (default = 10.0):'
     write(mystd,'(2X,a,$)') '>>> '
     read (*,*) beta
     write(mystd,*)

! check the parameters
     call s_assert2( norbs > 0 .and. norbs < 15, 'wrong number of orbitals' )
     call s_assert2( ntime > 0, 'wrong number of time slices' )
     call s_assert2( nbins > 0, 'wrong number of data bins' )
     call s_assert2( ctqmc > 0 .and. ctqmc < 5, 'wrong file type' )
     call s_assert2( nskip >= 1, 'wrong number of skipped points' )
     call s_assert2( beta > zero, 'wrong inversion of temperature' )

! allocate memory
     allocate(tau(ntime),       stat=istat)

     allocate(grn(ntime,norbs), stat=istat)
     allocate(err(ntime,norbs), stat=istat)

     allocate(grn_bin(ntime,norbs,nbins), stat=istat)
     if ( istat /= 0 ) then
         call s_print_error('maketau','can not allocate enough memory')
     endif ! back if ( istat /= 0 ) block

! initialize variables
     exists = .false.

! initialize arrays
     tau = zero

     grn = zero
     err = zero

     grn_bin = zero

! build time slice mesh
     call s_linspace_d(zero, beta, ntime, tau)

! A: std mode
     if ( ctqmc == 1 .or. ctqmc == 2 ) then

! inquire file status: solver.green.dat
         inquire (file = 'solver.green.dat', exist = exists)
         if ( exists .eqv. .false. ) then
             call s_print_error('maketau','file solver.green.dat does not exist')
         endif ! back if ( exists .eqv. .false. ) block

! open solver.green.dat file
         write(mystd,'(2X,a)') 'Reading solver.green.dat ...'
         open(mytmp, file='solver.green.dat', form='formatted', status='unknown')

! read green's function data
         do i=1,norbs
             do j=1,ntime
                 read(mytmp,*) itmp, jtmp, rtmp, grn(j,i), err(j,i)
             enddo ! over j={1,ntime} loop
             read(mytmp,*) ! skip two lines
             read(mytmp,*)
         enddo ! over i={1,norbs} loop

! close solver.green.dat file
         close(mytmp)

         write(mystd,'(2X,a)') '>>> status: OK'
         write(mystd,*)

! ensure grn is positive
         grn = abs(grn)

     endif ! back if ( ctqmc == 1 .or. ctqmc == 2 ) block

! B: bin mode
     if ( ctqmc == 3 .or. ctqmc == 4 ) then
         do ibin=1,nbins

! convert ibin to sbin
             write(sbin,'(i10)') ibin

! inquire file status: solver.green.dat.ibin
             inquire (file = 'solver.green.dat.'//trim(adjustl(sbin)), exist = exists)
             if ( exists .eqv. .false. ) then
                 call s_print_error('maketau','file solver.green.dat.* does not exist')
             endif ! back if ( exists .eqv. .false. ) block

! open solver.green.bin file
             write(mystd,'(2X,a,i4)') 'Reading solver.green.dat ...', ibin
             open(mytmp, file='solver.green.dat.'//trim(adjustl(sbin)), form='formatted', status='unknown')

! read green's function data
             do i=1,norbs
                 do j=1,ntime
                     read(mytmp,*) itmp, jtmp, rtmp, grn_bin(j,i,ibin)
                 enddo ! over j={1,ntime} loop
                 read(mytmp,*) ! skip two lines
                 read(mytmp,*)
             enddo ! over i={1,norbs} loop

! close solver.green.bin file
             close(mytmp)

             write(mystd,'(2X,a)') '>>> status: OK'
             write(mystd,*)

         enddo ! over ibin={1,nbins} loop

! ensure grn_bin is positive
         grn_bin = abs(grn_bin)

! calculate averaged green's function
         grn = zero
         do ibin=1,nbins
             grn = grn + grn_bin(:,:,ibin)
         enddo ! over ibin={1,nbins} loop
         grn = grn / real(nbins)

! calculate standard variance
         do i=1,norbs
             do j=1,ntime
                 err(j,i) = zero
                 do ibin=1,nbins
                     err(j,i) = err(j,i) + ( grn_bin(j,i,ibin) - grn(j,i) )**2
                 enddo ! over ibin={1,nbins} loop
                 err(j,i) = sqrt( err(j,i) / real( nbins - 1 ) )
             enddo ! over j={1,ntime} loop
         enddo ! over i={1,norbs} loop
     endif ! back if ( ctqmc == 3 .or. ctqmc == 4 ) block

! open tau.grn.dat file, which is used as the input for the
! hibiscus/entropy or hibiscus/stoch codes
     write(mystd,'(2X,a)') 'Writing tau.grn.dat ...'

     open(mytmp, file='tau.grn.dat', form='formatted', status='unknown')

! write out data
     do i=1,norbs
         do j=1,ntime,nskip
             write(mytmp,'(3f16.8)') tau(j), grn(j,i), err(j,i)
         enddo ! over j={1,ntime} loop

! write two blank lines
         write(mytmp,*)
         write(mytmp,*)
     enddo ! over i={1,norbs} loop

! close tau.grn.dat file
     close(mytmp)

     write(mystd,'(2X,a)') '>>> status: OK'
     write(mystd,*)

! deallocate memory
     deallocate(tau)

     deallocate(grn)
     deallocate(err)

     deallocate(grn_bin)

  end program maketau
