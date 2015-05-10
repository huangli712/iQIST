!!!=========+=========+=========+=========+=========+=========+=========+!
!!! HIBISCUS/toolbox/makestd @ iQIST                                     !
!!!                                                                      !
!!! This tool is used to convert solver.sgm.dat.* to std.sgm.dat, which  !
!!! is necessary for the self-energy function analytical continuation    !
!!! code (hibiscus/swing)                                                !
!!! author  : Li Huang (at IOP/CAS & SPCLab/CAEP & UNIFR)                !
!!! version : v2015.01.06T                                               !
!!! status  : WARNING: IN TESTING STAGE, USE IT IN YOUR RISK             !
!!! comment : any question, please contact with huangli712@gmail.com     !
!!!=========+=========+=========+=========+=========+=========+=========+!

!!
!!
!! Introduction
!! ============
!!
!! The makestd code is often used to postprocess the self-energy function
!! data to generate suitable input files for the hibiscus/swing code.
!!
!! At first, you have to run the ctqmc impurity solver for many times to
!! generate a series of solver.sgm.dat file. And then this code is used
!! to deal with these solver.sgm.dat files to generate std.sgm.dat which
!! is just what the hibiscus/swing code needs.
!!
!! Usage
!! =====
!!
!! # ./mstd or bin/mstd.x
!!
!! Input
!! =====
!!
!! solver.sgm.dat.* (necessary)
!!
!! Output
!! ======
!!
!! std.sgm.dat
!!
!! Documents
!! =========
!!
!! For more details, please go to iqist/doc/manual directory.
!!
!!

  program makestd
     use constants, only : dp, zero, mystd, mytmp

     implicit none

! local control parameters
! number of bands
     integer  :: nband = 1

! number of spin orientation
! note: do not modify it
     integer  :: nspin = 2

! number of orbitals, norbs = nspin * nband
     integer  :: norbs = 2

! number of matsubara frequency points
     integer  :: nfreq = 8193

! number of data bins
     integer  :: nbins = 1

! local variables
! loop index
     integer  :: i
     integer  :: j

! loop index over QMC data bins
     integer  :: ibin

! dummy integer variables
     integer  :: itmp

! status flag
     integer  :: istat

! used to check whether the input file exists
     logical  :: exists

! real(dp) dummy variables
     real(dp) :: rtmp

! current bin index, string representation
     character(len=10) :: sbin

! matsubara axis
     real(dp), allocatable :: msh(:)

! averaged self-energy function data
     real(dp), allocatable :: sre(:,:)
     real(dp), allocatable :: sim(:,:)

! standard deviation for self-energy function data
     real(dp), allocatable :: sre_std(:,:)
     real(dp), allocatable :: sim_std(:,:)

! binning data for self-energy function
     real(dp), allocatable :: sre_bin(:,:,:)
     real(dp), allocatable :: sim_bin(:,:,:)

! print program header
     write(mystd,'(2X,a)') 'HIBISCUS/toolbox/makestd'
     write(mystd,'(2X,a)') '>>> Making mean value and standard deviation for self-energy function'
     write(mystd,*) ! print blank line

     write(mystd,'(2X,a)') 'Version: 2015.01.06T '//'(built at '//__TIME__//" "//__DATE__//')'
     write(mystd,'(2X,a)') 'Develop: by li huang (at IOP/CAS & SPCLab/CAEP & UNIFR)'
     write(mystd,'(2X,a)') 'Support: huangli712@gmail.com'
     write(mystd,'(2X,a)') 'License: GNU General Public License version 3'
     write(mystd,*) ! print blank line

! setup necessary parameters
     write(mystd,'(2X,a)')   'Number of bands (default = 1):'
     write(mystd,'(2X,a,$)') '>>> '
     read (*,*) nband
     write(mystd,*)
     norbs = nband * nspin

     write(mystd,'(2X,a)')   'Number of matsubara frequency points (default = 8193):'
     write(mystd,'(2X,a,$)') '>>> '
     read (*,*) nfreq
     write(mystd,*)

     write(mystd,'(2X,a)')   'Number of data bins (default = 1):'
     write(mystd,'(2X,a,$)') '>>> '
     read (*,*) nbins
     write(mystd,*)

! check the parameters
     call s_assert2( nband > 0 .and. nband < 8, 'wrong number of bands' )
     call s_assert2( nfreq > 0, 'wrong number of matsubara frequency points' )
     call s_assert2( nbins > 0, 'wrong number of data bins' )

! allocate memory
     allocate(msh(nfreq),                 stat=istat)

     allocate(sre(nfreq,norbs),           stat=istat)
     allocate(sim(nfreq,norbs),           stat=istat)

     allocate(sre_std(nfreq,norbs),       stat=istat)
     allocate(sim_std(nfreq,norbs),       stat=istat)

     allocate(sre_bin(nfreq,norbs,nbins), stat=istat)
     allocate(sim_bin(nfreq,norbs,nbins), stat=istat)
     if ( istat /= 0 ) then
         call s_print_error('makestd','can not allocate enough memory')
     endif ! back if ( istat / = 0 ) block

! initialize variables
     exists = .false.

! initialize arrays
     msh = zero

     sre = zero
     sim = zero

     sre_std = zero
     sim_std = zero

     sre_bin = zero
     sim_bin = zero

! collect the self-energy function data
     std_loop: do ibin=1,nbins

! convert ibin to sbin
         write(sbin,'(i10)') ibin

! inquire file status: solver.sgm.dat.ibin
         inquire (file = 'solver.sgm.dat.'//trim(adjustl(sbin)), exist = exists)
         if ( exists .eqv. .false. ) then
             call s_print_error('makestd','file solver.sgm.dat does not exist')
         endif ! back if ( exists .eqv. .false. ) block

! open solver.sgm.dat file
         write(mystd,'(2X,a,i4)') 'Reading solver.sgm.dat ...', ibin
         open(mytmp, file='solver.sgm.dat.'//trim(adjustl(sbin)), form='formatted', status='unknown')

! read self-energy function data
         do i=1,nband
             do j=1,nfreq
                 read(mytmp,*) itmp, msh(j), sre_bin(j,i,ibin), &
                                             sim_bin(j,i,ibin), &
                                       sre_bin(j,i+nband,ibin), &
                                       sim_bin(j,i+nband,ibin)
             enddo ! over j={1,nfreq} loop
             read(mytmp,*) ! skip two lines
             read(mytmp,*)
         enddo ! over i={1,nband} loop

! close solver.sgm.dat file
         close(mytmp)

         write(mystd,'(2X,a)') '>>> status: OK'
         write(mystd,*)

     enddo std_loop ! over ibin={1,nbins} loop

! process the self-energy function data
     write(mystd,'(2X,a)') 'Postprocessing self-energy function data ...'

! calculate the averaged self-energy function
     do ibin=1,nbins
         sre = sre + sre_bin(:,:,ibin)
         sim = sim + sim_bin(:,:,ibin)
     enddo ! over ibin={1,nbins} loop
     sre = sre / real(nbins)
     sim = sim / real(nbins)

! calculate the standard deviation for self-energy function
     do i=1,norbs
         do j=1,nfreq

! just for the real part
             rtmp = zero
             do ibin=1,nbins
                 rtmp = rtmp + sre_bin(j,i,ibin)**2
             enddo ! over ibin={1,nbins} loop
             rtmp = rtmp / real(nbins)
             sre_std(j,i) = sqrt( rtmp - sre(j,i)**2 )

! just for the imag part
             rtmp = zero
             do ibin=1,nbins
                 rtmp = rtmp + sim_bin(j,i,ibin)**2
             enddo ! over ibin={1,nbins} loop
             rtmp = rtmp / real(nbins)
             sim_std(j,i) = sqrt( rtmp - sim(j,i)**2 )

         enddo ! over j={1,nfreq} loop
     enddo ! over i={1,norbs} loop
     write(mystd,'(2X,a)') '>>> status: OK'
     write(mystd,*)

! open std.sgm.dat file, which is used as the input for swing code
     write(mystd,'(2X,a)') 'Writing std.sgm.dat ...'

! open data file
     open(mytmp, file='std.sgm.dat', form='formatted', status='unknown')

! write out data
     do i=1,norbs
         do j=1,nfreq
             write(mytmp,'(5f16.8)') msh(j), sre(j,i), sim(j,i), sre_std(j,i), sim_std(j,i)
         enddo ! over j={1,nfreq} loop

! write two blank lines
         write(mytmp,*)
         write(mytmp,*)
     enddo ! over i={1,norbs} loop

! close std.sgm.dat file
     close(mytmp)

     write(mystd,'(2X,a)') '>>> status: OK'
     write(mystd,*)

! deallocate memory
     deallocate(msh)

     deallocate(sre)
     deallocate(sim)

     deallocate(sre_std)
     deallocate(sim_std)

     deallocate(sre_bin)
     deallocate(sim_bin)

  end program makestd
