!!!=========+=========+=========+=========+=========+=========+=========+!
!!! HIBISCUS/toolbox/makescr @ iQIST                                     !
!!!                                                                      !
!!! This tool is used to general the kernel function, i.e., K(\tau), for !
!!! the narcissus code                                                   !
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
!! The makescr code is often used to calculate the K(\tau) from classic
!! models or from the screening spectral function W(\omega). The results
!! are necessary input for the narcissus code. In this code and narcissus
!! code, we assume that K(\tau) is degenerated for multi-orbital system.
!!
!! The W(\omega) data are often obtained by the RPA calculations which
!! are not included in the iQIST software package.
!!
!! Usage
!! =====
!!
!! # ./mscr or bin/mscr.x
!!
!! Input
!! =====
!!
!! src.frq.dat
!!
!! Output
!! ======
!!
!! src.tau.dat
!!
!! Documents
!! =========
!!
!! For more details, please go to iqist/doc/manual directory.
!!
!!

  program makescr
     use constants, only : dp, zero, one, two, pi, epss, mystd, mytmp

     implicit none

! local control parameters
! number of time slices in [0, \beta]
     integer  :: ntime = 1024

! number of frequency points
     integer  :: nfreq = 400

! number of cutoff frequency points, ncut <= nfreq
     integer  :: ncut  = 400

! inversion of temperature
     real(dp) :: beta  = 10.00_dp

! experienced parameters to build default screening spectral function
     real(dp) :: lc    = 0.625_dp
     real(dp) :: wc    = 4.000_dp
     real(dp) :: step  = 0.020_dp

! local variables
! loop index
     integer  :: i
     integer  :: j

! status flag
     integer  :: istat

! used to check whether the input file exists
     logical  :: exists

! real(dp) dummy variables
     real(dp) :: f1, f2, f3

! imaginary time mesh for kernel function
     real(dp), allocatable :: kmsh(:)

! kernel function K(\tau)
     real(dp), allocatable :: ktau(:)

! frequency mesh for screening spectral function
     real(dp), allocatable :: wmsh(:)

! screening spectral function, real part
     real(dp), allocatable :: wref(:)

! screening spectral function, imaginary part
     real(dp), allocatable :: wimf(:)

! print program header
     write(mystd,'(2X,a)') 'HIBISCUS/toolbox/makescr'
     write(mystd,'(2X,a)') '>>> Making kernel function from screening spectral function'
     write(mystd,*) ! print blank line

     write(mystd,'(2X,a)') 'Version: 2014.10.11T '//'(built at '//__TIME__//" "//__DATE__//')'
     write(mystd,'(2X,a)') 'Develop: by li huang (at IOP/CAS & SPCLab/CAEP & UNIFR)'
     write(mystd,'(2X,a)') 'Support: huangli712@gmail.com'
     write(mystd,'(2X,a)') 'License: GNU General Public License version 3'
     write(mystd,*) ! print blank line

! setup necessary parameters
     write(mystd,'(2X,a)')   'Number of time slices (default = 1024):'
     write(mystd,'(2X,a,$)') '>>> '
     read (mystd,'(i)') ntime
     write(mystd,*)

     write(mystd,'(2X,a)')   'Number of frequency points (default = 400):'
     write(mystd,'(2X,a,$)') '>>> '
     read (mystd,'(i)') nfreq
     write(mystd,*)

     write(mystd,'(2X,a)')   'Number of cutoff frequency points (default = 400):'
     write(mystd,'(2X,a,$)') '>>> '
     read (mystd,'(i)') ncut
     write(mystd,*)

     write(mystd,'(2X,a)')   'Inversion of temperature (default = 10.0):'
     write(mystd,'(2X,a,$)') '>>> '
     read (mystd,  *  ) beta
     write(mystd,*)

! check the parameters
     call s_assert2( ntime > 0, 'wrong number of time slices' )
     call s_assert2( nfreq > 0, 'wrong number of frequency points' )
     call s_assert2( ncut > 0, 'wrong number of cutoff frequency points' )
     call s_assert2( ncut <= nfreq, 'wrong number of cutoff frequency points' )
     call s_assert2( beta > zero, 'wrong inversion of temperature' )

! allocate memory
     allocate(kmsh(ntime), stat=istat)
     allocate(ktau(ntime), stat=istat)

     allocate(wmsh(nfreq), stat=istat)
     allocate(wref(nfreq), stat=istat)
     allocate(wimf(nfreq), stat=istat)
     if ( istat /= 0 ) then
         call s_print_error('makescr','can not allocate enough memory')
     endif ! back if ( istat / = 0 ) block

! initialize variables
     exists = .false.

! initialize arrays
     kmsh = zero
     ktau = zero

     wmsh = zero
     wref = zero
     wimf = zero

! build imaginary time slice mesh
     call s_linspace_d(zero, beta, ntime, kmsh)

! inquire file status: scr.frq.dat
     inquire (file = 'scr.frq.dat', exist = exists)

! read screening spectral function from scr.frq.dat file
     if ( exists .eqv. .true. ) then
         write(mystd,'(2X,a)') 'Reading scr.frq.dat ...'

! open data file: scr.frq.dat
         open(mytmp, file='scr.frq.dat', form='formatted', status='unknown')

! read in screeing spectral function
         do i=1,nfreq
             read(mytmp,*) wmsh(i), wref(i), wimf(i)
         enddo ! over i={1,nfreq} loop

! close data file
         close(mytmp)

         write(mystd,'(2X,a)') '>>> status: OK'
         write(mystd,*)
! build screening spectral function using default ohmic model, only for debug
     else
         write(mystd,'(2X,a)') 'WARNING: scr.frq.dat file do not exist! using ohmic model'
         do i=1,nfreq
             wmsh(i) = real(i - one) * step
             wref(i) = lc * wmsh(i) * log(abs(( wc + wmsh(i) ) / ( wc - wmsh(i) ))) - two * lc * wc
             wimf(i) = - lc * pi * wmsh(i)
         enddo ! over i={1,nfreq} loop
     endif ! back if ( exists .eqv. .true. ) block

! shift the first point for wmsh
     if ( abs( wmsh(1) ) < epss ) then
         wmsh(1) = 0.0001_dp
     endif ! back if ( abs( wmsh(1) ) < epss ) block

! calculate kernel function and energy shift for U and mu
     write(mystd,'(2X,a)') 'Calculating kernel function ...'

     do i=1,ntime
         do j=1,nfreq-1
             step = wmsh(j+1) - wmsh(j)
             f1 = cosh( (kmsh(i) - beta/two) * wmsh(j) ) - cosh( - beta/two * wmsh(j) )
             f1 = f1 * wimf(j) / ( wmsh(j) ** 2 ) / sinh( wmsh(j) * beta/two )
             f2 = cosh( (kmsh(i) - beta/two) * wmsh(j+1) ) - cosh( - beta/two * wmsh(j+1) )
             f2 = f2 * wimf(j+1) / ( wmsh(j+1) ** 2 ) / sinh( wmsh(j+1) * beta/two )
             ktau(i) = ktau(i) + (f1 + f2) * step / (two * pi)
         enddo ! over j={1,nfreq-1} loop
     enddo ! over i={1,ntime} loop

     f3 = zero
     do j=1,nfreq-1
         step = wmsh(j+1) - wmsh(j)
         f1 = wimf(j) / wmsh(j)
         f2 = wimf(j+1) / wmsh(j+1)
         f3 = f3 + (f1 + f2) * step / (two * pi)
     enddo ! over j={1,nfreq-1} loop

     write(mystd,'(2X,a)') '>>> status: OK'
     write(mystd,*)

! dump data to scr.tau.dat file, which is used as the input for narcissus program
     write(mystd,'(2X,a)') 'Writing scr.tau.dat ...'

! open data file: scr.tau.dat
     open(mytmp, file='scr.tau.dat', form='formatted', status='unknown')

! write ktau data to disk file
! note: the third column is the reference data, we use the ohmic model
     write(mytmp,'(a,f16.8,2X,a,f16.8)') '# u shift:', two * f3, 'mu shift:', f3
     do i=1,ntime
         f1 = beta * wc * sin( pi * kmsh(i) / beta ) / pi
         write(mytmp,'(3f16.8)') kmsh(i), ktau(i), lc * log(one + f1)
     enddo ! over i={1,ntime} loop

! close data file
     close(mytmp)

     write(mystd,'(2X,a)') '>>> status: OK'
     write(mystd,*)

! deallocate memory
     deallocate(kmsh)
     deallocate(ktau)

     deallocate(wmsh)
     deallocate(wref)
     deallocate(wimf)

  end program makescr
