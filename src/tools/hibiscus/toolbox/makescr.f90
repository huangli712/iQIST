!!!=========+=========+=========+=========+=========+=========+=========+!
!!! HIBISCUS/toolbox/makescr @ iQIST                                     !
!!!                                                                      !
!!! This tool is used to general the screening function, i.e., K(\tau),  !
!!! for the narcissus code                                               !
!!! author  : Li Huang (at IOP/CAS & SPCLab/CAEP & UNIFR)                !
!!! version : v2016.02.13T                                               !
!!! status  : WARNING: IN TESTING STAGE, USE IT IN YOUR RISK             !
!!! comment : any question, please contact with lihuang.dmft@gmail.com   !
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
!! The W(\omega) data are often obtained by the cRPA calculations. This
!! feature is unfortunately not included in the iQIST software package.
!!
!! In order to be compatible with the narcissus code, you have to rename
!! the output file (i.e., scr.tau.dat) to solver.ktau.in.
!!
!! Usage
!! =====
!!
!! # ./mscr or build/mscr.x
!!
!! Input
!! =====
!!
!! scr.frq.dat
!!
!! Output
!! ======
!!
!! scr.tau.dat
!!
!! Documents
!! =========
!!
!! For more details, please see the on line reference manual.
!!
!!

  program makescr
     use constants, only : dp, zero, one, two, pi, epss, mystd, mytmp

     implicit none

! local control parameters
! model for screening functions
! model == 1, dervied from cRPA calculation
! model == 2, plasmon pole model
! model == 3, ohmic model
     integer  :: model = 1

! number of time slices in [0, \beta]
     integer  :: ntime = 1024

! number of frequency points
     integer  :: nfreq = 400

! number of cutoff frequency points, ncut <= nfreq
! we use it to avoid numerical overflow in sinh function
     integer  :: ncut  = 400

! inversion of temperature
     real(dp) :: beta  = 10.00_dp

! experienced parameters to build default screening function
     real(dp) :: lc    = 0.625_dp
     real(dp) :: wc    = 4.000_dp

! local variables
! loop index
     integer  :: i
     integer  :: j

! status flag
     integer  :: istat

! used to check whether the input file exists
     logical  :: exists

! real(dp) dummy variables
     real(dp) :: f1, f2, step

! Bose factor zb
     real(dp) :: zb

! imaginary time mesh for screening function
     real(dp), allocatable :: kmsh(:)

! screening function K(\tau)
     real(dp), allocatable :: ktau(:)

! first derivates of screening function K'(\tau)
     real(dp), allocatable :: ptau(:)

! frequency mesh for screening spectral function
     real(dp), allocatable :: wmsh(:)

! screening spectral function, real part
     real(dp), allocatable :: wref(:)

! screening spectral function, imaginary part
     real(dp), allocatable :: wimf(:)

! print program header
     write(mystd,'(2X,a)') 'HIBISCUS/toolbox/makescr'
     write(mystd,'(2X,a)') '>>> Making screening function and its first derivates'
     write(mystd,*) ! print blank line

     write(mystd,'(2X,a)') 'Version: 2016.02.13T '//'(built at '//__TIME__//" "//__DATE__//')'
     write(mystd,'(2X,a)') 'Develop: by li huang (at IOP/CAS & SPCLab/CAEP & UNIFR)'
     write(mystd,'(2X,a)') 'Support: lihuang.dmft@gmail.com'
     write(mystd,'(2X,a)') 'License: GNU General Public License version 3'
     write(mystd,*) ! print blank line

! setup necessary parameters
     write(mystd,'(2X,a)')   'Model for screening functions (default = 1):'
     write(mystd,'(2X,a)')   '  model : 1, dervied from cRPA calculation'
     write(mystd,'(2X,a)')   '  model : 2, plasmon pole model'
     write(mystd,'(2X,a)')   '  model : 3, ohmic model'
     write(mystd,'(2X,a,$)') '>>> '
     read (*,*) model
     write(mystd,*)

     write(mystd,'(2X,a)')   'Number of time slices (default = 1024):'
     write(mystd,'(2X,a,$)') '>>> '
     read (*,*) ntime
     write(mystd,*)

     write(mystd,'(2X,a)')   'Number of frequency points (default = 400):'
     write(mystd,'(2X,a,$)') '>>> '
     read (*,*) nfreq
     write(mystd,*)

     write(mystd,'(2X,a)')   'Number of cutoff frequency points (default = 400):'
     write(mystd,'(2X,a,$)') '>>> '
     read (*,*) ncut
     write(mystd,*)

     write(mystd,'(2X,a)')   'Inversion of temperature (default = 10.0):'
     write(mystd,'(2X,a,$)') '>>> '
     read (*,*) beta
     write(mystd,*)

     write(mystd,'(2X,a)')   'Experienced parameter for screening function lc (default = 0.625):'
     write(mystd,'(2X,a,$)') '>>> '
     read (*,*) lc
     write(mystd,*)

     write(mystd,'(2X,a)')   'Experienced parameter for screening function wc (default = 4.000):'
     write(mystd,'(2X,a,$)') '>>> '
     read (*,*) wc
     write(mystd,*)

! check the parameters
     call s_assert2( model > 0 .and. model < 4, 'wrong model for screening functions' )
     call s_assert2( ntime > 0, 'wrong number of time slices' )
     call s_assert2( nfreq > 0, 'wrong number of frequency points' )
     call s_assert2( ncut > 0 .and. ncut <= nfreq, 'wrong number of cutoff frequency points' )
     call s_assert2( beta > zero, 'wrong inversion of temperature' )
     call s_assert2( lc > zero, 'wrong experienced parameter for screening function lc' )
     call s_assert2( wc > zero, 'wrong experienced parameter for screening function wc' )

! allocate memory
     allocate(kmsh(ntime), stat=istat)
     allocate(ktau(ntime), stat=istat)
     allocate(ptau(ntime), stat=istat)

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
     ptau = zero

     wmsh = zero
     wref = zero
     wimf = zero

! build imaginary time slice mesh
     call s_linspace_d(zero, beta, ntime, kmsh)

! case 1: calculate screening function using screening spectral function
     if ( model == 1 ) then

! inquire file status: scr.frq.dat
         inquire (file = 'scr.frq.dat', exist = exists)
         if ( exists .eqv. .false. ) then
             call s_print_error('makescr','file scr.frq.dat does not exist')
         endif ! back if ( exists .eqv. .false. ) block

! read screening spectral function from scr.frq.dat file
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

! shift the first point for wmsh if it is zero
         if ( abs( wmsh(1) ) < epss ) then
             wmsh(1) = epss
         endif ! back if ( abs( wmsh(1) ) < epss ) block

! calculate screening function and its first derivates
         write(mystd,'(2X,a)') 'Calculating screening function ...'
         do i=1,ntime
             do j=1,ncut-1
                 step = wmsh(j+1) - wmsh(j)
                 f1 = cosh( (kmsh(i) - beta/two) * wmsh(j) ) - cosh( - beta/two * wmsh(j) )
                 f1 = f1 * wimf(j) / ( wmsh(j) ** 2 ) / sinh( wmsh(j) * beta/two )
                 f2 = cosh( (kmsh(i) - beta/two) * wmsh(j+1) ) - cosh( - beta/two * wmsh(j+1) )
                 f2 = f2 * wimf(j+1) / ( wmsh(j+1) ** 2 ) / sinh( wmsh(j+1) * beta/two )
                 ktau(i) = ktau(i) + (f1 + f2) * step / (two * pi)

                 f1 = sinh( (kmsh(i) - beta/two) * wmsh(j) ) / sinh( wmsh(j) * beta/two )
                 f1 = f1 * wimf(j) / wmsh(j)
                 f2 = sinh( (kmsh(i) - beta/two) * wmsh(j+1) ) / sinh( wmsh(j+1) * beta/two )
                 f2 = f2 * wimf(j+1) / wmsh(j+1)
                 ptau(i) = ptau(i) + (f1 + f2) * step / (two * pi)
             enddo ! over j={1,ncut-1} loop
         enddo ! over i={1,ntime} loop
         write(mystd,'(2X,a)') '>>> status: OK'
         write(mystd,*)

! calculate Bose factor Zb
         write(mystd,'(2X,a)') 'Calculating Bose factor Zb ...'
         zb = zero
         do j=1,nfreq-1
             step = wmsh(j+1) - wmsh(j)
             f1 = wimf(j) / ( wmsh(j) ** 2 )
             f2 = wimf(j+1) / ( wmsh(j+1) ** 2 )
             zb = zb + (f1 + f2) * step / (two * pi)
         enddo ! over j={1,nfreq-1} loop
         write(mystd,'(2X,a,f12.6)') 'Bose factor Zb:', exp(zb)
         write(mystd,'(2X,a)') '>>> status: OK'
         write(mystd,*)

     endif ! back if ( model == 1 ) block

! case 2: build screening function using default plasmon pole model
     if ( model == 2 ) then
         do i=1,ntime
             ktau(i) = (lc / wc)**2 / sinh(beta * wc / two)
             ktau(i) = ktau(i) * ( cosh(beta * wc / two) - cosh(beta * wc / two - kmsh(i) * wc) )
             ptau(i) = (lc / wc)**2 / sinh(beta * wc / two)
             ptau(i) = ptau(i) * sinh(beta * wc / two - kmsh(i) * wc) * wc
         enddo ! over i={1,ntime} loop
     endif ! back if ( model == 2 ) block

! case 3: build screening function using default ohmic model
     if ( model == 3 ) then
         do i=1,ntime
             ktau(i) = lc * log(one + beta * wc * sin(pi * kmsh(i) / beta) / pi)
             ptau(i) = lc * wc * cos(pi * kmsh(i) / beta) / (one + beta * wc * sin(pi * kmsh(i) / beta) / pi)
         enddo ! over i={1,ntime} loop
     endif ! back if ( model == 3 ) block

! dump data to scr.tau.dat file, which is served as the input for the
! narcissus code
     write(mystd,'(2X,a)') 'Writing scr.tau.dat ...'

! open data file: scr.tau.dat
     open(mytmp, file='scr.tau.dat', form='formatted', status='unknown')

! write ktau and ptau data to disk file
     write(mytmp,'(a,f16.8,2X,a,f16.8)') '# u shift:', two * ptau(1), 'mu shift:', ptau(1)
     do i=1,ntime
         write(mytmp,'(3f16.8)') kmsh(i), ktau(i), ptau(i)
     enddo ! over i={1,ntime} loop

! close data file
     close(mytmp)

     write(mystd,'(2X,a)') '>>> status: OK'
     write(mystd,*)

! deallocate memory
     deallocate(kmsh)
     deallocate(ktau)
     deallocate(ptau)

     deallocate(wmsh)
     deallocate(wref)
     deallocate(wimf)

  end program makescr
