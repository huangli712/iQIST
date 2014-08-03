!=========+=========+=========+=========+=========+=========+=========+>>>
! from screening spectral function W(\omega) calculate kernel function   !
! (K(\tau)), which can be used to feed narcissus code.                   !
! author  : li huang                                                     !
! version : v2011.08.18T                                                 !
! status  : WARNING: IN TESTING STAGE, USE IT IN YOUR RISK               !
! comment : any question, please contact with huangli712@yahoo.com.cn    !
!=========+=========+=========+=========+=========+=========+=========+>>>

  program makescr
     use constants

     implicit none

!-------------------------------------------------------------------------
! local setting parameters
!-------------------------------------------------------------------------
! number of time slice in [0, \beta]
     integer  :: ntime = 1024

! number of frequency points
     integer  :: nfreq = 400

! experienced parameters to build default screening spectral function
     real(dp) :: lc    = 0.625_dp
     real(dp) :: wc    = 4.000_dp
     real(dp) :: step  = 0.020_dp
     real(dp) :: scal  = 1.500_dp

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
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

! print program header
     write(mystd,'(2X,a)') 'MSCR'
     write(mystd,'(2X,a)') 'making kernel function from screening spectral function'
     write(mystd,'(2X,a)') 'version: 2011.08.18T'
     write(mystd,*) ! print blank line

! setup necessary parameters
     write(mystd,'(2X,a)')   '>>> number of time slice (default = 1024):'
     write(mystd,'(2X,a,$)') '>>> '
     read (mystd,'(i)') ntime
     write(mystd,*)

     write(mystd,'(2X,a)')   '>>> number of frequency points (default = 400):'
     write(mystd,'(2X,a,$)') '>>> '
     read (mystd,'(i)') nfreq
     write(mystd,*)

     write(mystd,'(2X,a)')   '>>> inversion of temperature (default = 10.0):'
     write(mystd,'(2X,a,$)') '>>> '
     read (mystd,  *  ) beta
     write(mystd,*)

! allocate memory
     allocate(kmsh(ntime), stat=istat)
     allocate(ktau(ntime), stat=istat)

     allocate(wmsh(nfreq), stat=istat)
     allocate(wref(nfreq), stat=istat)
     allocate(wimf(nfreq), stat=istat)

! initialize variables
     exists = .false.

! initialize arrays
     kmsh = zero
     ktau = zero

     wmsh = zero
     wref = zero
     wimf = zero

! build imaginary time slice mesh
     do i=1,ntime
         kmsh(i) = beta * real( i - 1 ) / real(ntime - 1)
     enddo ! over i={1,ntime} loop

! inquire file status: scr.frq.dat
     inquire (file = 'scr.frq.dat', exist = exists)

! read screening spectral function from scr.frq.dat file
     if ( exists == .true. ) then
         write(mystd,'(2X,a)') '>>> reading scr.frq.dat ...'

! open data file: scr.frq.dat
         open(mytmp, file='scr.frq.dat', form='formatted', status='unknown')

! read in screeing spectral function
         do i=1,nfreq
             read(mytmp,*) wmsh(i), wref(i), wimf(i)
         enddo ! over i={1,nfreq} loop
         scal = one

! close data file
         close(mytmp)

         write(mystd,'(2X,a)') '>>> status: OK'
         write(mystd,*)
! build screening spectral function using default ohmic model, only for debug
     else
         write(mystd,'(2X,a)') 'WARNING: scr.frq.dat file do not exist! using default screeing function'
         do i=1,nfreq
             wmsh(i) = real(i - one) * step
             wref(i) = lc * wmsh(i) * log(abs(( wc + wmsh(i) ) / ( wc - wmsh(i) ))) - two * lc * wc
             if ( wmsh(i) <= wc / scal ) then
                 wimf(i) = - lc * pi * wmsh(i)
             else
                 wimf(i) = zero
             endif
         enddo ! over i={1,nfreq} loop
     endif ! back if ( exists == .true. ) block

! calculate kernel function and energy shift for U and mu
     write(mystd,'(2X,a)') '>>> calculating kernel function ...'

     do i=1,ntime
         do j=1,nfreq-1
             step = wmsh(j+1) - wmsh(j)
!<             if ( j == 1 ) then
!<                  f1 = cosh( (kmsh(i) - beta/two) * epss ) - cosh( - beta/two * epss )
!<                  f1 = f1 * wimf(j) / ( epss ** 2 ) / sinh( epss * beta/two )
!<             else
                  f1 = cosh( (kmsh(i) - beta/two) * wmsh(j) ) - cosh( - beta/two * wmsh(j) )
                  f1 = f1 * wimf(j) / ( wmsh(j) ** 2 ) / sinh( wmsh(j) * beta/two )
!<             endif ! back if ( j == 1 ) block
             f2 = cosh( (kmsh(i) - beta/two) * wmsh(j+1) ) - cosh( - beta/two * wmsh(j+1) )
             f2 = f2 * wimf(j+1) / ( wmsh(j+1) ** 2 ) / sinh( wmsh(j+1) * beta/two )
             ktau(i) = ktau(i) + (f1 + f2) * step / (two * pi)
         enddo ! over j={1,nfreq-1} loop
     enddo ! over i={1,ntime} loop

     f3 = zero
     do j=1,nfreq-1
         step = wmsh(j+1) - wmsh(j)
!<         if ( j == 1 ) then
!<              f1 = wimf(j) / epss
!<         else
              f1 = wimf(j) / wmsh(j)
!<         endif ! back if ( j == 1 ) block
         f2 = wimf(j+1) / wmsh(j+1)
         f3 = f3 + (f1 + f2) * step / (two * pi)
     enddo ! over j={1,nfreq-1} loop
     f3 = f3 * scal

     write(mystd,'(2X,a)') '>>> status: OK'
     write(mystd,*)

! dump data to scr.tau.dat file, which is used as the input for narcissus program
     write(mystd,'(2X,a)') '>>> writing scr.tau.dat ...'

! open data file: scr.tau.dat
     open(mytmp, file='scr.tau.dat', form='formatted', status='unknown')

! write ktau data to disk file
! note: the third column is the reference data
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
