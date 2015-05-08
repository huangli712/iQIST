!!!=========+=========+=========+=========+=========+=========+=========+!
!!! HIBISCUS/toolbox/makesig @ iQIST                                     !
!!!                                                                      !
!!! This tool is used to perform analytical continuation for the self-   !
!!! energy function using the Pade approximation                         !
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
!! The makesig code is often used to transform self-energy functions from
!! matsubara frequency representation to real frequency representation
!! via the Pade approximation. The results are very sensitive to the data
!! noises in the self-energy function. So we do not recommend to use this
!! code to perform analytical continuation for the self-energy function.
!! However, the hibiscus/swing code may be a better choice.
!!
!! Usage
!! =====
!!
!! # ./msig or bin/msig.x
!!
!! Input
!! =====
!!
!! solver.sgm.dat (necessary)
!!
!! Output
!! ======
!!
!! sig.sgm.dat
!!
!! Documents
!! =========
!!
!! For more details, please go to iqist/doc/manual directory.
!!
!!

  program makesig
     use constants, only : dp, zero, one, two, pi, czero, czi, mystd, mytmp

     implicit none

! local control parameters
! number of orbitals, include spin degree of freedom
     integer  :: nq    = 2

! number of selected frequency points in matsubara mesh
     integer  :: nmesh = 256

! number of frequency points in real axis
     integer  :: ngrid = 1000

! number of frequency points for original self-energy function
     integer  :: nfreq = 8193

! inversion of temperature
     real(dp) :: beta  = 10.0_dp

! energy step, used to build real axis
     real(dp) :: dw    = 0.01_dp

! local parameters to build z = e + i*delta
     real(dp) :: delta = 0.0001_dp

! local variables
! loop index
     integer  :: i
     integer  :: j
     integer  :: k

! status flag
     integer  :: istat

! logical file exist flag
     logical  :: exists

! dummy variables
     real(dp) :: r
     real(dp) :: r1, r2
     real(dp) :: r3, r4

! matsubara frequency grid
     complex(dp), allocatable :: cmesh(:)

! real frequency grid
     complex(dp), allocatable :: rgrid(:)

! complex(dp) dummy matrices for Pade approximation (original data)
     complex(dp), allocatable :: cdummy(:)

! complex(dp) dummy matrices for Pade approximation (transformed data)
     complex(dp), allocatable :: rdummy(:)

! original self-energy function on matsubara frequency representation
     complex(dp), allocatable :: sigmaw(:,:)

! self-energy function on real frequency representation
     complex(dp), allocatable :: sigmat(:,:)

! print program header
     write(mystd,'(2X,a)') 'HIBISCUS/toolbox/makesig'
     write(mystd,'(2X,a)') '>>> Making self-energy function in real frequency axis'
     write(mystd,*) ! print blank line

     write(mystd,'(2X,a)') 'Version: 2015.01.06T '//'(built at '//__TIME__//" "//__DATE__//')'
     write(mystd,'(2X,a)') 'Develop: by li huang (at IOP/CAS & SPCLab/CAEP & UNIFR)'
     write(mystd,'(2X,a)') 'Support: huangli712@gmail.com'
     write(mystd,'(2X,a)') 'License: GNU General Public License version 3'
     write(mystd,*) ! print blank line

! setup necessary parameters
     write(mystd,'(2X,a)')   'Number of orbitals (default = 2):'
     write(mystd,'(2X,a,$)') '>>> '
     read (*,*) nq
     write(mystd,*)

     write(mystd,'(2X,a)')   'Number of selected frequency points in matsubara mesh (default = 256):'
     write(mystd,'(2X,a,$)') '>>> '
     read (*,*) nmesh
     write(mystd,*)

     write(mystd,'(2X,a)')   'Number of frequency points in real axis (default = 1000):'
     write(mystd,'(2X,a,$)') '>>> '
     read (*,*) ngrid
     write(mystd,*)

     write(mystd,'(2X,a)')   'Number of frequency points for original self-energy (default = 8193):'
     write(mystd,'(2X,a,$)') '>>> '
     read (*,*) nfreq
     write(mystd,*)

     write(mystd,'(2X,a)')   'Inversion of temperature (default = 10.0):'
     write(mystd,'(2X,a,$)') '>>> '
     read (*,*) beta
     write(mystd,*)

! check the parameters
     call s_assert2( nq > 0 .and. nq < 15, 'wrong number of orbitals' )
     call s_assert2( nmesh > 0, 'wrong number of selected frequency points for matsubara mesh' )
     call s_assert2( ngrid > 0, 'wrong number of frequency points for real axis' )
     call s_assert2( nfreq > 0, 'wrong number of frequency points for original self-energy' )
     call s_assert2( beta > zero, 'wrong inversion of temperature' )

! allocate memory
     allocate(cmesh(nmesh),            stat=istat)
     allocate(rgrid(-ngrid:ngrid),     stat=istat)

     allocate(cdummy(nmesh),           stat=istat)
     allocate(rdummy(-ngrid:ngrid),    stat=istat)

     allocate(sigmaw(nfreq,nq),        stat=istat)
     allocate(sigmat(-ngrid:ngrid,nq), stat=istat)
     if ( istat /= 0 ) then
         call s_print_error('makesig','can not allocate enough memory')
     endif ! back if ( istat / = 0 ) block

! initialize arrays
     cmesh  = czero
     rgrid  = czero
     cdummy = czero
     rdummy = czero
     sigmaw = czero
     sigmat = czero

! build matsubara frequency grid
     call s_linspace_z(czi, czi * ( two * real(nmesh - 1) + one ), nmesh, cmesh)
     cmesh = cmesh * pi / beta

! build real frequency grid
     call s_linspace_z(-dcmplx(ngrid * dw), dcmplx(ngrid * dw), 2 * ngrid + 1, rgrid)
     rgrid = rgrid + czi * delta

! inquire data file: solver.sgm.dat
     inquire(file = 'solver.sgm.dat', exist = exists)
     if ( exists .eqv. .false. ) then
         call s_print_error('makesig','file solver.sgm.dat does not exist')
     endif ! back if ( exists .eqv. .false. ) block

! open data file: solver.sgm.dat
     write(mystd,'(2X,a)') 'Reading solver.sgm.dat ...'
     open(mytmp, file='solver.sgm.dat', form='formatted', status='unknown')

! read in self-energy function
     do i=1,nq/2
         do j=1,nfreq
             read(mytmp,*) k, r, r1, r2, r3, r4
             sigmaw(j,i)      = dcmplx(r1,r2)
             sigmaw(j,i+nq/2) = dcmplx(r3,r4)
         enddo ! over j={1,nfreq} loop
         read(mytmp,*) ! skip two blank lines
         read(mytmp,*)
     enddo ! over i={1,nq/2} loop

! close data file
     close(mytmp)
     write(mystd,'(2X,a)') '>>> status: OK'
     write(mystd,*)

! using Pade approximation to deal with self-energy function
     PA_LOOP: do i=1,nq
         write(mystd,'(2X,a,i2)') 'Doing Pade transformation for self-energy function # ', i

! initialize dummy arrays
         cdummy = czero
         rdummy = czero

! get the first nmesh points
         do j=1,nmesh
             cdummy(j) = sigmaw(j,i)
         enddo ! over j={1,nmesh} loop

! perform Pade approximation
         call pade(nmesh, 2*ngrid+1, cmesh, rgrid, cdummy, rdummy)

! save the results
         do k=-ngrid,ngrid
             sigmat(k,i) = rdummy(k)
         enddo ! over k={-ngrid,ngrid} loop

         write(mystd,'(2X,a)') '>>> status: OK'
         write(mystd,*)
     enddo PA_LOOP ! over i={1,nq} loop

! open data file: sig.sgm.dat
     write(mystd,'(2X,a)') 'Writing sig.sgm.dat ...'
     open(mytmp, file='sig.sgm.dat', form='formatted', status='unknown')

! write out self-energy function
     do i=1,nq
         do j=-ngrid,ngrid
             write(mytmp,'(3f16.8)') real( rgrid(j) ), real( sigmat(j,i) ), aimag( sigmat(j,i) )
         enddo ! over j={-ngrid,ngrid} loop
         write(mytmp,*) ! write two blank lines
         write(mytmp,*)
     enddo ! over i={1,nq} loop

! close data file
     close(mytmp)
     write(mystd,'(2X,a)') '>>> status: OK'
     write(mystd,*)

! deallocate memory
     deallocate(cmesh)
     deallocate(rgrid)

     deallocate(cdummy)
     deallocate(rdummy)

     deallocate(sigmaw)
     deallocate(sigmat)

  end program makesig

!!>>> pade: using Pade approximation to transform green's function from
!!>>> matsubara frequency representation to real axis
  subroutine pade(nmesh, ngrid, cmesh, rgrid, matsubara, transform)
     use constants, only : dp, czero, cone

     implicit none

! external arguments
! number of frequency point
     integer, intent(in) :: nmesh
     integer, intent(in) :: ngrid

! complex frequency grid
     complex(dp), intent(in) :: cmesh(nmesh)

! real frequency grid
     complex(dp), intent(in) :: rgrid(ngrid)

! original data in complex frequency grid
     complex(dp), intent(in) :: matsubara(nmesh)

! transformed data in real frequency grid
     complex(dp), intent(inout) :: transform(ngrid)

! local variables
! loop index
     integer :: i
     integer :: j

! dummy variables
     complex(dp) :: avec(0:nmesh)
     complex(dp) :: bvec(0:nmesh)

! pade coefficience matrix
     complex(dp) :: p(nmesh,nmesh)

! step 1: evaluate the Pade coefficience, stored in p matrix
     do i=1,nmesh
         p(1,i) = matsubara(i)
     enddo ! over i={1,nmesh} loop

     do j=2,nmesh
         do i=2,j
             p(i,j) = ( p(i-1,i-1) - p(i-1,j) ) / ( cmesh(j) - cmesh(i-1) ) / p(i-1,j)
         enddo ! over i={2,j} loop
     enddo ! over j={2,nmesh} loop

! step 2: evaluate the Pade approximation
     do i=1,ngrid
         avec(0) = czero
         avec(1) = p(1,1)

         bvec(0) = cone
         bvec(1) = cone

         do j=1,nmesh-1
             avec(j+1) = avec(j) + ( rgrid(i) - cmesh(j) ) * p(j+1,j+1) * avec(j-1)
             bvec(j+1) = bvec(j) + ( rgrid(i) - cmesh(j) ) * p(j+1,j+1) * bvec(j-1)
         enddo ! over j={1,nmesh-1} loop

         transform(i) = avec(nmesh) / bvec(nmesh)
     enddo ! over i={1,ngrid} loop

     return
  end subroutine pade
