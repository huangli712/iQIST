!!!=========+=========+=========+=========+=========+=========+=========+!
!!! HIBISCUS/toolbox/makesig @ iQIST                                     !
!!!                                                                      !
!!! This tool is used to perform analytical continuation for the self-   !
!!! energy function using the Pade approximation.                        !
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
!! The makesig code is often used to transform self-energy functions from
!! matsubara frequency representation to real frequency representation
!! via the Pade approximation. The results are very sensitive to the data
!! noises in the self-energy function.
!!
!! Usage
!! =====
!!
!! # ./msig or bin/msig.x
!!
!! Input
!! =====
!!
!! Output
!! ======
!!
!! Documents
!! =========
!!
!! For more details, please go to iqist/doc/manual directory.
!!
!!

  program makesig
     use constants

     implicit none

!-------------------------------------------------------------------------
! local setting parameters
!-------------------------------------------------------------------------
! number of orbitals, include spin degree of freedom
     integer  :: nq    = 2

! number of frequency points for matsubara mesh
     integer  :: nmesh = 256

! number of frequency points for real axis
     integer  :: ngrid = 1000

! number of frequency points for original self-energy function
     integer  :: nfreq = 8193

! inversion temperature
     real(dp) :: beta  = 10.0_dp

! energy step, used to build real axis
     real(dp) :: dw    = 0.01_dp

! local parameters to build z = e + i*delta
     real(dp) :: delta = 0.0001_dp
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

!-------------------------------------------------------------------------
! local variables
!-------------------------------------------------------------------------
! loop index
     integer  :: i
     integer  :: j
     integer  :: k

! status flag
     integer  :: istat

! logical file exist flag
     logical  :: fexist

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
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

! print program header
     write(mystd,'(2X,a)') 'MSIG'
     write(mystd,'(2X,a)') 'making sigma in real frequency axis'
     write(mystd,'(2X,a)') 'version: 2011.08.18T'
     write(mystd,*) ! print blank line

! setup necessary parameters
     write(mystd,'(2X,a)')   '>>> number of orbitals (default = 2):'
     write(mystd,'(2X,a,$)') '>>> '
     read (mystd,'(i)') nq
     write(mystd,*)

     write(mystd,'(2X,a)')   '>>> number of frequency points for matsubara mesh (default = 256):'
     write(mystd,'(2X,a,$)') '>>> '
     read (mystd,'(i)') nmesh
     write(mystd,*)

     write(mystd,'(2X,a)')   '>>> number of frequency points for real axis (default = 1000):'
     write(mystd,'(2X,a,$)') '>>> '
     read (mystd,'(i)') ngrid
     write(mystd,*)

     write(mystd,'(2X,a)')   '>>> number of frequency points for self-energy (default = 8193):'
     write(mystd,'(2X,a,$)') '>>> '
     read (mystd,'(i)') nfreq
     write(mystd,*)

     write(mystd,'(2X,a)')   '>>> inversion of temperature (default = 10.0):'
     write(mystd,'(2X,a,$)') '>>> '
     read (mystd,  *  ) beta
     write(mystd,*)

! allocate memory
     allocate(cmesh(nmesh),            stat=istat)
     allocate(rgrid(-ngrid:ngrid),     stat=istat)

     allocate(cdummy(nmesh),           stat=istat)
     allocate(rdummy(-ngrid:ngrid),    stat=istat)

     allocate(sigmaw(nfreq,nq),        stat=istat)
     allocate(sigmat(-ngrid:ngrid,nq), stat=istat)

!-------------------------------------------------------------------------

! inquire data file: solver.sgm.dat
     inquire(file = 'solver.sgm.dat', exist = fexist)
     if ( fexist == .false. ) then
         write(mystd,'(2X,a)') 'file solver.sgm.dat does not exist'
         STOP ! terminate this program
     endif

! open data file: solver.sgm.dat
     write(mystd,'(2X,a)') '>>> reading solver.sgm.dat ...'
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

!-------------------------------------------------------------------------

! build matsubara frequency grid
     do i=1,nmesh
         cmesh(i) = czi * ( two * real(i - 1) + one ) * pi / beta
     enddo ! over i={1,nmesh} loop

! build real frequency grid
     do i=-ngrid,ngrid
         rgrid(i) = real(i) * dw + czi * delta
     enddo ! over i={-ngrid,ngrid} loop

!-------------------------------------------------------------------------

! using Pade approximation to deal with self-energy function
     pade_loop: do i=1,nq
         write(mystd,'(2X,a,i2)') '>>> Pade transformation for sigma function # ', i

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
     enddo pade_loop ! over i={1,nq} loop

!-------------------------------------------------------------------------

! open data file: sig.sgm.dat
     write(mystd,'(2X,a)') '>>> writing sig.sgm.dat ...'
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

!-------------------------------------------------------------------------

! deallocate memory
     deallocate(cmesh)
     deallocate(rgrid)

     deallocate(cdummy)
     deallocate(rdummy)

     deallocate(sigmaw)
     deallocate(sigmat)

  end program makesig

!>>> using Pade approximation to transform green's function from
! matsubara frequency representation to real axis
  subroutine pade(nmesh, ngrid, cmesh, rgrid, matsubara, transform)
     use constants

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
