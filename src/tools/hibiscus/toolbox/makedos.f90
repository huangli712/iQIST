!!!=========+=========+=========+=========+=========+=========+=========+!
!!! HIBISCUS/toolbox/makedos @ iQIST                                     !
!!!                                                                      !
!!! This tool is used to build gaussian, cubic lattice, bethe lattice,   !
!!! and lorentzian density of states, which will be used by the other    !
!!! hilbert transformation program                                       !
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
!! The makedos code is often used to generate typical density of states
!! for general lattice models. The output files can be processed by the
!! ctqmc_dmft_anydos() subroutine in the ctqmc_dmft.f90 file.
!!
!! Usage
!! =====
!!
!! # ./mdos or bin/mdos.x
!!
!! Input
!! =====
!!
!! N/A
!!
!! Output
!! ======
!!
!! dos.gauss.dat
!! dos.cubic.dat
!! dos.bethe.dat
!! dos.loren.dat
!!
!! Documents
!! =========
!!
!! For more details, please go to iqist/doc/manual directory.
!!
!!

  program makedos
     use constants, only : dp, zero, one, two, pi, mystd, mytmp

     implicit none

! local parameters
! \eta parameter, used to build delta function
     real(dp), parameter :: eta1  = 0.01_dp

! \eta^{2} parameter, used to build delta function
     real(dp), parameter :: eta2  = eta1**2

! local control parameters
! number of frequency points, the mesh is [-nw:nw]
     integer  :: nw   = 400

! number of k-points in one dimension
     integer  :: nk   = 200

! total number of k-points
     integer  :: nk3  = (200+1)**3

! hopping parameter t, bandwidth = 4*t
     real(dp) :: part = 1.00_dp

! maximum value of energy window
     real(dp) :: emax = 10.0_dp

! minimum value of energy window
     real(dp) :: emin =-10.0_dp

! local variables
! loop index
     integer  :: i

! loop index for k-point
     integer  :: k1, k2, k3

! status flag
     integer  :: istat

! energy interval
     real(dp) :: dw

! scale factor
     real(dp) :: fa

! frequency grid
     real(dp), allocatable :: mesh(:)

! density of states
     real(dp), allocatable :: pdos(:)

! precalculated cos function
     real(dp), allocatable :: cosk(:)

! print program header
     write(mystd,'(2X,a)') 'HIBISCUS/toolbox/makedos'
     write(mystd,'(2X,a)') '>>> Making density of states for general lattices'
     write(mystd,*) ! print blank line

     write(mystd,'(2X,a)') 'Version: 2014.10.11T '//'(built at '//__TIME__//" "//__DATE__//')'
     write(mystd,'(2X,a)') 'Develop: by li huang (at IOP/CAS & SPCLab/CAEP & UNIFR)'
     write(mystd,'(2X,a)') 'Support: huangli712@gmail.com'
     write(mystd,'(2X,a)') 'License: GNU General Public License version 3'
     write(mystd,*) ! print blank line

! setup necessary parameters
     write(mystd,'(2X,a)')   'Number of frequency points (default = 400):'
     write(mystd,'(2X,a,$)') '>>> '
     read (mystd,'(i)') nw
     write(mystd,*)

     write(mystd,'(2X,a)')   'Number of k points (default = 200):'
     write(mystd,'(2X,a,$)') '>>> '
     read (mystd,'(i)') nk
     write(mystd,*)
     nk3 = (nk+1)**3

     write(mystd,'(2X,a)')   'Hopping parameter t (default = 1.0):'
     write(mystd,'(2X,a,$)') '>>> '
     read (mystd,  *  ) part
     write(mystd,*)

     write(mystd,'(2X,a)')   'Energy window, maximum value (default = 10.0):'
     write(mystd,'(2X,a,$)') '>>> '
     read (mystd,  *  ) emax
     write(mystd,*)

     write(mystd,'(2X,a)')   'Energy window, minimum value (default =-10.0):'
     write(mystd,'(2X,a,$)') '>>> '
     read (mystd,  *  ) emin
     write(mystd,*)

! check the parameters
     call s_assert2( nw > 0, 'wrong number of frequency points' )
     call s_assert2( nk > 0, 'wrong number of k points' )
     call s_assert2( part > zero, 'wrong hopping parameter' )
     call s_assert2( emax > emin, 'wrong energy window' )

! allocate memory
     allocate(mesh(-nw:nw), stat=istat)
     allocate(pdos(-nw:nw), stat=istat)
     allocate(cosk( 0 :nk), stat=istat)
     if ( istat /= 0 ) then
         call s_print_error('makedos','can not allocate enough memory')
     endif ! back if ( istat /= 0 ) block

! initialize arrays
     mesh = zero
     pdos = zero
     cosk = zero

! build frequency mesh
     call s_linspace_d(emin, emax, 2 * nw + 1, mesh)

! build cosk
     call s_linspace_d(zero, pi, nk + 1, cosk)
     cosk = cos(cosk)

! gaussian density of states: d=\infty cubic lattice
!-------------------------------------------------------------------------
! D(\epsilon) = \frac{1}{t\sqrt{2\pi}} \exp{(-\frac{\epsilon^{2}}{2t^{2}})}
! see Rev. Mods. Phys. 68, 13, 1996, eq(19)
     write(mystd,'(2X,a)') 'Make gaussian   density of states ...'

     fa = two * part**2
     dw = one / ( part * sqrt( two * pi ) )
     do i=-nw,nw
         pdos(i) = dw * exp( -mesh(i)**2 / fa )
     enddo ! over i={-nw,nw} loop

! write it to data file
     open(mytmp, file='dos.gauss.dat', form='formatted', status='unknown')
     do i=-nw,nw
         write(mytmp,'(2f16.8)') mesh(i), pdos(i)
     enddo ! over i={-nw,nw} loop
     close(mytmp)

     write(mystd,'(2X,a)') '>>> status: OK'
     write(mystd,*)

! cubic lattice model
!-------------------------------------------------------------------------
! see Dieter Vollhardt
! Investigation of correlated electron systems using the limit of high dimensions
! url: http://www.physik.uni-augsburg.de/theo3/Research/high_dimensions.pdf
     write(mystd,'(2X,a)') 'Make cubic      density of states ...'

! setup factor, due to the rescaling of part
! \epsilon_k = -2t* \sum_{i=1}^{d} cos k_{i}
! t* = t/sqrt{Z}, t = part, Z = 2d
     fa = two * part / sqrt(6.0_dp)

! reset pdos
     pdos = zero
     BZ_K1_LOOP: do k1=0,nk         ! loop over kx
         BZ_K2_LOOP: do k2=0,nk     ! loop over ky
             BZ_K3_LOOP: do k3=0,nk ! loop over kz
                 dw = fa * ( cosk(k1) + cosk(k2) + cosk(k3) )
                 do i=-nw,nw        ! loop over freqency mesh
                     pdos(i) = pdos(i) + eta1 / ( ( mesh(i) + dw )**2 + eta2 )
                 enddo ! over i={-nw,nw} loop
             enddo BZ_K3_LOOP ! over k3={0,nk} loop
         enddo BZ_K2_LOOP ! over k2={0,nk} loop
     enddo BZ_K1_LOOP ! over k1={0,nk} loop

! renormalize pdos
     pdos = pdos / ( pi * nk3 )

! write it to data file
     open(mytmp, file='dos.cubic.dat', form='formatted', status='unknown')
     do i=-nw,nw
         write(mytmp,'(2f16.8)') mesh(i), pdos(i)
     enddo ! over i={-nw,nw} loop
     close(mytmp)

     write(mystd,'(2X,a)') '>>> status: OK'
     write(mystd,*)

! semicircular density of states: bethe lattice
!-------------------------------------------------------------------------
! D(\epsilon) = \frac{1} {2\pi t^{2}} \sqrt{4t^{2}-\epsilon^{2}}
! with \abs{\epsilon} < 2t
! see Rev. Mods. Phys. 68, 13, 1996, eq(21)
     write(mystd,'(2X,a)') 'Make bethe      density of states ...'

     fa = two * pi * part**2
     do i=-nw,nw
         if ( abs( mesh(i) ) > two * part ) then
             pdos(i) = zero
         else
             pdos(i) = sqrt( 4.0_dp * part**2 - mesh(i)**2 ) / fa
         endif ! back if ( abs( mesh(i) ) > two * part ) block
     enddo ! over i={-nw,nw} loop

! write it to data file
     open(mytmp, file='dos.bethe.dat', form='formatted', status='unknown')
     do i=-nw,nw
         write(mytmp,'(2f16.8)') mesh(i), pdos(i)
     enddo ! over i={-nw,nw} loop
     close(mytmp)

     write(mystd,'(2X,a)') '>>> status: OK'
     write(mystd,*)

! Lorentzian density of states
!-------------------------------------------------------------------------
! D(\epsilon) = \frac{t}{\pi (\epsilon^{2} + t^{2})}
! see Rev. Mods. Phys. 68, 13, 1996, eq(24)
     write(mystd,'(2X,a)') 'Make lorentzian density of states ...'

     do i=-nw,nw
         pdos(i) = part / ( pi * ( mesh(i)**2 + part**2 ) )
     enddo ! over i={-nw,nw} loop

! write it to data file
     open(mytmp, file='dos.loren.dat', form='formatted', status='unknown')
     do i=-nw,nw
         write(mytmp,'(2f16.8)') mesh(i), pdos(i)
     enddo ! over i={-nw,nw} loop
     close(mytmp)

     write(mystd,'(2X,a)') '>>> status: OK'
     write(mystd,*)

! deallocate memory
     deallocate(mesh)
     deallocate(pdos)
     deallocate(cosk)

  end program makedos
