!-------------------------------------------------------------------------
! project : lavender
! program : ctqmc_dmft_selfer
!           ctqmc_dmft_conver
!           ctqmc_dmft_mixer
!           ctqmc_dmft_bethe
! source  : ctqmc_dmft.f90
! type    : subroutine
! author  : li huang (email:huangli712@yahoo.com.cn)
! history : 09/16/2009 by li huang
!           09/18/2009 by li huang
!           09/21/2009 by li huang
!           09/22/2009 by li huang
!           10/28/2009 by li huang
!           12/01/2009 by li huang
!           12/06/2009 by li huang
!           12/24/2009 by li huang
!           01/13/2010 by li huang
! purpose : self-consistent engine for dynamical mean field theory (DMFT)
!           simulation. it is only suitable for hybridization expansion
!           version continuous time quantum Monte Carlo (CTQMC) quantum
!           impurity solver plus bethe lattice model.
! input   :
! output  :
! status  : unstable
! comment :
!-------------------------------------------------------------------------

!>>> the self-consistent engine for continuous time quantum Monte Carlo
! quantum impurity solver plus dynamical mean field theory simulation
  subroutine ctqmc_dmft_selfer()
     use constants
     use control
     use context

     implicit none

! local variables
! loop index over flavors
     integer  :: i

! loop index over frequencies
     integer  :: k

! status flag
     integer  :: istat

! effective chemical potential
     real(dp) :: qmune

! dummy hybridization function, in matsubara frequency axis, matrix form
     complex(dp), allocatable :: htmp(:,:,:)

! allocate memory
     allocate(htmp(mfreq,norbs,norbs), stat=istat)
     if ( istat /= 0 ) then
         call ctqmc_print_error('ctqmc_dmft_selfer','can not allocate enough memory')
     endif

! initialize htmp
     htmp = czero

! calculate new hybridization function using self-consistent condition
     call ctqmc_dmft_bethe(htmp, grnf)

! mixing new and old hybridization function: htmp and hybf
     call ctqmc_dmft_mixer(hybf, htmp)

! update original hybridization function
     hybf = htmp

! \mu_{eff} = (N - 0.5)*U - (N - 1)*2.5*J
     qmune = ( real(nband) - half ) * Uc - ( real(nband) - one ) * 2.5_dp * Jz
     qmune = mune - qmune

! calculate new bath weiss's function
     do i=1,norbs
         do k=1,mfreq
             wssf(k,i,i) = cmesh(k) + qmune - eimp(i) - hybf(k,i,i)
         enddo ! over k={1,mfreq} loop
     enddo ! over i={1,norbs} loop

     do k=1,mfreq
         call ctqmc_zmat_inv(norbs, wssf(k,:,:))
     enddo ! over k={1,mfreq} loop

! fourier transformation bath weiss's function from matsubara frequency
! space to imaginary time space
     call ctqmc_fourier_hybf(wssf, wtau)

! write out the new bath weiss's function in matsubara frequency axis
     if ( myid == master ) then ! only master node can do it
         call ctqmc_dump_wssf(rmesh, wssf)
     endif

! write out the new bath weiss's function in imaginary time axis
     if ( myid == master ) then ! only master node can do it
         call ctqmc_dump_wtau(tmesh, wtau)
     endif

! write out the new hybridization function
     if ( myid == master ) then ! only master node can do it
         call ctqmc_dump_hybf(rmesh, hybf)
     endif

! print necessary self-consistent simulation information
     if ( myid == master ) then ! only master node can do it
         write(mystd,'(2X,a)') 'LAVENDER >>> DMFT hybridization function is updated'
         write(mystd,*)
     endif

! deallocate memory
     deallocate(htmp)

     return
  end subroutine ctqmc_dmft_selfer

!>>> check the convergence of self-energy function
  subroutine ctqmc_dmft_conver(iter, convergence)
     use constants
     use control
     use context

     implicit none

! external arguments
! current iteration number
     integer, intent(in)    :: iter

! convergence flag
     logical, intent(inout) :: convergence

! local parameters
! required minimum iteration number to achieive convergence
     integer, parameter :: minit = 16

! local variables
! loop index over orbitals
     integer  :: i

! dummy variables
     real(dp) :: diff
     real(dp) :: norm
     real(dp) :: seps

! calculate diff and norm
! why not using the whole matrix? since the off-diagonal elementes may be NaN!
     diff = zero
     do i=1,norbs
         diff = diff + abs( sum( sig2(:,i,i) - sig1(:,i,i) ) )
     enddo ! over i={1,norbs} loop

     norm = zero
     do i=1,norbs
         norm = norm + abs( sum( sig2(:,i,i) + sig1(:,i,i) ) )
     enddo ! over i={1,norbs} loop
     norm = norm / two

! calculate seps
     seps = (diff / norm) / real(mfreq * norbs)

! judge convergence status
     convergence = ( ( seps <= eps8 ) .and. ( iter >= minit ) )

! update sig1
     sig1 = sig1 * (one - alpha) + sig2 * alpha

! write convergence information to screen
     if ( myid == master ) then ! only master node can do it
         write(mystd,'(3(2X,a,i3))') 'LAVENDER >>> cur_iter:', iter, 'min_iter:', minit, 'max_iter:', niter
         write(mystd,'(2(2X,a,E10.4))') 'LAVENDER >>> sig_curr:', seps, 'eps_curr:', eps8
         write(mystd,'( (2X,a,L1))') 'LAVENDER >>> self-consistent iteration convergence is ', convergence
         write(mystd,*)
     endif

     return
  end subroutine ctqmc_dmft_conver

!>>> complex(dp) version, mixing two vectors using linear mixing algorithm
  subroutine ctqmc_dmft_mixer(vec1, vec2)
     use constants
     use control

     implicit none

! external arguments
! older green/weiss/sigma function in input
     complex(dp), intent(inout) :: vec1(mfreq,norbs,norbs)

! newer green/weiss/sigma function in input, newest in output
     complex(dp), intent(inout) :: vec2(mfreq,norbs,norbs)

! linear mixing scheme
     vec2 = vec1 * (one - alpha) + vec2 * alpha

     return
  end subroutine ctqmc_dmft_mixer

!>>> self-consistent conditions, bethe lattice, semicircular density of
! states, force a paramagnetic order, equal bandwidth
  subroutine ctqmc_dmft_bethe(hybf, grnf)
     use constants
     use control

     implicit none

! external arguments
! hybridization function
     complex(dp), intent(out) :: hybf(mfreq,norbs,norbs)

! impurity green's function
     complex(dp), intent(in)  :: grnf(mfreq,norbs,norbs)

! local variables
! loop index over orbitals
     integer :: i
     integer :: j

     do i=1,norbs
         do j=1,norbs
             hybf(:,j,i) = part * part * grnf(:,j,i)
         enddo ! over j={1,norbs} loop
     enddo ! over i={1,norbs} loop

     return
  end subroutine ctqmc_dmft_bethe
