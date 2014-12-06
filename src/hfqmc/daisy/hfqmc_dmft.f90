!!!-----------------------------------------------------------------------
!!! project : daisy
!!! program : hfqmc_dmft_selfer
!!!           hfqmc_dmft_conver
!!!           hfqmc_dmft_nbethe
!!!           hfqmc_dmft_abethe
!!! source  : hfqmc_dmft.f90
!!! type    : subroutines
!!! author  : li huang (email:huangli712@gmail.com)
!!! history : 12/27/2005 by li huang
!!!           03/26/2010 by li huang
!!!           12/06/2014 by li huang
!!! purpose : the self-consistent engine for dynamical mean field theory
!!!           (DMFT) simulation. it is only suitable for Hirsch-Fye
!!!           quantum Monte Carlo (HFQMC) quantum impurity solver plus
!!!           bethe lattice model.
!!! status  : unstable
!!! comment :
!!!-----------------------------------------------------------------------

!!>>> hfqmc_dmft_selfer: the self-consistent engine for Hirsch-Fye quantum
!!>>> Monte Carlo quantum impurity solver plus dynamical mean field theory
!!>>> simulation
  subroutine hfqmc_dmft_selfer()
     use constants, only : dp, zero, one, czero, mystd

     use control, only : norbs
     use control, only : mfreq
     use control, only : ntime
     use control, only : alpha
     use control, only : myid, master
     use context, only : symm, tmesh, rmesh
     use context, only : gtau, wtau, wssf, sig2

     implicit none

! local variables
! loop index over orbitals
     integer :: i

! loop index over frequencies
     integer :: k

! status flag
     integer :: istat

! impurity green's function in imaginary time
     real(dp), allocatable :: grnt(:,:)

! bath weiss's function in imaginary time
     real(dp), allocatable :: wsst(:,:)

! impurity green's function in matsubara frequency
     complex(dp), allocatable :: grnw(:,:)

! bath weiss's function in matsubara frequency
     complex(dp), allocatable :: wssw(:,:)

! allocate memory for impurity green's function arrays
     allocate(grnt(ntime,norbs), stat=istat)
     if ( istat /= 0 ) then
         call s_print_error('hfqmc_dmft_selfer','can not allocate enough memory')
     endif ! back if ( istat /= 0 ) block

     allocate(grnw(mfreq,norbs), stat=istat)
     if ( istat /= 0 ) then
         call s_print_error('hfqmc_dmft_selfer','can not allocate enough memory')
     endif ! back if ( istat /= 0 ) block

! allocate memory for bath weiss's function arrays
     allocate(wsst(ntime,norbs), stat=istat)
     if ( istat /= 0 ) then
         call s_print_error('hfqmc_dmft_selfer','can not allocate enough memory')
     endif ! back if ( istat /= 0 ) block

     allocate(wssw(mfreq,norbs), stat=istat)
     if ( istat /= 0 ) then
         call s_print_error('hfqmc_dmft_selfer','can not allocate enough memory')
     endif ! back if ( istat /= 0 ) block

! preparing local arrays and matrices
! initialize real variables and arrays
     grnt = -gtau
     wsst =  zero

! initialize complex variables and arrays
     grnw = czero
     wssw = czero

! fourier impurity green's function from G(\tau) to G(i\omega)
     call hfqmc_fourier_t2w(grnt, grnw)

! select dynamical mean field theory self-consistent engine to calculate
! the new bath weiss's function
!-------------------------------------------------------------------------
! bethe lattice, semicircular density of states,                         !
! paramagnetic state, equal bandwidth                                    !
!-------------------------------------------------------------------------
     call hfqmc_dmft_nbethe(grnw, wssw)

!-------------------------------------------------------------------------
! bethe lattice, semicircular density of states,                         !
! antiferromagnetic state, equal bandwidth                               !
!-------------------------------------------------------------------------
!<     call hfqmc_dmft_abethe(grnw, wssw)

! mixing two bath weiss's function to produce a newer one
! before mixing: wssf is the old bath weiss's function, and wssw is the new one.
! after mixing: wssw is the newer bath weiss's function, and wssf is untouched.
     call s_mix_z(size(wssw), wssf, wssw, alpha)

! invfourier bath weiss's function from G0(i\omega) to G0(\tau)
     call hfqmc_fourier_w2t(wssw, wsst)

! symmetrize bath weiss's function if necessary
     call hfqmc_make_symm(symm, wsst)

! update original bath weiss's function
! in imaginary-time space
     wtau = -wsst

! in matsubara frequency space
     wssf =  wssw

! use Dyson equation to calculate the new self-energy function
     do i=1,norbs
         do k=1,mfreq
             sig2(k,i) = one / wssw(k,i) - one / grnw(k,i)
         enddo ! over k={1,mfreq} loop
     enddo ! over i={1,norbs} loop

! write new bath weiss's function to the disk file
     if ( myid == master ) then ! only master node can do it
         call hfqmc_dump_wtau(tmesh, wtau)
         call hfqmc_dump_wssf(rmesh, wssf)
     endif ! back if ( myid == master ) block

! print necessary self-consistent simulation information
     if ( myid == master ) then ! only master node can do it
         write(mystd,'(2X,a)') 'DAISY >>> DMFT bath weiss function is updated'
         write(mystd,*)
     endif ! back if ( myid == master ) block

! deallocate memory
     deallocate(grnt)
     deallocate(wsst)

     deallocate(grnw)
     deallocate(wssw)

     return
  end subroutine hfqmc_dmft_selfer

!!>>> hfqmc_dmft_conver: check the convergence of self-energy function
  subroutine hfqmc_dmft_conver(iter, convergence)
     use constants, only : dp, one, half, eps8, mystd

     use control, only : niter
     use control, only : alpha
     use control, only : myid, master
     use context, only : sig1, sig2

     implicit none

! external arguments
! current iteration number
     integer, intent(in)    :: iter

! convergence flag
     logical, intent(inout) :: convergence

! local parameters
! minimum iteration number to achieive convergence
     integer, parameter :: minit = 16

! local variables
! dummy variables
     real(dp) :: diff
     real(dp) :: norm
     real(dp) :: seps

! calculate diff and norm
     diff = sum ( abs( sig2 - sig1 ) )
     norm = sum ( abs( sig2 + sig1 ) ) * half

! calculate seps
     seps = (diff / norm) / real( size(sig2) )

! judge convergence status
     convergence = ( ( seps <= eps8 ) .and. ( iter >= minit ) )

! update sig1
     call s_mix_z(size(sig1), sig2, sig1, one - alpha)

! write convergence information to screen
     if ( myid == master ) then ! only master node can do it
         write(mystd,'(3(2X,a,i3))') 'DAISY >>> cur_iter:', iter, 'min_iter:', minit, 'max_iter:', niter
         write(mystd,'(2(2X,a,E12.4))') 'DAISY >>> sig_curr:', seps, 'eps_curr:', eps8
         write(mystd,'( (2X,a,L1))') 'DAISY >>> self-consistent iteration convergence is ', convergence
         write(mystd,*)
     endif ! back if ( myid == master ) block

     return
  end subroutine hfqmc_dmft_conver

!!>>> hfqmc_dmft_nbethe: self-consistent conditions, bethe lattice,
!!>>> semicircular density of states, force a paramagnetic order
  subroutine hfqmc_dmft_nbethe(grnw, wssw)
     use constants
     use control
     use context

     implicit none

! external arguments
! impurity green's function in matsubara frequency
     complex(dp), intent(in)  :: grnw(mfreq,norbs)

! bath weiss's function in matsubara frequency
     complex(dp), intent(out) :: wssw(mfreq,norbs)

! local variables
! loop index
     integer  :: i
     integer  :: j

! dummy real variables, t^2
     real(dp) :: t2

! effective chemical potential
     real(dp) :: qmune

! calculate t2 at first
     t2 = part * part

! \mu_{eff} = (N - 0.5)*U - (N - 1)*2.5*J
     qmune = ( real(nband) - half ) * Uc - ( real(nband) - one ) * 2.5_dp * Jz
     qmune = mune - qmune

     HUBBARD_BETHE_PARA: do i=1,mfreq
         do j=1,norbs
             wssw(i,j) = one / ( czi * rmesh(i) + qmune - eimp(j) - t2 * grnw(i,j) )
         enddo ! over j={1,norbs} loop
     enddo HUBBARD_BETHE_PARA ! over i={1,mfreq} loop

     return
  end subroutine hfqmc_dmft_nbethe

!!>>> hfqmc_dmft_abethe: self-cosistent conditions, bethe lattice,
!!>>> semicircular density of states, force a antiferromagnetic long
!!>>> range order
  subroutine hfqmc_dmft_abethe(grnw, wssw)
     use constants
     use control
     use context

     implicit none

! external arguments
! impurity green's function in matsubara frequency
     complex(dp), intent(in)  :: grnw(mfreq,norbs)

! bath weiss's function in matsubara frequency
     complex(dp), intent(out) :: wssw(mfreq,norbs)

! local variables
! loop index
     integer  :: i
     integer  :: j

! dummy real variables, t^2
     real(dp) :: t2

! effective chemical potential
     real(dp) :: qmune

! calculate t2 at first
     t2 = part * part

! \mu_{eff} = (N - 0.5)*U - (N - 1)*2.5*J
     qmune = ( real(nband) - half ) * Uc - ( real(nband) - one ) * 2.5_dp * Jz
     qmune = mune - qmune

     HUBBARD_BETHE_ANTI: do i=1,mfreq
         do j=1,nband
             wssw(i,j) = one / ( czi * rmesh(i) + qmune - eimp(j+nband)- t2 * grnw(i,j+nband) )
             wssw(i,j+nband) = one / ( czi * rmesh(i) + qmune - eimp(j)- t2 * grnw(i,j) )
         enddo ! over j={1,nband} loop
     enddo HUBBARD_BETHE_ANTI ! over i={1,mfreq} loop

     return
  end subroutine hfqmc_dmft_abethe
