!!!-----------------------------------------------------------------------
!!! project : narcissus
!!! program : ctqmc_dmft_selfer
!!!           ctqmc_dmft_conver
!!!           ctqmc_dmft_bethe
!!! source  : ctqmc_dmft.f90
!!! type    : subroutines
!!! author  : li huang (email:lihuang.dmft@gmail.com)
!!! history : 09/16/2009 by li huang (created)
!!!           04/24/2017 by li huang (last modified)
!!! purpose : implement a self-consistent engine for dynamical mean field
!!!           theory (DMFT) simulation. it is designed for hybridization
!!!           expansion version continuous time quantum Monte Carlo (CTQMC)
!!!           quantum impurity solver and Hubbard model on bethe lattice
!!! status  : unstable
!!! comment :
!!!-----------------------------------------------------------------------

!!
!! @sub ctqmc_dmft_selfer
!!
!! the self-consistent engine for continuous time quantum Monte Carlo
!! quantum impurity solver plus dynamical mean field theory simulation
!!
  subroutine ctqmc_dmft_selfer()
     use constants, only : dp, one, half, czi, mystd

     use control, only : cname
     use control, only : nband, norbs
     use control, only : mfreq
     use control, only : Uc, Jz
     use control, only : mune, alpha
     use control, only : myid, master
     use context, only : tmesh, rmesh
     use context, only : eimp
     use context, only : grnf
     use context, only : wtau, wssf, hybf

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
         call s_print_error('ctqmc_dmft_selfer','can not allocate enough memory')
     endif ! back if ( istat /= 0 ) block

! initialize htmp
     htmp = hybf

! calculate new hybridization function using self-consistent condition
     call ctqmc_dmft_bethe(hybf, grnf)

! mixing new and old hybridization function: htmp and hybf
     call s_mix_z(size(hybf), htmp, hybf, alpha)

! \mu_{eff} = (N - 0.5)*U - (N - 1)*2.5*J
     qmune = ( real(nband) - half ) * Uc - ( real(nband) - one ) * 2.5_dp * Jz
     qmune = mune - qmune

! calculate new bath weiss's function
! G^{-1}_0 = i\omega + mu - E_{imp} - \Delta(i\omega)
     do i=1,norbs
         do k=1,mfreq
             wssf(k,i,i) = czi * rmesh(k) + qmune - eimp(i) - hybf(k,i,i)
         enddo ! over k={1,mfreq} loop
     enddo ! over i={1,norbs} loop

     do k=1,mfreq
         call s_inv_z(norbs, wssf(k,:,:))
     enddo ! over k={1,mfreq} loop

! fourier transformation bath weiss's function from matsubara frequency
! space to imaginary time space
     call ctqmc_four_hybf(wssf, wtau)

! write out the new bath weiss's function in matsubara frequency axis
     if ( myid == master ) then ! only master node can do it
         call ctqmc_dump_wssf(rmesh, wssf)
     endif ! back if ( myid == master ) block

! write out the new bath weiss's function in imaginary time axis
     if ( myid == master ) then ! only master node can do it
         call ctqmc_dump_wtau(tmesh, wtau)
     endif ! back if ( myid == master ) block

! print necessary self-consistent simulation information
     if ( myid == master ) then ! only master node can do it
         write(mystd,'(2X,a)') cname//' >>> DMFT hybridization function is updated'
         write(mystd,*)
     endif ! back if ( myid == master ) block

! deallocate memory
     deallocate(htmp)

     return
  end subroutine ctqmc_dmft_selfer

!!
!! @sub ctqmc_dmft_conver
!!
!! check the convergence of self-energy function
!!
  subroutine ctqmc_dmft_conver(iter, convergence)
     use constants, only : dp, zero, one, two, eps8, mystd

     use control, only : cname
     use control, only : norbs, niter
     use control, only : mfreq
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
     call s_mix_z(size(sig1), sig2, sig1, one - alpha)

! write convergence information to screen
     if ( myid == master ) then ! only master node can do it
         write(mystd,'(3(2X,a,i3))') cname//' >>> cur_iter:', iter, 'min_iter:', minit, 'max_iter:', niter
         write(mystd,'(2(2X,a,E12.4))') cname//' >>> sig_curr:', seps, 'eps_curr:', eps8
         write(mystd,'( (2X,a,L1))') cname//' >>> self-consistent iteration convergence is ', convergence
         write(mystd,*)
     endif ! back if ( myid == master ) block

     return
  end subroutine ctqmc_dmft_conver

!!
!! @sub ctqmc_dmft_bethe
!!
!! dmft self-consistent conditions, bethe lattice, semicircular density
!! of states, force a paramagnetic order, equal band width
!!
  subroutine ctqmc_dmft_bethe(hybf, grnf)
     use constants, only : dp

     use control, only : norbs
     use control, only : mfreq
     use control, only : part

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

!!
!! @sub ctqmc_dmft_anydos
!!
!! dmft self-consistent conditions, general density of states, calculate
!! the new hybridization function by using hilbert transformation and
!! numerical integration
!!
  subroutine ctqmc_dmft_anydos(hybf, grnf, sigf)
     use constants, only : dp, zero, czi, czero, mytmp
     use mmpi, only : mp_bcast, mp_barrier

     use control, only : norbs
     use control, only : mfreq
     use control, only : mune
     use control, only : myid, master
     use context, only : rmesh
     use context, only : eimp

     implicit none

! external arguments
! hybridization function
     complex(dp), intent(out) :: hybf(mfreq,norbs,norbs)

! impurity green's function
     complex(dp), intent(out) :: grnf(mfreq,norbs,norbs)

! self-energy function
     complex(dp), intent(in)  :: sigf(mfreq,norbs,norbs)

! local parameters
! size of real frequency grid
     integer, parameter :: ngrid = 801

! local variables
! loop index
     integer  :: i
     integer  :: j
     integer  :: k

! used to check whether the input file (solver.anydos.in) exists
     logical  :: exists

! energy interval for the frequency grid
     real(dp) :: step

! normalization factor for the orbital-resolved density of states
     real(dp) :: dsum

! complex(dp) dummy variables, to evalute the new impurity green's function
     complex(dp) :: caux

! real frequency grid, we assume it is equidistant
     real(dp) :: epsi(ngrid)

! orbital-resolved density of states, they can be non-degenerate
     real(dp) :: pdos(ngrid,norbs)

! initialize essential array
     epsi = zero
     pdos = zero

     hybf = czero
     grnf = czero

! read in frequency grid and density of states if available
!-------------------------------------------------------------------------
     if ( myid == master ) then ! only master node can do it
         exists = .false.

! inquire about file's existence
         inquire (file = 'solver.anydos.in', exist = exists)
         if ( exists .eqv. .false. ) then
             call s_print_error('ctqmc_dmft_anydos','file solver.anydos.in does not exist')
         endif ! back if ( exists .eqv. .false. ) block

! find input file: solver.anydos.in, read it
! read in frequency grid and density of states from solver.anydos.in
! note: we assume the density of states for different orbitals are
! degenerate. please make sure it again before using this subroutine.
         open(mytmp, file='solver.anydos.in', form='formatted', status='unknown')
         do i=1,ngrid
             read(mytmp,*) epsi(i), pdos(i,1)
             do j=2,norbs
                 pdos(:,j) = pdos(:,1)
             enddo ! over j={2,norbs} loop
         enddo ! over i={1,ngrid} loop
         close(mytmp)
     endif ! back if ( myid == master ) block

! broadcast epsi and pdos from master node to all children nodes
# if defined (MPI)

! broadcast data
     call mp_bcast(epsi, master)

! broadcast data
     call mp_bcast(pdos, master)

! block until all processes have reached here
     call mp_barrier()

# endif  /* MPI */

! evaluate energy interval
     step = epsi(2) - epsi(1)

! calculate new impurity green's function by using hilbert transformation
     do j=1,norbs

! evaluate dsum, which is used to normalize the density of states
         dsum = sum( pdos(:,j) ) * step

         do k=1,mfreq

! caux = i\omega + \mune - E_{imp} - Sigma(i\omega)
             caux = czi * rmesh(k) + mune - eimp(j) - sigf(k,j,j)

! perform numerical integration
             do i=1,ngrid
                 grnf(k,j,j) = grnf(k,j,j) + pdos(i,j) / ( caux - epsi(i) )
             enddo ! over i={1,ngrid} loop

! renormalize grnf
             grnf(k,j,j) = grnf(k,j,j) * ( step / dsum )

         enddo ! over k={1,mfreq} loop
     enddo ! over j={1,norbs} loop

! calculate G^{-1}, now grnf contains G^{-1}
     do k=1,mfreq
         call s_inv_z(norbs, grnf(k,:,:))
     enddo ! over k={1,mfreq} loop

! calculate final hybridization function using dyson's equation
     do k=1,mfreq
         do i=1,norbs
             hybf(k,i,i) = czi * rmesh(k) + mune - eimp(i) - sigf(k,i,i) - grnf(k,i,i)
         enddo ! over i={1,norbs} loop
     enddo ! over k={1,mfreq} loop

     return
  end subroutine ctqmc_dmft_anydos
