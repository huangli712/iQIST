!!!-----------------------------------------------------------------------
!!! project : iqist @ narcissus
!!! program : ctqmc_dmft_selfer
!!!           ctqmc_dmft_conver
!!!           ctqmc_dmft_bethe
!!! source  : ctqmc_dmft.f90
!!! type    : subroutines
!!! author  : li huang (email:huangli@caep.cn)
!!! history : 09/16/2009 by li huang (created)
!!!           05/01/2025 by li huang (last modified)
!!! purpose : implement a hybridization expansion version continuous time
!!!           quantum Monte Carlo (CTQMC) quantum impurity solver plus
!!!           dynamical mean field theory (DMFT) self-consistent engine.
!!!           it is designed for the Hubbard model on a Bethe lattice.
!!! status  : unstable
!!! comment :
!!!-----------------------------------------------------------------------

!!
!! @sub ctqmc_dmft_selfer
!!
!! mini self-consistent engine for continuous time quantum Monte Carlo
!! quantum impurity solver plus dynamical mean field theory simulation
!!
  subroutine ctqmc_dmft_selfer()
     use constants, only : dp
     use constants, only : one, half
     use constants, only : czi
     use constants, only : mystd

     use control, only : cname
     use control, only : nband, norbs
     use control, only : mfreq
     use control, only : Uc, Jz
     use control, only : mune
     use control, only : alpha
     use control, only : myid, master

     use context, only : rmesh
     use context, only : eimp
     use context, only : grnf
     use context, only : wtau, wssf
     use context, only : hybf

     implicit none

!! local variables
     ! loop index over flavors
     integer  :: i

     ! loop index over frequencies
     integer  :: k

     ! status flag
     integer  :: istat

     ! effective chemical potential
     real(dp) :: qmune

     ! dummy hybridization function in matsubara frequency axis
     complex(dp), allocatable :: htmp(:,:,:)

!! [body

     ! allocate memory
     allocate(htmp(mfreq,norbs,norbs), stat=istat)
     if ( istat /= 0 ) then
         call s_print_error('ctqmc_dmft_selfer','can not allocate enough memory')
     endif ! back if ( istat /= 0 ) block

!!========================================================================
!!>>> starting self-consistent engine                                  <<<
!!========================================================================

     ! print necessary self-consistent simulation information
     if ( myid == master ) then ! only master node can do it
         write(mystd,'(2X,a)') cname//' >>> DMFT self-consistent engine running'
         write(mystd,'(4X,2a)') 'interacting lattice model  / ', 'Hubbard model'
         write(mystd,'(4X,2a)') 'density of states          / ', 'semicircular'
     endif ! back if ( myid == master ) block

     ! task 1: calculate new hybridization function
     !--------------------------------------------------------------------
     ! initialize htmp
     htmp = hybf

     ! apply the self-consistent condition. here we consider a Hubbard model
     ! on a bethe lattice. of course you can replace it with your implements
     call ctqmc_dmft_bethe(hybf, grnf)

     ! task 2: mix old and new hybridization functions
     !--------------------------------------------------------------------
     ! mix htmp and hybf using linear mixer
     call s_mix_z(size(hybf), htmp, hybf, alpha)

     ! task 3: calculate new bath weiss's function
     !--------------------------------------------------------------------
     ! determine effective chemical potential using
     !
     !     \mu_{eff} = (N - 0.5)*U - (N - 1)*2.5*J
     !
     ! where N is the number of bands
     qmune = mune
     qmune = qmune - ( real(nband) - half ) * Uc
     qmune = qmune + ( real(nband) - one ) * 2.5_dp * Jz

     ! apply dyson equation to get G^{-1}_0
     !
     !     G^{-1}_0 = i\omega + mu - E_{imp} - \Delta(i\omega)
     !
     ! where mu is the effective chemical potential
     do i=1,norbs
         do k=1,mfreq
             wssf(k,i,i) = czi * rmesh(k) + qmune - eimp(i) - hybf(k,i,i)
         enddo ! over k={1,mfreq} loop
     enddo ! over i={1,norbs} loop

     ! calculate G_0 via matrix inversion
     do k=1,mfreq
         call s_inv_z(norbs, wssf(k,:,:))
     enddo ! over k={1,mfreq} loop

     ! fourier transformation bath weiss's function from matsubara frequency
     ! space to imaginary time space
     call ctqmc_four_hybf(wssf, wtau)

     ! task 4: dump the calculated results
     !--------------------------------------------------------------------
     ! write out the new bath weiss's function in matsubara frequency axis
     if ( myid == master ) then ! only master node can do it
         call ctqmc_dump_wssf(wssf)
     endif ! back if ( myid == master ) block

     ! write out the new bath weiss's function in imaginary time axis
     if ( myid == master ) then ! only master node can do it
         call ctqmc_dump_wtau(wtau)
     endif ! back if ( myid == master ) block

!!========================================================================
!!>>> finishing self-consistent engine                                 <<<
!!========================================================================

     ! print necessary self-consistent simulation information
     if ( myid == master ) then ! only master node can do it
         write(mystd,'(4X,2a)') 'new hybridization function / ', 'calculated'
         write(mystd,'(4X,2a)') 'new bath weiss'//"'s"//' function  / ', 'calculated'
         write(mystd,'(2X,a)') cname//' >>> DMFT self-consistent engine shutdown'
         write(mystd,*)
     endif ! back if ( myid == master ) block

     ! deallocate memory
     deallocate(htmp)

!! body]

     return
  end subroutine ctqmc_dmft_selfer

!!
!! @sub ctqmc_dmft_conver
!!
!! check the convergence of matsubara self-energy function
!!
  subroutine ctqmc_dmft_conver(iter, conv)
     use constants, only : dp
     use constants, only : zero, one, two
     use constants, only : eps8
     use constants, only : mystd

     use control, only : cname
     use control, only : norbs
     use control, only : niter
     use control, only : mfreq
     use control, only : alpha
     use control, only : myid, master

     use context, only : sig1, sig2

     implicit none

!! external arguments
     ! current iteration number
     integer, intent(in)    :: iter

     ! convergence flag
     logical, intent(inout) :: conv

!! local parameters
     ! required minimum iteration number to achieive convergence
     integer, parameter :: minit = 16

!! local variables
     ! loop index over orbitals
     integer  :: i

     ! dummy variables
     real(dp) :: diff
     real(dp) :: norm
     real(dp) :: seps

!! [body

     ! write convergence information to screen
     if ( myid == master ) then ! only master node can do it
         write(mystd,'(2X,a)') cname//' >>> self-consistent iteration checker running'
         write(mystd,'(4X,a,i03.2)') 'maximum iteration / ', niter
         write(mystd,'(4X,a,e10.3)') 'maximum epsilon   / ', eps8
     endif ! back if ( myid == master ) block

     ! try to calculate diff and norm
     !
     ! why not using the whole matrix?
     ! since sometimes the off-diagonal elementes may be NaN!
     diff = zero
     do i=1,norbs
         diff = diff + abs( sum( sig2(:,i,i) - sig1(:,i,i) ) )
     enddo ! over i={1,norbs} loop
     !
     norm = zero
     do i=1,norbs
         norm = norm + abs( sum( sig2(:,i,i) + sig1(:,i,i) ) )
     enddo ! over i={1,norbs} loop
     norm = norm / two

     ! calculate seps
     seps = (diff / norm) / real(mfreq * norbs)

     ! judge convergence status
     conv = ( ( seps <= eps8 ) .and. ( iter >= minit ) )

     ! update sig1
     call s_mix_z(size(sig1), sig2, sig1, one - alpha)

     ! write convergence information to screen
     if ( myid == master ) then ! only master node can do it
         write(mystd,'(4X,a,i03.2)') 'current iteration / ', iter
         write(mystd,'(4X,a,e10.3)') 'current epsilon   / ', seps
         write(mystd,'(4X,a,X,a11)') 'selected watchdog / ', 'self-energy'
         write(mystd,'(4X,a,l2)')    'reach convergence / ', conv
         write(mystd,'(2X,a)') cname//' >>> self-consistent iteration checker shutdown'
         write(mystd,*)
     endif ! back if ( myid == master ) block

!! body]

     return
  end subroutine ctqmc_dmft_conver

!!
!! @sub ctqmc_dmft_bethe
!!
!! implement self-consistent condition: bethe lattice, semicircular
!! density of states, force a paramagnetic order, equal band width
!!
  subroutine ctqmc_dmft_bethe(hybf, grnf)
     use constants, only : dp

     use control, only : norbs
     use control, only : mfreq
     use control, only : part

     implicit none

!! external arguments
     ! hybridization function
     complex(dp), intent(out) :: hybf(mfreq,norbs,norbs)

     ! impurity green's function
     complex(dp), intent(in)  :: grnf(mfreq,norbs,norbs)

!! [body

     ! self-consistent condition is
     !
     !    Delta = t^2 G
     !
     ! this equation is valid for bethe lattice only
     hybf = part * part * grnf

!! body]

     return
  end subroutine ctqmc_dmft_bethe
