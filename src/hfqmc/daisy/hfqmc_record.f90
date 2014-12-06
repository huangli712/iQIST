!!!-----------------------------------------------------------------------
!!! project : daisy
!!! program : hfqmc_symm_spin
!!!           hfqmc_symm_band
!!!           hfqmc_make_symm
!!!           hfqmc_make_smth
!!!           hfqmc_make_freq
!!!           hfqmc_make_quas
!!!           hfqmc_make_nmat
!!!           hfqmc_reduce_gtau
!!!           hfqmc_reduce_nmat
!!! source  : hfqmc_record.f90
!!! type    : subroutine
!!! author  : li huang (email:huangli712@gmail.com)
!!! history : 01/07/2006 by li huang
!!!           08/25/2010 by li huang
!!!           12/06/2014 by li huang
!!! purpose : To build impurity green's function, bath weiss's function
!!!           and self-energy function in matsubara frequency space. some
!!!           auxiliary subroutines, such as symmetrizing and smoothing
!!!           tools are provided as well.
!!! status  : unstable
!!! comment :
!!!-----------------------------------------------------------------------

!!>>> hfqmc_symm_spin: enforce symmetry to the green's functions over spin
  subroutine hfqmc_symm_spin(norbs, gtau)
     use constants, only : dp, two

     implicit none

! external arguments
! number of orbitals
     integer, intent(in)     :: norbs

! green's function
     real(dp), intent(inout) :: gtau(norbs)

! local variables
! loop index
     integer  :: i

! number of correlated bands
     integer  :: nbnd

! dummy variables
     real(dp) :: raux

     nbnd = norbs / 2
     do i=1,nbnd
         raux = ( gtau(i) + gtau(i+nbnd) ) / two
         gtau(i) = raux
         gtau(i+nbnd) = raux
     enddo ! over i={1,nbnd} loop

     return
  end subroutine hfqmc_symm_spin

!!>>> hfqmc_symm_band: enforce symmetry to the green's functions over band
  subroutine hfqmc_symm_band(norbs, symm, gtau)
     use constants, only : dp, zero

     implicit none

! external arguments
! number of orbitals
     integer, intent(in)     :: norbs

! symmetry vector for correlated band
     integer, intent(in)     :: symm(norbs)

! green's function
     real(dp), intent(inout) :: gtau(norbs)

! local variables
! loop index for orbitals
     integer  :: ibnd
     integer  :: jbnd

! real(dp) dummy variables
     real(dp) :: raux

! histogram vector
     integer  :: hist(norbs)

! build histogram
     hist = 0
     do ibnd=1,norbs
         hist(symm(ibnd)) = hist(symm(ibnd)) + 1
     enddo ! over ibnd={1,norbs} loop

! perform symmetrization for those orbitals which symm index are identity
     do ibnd=1,norbs
         if ( hist(ibnd) > 0 ) then         ! need to enforce symmetry
             raux = zero

             do jbnd=1,norbs                ! gather the data
                 if ( symm(jbnd) == ibnd ) then
                     raux = raux + gtau(jbnd)
                 endif ! back if ( symm(jbnd) == ibnd ) block
             enddo ! over jbnd={1,norbs} loop

             raux = raux / real(hist(ibnd)) ! calculate average value

             do jbnd=1,norbs                ! setup it
                 if ( symm(jbnd) == ibnd ) then
                     gtau(jbnd) = raux
                 endif ! back if ( symm(jbnd) == ibnd ) block
             enddo ! over jbnd={1,norbs} loop
         endif ! back if ( hist(ibnd) > 0 ) block
     enddo ! over ibnd={1,norbs} loop

     return
  end subroutine hfqmc_symm_band

!!>>> hfqmc_make_symm: to symmetrize the final green's function
  subroutine hfqmc_make_symm(symm, gtau)
     use constants, only : dp

     use control, only : issun, isspn
     use control, only : norbs
     use control, only : ntime

     implicit none

! external arguments
! symmetry vector for correlated band
     integer, intent(in)     :: symm(norbs)

! green's function or weiss's function
     real(dp), intent(inout) :: gtau(ntime,norbs)

! local variables
! loop index
     integer :: i

! symmetrize the green's function over spin
! let spin up part = spin dn part for PM case
     if ( isspn == 1 ) then
         do i=1,ntime
             call hfqmc_symm_spin(norbs, gtau(i,1:norbs))
         enddo ! over i={1,ntime} loop
     endif ! back if ( isspn == 1 ) block

! symmetrize green's function over band
! deal with degenerate case, SU(N)
     if ( issun == 2 ) then
         do i=1,ntime
             call hfqmc_symm_band(norbs, symm, gtau(i,1:norbs))
         enddo ! over i={1,ntime} loop
     endif ! back if ( issun == 2 ) block

     return
  end subroutine hfqmc_make_symm

!!>>> hfqmc_make_smth: smooth impurity self-energy function
  subroutine hfqmc_make_smth(sigf)
     use constants, only : dp, czero

     use control, only : mfreq
     use control, only : ntime

     implicit none

! external arguments
! self-energy function to be smoothen
     complex(dp), intent(inout) :: sigf(mfreq)

! local variables
! loop index
     integer  :: i, j, k

! smooth radius
     integer  :: lim

! imaginary part of self-energy function
     real(dp) :: ti

! real part of self-energy function
     real(dp) :: tr

! dummy variables for summation
     complex(dp) :: ssum

! dummy self-energy function
     complex(dp) :: stmp(mfreq)

! determine smooth radius
     lim = ntime / 4

! deal with [1,lim]
     do k=1,lim
         stmp(k) = sigf(k)
     enddo ! over k={1,lim} loop

! deal with [lim+1,2*lim]
     do i=1,lim
         k = lim + i
         ssum = czero
         do j=-lim,lim
             ssum = ssum + sigf(k+j)
         enddo ! over j={-lim,lim} loop
         stmp(k) = ssum / real(2 * lim + 1)
         stmp(k) = ( ( lim - i ) * sigf(k) + i * stmp(k) ) / real(lim)
     enddo ! over i={1,lim} loop

! deal with [2*lim+1,mfreq-lim]
     do k=2*lim+1,mfreq-lim
         ssum = czero
         do j=-lim,lim
             ssum = ssum + sigf(k+j)
         enddo ! over j={-lim,lim} loop
         stmp(k) = ssum / real(2 * lim + 1)
     enddo ! over k={2*lim+1,mfreq-lim} loop

! deal with [mfreq-lim+1,mfreq]
     do k=mfreq-lim+1,mfreq
         tr = real( stmp(mfreq-lim) )
         ti = aimag( stmp(mfreq-lim) ) * real(mfreq - lim) / real(k)
         stmp(k) = dcmplx(tr,ti)
     enddo ! over k={mfreq-lim+1,mfreq} loop

! copy stmp to sigf
     do k=1,mfreq
         sigf(k) = stmp(k)
     enddo ! over k={1,mfreq} loop

     return
  end subroutine hfqmc_make_smth

!!>>> hfqmc_make_freq: calculate the final impurity green's function,
!!>>> bath weiss's function, and self-energy function in matsubara space
  subroutine hfqmc_make_freq()
     use constants, only : dp, zero, one

     use control, only : norbs
     use control, only : mfreq
     use control, only : myid, master
     use context, only : umat
     use context, only : rmesh
     use context, only : gtau, wtau
     use context, only : grnf, wssf, sig2

     implicit none

! local variables
! loop index
     integer  :: i
     integer  :: j
     integer  :: k

! 0-order expansion coefficients for self-energy function
     real(dp) :: u0(norbs)

! full Coulomb interaction matrix
     real(dp) :: uumat(norbs,norbs)

! fourier impurity green's function from G(\tau) to G(i\omega)
     call hfqmc_fourier_t2w(-gtau, grnf)

! fourier bath weiss's function from G0(\tau) to G0(i\omega)
     call hfqmc_fourier_t2w(-wtau, wssf)

! using the Dyson equation to calculate the self-energy function
     do i=1,norbs
         do j=1,mfreq
             sig2(j,i) = one / wssf(j,i) - one / grnf(j,i)
         enddo ! over j={1,mfreq} loop
     enddo ! over i={1,norbs} loop

! perform smoothness on self-energy function
     do i=1,norbs
         call hfqmc_make_smth( sig2(1:mfreq,i) )
     enddo ! over i={1,norbs} loop

! calculate full Coulomb interaction matrix
     k = 0
     uumat = zero
     do i=1,norbs-1
         do j=i+1,norbs
             k = k + 1
             uumat(i,j) = umat(k)
             uumat(j,i) = umat(k)
         enddo ! over j={i+1,norbs} loop
     enddo ! over i={1,norbs-1} loop

! calculate 0-order coefficient for self-energy function
     u0 = zero
     do i=1,norbs
         do j=1,norbs
             if ( i /= j ) then
                 u0(i) = u0(i) + uumat(j,i) * ( one - gtau(1,j) )
             endif ! back if ( i /= j ) block
         enddo ! over j={1,norbs} loop
     enddo ! over i={1,norbs} loop

! adjust the real part of self-energy function to fullfil the high
! frequency behavior
     do i=1,norbs
         do j=1,mfreq
             sig2(j,i) = sig2(j,i) - real( sig2(mfreq,i) ) + u0(i)
         enddo ! over j={1,mfreq} loop
     enddo ! over i={1,norbs} loop

! write out impurity green's function to disk file
     if ( myid == master ) then ! only master node can do it
         call hfqmc_dump_grnf(rmesh, grnf)
     endif ! back if ( myid == master ) block

! write out bath weiss's function to disk file
     if ( myid == master ) then ! only master node can do it
         call hfqmc_dump_wssf(rmesh, wssf)
     endif ! back if ( myid == master ) block

! write out self-energy function to disk file
     if ( myid == master ) then ! only master node can do it
         call hfqmc_dump_sigf(rmesh, sig2)
     endif ! back if ( myid == master ) block 

     return
  end subroutine hfqmc_make_freq

!>>> to calculate the quasiparticle weight
! original equation:
!
! Z = \frac{m}{m*}
!   = \frac{1}{1-\frac{\partial}{\partial\omega} Re \Sigma(\omega) |_{\omega=0}}
!
! in the context of QMC simulations, one usually approximates this
! quantity by its discrete Eliashberg estimate
!
! Z = \frac{1}{1-\frac{Im \Sigma(i\omega_1)}{\pi T}}
!
  subroutine hfqmc_make_quas()
     use constants
     use control
     use context

     implicit none

! local variables
! loop index
     integer :: i

! initialize it
     quas = zero

! calculate quasiparticle weight Z using Eliashberg estimate
     do i=1,norbs
         quas(i) = one / ( one - aimag( sig2(1,i) ) * ( beta / pi ) )
     enddo ! over i={1,norbs} loop

! write out quasiparticle weight to disk file
     if ( myid == master ) then ! only master node can do it
         call hfqmc_dump_quas(quas)
     endif

     return
  end subroutine hfqmc_make_quas

!>>> to calculate the occupation number
  subroutine hfqmc_make_nmat()
     use constants
     use control
     use context

     implicit none

! local variables
! loop index
     integer :: i

! initialize it
     nmat = zero

! calculate occupation number
     do i=1,norbs
         nmat(i) = one - gtau(1,i)
     enddo ! over i={1,norbs} loop

! write out occupation number to disk file
     if ( myid == master ) then ! only master node can do it
         call hfqmc_dump_nmat(nmat, nnmat)
     endif

     return
  end subroutine hfqmc_make_nmat

!>>> reduce the gtau from all children processes
  subroutine hfqmc_reduce_gtau(gtau_mpi)
     use constants
     use context

     use mmpi

     implicit none

! external arguments
! impurity green's function
     real(dp), intent(out) :: gtau_mpi(ntime,norbs)

! initialize gtau_mpi
     gtau_mpi = zero

! build gtau_mpi, collect data from all children processes
# if defined (MPI)

! collect data
     call mp_allreduce(gtau, gtau_mpi)

! block until all processes have reached here
     call mp_barrier()

# else  /* MPI */

     gtau_mpi = gtau

# endif /* MPI */

! calculate the average
     gtau_mpi = gtau_mpi / real(nprocs)

     return
  end subroutine hfqmc_reduce_gtau

!>>> reduce the nnmat from all children processes
  subroutine hfqmc_reduce_nmat(nnmat_mpi)
     use constants
     use context

     use mmpi

     implicit none

! external arguments
! double occupation number matrix
     real(dp), intent(out) :: nnmat_mpi(norbs,norbs)

! initialize nnmat_mpi
     nnmat_mpi = zero

! build nnmat_mpi, collect data from all children processes
# if defined (MPI)

! collect data
     call mp_allreduce(nnmat, nnmat_mpi)

! block until all processes have reached here
     call mp_barrier()

# else  /* MPI */

     nnmat_mpi = nnmat

# endif /* MPI */

! calculate the average
     nnmat_mpi = nnmat_mpi / real(nprocs)

     return
  end subroutine hfqmc_reduce_nmat
