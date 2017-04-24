!!!-----------------------------------------------------------------------
!!! project : narcissus
!!! program : ctqmc_impurity_solver
!!!           ctqmc_impurity_tester
!!! source  : ctqmc_solver.f90
!!! type    : subroutines
!!! author  : li huang (email:lihuang.dmft@gmail.com)
!!! history : 09/16/2009 by li huang (created)
!!!           04/24/2017 by li huang (last modified)
!!! purpose : the main subroutines for the hybridization expansion version
!!!           continuous time quantum Monte Carlo (CTQMC) quantum impurity
!!!           solver. they are the core engines
!!! status  : unstable
!!! comment :
!!!-----------------------------------------------------------------------

!!========================================================================
!!>>> core subroutines for quantum impurity solver                     <<<
!!========================================================================

!!
!! @sub ctqmc_impurity_solver
!!
!! core engine for hybridization expansion version continuous time quantum
!! Monte Carlo quantum impurity solver
!!
  subroutine ctqmc_impurity_solver(iter)
     use constants, only : dp, zero, one, mystd

     use control, only : cname               ! code name
     use control, only : isbnd, isspn        ! control symmetry
     use control, only : isobs, issus, isvrt ! control physical observables
     use control, only : nband, nspin, norbs ! size of model hamiltonian
     use control, only : ncfgs               ! size of hilbert space
     use control, only : mkink, mfreq        ! perturbation expansion order
     use control, only : nffrq, nbfrq, ntime ! matsubara frequency and time
     use control, only : nsweep, nwrite      ! monte carlo sampling
     use control, only : nmonte, ncarlo      ! interval for sampling
     use control, only : Uc, Jz              ! coulomb and hund's interaction
     use control, only : myid, master        ! mpi

     use context, only : tmesh, rmesh        ! frequency and time meshes
     use context, only : hist                ! histogram
     use context, only : prob                ! probability
     use context, only : nmat, nnmat         ! occupation
     use context, only : kmat, kkmat         ! kinetic energy fluctuation
     use context, only : lmat, rmat, lrmat   ! fidelity susceptibility
     use context, only : szpow               ! binder cumulant
     use context, only : schi, sschi, ssfom  ! spin susceptibility
     use context, only : ochi, oochi, oofom  ! charge susceptibility
     use context, only : g2_re, g2_im        ! two-particle green's function
     use context, only : h2_re, h2_im        ! two-particle green's function
     use context, only : ps_re, ps_im        ! pairing susceptibility
     use context, only : symm                ! symmetry
     use context, only : gtau, ftau          ! imaginary time green's function
     use context, only : grnf                ! matsubara green's function
     use context, only : sig2                ! self-energy function

     implicit none

! external arguments
! current self-consistent iteration number
     integer, intent(in) :: iter

! local variables
! loop index
     integer  :: i
     integer  :: j

! status flag
     integer  :: istat

! current QMC sweeping steps
     integer  :: cstep

! control flag, whether the solver is checked periodically
! cflag = 0  , do not check the quantum impurity solver
! cflag = 1  , check the quantum impurity solver periodically
! cflag = 99 , the quantum impurity solver is out of control
! cflag = 100, the quantum impurity solver has reached convergence
     integer  :: cflag

! starting time
     real(dp) :: time_begin

! ending time
     real(dp) :: time_end

! time consuming by current iteration
     real(dp) :: time_iter

! time consuming by total iteration
     real(dp) :: time_niter

! histogram for perturbation expansion series, for mpi case
     real(dp), allocatable :: hist_mpi(:)
     real(dp), allocatable :: hist_err(:)

! probability of atomic states, for mpi case
     real(dp), allocatable :: prob_mpi(:)
     real(dp), allocatable :: prob_err(:)

! impurity occupation number matrix, for mpi case
     real(dp), allocatable :: nmat_mpi(:)
     real(dp), allocatable :: nmat_err(:)

! impurity double occupation number matrix, for mpi case
     real(dp), allocatable :: nnmat_mpi(:,:)
     real(dp), allocatable :: nnmat_err(:,:)

! number of operators < k >, for mpi case
     real(dp), allocatable :: kmat_mpi(:)
     real(dp), allocatable :: kmat_err(:)

! square of number of operators < k^2 >, for mpi case
     real(dp), allocatable :: kkmat_mpi(:,:)
     real(dp), allocatable :: kkmat_err(:,:)

! number of operators at left half axis < k_l >, for mpi case
     real(dp), allocatable :: lmat_mpi(:)
     real(dp), allocatable :: lmat_err(:)

! number of operators at right half axis < k_r >, for mpi case
     real(dp), allocatable :: rmat_mpi(:)
     real(dp), allocatable :: rmat_err(:)

! used to evaluate fidelity susceptibility < k_l k_r >, for mpi case
     real(dp), allocatable :: lrmat_mpi(:,:)
     real(dp), allocatable :: lrmat_err(:,:)

! powers of local magnetization, for mpi case
     real(dp), allocatable :: szpow_mpi(:,:)
     real(dp), allocatable :: szpow_err(:,:)

! spin-spin correlation function, totally-averaged, for mpi case
     real(dp), allocatable :: schi_mpi(:)
     real(dp), allocatable :: schi_err(:)

! spin-spin correlation function, orbital-resolved, for mpi case
     real(dp), allocatable :: sschi_mpi(:,:)
     real(dp), allocatable :: sschi_err(:,:)

! spin-spin correlation function, orbital-resolved, for mpi case
     real(dp), allocatable :: ssfom_mpi(:,:)
     real(dp), allocatable :: ssfom_err(:,:)

! orbital-orbital correlation function, totally-averaged, for mpi case
     real(dp), allocatable :: ochi_mpi(:)
     real(dp), allocatable :: ochi_err(:)

! orbital-orbital correlation function, orbital-resolved, for mpi case
     real(dp), allocatable :: oochi_mpi(:,:,:)
     real(dp), allocatable :: oochi_err(:,:,:)

! orbital-orbital correlation function, orbital-resolved, for mpi case
     real(dp), allocatable :: oofom_mpi(:,:,:)
     real(dp), allocatable :: oofom_err(:,:,:)

! note: for the two-particle quantities, we don't measure the error bars
! used to measure two-particle green's function, real part, for mpi case
     real(dp), allocatable :: g2_re_mpi(:,:,:,:,:)

! used to measure two-particle green's function, imaginary part, for mpi case
     real(dp), allocatable :: g2_im_mpi(:,:,:,:,:)

! used to measure two-particle green's function, real part, for mpi case
     real(dp), allocatable :: h2_re_mpi(:,:,:,:,:)

! used to measure two-particle green's function, imaginary part, for mpi case
     real(dp), allocatable :: h2_im_mpi(:,:,:,:,:)

! used to measure particle-particle pair susceptibility, real part, for mpi case
     real(dp), allocatable :: ps_re_mpi(:,:,:,:,:)

! used to measure particle-particle pair susceptibility, imaginary part, for mpi case
     real(dp), allocatable :: ps_im_mpi(:,:,:,:,:)

! impurity green's function, imaginary time axis, for mpi case
     real(dp), allocatable :: gtau_mpi(:,:,:)
     real(dp), allocatable :: gtau_err(:,:,:)

! auxiliary correlation function, imaginary time axis, for mpi case
     real(dp), allocatable :: ftau_mpi(:,:,:)
     real(dp), allocatable :: ftau_err(:,:,:)

! impurity green's function, matsubara frequency axis, for mpi case
     complex(dp), allocatable :: grnf_mpi(:,:,:)
     complex(dp), allocatable :: grnf_err(:,:,:)

! allocate memory
     allocate(hist_mpi(mkink),             stat=istat)
     allocate(hist_err(mkink),             stat=istat)
     allocate(prob_mpi(ncfgs),             stat=istat)
     allocate(prob_err(ncfgs),             stat=istat)
     allocate(nmat_mpi(norbs),             stat=istat)
     allocate(nmat_err(norbs),             stat=istat)
     allocate(nnmat_mpi(norbs,norbs),      stat=istat)
     allocate(nnmat_err(norbs,norbs),      stat=istat)

     allocate(kmat_mpi(norbs),             stat=istat)
     allocate(kmat_err(norbs),             stat=istat)
     allocate(kkmat_mpi(norbs,norbs),      stat=istat)
     allocate(kkmat_err(norbs,norbs),      stat=istat)
     allocate(lmat_mpi(norbs),             stat=istat)
     allocate(lmat_err(norbs),             stat=istat)
     allocate(rmat_mpi(norbs),             stat=istat)
     allocate(rmat_err(norbs),             stat=istat)
     allocate(lrmat_mpi(norbs,norbs),      stat=istat)
     allocate(lrmat_err(norbs,norbs),      stat=istat)
     allocate(szpow_mpi(  4  ,norbs),      stat=istat)
     allocate(szpow_err(  4  ,norbs),      stat=istat)

     allocate(schi_mpi(ntime),             stat=istat)
     allocate(schi_err(ntime),             stat=istat)
     allocate(sschi_mpi(ntime,nband),      stat=istat)
     allocate(sschi_err(ntime,nband),      stat=istat)
     allocate(ssfom_mpi(nbfrq,nband),      stat=istat)
     allocate(ssfom_err(nbfrq,nband),      stat=istat)
     allocate(ochi_mpi(ntime),             stat=istat)
     allocate(ochi_err(ntime),             stat=istat)
     allocate(oochi_mpi(ntime,norbs,norbs),stat=istat)
     allocate(oochi_err(ntime,norbs,norbs),stat=istat)
     allocate(oofom_mpi(nbfrq,norbs,norbs),stat=istat)
     allocate(oofom_err(nbfrq,norbs,norbs),stat=istat)

     allocate(g2_re_mpi(nffrq,nffrq,nbfrq,norbs,norbs), stat=istat)
     allocate(g2_im_mpi(nffrq,nffrq,nbfrq,norbs,norbs), stat=istat)
     allocate(h2_re_mpi(nffrq,nffrq,nbfrq,norbs,norbs), stat=istat)
     allocate(h2_im_mpi(nffrq,nffrq,nbfrq,norbs,norbs), stat=istat)
     allocate(ps_re_mpi(nffrq,nffrq,nbfrq,norbs,norbs), stat=istat)
     allocate(ps_im_mpi(nffrq,nffrq,nbfrq,norbs,norbs), stat=istat)

     allocate(gtau_mpi(ntime,norbs,norbs), stat=istat)
     allocate(gtau_err(ntime,norbs,norbs), stat=istat)
     allocate(ftau_mpi(ntime,norbs,norbs), stat=istat)
     allocate(ftau_err(ntime,norbs,norbs), stat=istat)
     allocate(grnf_mpi(mfreq,norbs,norbs), stat=istat)
     allocate(grnf_err(mfreq,norbs,norbs), stat=istat)

     if ( istat /= 0 ) then
         call s_print_error('ctqmc_impurity_solver','can not allocate enough memory')
     endif ! back if ( istat /= 0 ) block

! setup cstep
     cstep = 0

! setup cflag, check the status of quantum impurity solver periodically
     cflag = 1

! setup timer
     time_iter = zero
     time_niter = zero

!!========================================================================
!!>>> starting quantum impurity solver                                 <<<
!!========================================================================

! print the header of continuous time quantum Monte Carlo quantum impurity solver
     if ( myid == master ) then ! only master node can do it
         write(mystd,'(2X,a)') cname//' >>> CTQMC quantum impurity solver running'
         write(mystd,'(4X,a,i10,4X,a,f10.5)') 'nband :', nband, 'Uc    :', Uc
         write(mystd,'(4X,a,i10,4X,a,f10.5)') 'nspin :', nspin, 'Jz    :', Jz
         write(mystd,*)
     endif ! back if ( myid == master ) block

!!========================================================================
!!>>> initializing quantum impurity solver                             <<<
!!========================================================================

! init the continuous time quantum Monte Carlo quantum impurity solver
! setup the key variables
     if ( myid == master ) then ! only master node can do it
         write(mystd,'(4X,a)') 'quantum impurity solver initializing'
     endif ! back if ( myid == master ) block

     call cpu_time(time_begin) ! record starting time
     call ctqmc_solver_init()
     call cpu_time(time_end)   ! record ending   time

! print the time information
     if ( myid == master ) then ! only master node can do it
         write(mystd,'(4X,a,f10.3,a)') 'time:', time_end - time_begin, 's'
         write(mystd,*)
     endif ! back if ( myid == master ) block

!!========================================================================
!!>>> retrieving quantum impurity solver                               <<<
!!========================================================================

! init the continuous time quantum Monte Carlo quantum impurity solver further
! retrieving the time series information produced by previous running
     if ( myid == master ) then ! only master node can do it
         write(mystd,'(4X,a)') 'quantum impurity solver retrieving'
     endif ! back if ( myid == master ) block

     call cpu_time(time_begin) ! record starting time
     call ctqmc_retrieve_status()
     call cpu_time(time_end)   ! record ending   time

! print the time information
     if ( myid == master ) then ! only master node can do it
         write(mystd,'(4X,a,f10.3,a)') 'time:', time_end - time_begin, 's'
         write(mystd,*)
     endif ! back if ( myid == master ) block

!!========================================================================
!!>>> warmming quantum impurity solver                                 <<<
!!========================================================================

! warmup the continuous time quantum Monte Carlo quantum impurity solver,
! in order to achieve equilibrium state
     if ( myid == master ) then ! only master node can do it
         write(mystd,'(4X,a)') 'quantum impurity solver warmming'
     endif ! back if ( myid == master ) block

     call cpu_time(time_begin) ! record starting time
     call ctqmc_warmup_diag()
     call cpu_time(time_end)   ! record ending   time

! print the time information
     if ( myid == master ) then ! only master node can do it
         write(mystd,'(4X,a,f10.3,a)') 'time:', time_end - time_begin, 's'
         write(mystd,*)
     endif ! back if ( myid == master ) block

!!========================================================================
!!>>> beginning main iteration                                         <<<
!!========================================================================

! start simulation
     if ( myid == master ) then ! only master node can do it
         write(mystd,'(4X,a)') 'quantum impurity solver sampling'
         write(mystd,*)
     endif ! back if ( myid == master ) block

     CTQMC_MAIN_ITERATION: do i=1, nsweep, nwrite

! record start time
         call cpu_time(time_begin)

         CTQMC_DUMP_ITERATION: do j=1, nwrite

!!========================================================================
!!>>> sampling perturbation expansion series                           <<<
!!========================================================================

! increase cstep by 1
             cstep = cstep + 1

! sampling the perturbation expansion feynman diagrams randomly
             call ctqmc_sample_diag(cstep)

!!========================================================================
!!>>> sampling the physical observables                                <<<
!!========================================================================

! record the histogram for perturbation expansion series
             call ctqmc_record_hist()

! record the impurity (double) occupation number matrix and other
! auxiliary physical observables
             if ( mod(cstep, nmonte) == 0 ) then
                 call ctqmc_record_nmat()
             endif ! back if ( mod(cstep, nmonte) == 0 ) block

! record the impurity green's function in matsubara frequency space
             if ( mod(cstep, nmonte) == 0 ) then
                 call ctqmc_record_grnf()
             endif ! back if ( mod(cstep, nmonte) == 0 ) block

! record the probability of eigenstates
             if ( mod(cstep, nmonte) == 0 ) then
                 call ctqmc_record_prob()
             endif ! back if ( mod(cstep, nmonte) == 0 ) block

! record the impurity green's function in imaginary time space
             if ( mod(cstep, ncarlo) == 0 ) then
                 call ctqmc_record_gtau()
             endif ! back if ( mod(cstep, ncarlo) == 0 ) block

! record the auxiliary correlation function, F(\tau)
             if ( mod(cstep, ncarlo) == 0 .and. isort >= 4 ) then
                 call ctqmc_record_ftau()
             endif ! back if ( mod(cstep, ncarlo) == 0 .and. isort >= 4 ) block

! record nothing
             if ( mod(cstep, nmonte) == 0 .and. btest(issus, 0) ) then
                 CONTINUE
             endif ! back if ( mod(cstep, nmonte) == 0 .and. btest(issus, 0) ) block

! record the spin-spin correlation function
             if ( mod(cstep, nmonte) == 0 .and. btest(issus, 1) ) then
                 call ctqmc_record_schi()
             endif ! back if ( mod(cstep, nmonte) == 0 .and. btest(issus, 1) ) block

! record the orbital-orbital correlation function
             if ( mod(cstep, nmonte) == 0 .and. btest(issus, 2) ) then
                 call ctqmc_record_ochi()
             endif ! back if ( mod(cstep, nmonte) == 0 .and. btest(issus, 2) ) block

! record the spin-spin correlation function
             if ( mod(cstep, nmonte) == 0 .and. btest(issus, 3) ) then
                 call ctqmc_record_sfom()
             endif ! back if ( mod(cstep, nmonte) == 0 .and. btest(issus, 3) ) block

! record the orbital-orbital correlation function
             if ( mod(cstep, nmonte) == 0 .and. btest(issus, 4) ) then
                 call ctqmc_record_ofom()
             endif ! back if ( mod(cstep, nmonte) == 0 .and. btest(issus, 4) ) block

! record the < k^2 > - < k >^2
             if ( mod(cstep, nmonte) == 0 .and. btest(issus, 5) ) then
                 call ctqmc_record_kmat()
             endif ! back if ( mod(cstep, nmonte) == 0 .and. btest(issus, 5) ) block

! record the fidelity susceptibility
             if ( mod(cstep, nmonte) == 0 .and. btest(issus, 6) ) then
                 call ctqmc_record_lmat()
             endif ! back if ( mod(cstep, nmonte) == 0 .and. btest(issus, 6) ) block

! record the powers of local magnetization
             if ( mod(cstep, nmonte) == 0 .and. btest(issus, 7) ) then
                 call ctqmc_record_szpw()
             endif ! back if ( mod(cstep, nmonte) == 0 .and. btest(issus, 7) ) block

! record nothing
             if ( mod(cstep, nmonte) == 0 .and. btest(isvrt, 0) ) then
                 CONTINUE
             endif ! back if ( mod(cstep, nmonte) == 0 .and. btest(isvrt, 0) ) block

! record the two-particle green's function
             if ( mod(cstep, nmonte) == 0 .and. btest(isvrt, 1) ) then
                 call ctqmc_record_twop()
             endif ! back if ( mod(cstep, nmonte) == 0 .and. btest(isvrt, 1) ) block

! record the two-particle green's function
             if ( mod(cstep, nmonte) == 0 .and. btest(isvrt, 2) ) then
                 call ctqmc_record_vrtx()
             endif ! back if ( mod(cstep, nmonte) == 0 .and. btest(isvrt, 2) ) block

! record the particle-particle pair susceptibility
             if ( mod(cstep, nmonte) == 0 .and. btest(isvrt, 3) ) then
                 call ctqmc_record_pair()
             endif ! back if ( mod(cstep, nmonte) == 0 .and. btest(isvrt, 3) ) block

         enddo CTQMC_DUMP_ITERATION ! over j={1,nwrite} loop

!!========================================================================
!!>>> reporting quantum impurity solver                                <<<
!!========================================================================

! it is time to write out the statistics results
         if ( myid == master ) then ! only master node can do it
             call ctqmc_print_runtime(iter, cstep)
         endif ! back if ( myid == master ) block

!!========================================================================
!!>>> reducing immediate results                                       <<<
!!========================================================================

! collect the histogram data from hist to hist_mpi
         call ctqmc_reduce_hist(hist_mpi, hist_err)

! collect the impurity green's function data from gtau to gtau_mpi
         call ctqmc_reduce_gtau(gtau_mpi, gtau_err)

! gtau_mpi need to be scaled properly before written
         gtau_mpi = gtau_mpi * real(ncarlo) / real(cstep)
         gtau_err = gtau_err * real(ncarlo) / real(cstep)

!!========================================================================
!!>>> symmetrizing immediate results                                   <<<
!!========================================================================

! symmetrize the impurity green's function over spin or over bands
         if ( issun == 2 .or. isspn == 1 ) then
             call ctqmc_symm_gtau(symm, gtau_mpi)
             call ctqmc_symm_gtau(symm, gtau_err)
         endif ! back if ( issun == 2 .or. isspn == 1 ) block

!!========================================================================
!!>>> writing immediate results                                        <<<
!!========================================================================

! write out the histogram data, hist_mpi
         if ( myid == master ) then ! only master node can do it
             call ctqmc_dump_hist(hist_mpi, hist_err)
         endif ! back if ( myid == master ) block

! write out the impurity green's function, gtau_mpi
         if ( myid == master ) then ! only master node can do it
             call ctqmc_dump_gtau(tmesh, gtau_mpi, gtau_err)
         endif ! back if ( myid == master ) block

!!========================================================================
!!>>> checking quantum impurity solver                                 <<<
!!========================================================================

! check the status at first
         call ctqmc_verify_diag(cflag)

! write out the snapshot for the current diagram configuration
         if ( myid == master ) then
             call ctqmc_dump_diag(iter, cstep)
         endif ! back if ( myid == master ) block

!!========================================================================
!!>>> timing quantum impurity solver                                   <<<
!!========================================================================

! record ending time for this iteration
         call cpu_time(time_end)

! calculate timing information
         time_iter = time_end - time_begin
         time_niter = time_niter + time_iter
         time_begin = time_end

! print out the result
         if ( myid == master ) then ! only master node can do it
             call s_time_analyzer(time_iter, time_niter)
             write(mystd,*)
         endif ! back if ( myid == master ) block

!!========================================================================
!!>>> escaping quantum impurity solver                                 <<<
!!========================================================================

! if the quantum impurity solver is out of control or reaches convergence
         if ( cflag == 99 .or. cflag == 100 ) then
             EXIT CTQMC_MAIN_ITERATION ! jump out the iteration
         endif ! back if ( cflag == 99 .or. cflag == 100 ) block

     enddo CTQMC_MAIN_ITERATION ! over i={1,nsweep} loop

!!========================================================================
!!>>> ending main iteration                                            <<<
!!========================================================================

!!========================================================================
!!>>> reducing final results                                           <<<
!!========================================================================

! collect the histogram data from hist to hist_mpi
     call ctqmc_reduce_hist(hist_mpi, hist_err)

! collect the probability data from prob to prob_mpi
     call ctqmc_reduce_prob(prob_mpi, prob_err)

! collect the occupation matrix data from nmat to nmat_mpi
! collect the double occupation matrix data from nnmat to nnmat_mpi
     call ctqmc_reduce_nmat(nmat_mpi, nnmat_mpi, nmat_err, nnmat_err)

! collect the < k^2 > - < k >^2 data from kmat to kmat_mpi
! collect the < k^2 > - < k >^2 data from kkmat to kkmat_mpi
     call ctqmc_reduce_kmat(kmat_mpi, kkmat_mpi, kmat_err, kkmat_err)

! collect the fidelity susceptibility data from lmat to lmat_mpi
! collect the fidelity susceptibility data from rmat to rmat_mpi
! collect the fidelity susceptibility data from lrmat to lrmat_mpi
     call ctqmc_reduce_lmat(lmat_mpi, rmat_mpi, lrmat_mpi, lmat_err, rmat_err, lrmat_err)

! collect the powers of local magnetization data from szpow to szpow_mpi
     call ctqmc_reduce_szpw(szpow_mpi, szpow_err)

! collect the spin-spin correlation function data from schi to schi_mpi
! collect the spin-spin correlation function data from sschi to sschi_mpi
! collect the spin-spin correlation function data from ssfom to ssfom_mpi
     call ctqmc_reduce_schi(schi_mpi, sschi_mpi, schi_err, sschi_err)
     call ctqmc_reduce_sfom(ssfom_mpi, ssfom_err)

! collect the orbital-orbital correlation function data from ochi to ochi_mpi
! collect the orbital-orbital correlation function data from oochi to oochi_mpi
! collect the orbital-orbital correlation function data from oofom to oofom_mpi
     call ctqmc_reduce_ochi(ochi_mpi, oochi_mpi, ochi_err, oochi_err)
     call ctqmc_reduce_ofom(oofom_mpi, oofom_err)

! collect the two-particle green's function from g2_re to g2_re_mpi
! collect the two-particle green's function from g2_im to g2_im_mpi
     call ctqmc_reduce_twop(g2_re_mpi, g2_im_mpi)

! collect the two-particle green's function from h2_re to h2_re_mpi
! collect the two-particle green's function from h2_im to h2_im_mpi
     call ctqmc_reduce_vrtx(h2_re_mpi, h2_im_mpi)

! collect the pair susceptibility from ps_re to ps_re_mpi
! collect the pair susceptibility from ps_im to ps_im_mpi
     call ctqmc_reduce_pair(ps_re_mpi, ps_im_mpi)

! collect the impurity green's function data from gtau to gtau_mpi
     call ctqmc_reduce_gtau(gtau_mpi, gtau_err)

! collect the auxiliary correlation function from ftau to ftau_mpi
     call ctqmc_reduce_ftau(ftau_mpi, ftau_err)

! collect the impurity green's function data from grnf to grnf_mpi
     call ctqmc_reduce_grnf(grnf_mpi, grnf_err)

! update original data and calculate the averages simultaneously
! average value section
     hist  = hist_mpi  * one
     prob  = prob_mpi  * real(ncarlo) / real(nsweep)

     nmat  = nmat_mpi  * real(nmonte) / real(nsweep)
     nnmat = nnmat_mpi * real(nmonte) / real(nsweep)
     kmat  = kmat_mpi  * real(nmonte) / real(nsweep)
     kkmat = kkmat_mpi * real(nmonte) / real(nsweep)
     lmat  = lmat_mpi  * real(nmonte) / real(nsweep)
     rmat  = rmat_mpi  * real(nmonte) / real(nsweep)
     lrmat = lrmat_mpi * real(nmonte) / real(nsweep)
     szpow = szpow_mpi * real(nmonte) / real(nsweep)
     schi  = schi_mpi  * real(nmonte) / real(nsweep)
     sschi = sschi_mpi * real(nmonte) / real(nsweep)
     ssfom = ssfom_mpi * real(nmonte) / real(nsweep)
     ochi  = ochi_mpi  * real(nmonte) / real(nsweep)
     oochi = oochi_mpi * real(nmonte) / real(nsweep)
     oofom = oofom_mpi * real(nmonte) / real(nsweep)

     g2_re = g2_re_mpi * real(nmonte) / real(nsweep)
     g2_im = g2_im_mpi * real(nmonte) / real(nsweep)
     h2_re = h2_re_mpi * real(nmonte) / real(nsweep)
     h2_im = h2_im_mpi * real(nmonte) / real(nsweep)
     ps_re = ps_re_mpi * real(nmonte) / real(nsweep)
     ps_im = ps_im_mpi * real(nmonte) / real(nsweep)

     gtau  = gtau_mpi  * real(ncarlo) / real(nsweep)
     ftau  = ftau_mpi  * real(ncarlo) / real(nsweep)
     grnf  = grnf_mpi  * real(nmonte) / real(nsweep)

! update original data and calculate the averages simultaneously
! error bar section
     hist_err  = hist_err  * one
     prob_err  = prob_err  * real(ncarlo) / real(nsweep)

     nmat_err  = nmat_err  * real(nmonte) / real(nsweep)
     nnmat_err = nnmat_err * real(nmonte) / real(nsweep)
     kmat_err  = kmat_err  * real(nmonte) / real(nsweep)
     kkmat_err = kkmat_err * real(nmonte) / real(nsweep)
     lmat_err  = lmat_err  * real(nmonte) / real(nsweep)
     rmat_err  = rmat_err  * real(nmonte) / real(nsweep)
     lrmat_err = lrmat_err * real(nmonte) / real(nsweep)
     szpow_err = szpow_err * real(nmonte) / real(nsweep)
     schi_err  = schi_err  * real(nmonte) / real(nsweep)
     sschi_err = sschi_err * real(nmonte) / real(nsweep)
     ssfom_err = ssfom_err * real(nmonte) / real(nsweep)
     ochi_err  = ochi_err  * real(nmonte) / real(nsweep)
     oochi_err = oochi_err * real(nmonte) / real(nsweep)
     oofom_err = oofom_err * real(nmonte) / real(nsweep)

     gtau_err  = gtau_err  * real(ncarlo) / real(nsweep)
     ftau_err  = ftau_err  * real(ncarlo) / real(nsweep)
     grnf_err  = grnf_err  * real(nmonte) / real(nsweep)

! build atomic green's function and self-energy function using improved
! Hubbard-I approximation, and then make interpolation for self-energy
! function between low frequency QMC data and high frequency Hubbard-I
! approximation data, the impurity green's function can be obtained by
! using dyson's equation finally
     if ( isort <= 3 ) then
         call ctqmc_make_hub1()
! build atomic green's function and self-energy function using improved
! Hubbard-I approximation, and then make forward fourier transformation
! for impurity green's function and auxiliary correlation function. then
! the final self-energy function is obtained by analytical formula
     else
         call ctqmc_make_hub2()
     endif ! back if ( isort <= 3 ) block

!!========================================================================
!!>>> symmetrizing final results                                       <<<
!!========================================================================

! symmetrize the occupation number matrix (nmat) over spin or over bands
     if ( issun == 2 .or. isspn == 1 ) then
         call ctqmc_symm_nmat(symm, nmat)
         call ctqmc_symm_nmat(symm, nmat_err)
     endif ! back if ( issun == 2 .or. isspn == 1 ) block

! symmetrize the impurity green's function (gtau) over spin or over bands
     if ( issun == 2 .or. isspn == 1 ) then
         call ctqmc_symm_gtau(symm, gtau)
         call ctqmc_symm_gtau(symm, gtau_err)
     endif ! back if ( issun == 2 .or. isspn == 1 ) block

! symmetrize the impurity green's function (grnf) over spin or over bands
     if ( issun == 2 .or. isspn == 1 ) then
         call ctqmc_symm_grnf(symm, grnf)
         call ctqmc_symm_grnf(symm, grnf_err)
     endif ! back if ( issun == 2 .or. isspn == 1 ) block

! symmetrize the impurity self-energy function (sig2) over spin or over bands
     if ( issun == 2 .or. isspn == 1 ) then
         call ctqmc_symm_grnf(symm, sig2)
     endif ! back if ( issun == 2 .or. isspn == 1 ) block

!!========================================================================
!!>>> writing final results                                            <<<
!!========================================================================

! write out the final histogram data, hist
     if ( myid == master ) then ! only master node can do it
         call ctqmc_dump_hist(hist, hist_err)
     endif ! back if ( myid == master ) block

! write out the final probability data, prob
     if ( myid == master ) then ! only master node can do it
         call ctqmc_dump_prob(prob, prob_err)
     endif ! back if ( myid == master ) block

! write out the final (double) occupation matrix data, nmat and nnmat
     if ( myid == master ) then ! only master node can do it
         call ctqmc_dump_nmat(nmat, nnmat, nmat_err, nnmat_err)
     endif ! back if ( myid == master ) block

! write out the final < k^2 > - < k >^2 data, kmat and kkmat
     if ( myid == master ) then ! only master node can do it
         call ctqmc_dump_kmat(kmat, kkmat, kmat_err, kkmat_err)
     endif ! back if ( myid == master ) block

! write out the final fidelity susceptibility data, lmat, rmat, and lrmat
     if ( myid == master ) then ! only master node can do it
         call ctqmc_dump_lmat(lmat, rmat, lrmat, lmat_err, rmat_err, lrmat_err)
     endif ! back if ( myid == master ) block

! write out the final powers of local magnetization data, szpow
     if ( myid == master ) then ! only master node can do it
         call ctqmc_dump_szpw(szpow, szpow_err)
     endif ! back if ( myid == master ) block

! write out the final spin-spin correlation function data, schi, sschi, and ssfom
     if ( myid == master ) then ! only master node can do it
         call ctqmc_dump_schi(schi, sschi, schi_err, sschi_err)
         call ctqmc_dump_sfom(ssfom, ssfom_err)
     endif ! back if ( myid == master ) block

! write out the final orbital-orbital correlation function data, ochi, oochi, and oofom
     if ( myid == master ) then ! only master node can do it
         call ctqmc_dump_ochi(ochi, oochi, ochi_err, oochi_err)
         call ctqmc_dump_ofom(oofom, oofom_err)
     endif ! back if ( myid == master ) block

! write out the final two-particle green's function data, g2_re and g2_im
     if ( myid == master ) then ! only master node can do it
         call ctqmc_dump_twop(g2_re, g2_im)
     endif ! back if ( myid == master ) block

! write out the final two-particle green's function data, h2_re and h2_im
     if ( myid == master ) then ! only master node can do it
         call ctqmc_dump_vrtx(h2_re, h2_im)
     endif ! back if ( myid == master ) block

! write out the final particle-particle pair susceptibility data, ps_re and ps_im
     if ( myid == master ) then ! only master node can do it
         call ctqmc_dump_pair(ps_re, ps_im)
     endif ! back if ( myid == master ) block

! write out the final impurity green's function data, gtau
     if ( myid == master ) then ! only master node can do it
         call ctqmc_dump_gtau(tmesh, gtau, gtau_err)
     endif ! back if ( myid == master ) block

! write out the final impurity green's function data, grnf
     if ( myid == master ) then ! only master node can do it
         call ctqmc_dump_grnf(rmesh, grnf, grnf_err)
     endif ! back if ( myid == master ) block

! write out the final self-energy function data, sig2
     if ( myid == master ) then ! only master node can do it
         call ctqmc_dump_sigf(rmesh, sig2)
     endif ! back if ( myid == master ) block

!!========================================================================
!!>>> saving quantum impurity solver                                   <<<
!!========================================================================

! save the perturbation expansion series information to the disk file
     if ( myid == master ) then ! only master node can do it
         call ctqmc_save_status()
     endif ! back if ( myid == master ) block

!!========================================================================
!!>>> finishing quantum impurity solver                                <<<
!!========================================================================

! print the footer of continuous time quantum Monte Carlo quantum impurity solver
     if ( myid == master ) then ! only master node can do it
         write(mystd,'(2X,a)') cname//' >>> CTQMC quantum impurity solver shutdown'
         write(mystd,*)
     endif ! back if ( myid == master ) block

! deallocate memory
     deallocate(hist_mpi )
     deallocate(hist_err )
     deallocate(prob_mpi )
     deallocate(prob_err )
     deallocate(nmat_mpi )
     deallocate(nmat_err )
     deallocate(nnmat_mpi)
     deallocate(nnmat_err)
     deallocate(kmat_mpi )
     deallocate(kmat_err )
     deallocate(kkmat_mpi)
     deallocate(kkmat_err)
     deallocate(lmat_mpi )
     deallocate(lmat_err )
     deallocate(rmat_mpi )
     deallocate(rmat_err )
     deallocate(lrmat_mpi)
     deallocate(lrmat_err)
     deallocate(szpow_mpi)
     deallocate(szpow_err)
     deallocate(schi_mpi )
     deallocate(schi_err )
     deallocate(sschi_mpi)
     deallocate(sschi_err)
     deallocate(ssfom_mpi)
     deallocate(ssfom_err)
     deallocate(ochi_mpi )
     deallocate(ochi_err )
     deallocate(oochi_mpi)
     deallocate(oochi_err)
     deallocate(oofom_mpi)
     deallocate(oofom_err)
     deallocate(g2_re_mpi)
     deallocate(g2_im_mpi)
     deallocate(h2_re_mpi)
     deallocate(h2_im_mpi)
     deallocate(ps_re_mpi)
     deallocate(ps_im_mpi)
     deallocate(gtau_mpi )
     deallocate(gtau_err )
     deallocate(ftau_mpi )
     deallocate(ftau_err )
     deallocate(grnf_mpi )
     deallocate(grnf_err )

     return
  end subroutine ctqmc_impurity_solver

!!========================================================================
!!>>> debug layer                                                      <<<
!!========================================================================

!!
!! @sub ctqmc_impurity_tester
!!
!! testing subroutine, please try to active it on the ctqmc_sample_diag()
!! subroutine
!!
  subroutine ctqmc_impurity_tester()
     use constants ! ALL

     use control   ! ALL
     use context   ! ALL

     implicit none

!-------------------------------------------------------------------------
! please insert your debug code here
!-------------------------------------------------------------------------

     call cat_disp_segments(2)
     call s_print_error('ctqmc_impurity_tester','in debug mode')

     return
  end subroutine ctqmc_impurity_tester
