!!!-----------------------------------------------------------------------
!!! project : gardenia
!!! program : ctqmc_impurity_solver
!!!           ctqmc_impurity_tester
!!!           ctqmc_warmup_diag
!!!           ctqmc_sample_diag
!!!           ctqmc_temper_diag
!!!           ctqmc_verify_diag
!!!           ctqmc_insert_kink
!!!           ctqmc_remove_kink
!!!           ctqmc_lshift_kink
!!!           ctqmc_rshift_kink
!!!           ctqmc_reswap_kink
!!!           ctqmc_reflip_kink
!!!           ctqmc_reload_kink
!!! source  : ctqmc_solver.f90
!!! type    : subroutines
!!! author  : li huang (email:lihuang.dmft@gmail.com)
!!! history : 09/16/2009 by li huang (created)
!!!           03/28/2017 by li huang (last modified)
!!! purpose : the main subroutine for the hybridization expansion version
!!!           continuous time quantum Monte Carlo (CTQMC) quantum impurity
!!!           solver
!!! status  : unstable
!!! comment :
!!!-----------------------------------------------------------------------

!!========================================================================
!!>>> core subroutines for quantum impurity solver                     <<<
!!========================================================================

!!>>> ctqmc_impurity_solver: core engine for hybridization expansion version
!!>>> continuous time quantum Monte Carlo quantum impurity solver
  subroutine ctqmc_impurity_solver(iter)
     use constants, only : dp, zero, one, mystd

     use control, only : cname
     use control, only : issun, isspn, isort, issus, isvrt
     use control, only : nband, nspin, norbs, ncfgs
     use control, only : mkink, mfreq
     use control, only : nffrq, nbfrq, ntime, nsweep, nwrite, nmonte, ncarlo
     use control, only : Uc, Jz
     use control, only : beta
     use control, only : myid, master
     use context, only : tmesh, rmesh
     use context, only : hist, prob
     use context, only : nmat, nnmat
     use context, only : kmat, kkmat
     use context, only : lmat, rmat, lrmat
     use context, only : szpow
     use context, only : schi, sschi, ssfom, ochi, oochi, oofom
     use context, only : g2_re, g2_im, h2_re, h2_im, ps_re, ps_im
     use context, only : symm
     use context, only : gtau, ftau, grnf
     use context, only : sig2

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

! setup nsweep
! whether it is time to enter QMC data accumulating mode
     if ( iter == 999 ) then
         nsweep = nsweep * 10
         nwrite = nwrite * 10
     endif ! back if ( iter == 999 ) block

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
! ctqmc_sample_diag() is suitable for low temperature region, while
! ctqmc_temper_diag() is suitable for extreme high temperature region
             if ( beta > one ) then
                 call ctqmc_sample_diag(cstep)
             else
                 call ctqmc_temper_diag(cstep)
             endif ! back if ( beta > one ) block

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
             if ( mod(cstep, ncarlo) == 0 ) then
                 call ctqmc_record_prob()
             endif ! back if ( mod(cstep, ncarlo) == 0 ) block

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

!!>>> ctqmc_impurity_tester: testing subroutine, please try to active it
!!>>> on ctqmc_sample_diag() subroutine
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

!!========================================================================
!!>>> service layer: updating perturbation expansion series 1          <<<
!!========================================================================

!!>>> ctqmc_warmup_diag: perform thermalization or warmup on the
!!>>> perturbation expansion series to achieve thermodynamics stable
!!>>> equilibrium state
  subroutine ctqmc_warmup_diag()
     use constants, only : zero

     use control, only : ntherm
     use context, only : insert_tcount, insert_accept, insert_reject
     use context, only : remove_tcount, remove_accept, remove_reject
     use context, only : lshift_tcount, lshift_accept, lshift_reject
     use context, only : rshift_tcount, rshift_accept, rshift_reject
     use context, only : reswap_tcount, reswap_accept, reswap_reject
     use context, only : reflip_tcount, reflip_accept, reflip_reject

     implicit none

! local variables
! loop index
     integer :: i

! warm up the diagram series
     do i=1,ntherm
         call ctqmc_sample_diag(i)
     enddo ! over i={1,ntherm} loop

! reinit statistics variables
     insert_tcount = zero
     insert_accept = zero
     insert_reject = zero

     remove_tcount = zero
     remove_accept = zero
     remove_reject = zero

     lshift_tcount = zero
     lshift_accept = zero
     lshift_reject = zero

     rshift_tcount = zero
     rshift_accept = zero
     rshift_reject = zero

     reswap_tcount = zero
     reswap_accept = zero
     reswap_reject = zero

     reflip_tcount = zero
     reflip_accept = zero
     reflip_reject = zero

     return
  end subroutine ctqmc_warmup_diag

!!>>> ctqmc_sample_diag: visit the perturbation expansion diagrams randomly
  subroutine ctqmc_sample_diag(cstep)
     use constants, only : dp
     use spring, only : spring_sfmt_stream

     use control, only : nflip, nclean

     implicit none

! external arguments
! current QMC sweep steps
     integer, intent(in) :: cstep

! change the order of perturbation expansion series
     if ( spring_sfmt_stream() < 0.9_dp ) then
         if ( spring_sfmt_stream() > 0.5_dp ) then
             call ctqmc_insert_kink()  ! insert one new kink
         else
             call ctqmc_remove_kink()  ! remove one old kink
         endif ! back if ( spring_sfmt_stream() > 0.5_dp ) block
! do not change the order of perturbation expansion series
     else
         if ( spring_sfmt_stream() > 0.5_dp ) then
             call ctqmc_lshift_kink()  ! shift the left  endpoints
         else
             call ctqmc_rshift_kink()  ! shift the right endpoints
         endif ! back if ( spring_sfmt_stream() > 0.5_dp ) block
     endif ! back if ( spring_sfmt_stream() < 0.9_dp ) block

! numerical trick: perform global spin flip periodically
     if ( nflip > 0  .and. mod(cstep, +nflip) == 0 ) then
         if ( spring_sfmt_stream() < 0.8_dp ) then
             call ctqmc_reflip_kink(2) ! flip intra-orbital spins one by one
         else
             call ctqmc_reflip_kink(3) ! flip intra-orbital spins globally
         endif ! back if ( spring_sfmt_stream() < 0.8_dp ) block
     endif ! back if ( nflip > 0  .and. mod(cstep, +nflip) == 0 ) block

     if ( nflip < 0  .and. mod(cstep, -nflip) == 0 ) then
         if ( spring_sfmt_stream() < 0.8_dp ) then
             call ctqmc_reflip_kink(1) ! flip inter-orbital spins randomly
         else
             call ctqmc_reflip_kink(3) ! flip intra-orbital spins globally
         endif ! back if ( spring_sfmt_stream() < 0.8_dp ) block
     endif ! back if ( nflip < 0  .and. mod(cstep, -nflip) == 0 ) block

! numerical trick: perform global update periodically
     if ( nclean > 0 .and. mod(cstep, nclean) == 0 ) then
         call ctqmc_reload_kink()
     endif ! back if ( nclean > 0 .and. mod(cstep, nclean) == 0 ) block

     return
  end subroutine ctqmc_sample_diag

!!>>> ctqmc_temper_diag: visit the perturbation expansion diagrams randomly
!!>>> at very high temperature
  subroutine ctqmc_temper_diag(cstep)
     use constants, only : dp
     use spring, only : spring_sfmt_stream

     use control, only : nflip, nclean

     implicit none

! external arguments
! current QMC sweep steps
     integer, intent(in) :: cstep

! change the order of perturbation expansion series
     if ( spring_sfmt_stream() < 0.1_dp ) then
         if ( spring_sfmt_stream() > 0.5_dp ) then
             call ctqmc_insert_kink()  ! insert one new kink
         else
             call ctqmc_remove_kink()  ! remove one old kink
         endif ! back if ( spring_sfmt_stream() > 0.5_dp ) block
! do not change the order of perturbation expansion series
     else
         if ( spring_sfmt_stream() > 0.5_dp ) then
             call ctqmc_lshift_kink()  ! shift the left  endpoints
         else
             call ctqmc_reswap_kink()  ! swap creator and destroyer
         endif ! back if ( spring_sfmt_stream() > 0.5_dp ) block
     endif ! back if ( spring_sfmt_stream() < 0.1_dp ) block

! numerical trick: perform global spin flip periodically
     if ( nflip > 0  .and. mod(cstep, +nflip) == 0 ) then
         if ( spring_sfmt_stream() < 0.8_dp ) then
             call ctqmc_reflip_kink(2) ! flip intra-orbital spins one by one
         else
             call ctqmc_reflip_kink(3) ! flip intra-orbital spins globally
         endif ! back if ( spring_sfmt_stream() < 0.8_dp ) block
     endif ! back if ( nflip > 0  .and. mod(cstep, +nflip) == 0 ) block

     if ( nflip < 0  .and. mod(cstep, -nflip) == 0 ) then
         if ( spring_sfmt_stream() < 0.8_dp ) then
             call ctqmc_reflip_kink(1) ! flip inter-orbital spins randomly
         else
             call ctqmc_reflip_kink(3) ! flip intra-orbital spins globally
         endif ! back if ( spring_sfmt_stream() < 0.8_dp ) block
     endif ! back if ( nflip < 0  .and. mod(cstep, -nflip) == 0 ) block

! numerical trick: perform global update periodically
     if ( nclean > 0 .and. mod(cstep, nclean) == 0 ) then
         call ctqmc_reload_kink()
     endif ! back if ( nclean > 0 .and. mod(cstep, nclean) == 0 ) block

     return
  end subroutine ctqmc_temper_diag

!!>>> ctqmc_verify_diag: checking whether the quantum impurity solver is
!!>>> consistent internally
  subroutine ctqmc_verify_diag(cflag)
     use constants, only : mystd

     use control, only : norbs
     use control, only : myid, master
     use context, only : index_s, index_e, time_s, time_e
     use context, only : rank, stts

     implicit none

! external arguments
! control flag
     integer, intent(inout) :: cflag

! local variables
! loop index
     integer :: i
     integer :: j

     if ( cflag == 1 ) then

! check perturbation expansion order
         do i=1,norbs
             if ( stts(i) == 0 ) then
                 if ( rank(i) /= 0 ) cflag = 99
             endif ! back if ( stts(i) == 0 ) block

             if ( stts(i) == 1 ) then
                 if ( rank(i) == 0 ) cflag = 99
             endif ! back if ( stts(i) == 1 ) block

             if ( stts(i) == 2 ) then
                 if ( rank(i) == 0 ) cflag = 99
             endif ! back if ( stts(i) == 2 ) block

             if ( stts(i) == 3 ) then
                 if ( rank(i) /= 0 ) cflag = 99
             endif ! back if ( stts(i) == 3 ) block
         enddo ! over i={1,norbs} loop

! check time order of operators
         do i=1,norbs
             do j=1,rank(i)-1
                 if ( time_s( index_s(j, i), i ) > time_s( index_s(j+1, i), i ) ) then
                     cflag = 99
                 endif ! back if ( time_s( index_s(j, i), i ) > time_s( index_s(j+1, i), i ) ) block
                 if ( time_e( index_e(j, i), i ) > time_e( index_e(j+1, i), i ) ) then
                     cflag = 99
                 endif ! back if ( time_e( index_e(j, i), i ) > time_e( index_e(j+1, i), i ) ) block
             enddo ! over j={1,rank(i)-1} loop
         enddo ! over i={1,norbs} loop

! check segment and anti-segment
         do i=1,norbs
             if ( stts(i) == 1 ) then
                 if ( time_s( index_s(1, i), i ) > time_e( index_e(1, i), i ) ) then
                     cflag = 99
                 endif ! back if ( time_s( index_s(1, i), i ) > time_e( index_e(1, i), i ) ) block
             endif ! back if ( stts(i) == 1 ) block

             if ( stts(i) == 2 ) then
                 if ( time_s( index_s(1, i), i ) < time_e( index_e(1, i), i ) ) then
                     cflag = 99
                 endif ! back if ( time_s( index_s(1, i), i ) < time_e( index_e(1, i), i ) ) block
             endif ! back if ( stts(i) == 2 ) block
         enddo ! over i={1,norbs} loop

! write the results, only master node can do it
         if ( myid == master ) then
             if ( cflag == 99 ) then
                 write(mystd,'(4X,a)') '>>> quantum impurity solver status: error?'
                 write(mystd,'(4X,a)') '>>> please check the status file: solver.status.dat'
                 call ctqmc_save_status()
                 call s_print_error('ctqmc_verify_diag','unknown fatal error occur')
             else
                 write(mystd,'(4X,a)') '>>> quantum impurity solver status: normal'
             endif ! back if ( cflag == 99 ) block
         endif ! back if ( myid == master ) block

     endif ! back if ( cflag == 1 ) block

     return
  end subroutine ctqmc_verify_diag

!!========================================================================
!!>>> service layer: updating perturbation expansion series 2          <<<
!!========================================================================

!!>>> ctqmc_insert_kink: insert new segment or anti-segment in the
!!>>> perturbation expansion series
  subroutine ctqmc_insert_kink()
     use constants, only : dp, zero, one
     use spring, only : spring_sfmt_stream

     use control, only : norbs
     use control, only : mkink
     use control, only : beta
     use context, only : ckink, cstat
     use context, only : insert_tcount, insert_accept, insert_reject
     use context, only : rank, stts

     implicit none

! local variables
! whether the new segment or anti-segment can be inserted diagrammatically
     logical  :: ladd

! whether it is an anti-segment
! anti = .true., anti-segment
! anti = .false., segment
     logical  :: anti

! whether the update operation is accepted
     logical  :: pass

! current flavor channel for both band and spin
     integer  :: flvr

! index address to insert new segment or anti-segment
! is and ie are for start and end points, respectively
     integer  :: is, ie

! transition probability for insert new segment or anti-segment
     real(dp) :: p

! \tau_s, start point of the new segment
     real(dp) :: tau_start

! \tau_e, end   point of the new segment
     real(dp) :: tau_end

! the possible maximum length of the new segment
     real(dp) :: tau_max

! ratio between old and new configurations, the local trace part
     real(dp) :: trace_ratio

! ratio between old and new configurations, the determinant part
     real(dp) :: deter_ratio

! initialize logical variables
     ladd = .false.
     anti = .false.
     pass = .false.

! select the flavor channel randomly among 1 ~ norbs
     flvr = ceiling( spring_sfmt_stream() * norbs )

! get the perturbation expansion order ( number of existing segments or
! anti-segments ) for current flavor channel
     ckink = rank(flvr)
     if ( ckink == mkink ) then
!<         call s_print_exception('ctqmc_insert_kink','can not insert any segments')
         insert_tcount = insert_tcount + one
         insert_reject = insert_reject + one
         RETURN
     endif ! back if ( ckink == mkink ) block

! randomly choose anti and tau_start, and then check whether tau_start is
! valid, if tau_start is valid, then determine tau_end, tau_max, is, and
! ie consistently, and set ladd to .true., if tau_start is not valid, then
! set ladd to .false.
     call cat_insert_flavor(flvr, is, ie, anti, ladd, tau_start, tau_end, tau_max)

! calculate the transition ratio between old and new configurations,
! for the local trace part
     if ( ladd .eqv. .true. ) then
         call cat_insert_ztrace(flvr, anti, tau_start, tau_end, trace_ratio)
     else
         trace_ratio = zero
     endif ! back if ( ladd .eqv. .true. ) block

! calculate the transition ratio between old and new configurations,
! for the determinant part
     if ( ladd .eqv. .true. ) then
         call cat_insert_detrat(flvr, tau_start, tau_end, deter_ratio)
     else
         deter_ratio = zero
     endif ! back if ( ladd .eqv. .true. ) block

! calculate the transition probability for insert new segment or anti-segment
     p = deter_ratio * trace_ratio * ( beta * tau_max / real( ckink + 1 ) )

! determine pass, using important sampling algorithm (metropolis algorithm)
     pass = ( min( one, abs(p) ) > spring_sfmt_stream() )

! if the update action is accepted
     if ( pass .eqv. .true. ) then

! update the mmat matrix and gmat matrix, respectively,
! cat_insert_action() subroutine is invoked internally to update the
! perturbation expansion series
         call cat_insert_matrix(flvr, is, ie, tau_start, tau_end, deter_ratio)

! update ckink for current flavor channel
         ckink = ckink + 1

! update stts for current flavor channel
         stts(flvr) = cstat

! update rank for current flavor channel
         rank(flvr) = rank(flvr) + 1

     endif ! back if ( pass .eqv. .true. ) block

! update the insert statistics
     insert_tcount = insert_tcount + one
     if ( pass .eqv. .true. ) then
         insert_accept = insert_accept + one
     else
         insert_reject = insert_reject + one
     endif ! back if ( pass .eqv. .true. ) block

     return
  end subroutine ctqmc_insert_kink

!!>>> ctqmc_remove_kink: remove old segment or anti-segment in the
!!>>> perturbation expansion series
  subroutine ctqmc_remove_kink()
     use constants, only : dp, one
     use spring, only : spring_sfmt_stream

     use control, only : norbs
     use control, only : beta
     use context, only : ckink, cstat
     use context, only : remove_tcount, remove_accept, remove_reject
     use context, only : rank, stts

     implicit none

! local variables
! whether it is an anti-segment
! anti = .true., anti-segment
! anti = .false., segment
     logical  :: anti

! whether the update operation is accepted
     logical  :: pass

! current flavor channel for both band and spin
     integer  :: flvr

! index address to remove old segment or anti-segment
! is and ie are for start and end points, respectively
     integer  :: is, ie

! transition probability for remove old segment or anti-segment
     real(dp) :: p

! \tau_s, start point of the old segment
     real(dp) :: tau_start

! \tau_e, end   point of the old segment
     real(dp) :: tau_end

! the possible maximum length of the old segment
     real(dp) :: tau_max

! ratio between old and new configurations, the local trace part
     real(dp) :: trace_ratio

! ratio between old and new configurations, the determinant part
     real(dp) :: deter_ratio

! initialize logical variables
     anti = .false.
     pass = .false.

! select the flavor channel randomly among 1 ~ norbs
     flvr = ceiling( spring_sfmt_stream() * norbs )

! get the perturbation expansion order ( number of existing segments or
! anti-segments ) for current flavor channel
     ckink = rank(flvr)
     if ( ckink == 0 ) then
!<         call s_print_exception('ctqmc_remove_kink','can not remove any segments')
         remove_tcount = remove_tcount + one
         remove_reject = remove_reject + one
         RETURN
     endif ! back if ( ckink == 0 ) block

! at first determine anti and is randomly, then tau_start is obtained by
! is. and then ie, tau_end, and tau_max are evaluated carefully according
! to is and ie
     call cat_remove_flavor(flvr, is, ie, anti, tau_start, tau_end, tau_max)

! calculate the transition ratio between old and new configurations,
! for the local trace part
     call cat_remove_ztrace(flvr, anti, tau_start, tau_end, trace_ratio)

! calculate the transition ratio between old and new configurations,
! for the determinant part
     call cat_remove_detrat(flvr, is, ie, deter_ratio)

! calculate the transition probability for remove old segment or anti-segment
     p = deter_ratio * trace_ratio * real( ckink ) / ( beta * tau_max )

! determine pass, using important sampling algorithm (metropolis algorithm)
     pass = ( min( one, abs(p) ) > spring_sfmt_stream() )

! if update action is accepted
     if ( pass .eqv. .true. ) then

! update the mmat matrix and gmat matrix, respectively,
! cat_remove_action() subroutine is invoked internally to update the
! perturbation expansion series
         call cat_remove_matrix(flvr, is, ie)

! decrease ckink for current flavor channel
         ckink = ckink - 1

! update stts for current flavor channel
         stts(flvr) = cstat

! update rank for current flavor channel
         rank(flvr) = rank(flvr) - 1

     endif ! back if ( pass .eqv. .true. ) block

! update the remove statistics
     remove_tcount = remove_tcount + one
     if ( pass .eqv. .true. ) then
         remove_accept = remove_accept + one
     else
         remove_reject = remove_reject + one
     endif ! back if ( pass .eqv. .true. ) block

     return
  end subroutine ctqmc_remove_kink

!!>>> ctqmc_lshift_kink: left shift old segment or anti-segment in the
!!>>> perturbation expansion series
  subroutine ctqmc_lshift_kink()
     use constants, only : dp, one
     use spring, only : spring_sfmt_stream

     use control, only : norbs
     use context, only : ckink, cstat
     use context, only : lshift_tcount, lshift_accept, lshift_reject
     use context, only : rank, stts

     implicit none

! local variables
! whether the update operation winds around the circle
     logical  :: ring

! whether the update operation is accepted
     logical  :: pass

! current flavor channel for both band and spin
     integer  :: flvr

! index address to left shift old segment or anti-segment
! iso and isn are for old and new indices, respectively
     integer  :: iso, isn

! transition probability for left shift old segment or anti-segment
     real(dp) :: p

! \tau_s, start point of the old segment (old point)
     real(dp) :: tau_start1

! \tau_s, start point of the old segment (new point)
     real(dp) :: tau_start2

! ratio between old and new configurations, the local trace part
     real(dp) :: trace_ratio

! ratio between old and new configurations, the determinant part
     real(dp) :: deter_ratio

! initialize logical variables
     ring = .false.
     pass = .false.

! select the flavor channel randomly among 1 ~ norbs
     flvr = ceiling( spring_sfmt_stream() * norbs )

! get the perturbation expansion order ( number of existing segments or
! anti-segments ) for current flavor channel
     ckink = rank(flvr)
     if ( ckink == 0 ) then
!<         call s_print_exception('ctqmc_lshift_kink','can not lshift any segments')
         lshift_tcount = lshift_tcount + one
         lshift_reject = lshift_reject + one
         RETURN
     endif ! back if ( ckink == 0 ) block

! at first, we select iso randomly, and then obtain tau_start1. according
! to the existing segments, we determine tau_start2 and related index isn,
! finally ring is evaluated.
     call cat_lshift_flavor(flvr, iso, isn, ring, tau_start1, tau_start2)

! calculate the transition ratio between old and new configurations,
! for the local trace part
     call cat_lshift_ztrace(flvr, ring, tau_start1, tau_start2, trace_ratio)

! calculate the transition ratio between old and new configurations,
! for the determinant part
     call cat_lshift_detrat(flvr, iso, tau_start1, tau_start2, deter_ratio)

! calculate the transition probability for left shift old segment or anti-segment
     p = deter_ratio * trace_ratio

! determine pass, using important sampling algorithm (metropolis algorithm)
     pass = ( min( one, abs(p) ) > spring_sfmt_stream() )

! if update action is accepted
     if ( pass .eqv. .true. ) then

! update the mmat matrix and gmat matrix, respectively,
! cat_lshift_action() subroutine is invoked internally to update the
! perturbation expansion series
         call cat_lshift_matrix(flvr, iso, isn, tau_start1, tau_start2, deter_ratio)

! update stts for current flavor channel
         stts(flvr) = cstat

     endif ! back if ( pass .eqv. .true. ) block

! update the lshift statistics
     lshift_tcount = lshift_tcount + one
     if ( pass .eqv. .true. ) then
         lshift_accept = lshift_accept + one
     else
         lshift_reject = lshift_reject + one
     endif ! back if ( pass .eqv. .true. ) block

     return
  end subroutine ctqmc_lshift_kink

!!>>> ctqmc_rshift_kink: right shift old segment or anti-segment in the
!!>>> perturbation expansion series
  subroutine ctqmc_rshift_kink()
     use constants, only : dp, one
     use spring, only : spring_sfmt_stream

     use control, only : norbs
     use context, only : ckink, cstat
     use context, only : rshift_tcount, rshift_accept, rshift_reject
     use context, only : rank, stts

     implicit none

! local variables
! whether the update operation winds around the circle
     logical  :: ring

! whether the update operation is accepted
     logical  :: pass

! current flavor channel for both band and spin
     integer  :: flvr

! index address to right shift old segment or anti-segment
! ieo and ien are for old and new indices, respectively
     integer  :: ieo, ien

! transition probability for right shift old segment or anti-segment
     real(dp) :: p

! \tau_e, end   point of the old segment (old point)
     real(dp) :: tau_end1

! \tau_e, end   point of the old segment (new point)
     real(dp) :: tau_end2

! ratio between old and new configurations, the local trace part
     real(dp) :: trace_ratio

! ratio between old and new configurations, the determinant part
     real(dp) :: deter_ratio

! initialize logical variables
     ring = .false.
     pass = .false.

! select the flavor channel randomly among 1 ~ norbs
     flvr = ceiling( spring_sfmt_stream() * norbs )

! get the perturbation expansion order ( number of existing segments or
! anti-segments ) for current flavor channel
     ckink = rank(flvr)
     if ( ckink == 0 ) then
!<         call s_print_exception('ctqmc_rshift_kink','can not rshift any segments')
         rshift_tcount = rshift_tcount + one
         rshift_reject = rshift_reject + one
         RETURN
     endif ! back if ( ckink == 0 ) block

! at first, we select ieo randomly, and then obtain tau_end1. according
! to the existing segments, we determine tau_end2 and related index ien,
! finally ring is evaluated.
     call cat_rshift_flavor(flvr, ieo, ien, ring, tau_end1, tau_end2)

! calculate the transition ratio between old and new configurations,
! for the local trace part
     call cat_rshift_ztrace(flvr, ring, tau_end1, tau_end2, trace_ratio)

! calculate the transition ratio between old and new configurations,
! for the determinant part
     call cat_rshift_detrat(flvr, ieo, tau_end1, tau_end2, deter_ratio)

! calculate the transition probability for right shift old segment or anti-segment
     p = deter_ratio * trace_ratio

! determine pass, using important sampling algorithm (metropolis algorithm)
     pass = ( min( one, abs(p) ) > spring_sfmt_stream() )

! if update action is accepted
     if ( pass .eqv. .true. ) then

! update the mmat matrix and gmat matrix, respectively,
! cat_rshift_action() subroutine is invoked internally to update the
! perturbation expansion series
         call cat_rshift_matrix(flvr, ieo, ien, tau_end1, tau_end2, deter_ratio)

! update stts for current flavor channel
         stts(flvr) = cstat

     endif ! back if ( pass .eqv. .true. ) block

! update the rshift statistics
     rshift_tcount = rshift_tcount + one
     if ( pass .eqv. .true. ) then
         rshift_accept = rshift_accept + one
     else
         rshift_reject = rshift_reject + one
     endif ! back if ( pass .eqv. .true. ) block

     return
  end subroutine ctqmc_rshift_kink

!!>>> ctqmc_reswap_kink: perform a global update, exchange the segments
!!>>> and anti-segments
  subroutine ctqmc_reswap_kink()
     use constants, only : dp, one
     use spring, only : spring_sfmt_stream

     use control, only : norbs
     use context, only : ckink, cstat
     use context, only : reswap_tcount, reswap_accept, reswap_reject
     use context, only : rank, stts

     implicit none

! local variables
! whether the update operation is accepted
     logical  :: pass

! current flavor channel for both band and spin
     integer  :: flvr

! transition probability for global swap
     real(dp) :: p

! ratio between old and new configurations, the local trace part
     real(dp) :: trace_ratio

! ratio between old and new configurations, the determinant part
     real(dp) :: deter_ratio

! initialize logical variables
     pass = .false.

! select the flavor channel randomly among 1 ~ norbs
     flvr = ceiling( spring_sfmt_stream() * norbs )

! get the perturbation expansion order ( number of existing segments or
! anti-segments ) for current flavor channel
     ckink = rank(flvr)

! calculate the transition ratio between old and new configurations,
! for the local trace part
     call cat_reswap_ztrace(flvr, trace_ratio)

! calculate the transition ratio between old and new configurations,
! for the determinant part
! note: here deter_ratio is -1 or 1 in fact
     call cat_reswap_detrat(flvr, deter_ratio)

! calculate the transition probability for global swap
     p = deter_ratio * trace_ratio

! determine pass, using important sampling algorithm (metropolis algorithm)
     pass = ( min( one, abs(p) ) > spring_sfmt_stream() )

! if update action is accepted
     if ( pass .eqv. .true. ) then

! update the mmat matrix and gmat matrix, respectively,
! cat_reswap_action() subroutine is invoked internally to update the
! perturbation expansion series
         call cat_reswap_matrix(flvr)

! update stts for current flavor channel
         stts(flvr) = cstat

     endif ! back if ( pass .eqv. .true. ) block

! update the reswap statistics
     reswap_tcount = reswap_tcount + one
     if ( pass .eqv. .true. ) then
         reswap_accept = reswap_accept + one
     else
         reswap_reject = reswap_reject + one
     endif ! back if ( pass .eqv. .true. ) block

     return
  end subroutine ctqmc_reswap_kink

!!>>> ctqmc_reflip_kink: perform a global update, exchange the states
!!>>> between spin up and spin down, it maybe useful for magnetic systems
  subroutine ctqmc_reflip_kink(cflip)
     use constants, only : dp, one
     use spring, only : spring_sfmt_stream

     use control, only : nband, norbs
     use context, only : reflip_tcount, reflip_accept, reflip_reject
     use context, only : rank, symm

     implicit none

! external arguments
! control flag
! if cflip = 1, flip inter-orbital spins randomly;
! if cflip = 2, flip intra-orbital spins one by one;
! if cflip = 3, flip intra-orbital spins globally
     integer, intent(in) :: cflip

! local variables
! whether the update operation is accepted
     logical  :: pass

! selected flavor pairs
     integer  :: fup, fdn

! loop index for flavor channel
     integer  :: flvr

! maximum rank order
     integer  :: kmax

! transition probability for global spin flip
     real(dp) :: p

! global flip determinant ratio
     real(dp) :: ratup
     real(dp) :: ratdn

! initialize logical variables
     pass = .false.

! initialize transition probability
     p = one

     if ( cflip == 1 ) then
! determine fup and fdn, and fup /= fdn
         fup = ceiling( spring_sfmt_stream() * norbs )
         do while ( .true. )
             fdn = ceiling( spring_sfmt_stream() * norbs )
             if ( fup /= fdn .and. symm(fup) == symm(fdn) ) EXIT
         enddo ! over do while loop

! calculate the transition ratio between old and new configurations,
! for the determinant part, spin up case
         call cat_reflip_detrat(fup, fdn, ratup)

! calculate the transition ratio between old and new configurations,
! for the determinant part, spin dn case
         call cat_reflip_detrat(fdn, fup, ratdn)

! calculate the transition probability for global spin flip
         p = ratup * ratdn

! determine pass, using important sampling algorithm (metropolis algorithm)
         pass = ( min( one, abs(p) ) > spring_sfmt_stream() )

! if update action is accepted
         if ( pass .eqv. .true. ) then

! get maximum rank order in spin up and spin down states
             kmax = max( rank(fup), rank(fdn) )

! swap global variables between spin up and spin down states
             call cat_reflip_matrix(fup, fdn, kmax)

         endif ! back if ( pass .eqv. .true. ) block

! update the reflip statistics
         reflip_tcount = reflip_tcount + one
         if ( pass .eqv. .true. ) then
             reflip_accept = reflip_accept + one
         else
             reflip_reject = reflip_reject + one
         endif ! back if ( pass .eqv. .true. ) block

     else if ( cflip == 2 ) then ! cflip = 2, local flip
         do flvr=1,nband

! get fup and fdn
             fup = flvr; fdn = flvr + nband

! calculate the transition ratio between old and new configurations,
! for the determinant part, spin up case
             call cat_reflip_detrat(fup, fdn, ratup)

! calculate the transition ratio between old and new configurations,
! for the determinant part, spin dn case
             call cat_reflip_detrat(fdn, fup, ratdn)

! calculate the transition probability for global spin flip
             p = ratup * ratdn

! determine pass, using important sampling algorithm (metropolis algorithm)
             pass = ( min( one, abs(p) ) > spring_sfmt_stream() )

! if update action is accepted
             if ( pass .eqv. .true. ) then

! get maximum rank order in spin up and spin down states
                 kmax = max( rank(fup), rank(fdn) )

! swap global variables between spin up and spin down states
                 call cat_reflip_matrix(fup, fdn, kmax)

             endif ! back if ( pass .eqv. .true. ) block

! update the reflip statistics
             reflip_tcount = reflip_tcount + one
             if ( pass .eqv. .true. ) then
                 reflip_accept = reflip_accept + one
             else
                 reflip_reject = reflip_reject + one
             endif ! back if ( pass .eqv. .true. ) block

         enddo ! over flvr={1,nband} loop

     else if ( cflip == 3 ) then ! cflip = 3, global flip
         do flvr=1,nband

! get fup and fdn
             fup = flvr; fdn = flvr + nband

! calculate the transition ratio between old and new configurations,
! for the determinant part, spin up case
             call cat_reflip_detrat(fup, fdn, ratup)

! calculate the transition ratio between old and new configurations,
! for the determinant part, spin dn case
             call cat_reflip_detrat(fdn, fup, ratdn)

! calculate the transition probability for global spin flip
             p = p * ( ratup * ratdn )

         enddo ! over flvr={1,nband} loop

! determine pass, using important sampling algorithm (metropolis algorithm)
         pass = ( min( one, abs(p) ) > spring_sfmt_stream() )

! if update action is accepted
         if ( pass .eqv. .true. ) then

             do flvr=1,nband

! get fup and fdn
                 fup = flvr; fdn = flvr + nband

! get maximum rank order in spin up and spin down states
                 kmax = max( rank(fup), rank(fdn) )

! swap global variables between spin up and spin down states
                 call cat_reflip_matrix(fup, fdn, kmax)

             enddo ! over flvr={1,nband} loop

         endif ! back if ( pass .eqv. .true. ) block

! update the reflip statistics
         reflip_tcount = reflip_tcount + one
         if ( pass .eqv. .true. ) then
             reflip_accept = reflip_accept + one
         else
             reflip_reject = reflip_reject + one
         endif ! back if ( pass .eqv. .true. ) block

     endif ! back if ( cflip == 1 ) block

     return
  end subroutine ctqmc_reflip_kink

!!>>> ctqmc_reload_kink: global update all segments or anti-segments in
!!>>> the perturbation expansion series
  subroutine ctqmc_reload_kink()
     use control, only : norbs
     use context, only : rank

     implicit none

! local variables
! loop index for flavor channel
     integer :: flvr

     do flvr=1,norbs

! check the perturbation expansion order ( number of existing segments or
! anti-segments ) for current flavor channel
         if ( rank(flvr) == 0 ) CYCLE

! regenerate the mmat matrix and gmat matrix from scratch
         call cat_reload_matrix(flvr)

     enddo ! over flvr={1,norbs} loop

     return
  end subroutine ctqmc_reload_kink
