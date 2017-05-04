!!!-----------------------------------------------------------------------
!!! project : narcissus
!!! program : ctqmc_impurity_solver
!!!           ctqmc_impurity_tester
!!! source  : ctqmc_solver.f90
!!! type    : subroutines
!!! author  : li huang (email:lihuang.dmft@gmail.com)
!!! history : 09/16/2009 by li huang (created)
!!!           05/05/2017 by li huang (last modified)
!!! purpose : the main subroutines for the hybridization expansion version
!!!           continuous time quantum Monte Carlo (CTQMC) quantum impurity
!!!           solver. they are the most important subroutines
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
                                             !
     use control, only : isobs               ! control physical observables
     use control, only : issus               ! control sp/ch susceptibility
     use control, only : isvrt               ! control two-particle quantities
                                             !
     use control, only : nband, norbs, ncfgs ! size of model Hamiltonian
     use control, only : mkink               ! perturbation expansion order
     use control, only : mfreq               ! matsubara frequency
     use control, only : nffrq, nbfrq        ! fermionic and bosonic frequencies
     use control, only : ntime               ! imaginary time
     use control, only : nsweep, nwrite      ! monte carlo sampling
     use control, only : nmonte, ncarlo      ! interval for monte carlo sampling
     use control, only : myid, master        ! mpi environment
                                             !
     use context, only : hist                ! histogram
     use context, only : prob                ! atomic eigenstate probability
     use context, only : paux                ! auxiliary physical observables
     use context, only : nimp, nmat          ! occupation and double occupation
                                             !
     use context, only : knop, kmat          ! kinetic energy fluctuation
     use context, only : lnop, rnop, lrmm    ! fidelity susceptibility
     use context, only : szpw                ! binder cumulant
                                             !
     use context, only : schi, sp_t, sp_w    ! spin susceptibility
     use context, only : cchi, ch_t, ch_w    ! charge susceptibility
                                             !
     use context, only : g2pw                ! two-particle green's function
     use context, only : h2pw                ! irreducible vertex function
     use context, only : p2pw                ! pairing susceptibility
                                             !
     use context, only : symm                ! symmetry
     use context, only : gtau, ftau          ! imaginary time green's function
     use context, only : grnf, frnf          ! matsubara green's function
     use context, only : sig2                ! matsubara self-energy function

     implicit none

! external arguments
! current iteration number for self-consistent cycle
     integer, intent(in) :: iter

! local variables
! loop index
     integer  :: i
     integer  :: j

! status flag
     integer  :: istat

! current QMC sweeping steps
     integer  :: cstep

! control flag, whether the solver should be checked periodically
! cflag = 0 , do not check the quantum impurity solver
! cflag = 1 , check the quantum impurity solver periodically
! cflag = 99, the quantum impurity solver is out of control
     integer  :: cflag

! starting time
     real(dp) :: time_begin

! ending time
     real(dp) :: time_end

! time consuming by current iteration
     real(dp) :: time_cur

! time consuming by total iteration
     real(dp) :: time_sum

! histogram for perturbation expansion series
     real(dp), allocatable :: hist_mpi(:)
     real(dp), allocatable :: hist_err(:)

! probability of atomic eigenstates
     real(dp), allocatable :: prob_mpi(:)
     real(dp), allocatable :: prob_err(:)

! auxiliary physical observables
     real(dp), allocatable :: paux_mpi(:)
     real(dp), allocatable :: paux_err(:)

! impurity occupation number matrix
     real(dp), allocatable :: nimp_mpi(:)
     real(dp), allocatable :: nimp_err(:)

! impurity double occupation number matrix
     real(dp), allocatable :: nmat_mpi(:,:)
     real(dp), allocatable :: nmat_err(:,:)

! impurity green's function in imaginary time axis
     real(dp), allocatable :: gtau_mpi(:,:,:)
     real(dp), allocatable :: gtau_err(:,:,:)

! auxiliary correlation function in imaginary time axis
     real(dp), allocatable :: ftau_mpi(:,:,:)
     real(dp), allocatable :: ftau_err(:,:,:)

! number of operators < k >
     real(dp), allocatable :: knop_mpi(:)
     real(dp), allocatable :: knop_err(:)

! square of number of operators < k^2 >
     real(dp), allocatable :: kmat_mpi(:,:)
     real(dp), allocatable :: kmat_err(:,:)

! number of operators at left half axis < k_l >
     real(dp), allocatable :: lnop_mpi(:)
     real(dp), allocatable :: lnop_err(:)

! number of operators at right half axis < k_r >
     real(dp), allocatable :: rnop_mpi(:)
     real(dp), allocatable :: rnop_err(:)

! used to evaluate fidelity susceptibility < k_l k_r >
     real(dp), allocatable :: lrmm_mpi(:,:)
     real(dp), allocatable :: lrmm_err(:,:)

! powers of local magnetization
     real(dp), allocatable :: szpw_mpi(:,:)
     real(dp), allocatable :: szpw_err(:,:)

! totally-averaged spin-spin correlation function
     real(dp), allocatable :: schi_mpi(:)
     real(dp), allocatable :: schi_err(:)

! orbital-resolved spin-spin correlation function
     real(dp), allocatable :: sp_t_mpi(:,:)
     real(dp), allocatable :: sp_t_err(:,:)

! orbital-resolved spin-spin correlation function
     real(dp), allocatable :: sp_w_mpi(:,:)
     real(dp), allocatable :: sp_w_err(:,:)

! totally-averaged orbital-orbital correlation function
     real(dp), allocatable :: cchi_mpi(:)
     real(dp), allocatable :: cchi_err(:)

! orbital-resolved orbital-orbital correlation function
     real(dp), allocatable :: ch_t_mpi(:,:,:)
     real(dp), allocatable :: ch_t_err(:,:,:)

! orbital-resolved orbital-orbital correlation function
     real(dp), allocatable :: ch_w_mpi(:,:,:)
     real(dp), allocatable :: ch_w_err(:,:,:)

! two-particle green's function
     complex(dp), allocatable :: g2pw_mpi(:,:,:,:,:)
     complex(dp), allocatable :: g2pw_err(:,:,:,:,:)

! irreducible vertex function
     complex(dp), allocatable :: h2pw_mpi(:,:,:,:,:)
     complex(dp), allocatable :: h2pw_err(:,:,:,:,:)

! particle-particle pairing susceptibility
     complex(dp), allocatable :: p2pw_mpi(:,:,:,:,:)
     complex(dp), allocatable :: p2pw_err(:,:,:,:,:)

! impurity green's function in matsubara frequency axis
     complex(dp), allocatable :: grnf_mpi(:,:,:)
     complex(dp), allocatable :: grnf_err(:,:,:)

! allocate memory
     allocate(hist_mpi(mkink),             stat=istat)
     allocate(hist_err(mkink),             stat=istat)
     allocate(prob_mpi(ncfgs),             stat=istat)
     allocate(prob_err(ncfgs),             stat=istat)
     allocate(paux_mpi(  9  ),             stat=istat)
     allocate(paux_err(  9  ),             stat=istat)
     allocate(nimp_mpi(norbs),             stat=istat)
     allocate(nimp_err(norbs),             stat=istat)
     allocate(nmat_mpi(norbs,norbs),       stat=istat)
     allocate(nmat_err(norbs,norbs),       stat=istat)

     allocate(knop_mpi(norbs),             stat=istat)
     allocate(knop_err(norbs),             stat=istat)
     allocate(kmat_mpi(norbs,norbs),       stat=istat)
     allocate(kmat_err(norbs,norbs),       stat=istat)
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
     time_cur = zero
     time_sum = zero

!!========================================================================
!!>>> starting quantum impurity solver                                 <<<
!!========================================================================

! print the header of continuous time quantum Monte Carlo quantum impurity
! solver. it contains important information about the control parameters
     if ( myid == master ) then ! only master node can do it
         call ctqmc_print_control()
     endif ! back if ( myid == master ) block

!!========================================================================
!!>>> initializing quantum impurity solver                             <<<
!!========================================================================

! init the continuous time quantum Monte Carlo quantum impurity solver,
! the key variables and arrays should be prepared here
     if ( myid == master ) then ! only master node can do it
         write(mystd,'(4X,a)') 'quantum impurity solver initializing'
     endif ! back if ( myid == master ) block

     call cpu_time(time_begin) ! record starting time
     call ctqmc_reset_array()
     call cpu_time(time_end)   ! record ending   time

! print the time information
     if ( myid == master ) then ! only master node can do it
         write(mystd,'(4X,a,f10.3,a)') 'time:', time_end - time_begin, 's'
         write(mystd,*)
     endif ! back if ( myid == master ) block

!!========================================================================
!!>>> retrieving quantum impurity solver                               <<<
!!========================================================================

! init the continuous time quantum Monte Carlo quantum impurity solver
! further, retrieving the time series information produced by previous run
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
!!>>> warming quantum impurity solver                                  <<<
!!========================================================================

! warmup the continuous time quantum Monte Carlo quantum impurity solver
! in order to achieve equilibrium state
     if ( myid == master ) then ! only master node can do it
         write(mystd,'(4X,a)') 'quantum impurity solver warmming'
     endif ! back if ( myid == master ) block

     call cpu_time(time_begin) ! record starting time
     call ctqmc_warming()
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

     MC_SWEEP: do i=1,nsweep,nwrite

! record start time
         call cpu_time(time_begin)

         MC_BLOCK: do j=1,nwrite

!!========================================================================
!!>>> visiting perturbation expansion series                           <<<
!!========================================================================

! increase cstep by 1
             cstep = cstep + 1

! sampling the perturbation expansion feynman diagrams randomly
             call ctqmc_walking(cstep)

!!========================================================================
!!>>> sampling the physical observables 1 (always)                     <<<
!!========================================================================

! record the histogram for perturbation expansion series
             call ctqmc_record_hist()

! record the probability of eigenstates
             if ( mod(cstep, nmonte) == 0 ) then
                 call ctqmc_record_prob()
             endif ! back if ( mod(cstep, nmonte) == 0 ) block

! record the auxiliary physical observables
             if ( mod(cstep, nmonte) == 0 ) then
                 call ctqmc_record_paux()
             endif ! back if ( mod(cstep, nmonte) == 0 ) block

! record the impurity (double) occupation number matrix
             if ( mod(cstep, nmonte) == 0 ) then
                 call ctqmc_record_nmat()
             endif ! back if ( mod(cstep, nmonte) == 0 ) block

!!========================================================================
!!>>> sampling the physical observables 2 (always)                     <<<
!!========================================================================

! record the impurity green's function in imaginary time space
             if ( mod(cstep, nmonte) == 0 ) then
                 call ctqmc_record_gtau()
             endif ! back if ( mod(cstep, nmonte) == 0 ) block

! record the auxiliary correlation function to calculate self-energy function
             if ( mod(cstep, nmonte) == 0 ) then
                 call ctqmc_record_ftau()
             endif ! back if ( mod(cstep, nmonte) == 0 ) block

! record the impurity green's function in matsubara frequency space
             if ( mod(cstep, nmonte) == 0 ) then
                 call ctqmc_record_grnf()
             endif ! back if ( mod(cstep, nmonte) == 0 ) block

!!========================================================================
!!>>> sampling the physical observables 3 (optional)                   <<<
!!========================================================================

! record the < k^2 > - < k >^2
             if ( mod(cstep, nmonte) == 0 .and. btest(isobs, 1) ) then
                 call ctqmc_record_kmat()
             endif ! back if ( mod(cstep, nmonte) == 0 .and. btest(isobs, 1) ) block

! record the fidelity susceptibility
             if ( mod(cstep, nmonte) == 0 .and. btest(isobs, 2) ) then
                 call ctqmc_record_lmat()
             endif ! back if ( mod(cstep, nmonte) == 0 .and. btest(isobs, 2) ) block

! record the powers of local magnetization
             if ( mod(cstep, nmonte) == 0 .and. btest(isobs, 3) ) then
                 call ctqmc_record_szpw()
             endif ! back if ( mod(cstep, nmonte) == 0 .and. btest(isobs, 3) ) block

!!========================================================================
!!>>> sampling the physical observables 4 (optional)                   <<<
!!========================================================================

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

!!========================================================================
!!>>> sampling the physical observables 5 (optional)                   <<<
!!========================================================================

! record the two-particle green's function
             if ( mod(cstep, nmonte) == 0 .and. btest(isvrt, 1) ) then
                 call ctqmc_record_twop()
             endif ! back if ( mod(cstep, nmonte) == 0 .and. btest(isvrt, 1) ) block

! record the particle-particle pairing susceptibility
             if ( mod(cstep, nmonte) == 0 .and. btest(isvrt, 2) ) then
                 call ctqmc_record_pair()
             endif ! back if ( mod(cstep, nmonte) == 0 .and. btest(isvrt, 2) ) block

         enddo MC_BLOCK ! over j={1,nwrite} loop

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

! the data need to be scaled properly before written
         hist_mpi = hist_mpi * one
         hist_err = hist_err * one
         gtau_mpi = gtau_mpi * real(nmonte) / real(cstep)
         gtau_err = gtau_err * real(nmonte) / real(cstep)

!!========================================================================
!!>>> symmetrizing immediate results                                   <<<
!!========================================================================

! symmetrize the impurity green's function over spin or over bands
         call ctqmc_symm_gtau(symm, gtau_mpi)
         call ctqmc_symm_gtau(symm, gtau_err)

!!========================================================================
!!>>> writing immediate results                                        <<<
!!========================================================================

! write out the histogram data, hist_mpi
         if ( myid == master ) then ! only master node can do it
             call ctqmc_dump_hist(hist_mpi, hist_err)
         endif ! back if ( myid == master ) block

! write out the impurity green's function, gtau_mpi
         if ( myid == master ) then ! only master node can do it
             call ctqmc_dump_gtau(gtau_mpi, gtau_err)
         endif ! back if ( myid == master ) block

!!========================================================================
!!>>> checking quantum impurity solver                                 <<<
!!========================================================================

! check the status at first
         call ctqmc_warning(cflag)

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
         time_cur = time_end - time_begin
         time_sum = time_sum + time_cur

! reset timer
         time_begin = time_end

! print out the result
         if ( myid == master ) then ! only master node can do it
             call s_time_analyzer(time_cur, time_sum)
             write(mystd,*)
         endif ! back if ( myid == master ) block

!!========================================================================
!!>>> escaping quantum impurity solver                                 <<<
!!========================================================================

! if the quantum impurity solver is out of control
         if ( cflag == 99 ) then
             EXIT MC_SWEEP ! jump out the iteration
         endif ! back if ( cflag == 99 ) block

     enddo MC_SWEEP ! over i={1,nsweep} loop

!!========================================================================
!!>>> ending main iteration                                            <<<
!!========================================================================

!!========================================================================
!!>>> reducing final results                                           <<<
!!========================================================================

! collect data from all children processes
     call ctqmc_reduce_hist(hist_mpi, hist_err)
     call ctqmc_reduce_prob(prob_mpi, prob_err)
     call ctqmc_reduce_paux(paux_mpi, paux_err)
     call ctqmc_reduce_nmat(nmat_mpi, nnmat_mpi, nmat_err, nnmat_err)

     call ctqmc_reduce_gtau(gtau_mpi, gtau_err)
     call ctqmc_reduce_ftau(ftau_mpi, ftau_err)
     call ctqmc_reduce_grnf(grnf_mpi, grnf_err)

     call ctqmc_reduce_kmat(kmat_mpi, kkmat_mpi, kmat_err, kkmat_err)
     call ctqmc_reduce_lmat(lmat_mpi, rmat_mpi, lrmat_mpi, lmat_err, rmat_err, lrmat_err)
     call ctqmc_reduce_szpw(szpow_mpi, szpow_err)

     call ctqmc_reduce_schi(schi_mpi, sschi_mpi, schi_err, sschi_err)
     call ctqmc_reduce_sfom(ssfom_mpi, ssfom_err)
     call ctqmc_reduce_ochi(ochi_mpi, oochi_mpi, ochi_err, oochi_err)
     call ctqmc_reduce_ofom(oofom_mpi, oofom_err)

     call ctqmc_reduce_twop(g2_re_mpi, g2_im_mpi, h2_re_mpi, h2_im_mpi)
     call ctqmc_reduce_pair(ps_re_mpi, ps_im_mpi)

! update original data and calculate the averages simultaneously
! average value section
     hist  = hist_mpi  * one
     prob  = prob_mpi  * real(nmonte) / real(nsweep)
     paux  = paux_mpi  * real(nmonte) / real(nsweep)
     nmat  = nmat_mpi  * real(nmonte) / real(nsweep)
     nnmat = nnmat_mpi * real(nmonte) / real(nsweep)

     gtau  = gtau_mpi  * real(nmonte) / real(nsweep)
     ftau  = ftau_mpi  * real(nmonte) / real(nsweep)
     grnf  = grnf_mpi  * real(nmonte) / real(nsweep)

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

! update original data and calculate the averages simultaneously
! error bar section
     hist_err  = hist_err  * one
     prob_err  = prob_err  * real(nmonte) / real(nsweep)
     paux_err  = paux_err  * real(nmonte) / real(nsweep)
     nmat_err  = nmat_err  * real(nmonte) / real(nsweep)
     nnmat_err = nnmat_err * real(nmonte) / real(nsweep)

     gtau_err  = gtau_err  * real(nmonte) / real(nsweep)
     ftau_err  = ftau_err  * real(nmonte) / real(nsweep)
     grnf_err  = grnf_err  * real(nmonte) / real(nsweep)

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

! try to evaluate the impurity green's function and self-energy function
     call ctqmc_make_hub2()

!!========================================================================
!!>>> symmetrizing final results                                       <<<
!!========================================================================

! symmetrize the occupation number matrix over spin or over bands
     call ctqmc_symm_nmat(symm, nmat)
     call ctqmc_symm_nmat(symm, nmat_err)

! symmetrize the impurity green's function over spin or over bands
     call ctqmc_symm_gtau(symm, gtau)
     call ctqmc_symm_gtau(symm, gtau_err)
     call ctqmc_symm_grnf(symm, grnf)

! symmetrize the auxiliary correlation function over spin or over bands
     call ctqmc_symm_gtau(symm, ftau)
     call ctqmc_symm_gtau(symm, ftau_err)
     call ctqmc_symm_grnf(symm, frnf)

! symmetrize the self-energy function over spin or over bands
     call ctqmc_symm_grnf(symm, sig2)

!!========================================================================
!!>>> writing final results                                            <<<
!!========================================================================

! write out the final data to external files
     if ( myid == master ) then ! only master node can do it
         call ctqmc_dump_hist(hist, hist_err)
         call ctqmc_dump_prob(prob, prob_err)
         call ctqmc_dump_paux(paux, paux_err)
         call ctqmc_dump_nmat(nmat, nnmat, nmat_err, nnmat_err)

         call ctqmc_dump_gtau(gtau, gtau_err)
         call ctqmc_dump_grnf(grnf, grnf_err)
         call ctqmc_dump_ftau(ftau, ftau_err)
         call ctqmc_dump_sigf(sig2)

         call ctqmc_dump_kmat(kmat, kkmat, kmat_err, kkmat_err)
         call ctqmc_dump_lmat(lmat, rmat, lrmat, lmat_err, rmat_err, lrmat_err)
         call ctqmc_dump_szpw(szpow, szpow_err)

         call ctqmc_dump_schi(schi, sschi, schi_err, sschi_err)
         call ctqmc_dump_sfom(ssfom, ssfom_err)
         call ctqmc_dump_ochi(ochi, oochi, ochi_err, oochi_err)
         call ctqmc_dump_ofom(oofom, oofom_err)

         call ctqmc_dump_twop(g2_re, g2_im, h2_re, h2_im)
         call ctqmc_dump_pair(ps_re, ps_im)
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

! print the footer of continuous time quantum Monte Carlo quantum impurity
! solver. to tell the user it is done
     if ( myid == master ) then ! only master node can do it
         write(mystd,'(2X,a)') cname//' >>> CTQMC quantum impurity solver shutdown'
         write(mystd,*)
     endif ! back if ( myid == master ) block

! deallocate memory
     deallocate(hist_mpi )
     deallocate(hist_err )
     deallocate(prob_mpi )
     deallocate(prob_err )
     deallocate(paux_mpi )
     deallocate(paux_err )
     deallocate(nmat_mpi )
     deallocate(nmat_err )
     deallocate(nnmat_mpi)
     deallocate(nnmat_err)

     deallocate(gtau_mpi )
     deallocate(gtau_err )
     deallocate(ftau_mpi )
     deallocate(ftau_err )
     deallocate(grnf_mpi )
     deallocate(grnf_err )

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

     return
  end subroutine ctqmc_impurity_solver

!!========================================================================
!!>>> debug subroutines for quantum impurity solver                    <<<
!!========================================================================

!!
!! @sub ctqmc_impurity_tester
!!
!! debug subroutine, please try to active it on the ctqmc_walking()
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
