!!!-----------------------------------------------------------------------
!!! project : narcissus
!!! program : ctqmc_impurity_solver
!!!           ctqmc_impurity_tester
!!! source  : ctqmc_solver.f90
!!! type    : subroutines
!!! author  : li huang (email:lihuang.dmft@gmail.com)
!!! history : 09/16/2009 by li huang (created)
!!!           04/16/2019 by li huang (last modified)
!!! purpose : the main subroutines for the hybridization expansion version
!!!           continuous time quantum Monte Carlo (CTQMC) quantum impurity
!!!           solver. they implement the initialization, thermalization,
!!!           random walk, measurement, and finalization algorithms.
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
     use constants, only : dp
     use constants, only : zero, one
     use constants, only : mystd

     use control, only : cname               ! code name
                                             !
     use control, only : iswor               ! worm algorithm
     use control, only : isobs               ! control physical observables
     use control, only : issus               ! control spin and charge susceptibility
     use control, only : isvrt               ! control two-particle quantities
                                             !
     use control, only : nband, norbs, ncfgs ! size of model hamiltonian
     use control, only : mkink               ! perturbation expansion order
     use control, only : mfreq               ! matsubara frequency
     use control, only : nffrq, nbfrq        ! fermionic and bosonic frequencies
     use control, only : ntime               ! imaginary time slice
     use control, only : nsweep, nwrite      ! monte carlo sampling
     use control, only : nmonte              ! interval between successive measurements
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
     use context, only : g2ph                ! two-particle green's function (ph)
     use context, only : h2ph                ! two-particle vertex function (ph)
     use context, only : g2pp                ! two-particle green's function (pp)
     use context, only : h2pp                ! two-particle vertex function (pp)
                                             !
     use context, only : symm                ! symmetry indicator
                                             !
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

! the following observables are always measured: real(dp)
!-------------------------------------------------------------------------
! histogram for perturbation expansion series
     real(dp), allocatable :: hist_mpi(:)
     real(dp), allocatable :: hist_err(:)

! probability of atomic eigenstates
     real(dp), allocatable :: prob_mpi(:)
     real(dp), allocatable :: prob_err(:)

! auxiliary physical observables
     real(dp), allocatable :: paux_mpi(:)
     real(dp), allocatable :: paux_err(:)

! impurity occupation number, < n_i >
     real(dp), allocatable :: nimp_mpi(:)
     real(dp), allocatable :: nimp_err(:)

! impurity double occupation number matrix, < n_i n_j >
     real(dp), allocatable :: nmat_mpi(:,:)
     real(dp), allocatable :: nmat_err(:,:)

! impurity green's function in imaginary time axis
     real(dp), allocatable :: gtau_mpi(:,:,:)
     real(dp), allocatable :: gtau_err(:,:,:)

! auxiliary correlation function in imaginary time axis
     real(dp), allocatable :: ftau_mpi(:,:,:)
     real(dp), allocatable :: ftau_err(:,:,:)
!-------------------------------------------------------------------------

! the following observables are always measured: complex(dp)
!-------------------------------------------------------------------------
! impurity green's function in matsubara frequency axis
     complex(dp), allocatable :: grnf_mpi(:,:,:)
     complex(dp), allocatable :: grnf_err(:,:,:)

! auxiliary correlation function in matsubara frequency axis
     complex(dp), allocatable :: frnf_mpi(:,:,:)
     complex(dp), allocatable :: frnf_err(:,:,:)

! self-energy function in matsubara frequency axis
     complex(dp), allocatable :: sig2_mpi(:,:,:)
     complex(dp), allocatable :: sig2_err(:,:,:)
!-------------------------------------------------------------------------

! the following observables are measured optionally: real(dp)
!-------------------------------------------------------------------------
! number of operators, < k >
     real(dp), allocatable :: knop_mpi(:)
     real(dp), allocatable :: knop_err(:)

! crossing product of k_i and k_j, < k_i k_j >
     real(dp), allocatable :: kmat_mpi(:,:)
     real(dp), allocatable :: kmat_err(:,:)

! number of operators at left half axis, < k_l >
     real(dp), allocatable :: lnop_mpi(:)
     real(dp), allocatable :: lnop_err(:)

! number of operators at right half axis, < k_r >
     real(dp), allocatable :: rnop_mpi(:)
     real(dp), allocatable :: rnop_err(:)

! crossing product of k_l and k_r, < k_l k_r >
     real(dp), allocatable :: lrmm_mpi(:,:)
     real(dp), allocatable :: lrmm_err(:,:)

! powers of local magnetization, < S^n_z>
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

! totally-averaged charge-charge correlation function
     real(dp), allocatable :: cchi_mpi(:)
     real(dp), allocatable :: cchi_err(:)

! orbital-resolved charge-charge correlation function
     real(dp), allocatable :: ch_t_mpi(:,:,:)
     real(dp), allocatable :: ch_t_err(:,:,:)

! orbital-resolved charge-charge correlation function
     real(dp), allocatable :: ch_w_mpi(:,:,:)
     real(dp), allocatable :: ch_w_err(:,:,:)
!-------------------------------------------------------------------------

! the following observables are measured optionally: complex(dp)
!-------------------------------------------------------------------------
! two-particle green's function (ph channel)
     complex(dp), allocatable :: g2ph_mpi(:,:,:,:,:)
     complex(dp), allocatable :: g2ph_err(:,:,:,:,:)

! two-particle vertex function (ph channel)
     complex(dp), allocatable :: h2ph_mpi(:,:,:,:,:)
     complex(dp), allocatable :: h2ph_err(:,:,:,:,:)

! two-particle green's function (pp channel)
     complex(dp), allocatable :: g2pp_mpi(:,:,:,:,:)
     complex(dp), allocatable :: g2pp_err(:,:,:,:,:)

! two-particle vertex function (pp channel)
     complex(dp), allocatable :: h2pp_mpi(:,:,:,:,:)
     complex(dp), allocatable :: h2pp_err(:,:,:,:,:)
!-------------------------------------------------------------------------

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

     allocate(gtau_mpi(ntime,norbs,norbs), stat=istat)
     allocate(gtau_err(ntime,norbs,norbs), stat=istat)
     allocate(ftau_mpi(ntime,norbs,norbs), stat=istat)
     allocate(ftau_err(ntime,norbs,norbs), stat=istat)
     allocate(grnf_mpi(mfreq,norbs,norbs), stat=istat)
     allocate(grnf_err(mfreq,norbs,norbs), stat=istat)
     allocate(frnf_mpi(mfreq,norbs,norbs), stat=istat)
     allocate(frnf_err(mfreq,norbs,norbs), stat=istat)
     allocate(sig2_mpi(mfreq,norbs,norbs), stat=istat)
     allocate(sig2_err(mfreq,norbs,norbs), stat=istat)

     allocate(knop_mpi(norbs),             stat=istat)
     allocate(knop_err(norbs),             stat=istat)
     allocate(kmat_mpi(norbs,norbs),       stat=istat)
     allocate(kmat_err(norbs,norbs),       stat=istat)
     allocate(lnop_mpi(norbs),             stat=istat)
     allocate(lnop_err(norbs),             stat=istat)
     allocate(rnop_mpi(norbs),             stat=istat)
     allocate(rnop_err(norbs),             stat=istat)
     allocate(lrmm_mpi(norbs,norbs),       stat=istat)
     allocate(lrmm_err(norbs,norbs),       stat=istat)
     allocate(szpw_mpi(  4  ,norbs),       stat=istat)
     allocate(szpw_err(  4  ,norbs),       stat=istat)

     allocate(schi_mpi(ntime),             stat=istat)
     allocate(schi_err(ntime),             stat=istat)
     allocate(sp_t_mpi(ntime,nband),       stat=istat)
     allocate(sp_t_err(ntime,nband),       stat=istat)
     allocate(sp_w_mpi(nbfrq,nband),       stat=istat)
     allocate(sp_w_err(nbfrq,nband),       stat=istat)
     allocate(cchi_mpi(ntime),             stat=istat)
     allocate(cchi_err(ntime),             stat=istat)
     allocate(ch_t_mpi(ntime,norbs,norbs), stat=istat)
     allocate(ch_t_err(ntime,norbs,norbs), stat=istat)
     allocate(ch_w_mpi(nbfrq,norbs,norbs), stat=istat)
     allocate(ch_w_err(nbfrq,norbs,norbs), stat=istat)

     allocate(g2ph_mpi(nffrq,nffrq,nbfrq,norbs,norbs), stat=istat)
     allocate(g2ph_err(nffrq,nffrq,nbfrq,norbs,norbs), stat=istat)
     allocate(h2ph_mpi(nffrq,nffrq,nbfrq,norbs,norbs), stat=istat)
     allocate(h2ph_err(nffrq,nffrq,nbfrq,norbs,norbs), stat=istat)
     allocate(g2pp_mpi(nffrq,nffrq,nbfrq,norbs,norbs), stat=istat)
     allocate(g2pp_err(nffrq,nffrq,nbfrq,norbs,norbs), stat=istat)
     allocate(h2pp_mpi(nffrq,nffrq,nbfrq,norbs,norbs), stat=istat)
     allocate(h2pp_err(nffrq,nffrq,nbfrq,norbs,norbs), stat=istat)

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

! init or reset the continuous time quantum Monte Carlo quantum impurity
! solver, the key variables and arrays should be prepared here
     if ( myid == master ) then ! only master node can do it
         write(mystd,'(4X,a)') 'quantum impurity solver initializing'
     endif ! back if ( myid == master ) block

     call cpu_time(time_begin) ! record starting time
     call ctqmc_reset_array()
     call cpu_time(time_end) ! record ending time

! print the time information
     if ( myid == master ) then ! only master node can do it
         write(mystd,'(4X,a,f10.3,a)') 'time:', time_end - time_begin, 's'
         write(mystd,*)
     endif ! back if ( myid == master ) block

!!========================================================================
!!>>> retrieving quantum impurity solver                               <<<
!!========================================================================

! init the continuous time quantum Monte Carlo quantum impurity solver
! further, retrieving the diagrammatic series produced by previous run
     if ( myid == master ) then ! only master node can do it
         write(mystd,'(4X,a)') 'quantum impurity solver retrieving'
     endif ! back if ( myid == master ) block

     call cpu_time(time_begin) ! record starting time
     call ctqmc_retrieve_status()
     call cpu_time(time_end) ! record ending time

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
     call ctqmc_try_warming()
     call cpu_time(time_end) ! record ending time

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
         write(mystd,'(4X,a)',advance='no') 'RANDOM WALKING BEGINS'
         select case ( iswor )

             case (1)
                 write(mystd,'(1X,a)') '(CONFIGURATION SPACE: C_Z)'

             case (2)
                 write(mystd,'(1X,a)') '(CONFIGURATION SPACE: C_Z U C_G)'

         end select
         write(mystd,*)
     endif ! back if ( myid == master ) block

!!========================================================================
!!>>> visiting perturbation expansion series                           <<<
!!========================================================================

     MC_SWEEP: do i=1,nsweep,nwrite

! record start time
         call cpu_time(time_begin)

         MC_WRITE: do j=1,nwrite

! increase cstep by 1
             cstep = cstep + 1

! sampling the perturbation expansion feynman diagrams randomly
             call ctqmc_try_walking(cstep)

!!========================================================================
!!>>> sampling the physical observables 1 (always)                     <<<
!!========================================================================

! the following physical observables are always measured
! record the histogram for perturbation expansion series
             call ctqmc_record_hist()

! record the probability of atomic eigenstates
             if ( mod(cstep, nmonte) == 0 ) then
                 call ctqmc_record_prob()
             endif ! back if ( mod(cstep, nmonte) == 0 ) block

! record the auxiliary physical observables
             if ( mod(cstep, nmonte) == 0 ) then
                 call ctqmc_record_paux()
             endif ! back if ( mod(cstep, nmonte) == 0 ) block

! record the impurity (double) occupation number (matrix)
             if ( mod(cstep, nmonte) == 0 ) then
                 call ctqmc_record_nmat() ! AN EMPTY CALL
             endif ! back if ( mod(cstep, nmonte) == 0 ) block

!!========================================================================
!!>>> sampling the physical observables 2 (always)                     <<<
!!========================================================================

! the following physical observables are always measured
! record the impurity green's function in imaginary time space
             if ( mod(cstep, nmonte) == 0 ) then
                 call ctqmc_record_gtau()
             endif ! back if ( mod(cstep, nmonte) == 0 ) block

! record the auxiliary correlation function in imaginary time space
             if ( mod(cstep, nmonte) == 0 ) then
                 call ctqmc_record_ftau()
             endif ! back if ( mod(cstep, nmonte) == 0 ) block

! record the impurity green's function in matsubara frequency space
             if ( mod(cstep, nmonte) == 0 ) then
                 call ctqmc_record_grnf()
             endif ! back if ( mod(cstep, nmonte) == 0 ) block

! record the auxiliary correlation function in matsubara frequency space
             if ( mod(cstep, nmonte) == 0 ) then
                 call ctqmc_record_frnf() ! AN EMPTY CALL
             endif ! back if ( mod(cstep, nmonte) == 0 ) block

! record the self-energy function in matsubara frequency space
             if ( mod(cstep, nmonte) == 0 ) then
                 call ctqmc_record_sig2() ! AN EMPTY CALL
             endif ! back if ( mod(cstep, nmonte) == 0 ) block

!!========================================================================
!!>>> sampling the physical observables 3 (optional)                   <<<
!!========================================================================

! the following physical observables are measured optionally (by isobs)
! record the kinetic energy fluctuation
             if ( mod(cstep, nmonte) == 0 .and. btest(isobs, 1) ) then
                 call ctqmc_record_kmat()
             endif ! back if ( mod(cstep, nmonte) == 0 .and. btest(isobs, 1) ) block

! record the fidelity susceptibility
             if ( mod(cstep, nmonte) == 0 .and. btest(isobs, 2) ) then
                 call ctqmc_record_lrmm()
             endif ! back if ( mod(cstep, nmonte) == 0 .and. btest(isobs, 2) ) block

! record the powers of local magnetization
             if ( mod(cstep, nmonte) == 0 .and. btest(isobs, 3) ) then
                 call ctqmc_record_szpw()
             endif ! back if ( mod(cstep, nmonte) == 0 .and. btest(isobs, 3) ) block

!!========================================================================
!!>>> sampling the physical observables 4 (optional)                   <<<
!!========================================================================

! the following physical observables are measured optionally (by issus)
! record the spin-spin correlation function
             if ( mod(cstep, nmonte) == 0 .and. btest(issus, 1) ) then
                 call ctqmc_record_sp_t()
             endif ! back if ( mod(cstep, nmonte) == 0 .and. btest(issus, 1) ) block

! record the charge-charge correlation function
             if ( mod(cstep, nmonte) == 0 .and. btest(issus, 2) ) then
                 call ctqmc_record_ch_t()
             endif ! back if ( mod(cstep, nmonte) == 0 .and. btest(issus, 2) ) block

! record the spin-spin correlation function
             if ( mod(cstep, nmonte) == 0 .and. btest(issus, 3) ) then
                 call ctqmc_record_sp_w()
             endif ! back if ( mod(cstep, nmonte) == 0 .and. btest(issus, 3) ) block

! record the charge-charge correlation function
             if ( mod(cstep, nmonte) == 0 .and. btest(issus, 4) ) then
                 call ctqmc_record_ch_w()
             endif ! back if ( mod(cstep, nmonte) == 0 .and. btest(issus, 4) ) block

!!========================================================================
!!>>> sampling the physical observables 5 (optional)                   <<<
!!========================================================================

! the following physical observables are measured optionally (by isvrt)
! record the two-particle green's function, particle-hole channel, AABB
             if ( mod(cstep, nmonte) == 0 .and. btest(isvrt, 1) ) then
                 call ctqmc_record_g2ph()
             endif ! back if ( mod(cstep, nmonte) == 0 .and. btest(isvrt, 1) ) block

! record the two-particle green's function, particle-hole channel, ABBA
             if ( mod(cstep, nmonte) == 0 .and. btest(isvrt, 2) ) then
                 call ctqmc_record_g2ph()
             endif ! back if ( mod(cstep, nmonte) == 0 .and. btest(isvrt, 2) ) block

! record the two-particle green's function, particle-particle channel, AABB
             if ( mod(cstep, nmonte) == 0 .and. btest(isvrt, 3) ) then
                 call ctqmc_record_g2pp()
             endif ! back if ( mod(cstep, nmonte) == 0 .and. btest(isvrt, 3) ) block

! record the two-particle green's function, particle-particle channel, ABBA
             if ( mod(cstep, nmonte) == 0 .and. btest(isvrt, 4) ) then
                 call ctqmc_record_g2pp()
             endif ! back if ( mod(cstep, nmonte) == 0 .and. btest(isvrt, 4) ) block

         enddo MC_WRITE ! over j={1,nwrite} loop

!!========================================================================
!!>>> reporting quantum impurity solver                                <<<
!!========================================================================

! it is time to write out the statistics results
         if ( myid == master ) then ! only master node can do it
             call ctqmc_print_runtime(iter, cstep)
         endif ! back if ( myid == master ) block

! write out the snapshot for the current configuration if necessary
         if ( myid == master ) then ! only master node can do it
             call ctqmc_dump_diag(iter, cstep)
         endif ! back if ( myid == master ) block

!!========================================================================
!!>>> reducing immediate results                                       <<<
!!========================================================================

! collect the histogram data from hist to hist_mpi
         call ctqmc_reduce_hist(hist_mpi, hist_err)

!!========================================================================
!!>>> writing immediate results                                        <<<
!!========================================================================

! write out the histogram data, hist_mpi
         if ( myid == master ) then ! only master node can do it
             call ctqmc_dump_hist(hist_mpi, hist_err)
         endif ! back if ( myid == master ) block

!!========================================================================
!!>>> checking quantum impurity solver                                 <<<
!!========================================================================

! check the status at first
         call ctqmc_try_warning(cflag)

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

! print out the timing result
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

! end simulation
     if ( myid == master ) then ! only master node can do it
         write(mystd,'(4X,a)',advance='no') 'RANDOM WALKING STOPS'
         select case ( iswor )

             case (1)
                 write(mystd,'(1X,a)') '(CONFIGURATION SPACE: C_Z)'

             case (2)
                 write(mystd,'(1X,a)') '(CONFIGURATION SPACE: C_Z U C_G)'

         end select
         write(mystd,*)
     endif ! back if ( myid == master ) block

!!========================================================================
!!>>> reducing final results                                           <<<
!!========================================================================

! start to collect data
     if ( myid == master ) then ! only master node can do it
         write(mystd,'(4X,a)') 'quantum impurity solver reducing'
     endif ! back if ( myid == master ) block

     call cpu_time(time_begin) ! record starting time

! calculate the averaged values
     AVERAGE_DATA: BLOCK

         hist = hist * one
         prob = prob * real(nmonte) / real(nsweep)
         paux = paux * real(nmonte) / real(nsweep)
         nimp = nimp * real(nmonte) / real(nsweep)
         nmat = nmat * real(nmonte) / real(nsweep)

         gtau = gtau * real(nmonte) / real(nsweep)
         ftau = ftau * real(nmonte) / real(nsweep)
         grnf = grnf * real(nmonte) / real(nsweep)
         frnf = frnf * real(nmonte) / real(nsweep)
         sig2 = sig2 * real(nmonte) / real(nsweep)

         knop = knop * real(nmonte) / real(nsweep)
         kmat = kmat * real(nmonte) / real(nsweep)
         lnop = lnop * real(nmonte) / real(nsweep)
         rnop = rnop * real(nmonte) / real(nsweep)
         lrmm = lrmm * real(nmonte) / real(nsweep)
         szpw = szpw * real(nmonte) / real(nsweep)

         schi = schi * real(nmonte) / real(nsweep)
         sp_t = sp_t * real(nmonte) / real(nsweep)
         sp_w = sp_w * real(nmonte) / real(nsweep)
         cchi = cchi * real(nmonte) / real(nsweep)
         ch_t = ch_t * real(nmonte) / real(nsweep)
         ch_w = ch_w * real(nmonte) / real(nsweep)

         g2ph = g2ph * real(nmonte) / real(nsweep)
         h2ph = h2ph * real(nmonte) / real(nsweep)
         g2pp = g2pp * real(nmonte) / real(nsweep)
         h2pp = h2pp * real(nmonte) / real(nsweep)

     END BLOCK AVERAGE_DATA

! calculate some essential observables which are not measured directly
     COMPUTE_DATA: BLOCK

! try to evaluate the self-energy function
         call ctqmc_make_hub2()

! try to evaluate the imaginary time green's function
         call ctqmc_tran_gtau(gtau, gtau_mpi); gtau = gtau_mpi
         call ctqmc_tran_gtau(ftau, ftau_mpi); ftau = ftau_mpi

! try to evaluate the two-particle green's function (ph channel)
         if ( btest(isvrt, 1) .or. btest(isvrt, 2) ) then
             call ctqmc_tran_twop(g2ph, g2ph_mpi); g2ph = g2ph_mpi
             call ctqmc_tran_twop(h2ph, h2ph_mpi); h2ph = h2ph_mpi
         endif ! back if ( btest(isvrt, 1) .or. btest(isvrt, 2) ) block

! try to evaluate the two-particle green's function (pp channel)
         if ( btest(isvrt, 3) .or. btest(isvrt, 4) ) then
             call ctqmc_tran_twop(g2pp, g2pp_mpi); g2pp = g2pp_mpi
             call ctqmc_tran_twop(h2pp, h2pp_mpi); h2pp = h2pp_mpi
         endif ! back if ( btest(isvrt, 3) .or. btest(isvrt, 4) ) block

     END BLOCK COMPUTE_DATA

! collect data from all children processes
     COLLECT_DATA: BLOCK

         call ctqmc_reduce_hist(hist_mpi, hist_err)
         call ctqmc_reduce_prob(prob_mpi, prob_err)
         call ctqmc_reduce_paux(paux_mpi, paux_err)
         call ctqmc_reduce_nmat(nimp_mpi, nmat_mpi, nimp_err, nmat_err)

         call ctqmc_reduce_gtau(gtau_mpi, gtau_err)
         call ctqmc_reduce_ftau(ftau_mpi, ftau_err)
         call ctqmc_reduce_grnf(grnf_mpi, grnf_err)
         call ctqmc_reduce_frnf(frnf_mpi, frnf_err)
         call ctqmc_reduce_sig2(sig2_mpi, sig2_err)

         call ctqmc_reduce_kmat(knop_mpi, kmat_mpi, knop_err, kmat_err)
         call ctqmc_reduce_lrmm(lnop_mpi, rnop_mpi, lrmm_mpi, lnop_err, rnop_err, lrmm_err)
         call ctqmc_reduce_szpw(szpw_mpi, szpw_err)

         call ctqmc_reduce_sp_t(schi_mpi, sp_t_mpi, schi_err, sp_t_err)
         call ctqmc_reduce_sp_w(sp_w_mpi, sp_w_err)
         call ctqmc_reduce_ch_t(cchi_mpi, ch_t_mpi, cchi_err, ch_t_err)
         call ctqmc_reduce_ch_w(ch_w_mpi, ch_w_err)

         call ctqmc_reduce_g2ph(g2ph_mpi, h2ph_mpi, g2ph_err, h2ph_err)
         call ctqmc_reduce_g2pp(g2pp_mpi, h2pp_mpi, g2pp_err, h2pp_err)

     END BLOCK COLLECT_DATA

! update original data
     REPLACE_DATA: BLOCK

         hist = hist_mpi
         prob = prob_mpi
         paux = paux_mpi
         nimp = nimp_mpi
         nmat = nmat_mpi

         gtau = gtau_mpi
         ftau = ftau_mpi
         grnf = grnf_mpi
         frnf = frnf_mpi
         sig2 = sig2_mpi

         knop = knop_mpi
         kmat = kmat_mpi
         lnop = lnop_mpi
         rnop = rnop_mpi
         lrmm = lrmm_mpi
         szpw = szpw_mpi

         schi = schi_mpi
         sp_t = sp_t_mpi
         sp_w = sp_w_mpi
         cchi = cchi_mpi
         ch_t = ch_t_mpi
         ch_w = ch_w_mpi

         g2ph = g2ph_mpi
         h2ph = h2ph_mpi
         g2pp = g2pp_mpi
         h2pp = h2pp_mpi

     END BLOCK REPLACE_DATA

     call cpu_time(time_end) ! record ending time

! print the time information
     if ( myid == master ) then ! only master node can do it
         write(mystd,'(4X,a,f10.3,a)') 'time:', time_end - time_begin, 's'
         write(mystd,*)
     endif ! back if ( myid == master ) block

!!========================================================================
!!>>> symmetrizing final results                                       <<<
!!========================================================================

! start to symmetrize data
     if ( myid == master ) then ! only master node can do it
         write(mystd,'(4X,a)') 'quantum impurity solver symmetrizing'
     endif ! back if ( myid == master ) block

     call cpu_time(time_begin) ! record starting time

! symmetrize the occupation number matrix over spin or over bands
     call ctqmc_symm_nimp(symm, nimp)
     call ctqmc_symm_nimp(symm, nimp_err)

! symmetrize the impurity green's function over spin or over bands
     call ctqmc_symm_gtau(symm, gtau)
     call ctqmc_symm_gtau(symm, gtau_err)
     call ctqmc_symm_grnf(symm, grnf)
     call ctqmc_symm_grnf(symm, grnf_err)

! symmetrize the auxiliary correlation function over spin or over bands
     call ctqmc_symm_gtau(symm, ftau)
     call ctqmc_symm_gtau(symm, ftau_err)
     call ctqmc_symm_grnf(symm, frnf)
     call ctqmc_symm_grnf(symm, frnf_err)

! symmetrize the self-energy function over spin or over bands
     call ctqmc_symm_grnf(symm, sig2)
     call ctqmc_symm_grnf(symm, sig2_err)

     call cpu_time(time_end) ! record ending time

! print the time information
     if ( myid == master ) then ! only master node can do it
         write(mystd,'(4X,a,f10.3,a)') 'time:', time_end - time_begin, 's'
         write(mystd,*)
     endif ! back if ( myid == master ) block

!!========================================================================
!!>>> writing final results                                            <<<
!!========================================================================

! start to write data
     if ( myid == master ) then ! only master node can do it
         write(mystd,'(4X,a)') 'quantum impurity solver writing'
     endif ! back if ( myid == master ) block

     call cpu_time(time_begin) ! record starting time

! write out the final data to external files
     if ( myid == master ) then ! only master node can do it
         call ctqmc_dump_hist(hist, hist_err)
         call ctqmc_dump_prob(prob, prob_err)
         call ctqmc_dump_paux(paux, paux_err)
         call ctqmc_dump_nmat(nimp, nmat, nimp_err, nmat_err)

         call ctqmc_dump_gtau(gtau, gtau_err)
         call ctqmc_dump_grnf(grnf, grnf_err)
         call ctqmc_dump_ftau(ftau, ftau_err)
         call ctqmc_dump_frnf(frnf, frnf_err)
         call ctqmc_dump_sig2(sig2, sig2_err)

         call ctqmc_dump_kmat(knop, kmat, knop_err, kmat_err)
         call ctqmc_dump_lrmm(lnop, rnop, lrmm, lnop_err, rnop_err, lrmm_err)
         call ctqmc_dump_szpw(szpw, szpw_err)

         call ctqmc_dump_sp_t(schi, sp_t, schi_err, sp_t_err)
         call ctqmc_dump_sp_w(sp_w, sp_w_err)
         call ctqmc_dump_ch_t(cchi, ch_t, cchi_err, ch_t_err)
         call ctqmc_dump_ch_w(ch_w, ch_w_err)

         call ctqmc_dump_g2ph(g2ph, h2ph, g2ph_err, h2ph_err)
         call ctqmc_dump_g2pp(g2pp, h2pp, g2pp_err, h2pp_err)
     endif ! back if ( myid == master ) block

     call cpu_time(time_end) ! record ending time

! print the time information
     if ( myid == master ) then ! only master node can do it
         write(mystd,'(4X,a,f10.3,a)') 'time:', time_end - time_begin, 's'
         write(mystd,*)
     endif ! back if ( myid == master ) block

!!========================================================================
!!>>> saving quantum impurity solver                                   <<<
!!========================================================================

! start to save the diagrammatic information
     if ( myid == master ) then ! only master node can do it
         write(mystd,'(4X,a)') 'quantum impurity solver saving'
     endif ! back if ( myid == master ) block

     call cpu_time(time_begin) ! record starting time

! save the perturbation expansion series information to the disk file
     if ( myid == master ) then ! only master node can do it
         call ctqmc_save_status()
     endif ! back if ( myid == master ) block

     call cpu_time(time_end) ! record ending time

! print the time information
     if ( myid == master ) then ! only master node can do it
         write(mystd,'(4X,a,f10.3,a)') 'time:', time_end - time_begin, 's'
         write(mystd,*)
     endif ! back if ( myid == master ) block

!!========================================================================
!!>>> finishing quantum impurity solver                                <<<
!!========================================================================

! print the footer of continuous time quantum Monte Carlo quantum impurity
! solver. to tell the user it is over
     if ( myid == master ) then ! only master node can do it
         write(mystd,'(2X,a)') cname//' >>> CTQMC quantum impurity solver shutdown'
         write(mystd,*)
     endif ! back if ( myid == master ) block

! deallocate memory
     deallocate(hist_mpi)
     deallocate(hist_err)
     deallocate(prob_mpi)
     deallocate(prob_err)
     deallocate(paux_mpi)
     deallocate(paux_err)
     deallocate(nimp_mpi)
     deallocate(nimp_err)
     deallocate(nmat_mpi)
     deallocate(nmat_err)

     deallocate(gtau_mpi)
     deallocate(gtau_err)
     deallocate(ftau_mpi)
     deallocate(ftau_err)
     deallocate(grnf_mpi)
     deallocate(grnf_err)
     deallocate(frnf_mpi)
     deallocate(frnf_err)
     deallocate(sig2_mpi)
     deallocate(sig2_err)

     deallocate(knop_mpi)
     deallocate(knop_err)
     deallocate(kmat_mpi)
     deallocate(kmat_err)
     deallocate(lnop_mpi)
     deallocate(lnop_err)
     deallocate(rnop_mpi)
     deallocate(rnop_err)
     deallocate(lrmm_mpi)
     deallocate(lrmm_err)
     deallocate(szpw_mpi)
     deallocate(szpw_err)

     deallocate(schi_mpi)
     deallocate(schi_err)
     deallocate(sp_t_mpi)
     deallocate(sp_t_err)
     deallocate(sp_w_mpi)
     deallocate(sp_w_err)
     deallocate(cchi_mpi)
     deallocate(cchi_err)
     deallocate(ch_t_mpi)
     deallocate(ch_t_err)
     deallocate(ch_w_mpi)
     deallocate(ch_w_err)

     deallocate(g2ph_mpi)
     deallocate(g2ph_err)
     deallocate(h2ph_mpi)
     deallocate(h2ph_err)
     deallocate(g2pp_mpi)
     deallocate(g2pp_err)
     deallocate(h2pp_mpi)
     deallocate(h2pp_err)

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

     call cat_disp_diagrams(2)
     call s_print_error('ctqmc_impurity_tester','in debug mode')

     return
  end subroutine ctqmc_impurity_tester
