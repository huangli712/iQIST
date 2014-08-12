!-------------------------------------------------------------------------
! project : lavender
! program : ctqmc_impurity_solver
!           ctqmc_diagram_warmming
!           ctqmc_diagram_sampling
!           ctqmc_diagram_checking
!           ctqmc_impurity_tester
! source  : ctqmc_solver.f90
! type    : subroutines
! author  : li huang (email:huangli712@yahoo.com.cn)
! history : 09/16/2009 by li huang
!           09/20/2009 by li huang
!           09/24/2009 by li huang
!           09/26/2009 by li huang
!           10/20/2009 by li huang
!           10/29/2009 by li huang
!           11/01/2009 by li huang
!           11/17/2009 by li huang
!           11/22/2009 by li huang
!           12/02/2009 by li huang
!           12/04/2009 by li huang
!           12/06/2009 by li huang
!           12/17/2009 by li huang
!           12/22/2009 by li huang
!           12/26/2009 by li huang
!           12/30/2009 by li huang
!           01/13/2010 by li huang
!           02/27/2010 by li huang
!           06/09/2010 by li huang
!           06/21/2010 by li huang
! purpose : the main subroutine for the hybridization expansion version
!           continuous time quantum Monte Carlo (CTQMC) quantum impurity
!           solver
! input   :
! output  :
! status  : unstable
! comment :
!-------------------------------------------------------------------------

!>>> core engine for hybridization expansion version continuous time
! quantum Monte Carlo quantum impurity solver
  subroutine ctqmc_impurity_solver(iter)
     use constants
     use control
     use context

     use m_sector
     use m_npart

     implicit none

! external arguments
! current self-consistent iteration number
     integer, intent(in) :: iter

! local variables
! loop index
     integer  :: i
     integer  :: j
     integer  :: m
     integer  :: n

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
     integer, allocatable  :: hist_mpi(:)

! spin-spin correlation function, totally-averaged, for mpi case
     real(dp), allocatable :: schi_mpi(:)

! orbital-orbital correlation function, totally-averaged, for mpi case
     real(dp), allocatable :: ochi_mpi(:)

! impurity occupation number matrix, for mpi case
     real(dp), allocatable :: nmat_mpi(:)

! probability of atomic states, for mpi case
     real(dp), allocatable :: prob_mpi(:)

! spin-spin correlation function, orbital-resolved, for mpi case
     real(dp), allocatable :: sschi_mpi(:,:)

! orbital-orbital correlation function, orbital-resolved, for mpi case
     real(dp), allocatable :: oochi_mpi(:,:)

! impurity double occupation number matrix, for mpi case
     real(dp), allocatable :: nnmat_mpi(:,:)

! impurity green's function, imaginary time axis, for mpi case
     real(dp), allocatable :: gtau_mpi(:,:,:)

! auxiliary correlation function, imaginary time axis, for mpi case
     real(dp), allocatable :: ftau_mpi(:,:,:)

! used to measure two-particle green's function, real part, for mpi case
     real(dp), allocatable :: g2_re_mpi(:,:,:,:,:)

! used to measure two-particle green's function, imaginary part, for mpi case
     real(dp), allocatable :: g2_im_mpi(:,:,:,:,:)

! used to measure vertex function, real part, for mpi case
     real(dp), allocatable :: h2_re_mpi(:,:,:,:,:)

! used to measure vertex function, imaginary part, for mpi case
     real(dp), allocatable :: h2_im_mpi(:,:,:,:,:)

! impurity green's function, matsubara frequency axis, for mpi case
     complex(dp), allocatable :: grnf_mpi(:,:,:)

! allocate memory
     allocate(hist_mpi(mkink),             stat=istat)
     if ( istat /= 0 ) then
         call s_print_error('ctqmc_impurity_solver','can not allocate enough memory')
     endif

     allocate(schi_mpi(ntime),             stat=istat)
     if ( istat /= 0 ) then
         call s_print_error('ctqmc_impurity_solver','can not allocate enough memory')
     endif

     allocate(ochi_mpi(ntime),             stat=istat)
     if ( istat /= 0 ) then
         call s_print_error('ctqmc_impurity_solver','can not allocate enough memory')
     endif

     allocate(nmat_mpi(norbs),             stat=istat)
     if ( istat /= 0 ) then
         call s_print_error('ctqmc_impurity_solver','can not allocate enough memory')
     endif

     allocate(prob_mpi(ncfgs),             stat=istat)
     if ( istat /= 0 ) then
         call s_print_error('ctqmc_impurity_solver','can not allocate enough memory')
     endif

     allocate(sschi_mpi(ntime,nband),      stat=istat)
     if ( istat /= 0 ) then
         call s_print_error('ctqmc_impurity_solver','can not allocate enough memory')
     endif

     allocate(oochi_mpi(ntime,norbs),      stat=istat)
     if ( istat /= 0 ) then
         call s_print_error('ctqmc_impurity_solver','can not allocate enough memory')
     endif

     allocate(nnmat_mpi(norbs,norbs),      stat=istat)
     if ( istat /= 0 ) then
         call s_print_error('ctqmc_impurity_solver','can not allocate enough memory')
     endif

     allocate(gtau_mpi(ntime,norbs,norbs), stat=istat)
     if ( istat /= 0 ) then
         call s_print_error('ctqmc_impurity_solver','can not allocate enough memory')
     endif

     allocate(ftau_mpi(ntime,norbs,norbs), stat=istat)
     if ( istat /= 0 ) then
         call s_print_error('ctqmc_impurity_solver','can not allocate enough memory')
     endif

     allocate(g2_re_mpi(norbs,norbs,nffrq,nffrq,nbfrq), stat=istat)
     if ( istat /= 0 ) then
         call s_print_error('ctqmc_impurity_solver','can not allocate enough memory')
     endif

     allocate(g2_im_mpi(norbs,norbs,nffrq,nffrq,nbfrq), stat=istat)
     if ( istat /= 0 ) then
         call s_print_error('ctqmc_impurity_solver','can not allocate enough memory')
     endif

     allocate(h2_re_mpi(norbs,norbs,nffrq,nffrq,nbfrq), stat=istat)
     if ( istat /= 0 ) then
         call s_print_error('ctqmc_impurity_solver','can not allocate enough memory')
     endif

     allocate(h2_im_mpi(norbs,norbs,nffrq,nffrq,nbfrq), stat=istat)
     if ( istat /= 0 ) then
         call s_print_error('ctqmc_impurity_solver','can not allocate enough memory')
     endif

     allocate(grnf_mpi(mfreq,norbs,norbs), stat=istat)
     if ( istat /= 0 ) then
         call s_print_error('ctqmc_impurity_solver','can not allocate enough memory')
     endif

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
     endif

!=========================================================================
!>>> starting quantum impurity solver                                  <<<
!=========================================================================

! print the header of continuous time quantum Monte Carlo quantum impurity solver
     if ( myid == master ) then ! only master node can do it
         write(mystd,'(2X,a)') 'MANJUSHAKA >>> CTQMC quantum impurity solver running'
         write(mystd,'(4X,a,i10,4X,a,f10.5)') 'nband :', nband, 'Uc    :', Uc
         write(mystd,'(4X,a,i10,4X,a,f10.5)') 'nspin :', nspin, 'Jz    :', Jz
         write(mystd,*)
     endif

!=========================================================================
!>>> initializing quantum impurity solver                              <<<
!=========================================================================

! init the continuous time quantum Monte Carlo quantum impurity solver
! setup the key variables
     if ( myid == master ) then ! only master node can do it
         write(mystd,'(4X,a)') 'quantum impurity solver initializing'
     endif

     call cpu_time(time_begin) ! record starting time
     call ctqmc_solver_init()
     call cpu_time(time_end)   ! record ending   time

! print the time information
     if ( myid == master ) then ! only master node can do it
         write(mystd,'(4X,a,f10.3,a)') 'time:', time_end - time_begin, 's'
         write(mystd,*)
     endif

!=========================================================================
!>>> retrieving quantum impurity solver                                <<<
!=========================================================================

! init the continuous time quantum Monte Carlo quantum impurity solver further
! retrieving the time series information produced by previous running
     if ( myid == master ) then ! only master node can do it
         write(mystd,'(4X,a)') 'quantum impurity solver retrieving'
     endif

     call cpu_time(time_begin) ! record starting time
! for dynamically truncate high energy states, the trace of saved diagramm 
! may be zero, so we don't retrieve it for itrun == 3
     if (itrun == 1 .or. itrun == 2) then
         call ctqmc_retrieve_status()
     endif
     call cpu_time(time_end)   ! record ending   time

! print the time information
     if ( myid == master ) then ! only master node can do it
         write(mystd,'(4X,a,f10.3,a)') 'time:', time_end - time_begin, 's'
         write(mystd,*)
     endif

!=========================================================================
!>>> warmming quantum impurity solver                                  <<<
!=========================================================================

! warmup the continuous time quantum Monte Carlo quantum impurity solver,
! in order to achieve equilibrium state
     if ( myid == master ) then ! only master node can do it
         write(mystd,'(4X,a)') 'quantum impurity solver warmming'
     endif

     call cpu_time(time_begin) ! record starting time
     call ctqmc_diagram_warmming()
     call cpu_time(time_end)   ! record ending   time

! print the time information
     if ( myid == master ) then ! only master node can do it
         write(mystd,'(4X,a,f10.3,a)') 'time:', time_end - time_begin, 's'
         write(mystd,*)
     endif

!=========================================================================
!>>> beginning main iteration                                          <<<
!=========================================================================

! start simulation
     if ( myid == master ) then ! only master node can do it
         write(mystd,'(4X,a)') 'quantum impurity solver sampling'
         write(mystd,*)
     endif

     CTQMC_MAIN_ITERATION: do i=1, nsweep, nwrite

! record start time
         call cpu_time(time_begin)

         CTQMC_DUMP_ITERATION: do j=1, nwrite

!=========================================================================
!>>> sampling perturbation expansion series                            <<<
!=========================================================================

! increase cstep by 1
             cstep = cstep + 1

! sampling the perturbation expansion feynman diagrams randomly
             call ctqmc_diagram_sampling(cstep)

!=========================================================================
!>>> sampling the physical observables                                 <<<
!=========================================================================

! record the histogram for perturbation expansion series
             call ctqmc_record_hist()

! record nothing
             if ( mod(cstep, nmonte) == 0 .and. isvrt == 1 ) then
                 CONTINUE
             endif

! record the spin-spin correlation function
             if ( mod(cstep, nmonte) == 0 .and. isvrt == 2 ) then
                 call ctqmc_record_schi()
             endif

! record the orbital-orbital correlation function
             if ( mod(cstep, nmonte) == 0 .and. isvrt == 3 ) then
                 call ctqmc_record_ochi()
             endif

! record the two-particle green's function
             if ( mod(cstep, nmonte) == 0 .and. isvrt == 4 ) then
                 call ctqmc_record_twop()
             endif

! record the vertex function
             if ( mod(cstep, nmonte) == 0 .and. isvrt == 5 ) then
                 call ctqmc_record_vrtx()
             endif

! record the impurity (double) occupation number matrix and other
! auxiliary physical observables
             if ( mod(cstep, nmonte) == 0 ) then
                 call ctqmc_record_nmat()
             endif

! record the impurity green's function in matsubara frequency space
             if ( mod(cstep, nmonte) == 0 ) then
                 call ctqmc_record_grnf()
             endif

! record the auxiliary correlation function, F^{j}(\tau)
             if ( mod(cstep, ncarlo) == 0 .and. isort >= 4 ) then
                 call ctqmc_record_ftau()
             endif

! record the probability of eigenstates
             if ( mod(cstep, ncarlo) == 0 ) then
                 call ctqmc_record_prob()
             endif

! record the impurity green's function in imaginary time space
             if ( mod(cstep, ncarlo) == 0 ) then
                 call ctqmc_record_gtau()
             endif

         enddo CTQMC_DUMP_ITERATION ! over j={1,nwrite} loop

!=========================================================================
!>>> reporting quantum impurity solver                                 <<<
!=========================================================================

! it is time to write out the statistics results
         if ( myid == master ) then ! only master node can do it
             call ctqmc_print_runtime(iter, cstep)
         endif ! back if ( myid == master ) block

!=========================================================================
!>>> reducing immediate results                                        <<<
!=========================================================================

! collect the histogram data from hist to hist_mpi
         call ctqmc_reduce_hist(hist_mpi)

! collect the impurity green's function data from gtau to gtau_mpi
         call ctqmc_reduce_gtau(gtau_mpi)

! gtau_mpi need to be scaled properly before written
         do m=1,norbs
             do n=1,ntime
                 gtau_mpi(n,m,m) = gtau_mpi(n,m,m) * real(ncarlo) / real(cstep)
             enddo ! over n={1,ntime} loop
         enddo ! over m={1,norbs} loop

!=========================================================================
!>>> symmetrizing immediate results                                    <<<
!=========================================================================

! symmetrize the impurity green's function over spin or over bands
         if ( issun == 2 .or. isspn == 1 ) then
             call ctqmc_symm_gtau(symm, gtau_mpi)
         endif

!=========================================================================
!>>> writing immediate results                                         <<<
!=========================================================================

! write out the histogram data, hist_mpi
         if ( myid == master ) then ! only master node can do it
             call ctqmc_dump_hist(hist_mpi)
         endif

! write out the impurity green's function, gtau_mpi
         if ( myid == master ) then ! only master node can do it
             if ( iter /= 999 ) then
                 call ctqmc_dump_gtau(tmesh, gtau_mpi)
             else
                 call ctqmc_dump_gbin(cstep / nwrite, tmesh, gtau_mpi)
                 write(mystd,'(4X,a)') '>>> quantum impurity solver status: binned'
             endif ! back if ( iter /= 999 ) block
         endif ! back if ( myid == master ) block

!=========================================================================
!>>> checking quantum impurity solver                                  <<<
!=========================================================================

         call ctqmc_diagram_checking(cflag)

!=========================================================================
!>>> timing quantum impurity solver                                    <<<
!=========================================================================

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
         endif

!=========================================================================
!>>> escaping quantum impurity solver                                  <<<
!=========================================================================

! if the quantum impurity solver is out of control or reaches convergence
         if ( cflag == 99 .or. cflag == 100 ) then
             EXIT CTQMC_MAIN_ITERATION ! jump out the iteration
         endif

     enddo CTQMC_MAIN_ITERATION ! over i={1,nsweep} loop

!=========================================================================
!>>> ending main iteration                                             <<<
!=========================================================================

!=========================================================================
!>>> reducing final results                                            <<<
!=========================================================================

! special considerations for prob and grnf, since for children processes,
! caves may be different
     prob = prob / real(caves); grnf = grnf / real(caves)

! collect the histogram data from hist to hist_mpi
     call ctqmc_reduce_hist(hist_mpi)

! collect the spin-spin correlation function data from schi to schi_mpi,
! from sschi to sschi_mpi
     call ctqmc_reduce_schi(schi_mpi, sschi_mpi)

! collect the orbital-orbital correlation function data from ochi to
! ochi_mpi, from oochi to oochi_mpi
     call ctqmc_reduce_ochi(ochi_mpi, oochi_mpi)

! collect the (double) occupation matrix data from nmat to nmat_mpi, from
! nnmat to nnmat_mpi
     call ctqmc_reduce_nmat(nmat_mpi, nnmat_mpi)

! collect the probability data from prob to prob_mpi
     call ctqmc_reduce_prob(prob_mpi)

! collect the impurity green's function data from gtau to gtau_mpi
     call ctqmc_reduce_gtau(gtau_mpi)

! collect the impurity green's function data from grnf to grnf_mpi
     call ctqmc_reduce_grnf(grnf_mpi)

! collect the auxiliary correlation function from ftau to ftau_mpi
     call ctqmc_reduce_ftau(ftau_mpi)

! collect the two-particle green's function from g2_re to g2_re_mpi, etc.
     call ctqmc_reduce_twop(g2_re_mpi, g2_im_mpi)

! collect the vertex function from h2_re to h2_re_mpi, etc.
     call ctqmc_reduce_vrtx(h2_re_mpi, h2_im_mpi)

! update original data and calculate the averages simultaneously
     hist  = hist_mpi

     schi  = schi_mpi  * real(nmonte) / real(nsweep)
     ochi  = ochi_mpi  * real(nmonte) / real(nsweep)
     nmat  = nmat_mpi  * real(nmonte) / real(nsweep)
     prob  = prob_mpi  * real(ncarlo)

     do m=1,nband
         do n=1,ntime
             sschi(n,m) = sschi_mpi(n,m)   * real(nmonte) / real(nsweep)
         enddo ! over n={1,ntime} loop
     enddo ! over m={1,nband} loop

     do m=1,norbs
         do n=1,ntime
             oochi(n,m) = oochi_mpi(n,m)   * real(nmonte) / real(nsweep)
         enddo ! over n={1,ntime} loop
     enddo ! over m={1,norbs} loop

     do m=1,norbs
         do n=1,norbs
             nnmat(n,m) = nnmat_mpi(n,m)   * real(nmonte) / real(nsweep)
         enddo ! over n={1,norbs} loop
     enddo ! over m={1,norbs} loop

     do m=1,norbs
         do n=1,ntime
             gtau(n,m,m) = gtau_mpi(n,m,m) * real(ncarlo) / real(nsweep)
         enddo ! over n={1,ntime} loop
     enddo ! over m={1,norbs} loop

     do m=1,norbs
         do n=1,nfreq
             grnf(n,m,m) = grnf_mpi(n,m,m) * real(nmonte)
         enddo ! over n={1,nfreq} loop
     enddo ! over m={1,norbs} loop

     do m=1,norbs
         do n=1,ntime
             ftau(n,m,:) = ftau_mpi(n,m,:) * real(ncarlo) / real(nsweep)
         enddo ! over n={1,ntime} loop
     enddo ! over m={1,norbs} loop

     do m=1,norbs
         do n=1,norbs
             g2_re(n,m,:,:,:) = g2_re_mpi(n,m,:,:,:) * real(nmonte) / real(nsweep)
             g2_im(n,m,:,:,:) = g2_im_mpi(n,m,:,:,:) * real(nmonte) / real(nsweep)
         enddo ! over n={1,norbs} loop
     enddo ! over m={1,norbs} loop

     do m=1,norbs
         do n=1,norbs
             h2_re(n,m,:,:,:) = h2_re_mpi(n,m,:,:,:) * real(nmonte) / real(nsweep)
             h2_im(n,m,:,:,:) = h2_im_mpi(n,m,:,:,:) * real(nmonte) / real(nsweep)
         enddo ! over n={1,norbs} loop
     enddo ! over m={1,norbs} loop

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
         call ctqmc_make_hub1() ! call ctqmc_make_hub2() ! not implemented
     endif ! back if ( isort <= 3 ) block

!=========================================================================
!>>> symmetrizing final results                                        <<<
!=========================================================================

! symmetrize the occupation number matrix (nmat) over spin or over bands
     if ( issun == 2 .or. isspn == 1 ) then
         call ctqmc_symm_nmat(symm, nmat)
     endif

! symmetrize the impurity green's function (gtau) over spin or over bands
     if ( issun == 2 .or. isspn == 1 ) then
         call ctqmc_symm_gtau(symm, gtau)
     endif

! symmetrize the impurity green's function (grnf) over spin or over bands
     if ( issun == 2 .or. isspn == 1 ) then
         call ctqmc_symm_grnf(symm, grnf)
     endif

! symmetrize the impurity self-energy function (sig2) over spin or over bands
     if ( issun == 2 .or. isspn == 1 ) then
         call ctqmc_symm_grnf(symm, sig2)
     endif

!=========================================================================
!>>> writing final results                                             <<<
!=========================================================================

! write out the final histogram data, hist
     if ( myid == master ) then ! only master node can do it
         call ctqmc_dump_hist(hist)
     endif

! write out the final spin-spin correlation function data, schi and sschi
     if ( myid == master ) then ! only master node can do it
         call ctqmc_dump_schi(schi, sschi)
     endif

! write out the final orbital-orbital correlation function data, ochi and oochi
     if ( myid == master ) then ! only master node can do it
         call ctqmc_dump_ochi(ochi, oochi)
     endif

! write out the final (double) occupation matrix data, nmat and nnmat
     if ( myid == master ) then ! only master node can do it
         call ctqmc_dump_nmat(nmat, nnmat)
     endif

! write out the final probability data, prob
     if ( myid == master ) then ! only master node can do it
         call ctqmc_dump_prob(prob, naux, saux)
     endif

! write out the final probability data of sectors
     if ( myid == master ) then ! only master node can do it
         call ctqmc_dump_psect()
     endif

! write out the final impurity green's function data, gtau
     if ( myid == master ) then ! only master node can do it
         call ctqmc_dump_gtau(tmesh, gtau)
     endif

! write out the final impurity green's function data, grnf
     if ( myid == master ) then ! only master node can do it
         call ctqmc_dump_grnf(rmesh, grnf)
     endif

! write out the final self-energy function data, sig2
     if ( myid == master ) then ! only master node can do it
         call ctqmc_dump_sigf(rmesh, sig2)
     endif

! write out the final two-particle green's function data, g2_re and g2_im
     if ( myid == master ) then ! only master node can do it
         call ctqmc_dump_twop(g2_re, g2_im)
     endif

! write out the final vertex function data, h2_re and h2_im
     if ( myid == master ) then ! only master node can do it
         call ctqmc_dump_vrtx(h2_re, h2_im)
     endif

!=========================================================================
!>>> saving quantum impurity solver                                    <<<
!=========================================================================

! save the perturbation expansion series information to the disk file
     if ( myid == master ) then ! only master node can do it
         call ctqmc_save_status()
     endif

!=========================================================================
!>>> finishing quantum impurity solver                                 <<<
!=========================================================================

! print the footer of continuous time quantum Monte Carlo quantum impurity solver
     if ( myid == master ) then ! only master node can do it
         write(mystd,'(2X,a)') 'MANJUSHAKA >>> CTQMC quantum impurity solver shutdown'
         write(mystd,*)
     endif

! deallocate memory for occu and npart
     call ctqmc_deallocate_memory_occu()
     call ctqmc_deallocate_memory_part()

! deallocate memory
     deallocate(hist_mpi)
     deallocate(schi_mpi)
     deallocate(ochi_mpi)
     deallocate(nmat_mpi)
     deallocate(prob_mpi)
     deallocate(gtau_mpi)
     deallocate(ftau_mpi)
     deallocate(grnf_mpi)
     deallocate(sschi_mpi)
     deallocate(oochi_mpi)
     deallocate(nnmat_mpi)
     deallocate(g2_re_mpi)
     deallocate(g2_im_mpi)
     deallocate(h2_re_mpi)
     deallocate(h2_im_mpi)

     return
  end subroutine ctqmc_impurity_solver

!>>> perform thermalization on perturbation expansion series to achieve
! thermodynamics equilibrium state
  subroutine ctqmc_diagram_warmming()
     use constants, only : zero
     use control, only : ntherm
     use context

     implicit none

! local variables
! loop index
     integer :: i

! warm up the diagram series
     do i=1,ntherm
         call ctqmc_diagram_sampling(i)
     enddo ! over i={1,ntherm} loop

! reset cnegs
     cnegs = 0

! reset caves
     caves = 0

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

     reflip_tcount = zero
     reflip_accept = zero
     reflip_reject = zero

     return
  end subroutine ctqmc_diagram_warmming

!>>> visit the perturbation expansion diagrams randomly
  subroutine ctqmc_diagram_sampling(cstep)
     use constants, only : dp
     use control, only : nflip, nclean

     use spring

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
         endif
! do not change the order of perturbation expansion series
     else
         if ( spring_sfmt_stream() > 0.5_dp ) then
             call ctqmc_lshift_kink()  ! shift the create  operators
         else
             call ctqmc_rshift_kink()  ! shift the destroy operators
         endif
     endif ! back if ( spring_sfmt_stream() < 0.9_dp ) block

! numerical trick: perform global spin flip periodically
     if ( nflip > 0  .and. mod(cstep, +nflip) == 0 ) then
         if ( spring_sfmt_stream() < 0.8_dp ) then
             call ctqmc_reflip_kink(2) ! flip intra-orbital spins one by one
         else
             call ctqmc_reflip_kink(3) ! flip intra-orbital spins globally
         endif
     endif

     if ( nflip < 0  .and. mod(cstep, -nflip) == 0 ) then
         if ( spring_sfmt_stream() < 0.8_dp ) then
             call ctqmc_reflip_kink(1) ! flip inter-orbital spins randomly
         else
             call ctqmc_reflip_kink(3) ! flip intra-orbital spins globally
         endif
     endif

! numerical trick: perform global update periodically
     if ( nclean > 0 .and. mod(cstep, nclean) == 0 ) then
         call ctqmc_reload_kink()
     endif

     return
  end subroutine ctqmc_diagram_sampling

!>>> checking whether the quantum impurity solver is consistent internally
  subroutine ctqmc_diagram_checking(cflag)
     use constants
     use control
     use context

     implicit none

! external arguments
! control flag
     integer, intent(inout) :: cflag

! local variables
! loop index
     integer :: i
     integer :: j

     if ( cflag == 1 ) then

! check time order of operators in colour part
         do i=1,norbs
             do j=1,rank(i)-1
                 if ( time_s( index_s(j, i), i ) > time_s( index_s(j+1, i), i ) ) then
                     cflag = 99
                 endif
                 if ( time_e( index_e(j, i), i ) > time_e( index_e(j+1, i), i ) ) then
                     cflag = 99
                 endif
             enddo ! over j={1,rank(i)-1} loop
         enddo ! over i={1,norbs} loop

! check time order of operators in flavor part
         do j=1,2*sum(rank)-1
             if ( index_v(j) <= 0 .or. index_v(j+1) <= 0 ) then
                 cflag = 99
             endif
         enddo ! over j={1,2*sum(rank)-1} loop

         do j=1,2*sum(rank)-1
             if ( time_v( index_v(j) ) > time_v( index_v(j+1) ) ) then
                 cflag = 99
             endif
         enddo ! over j={1,2*sum(rank)-1} loop

! write the results, only master node can do it
         if ( myid == master ) then
             if ( cflag == 99 ) then
                 write(mystd,'(4X,a)') '>>> quantum impurity solver status: error?'
                 write(mystd,'(4X,a)') '>>> please check the status file: solver.status.dat'
                 call ctqmc_save_status()
                 call s_print_error('ctqmc_diagram_checking','unknown fatal error occur')
             else
                 write(mystd,'(4X,a)') '>>> quantum impurity solver status: normal'
             endif
         endif ! back if ( myid == master ) block

     endif ! back if ( cflag == 1 ) block

     return
  end subroutine ctqmc_diagram_checking

!>>> testing subroutine, please active it on ctqmc_diagram_sampling()
  subroutine ctqmc_impurity_tester()
     use constants
     use control
     use context

     implicit none

!-------------------------------------------------------------------------
! insert your debug code here
!-------------------------------------------------------------------------

     call ctqmc_make_display(1)
     call ctqmc_make_display(2)
     call s_print_error('ctqmc_impurity_tester','in debug mode')

     return
  end subroutine ctqmc_impurity_tester
