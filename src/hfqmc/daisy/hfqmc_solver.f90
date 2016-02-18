!!!-----------------------------------------------------------------------
!!! project : daisy
!!! program : hfqmc_impurity_solver
!!! source  : hfqmc_solver.f90
!!! type    : subroutine
!!! author  : li huang (email:lihuang.dmft@gmail.com)
!!! history : 01/05/2006 by li huang (created)
!!!           08/17/2015 by li huang (last modified)
!!! purpose : it is a mature parallel implemention of the famous Hirsch
!!!           -Fye quantum Monte Carlo algorithm, which is used to solve
!!!           the multi-orbital Anderson impurity model as usually. we
!!!           denote it as HFQMC (Hirsch-Fye quantum Monte Carlo) quantum
!!!           impurity solver.
!!! status  : unstable
!!! comment : this subroutine is fully parallelism by mpi
!!!-----------------------------------------------------------------------

!!>>> hfqmc_impurity_solver: solve the Anderson impurity model using the
!!>>> Hirsch-Fye quantum Monte Carlo algorithm
  subroutine hfqmc_impurity_solver(iter)
     use constants, only : dp, zero, one, mystd
     use spring, only : spring_sfmt_stream

     use control, only : cname
     use control, only : isscf
     use control, only : norbs
     use control, only : mstep
     use control, only : nsing, ntime, ntherm, nsweep, ncarlo
     use control, only : myid, master
     use context, only : ktep, diag, atep, btep
     use context, only : gmat, wmat
     use context, only : symm, nmat, tmesh, rmesh, nnmat
     use context, only : gtau, wtau, grnf, wssf, sig2

     implicit none

! external arguments
! current iteration number
     integer, intent(in) :: iter

! local variables
! loop index
     integer  :: i
     integer  :: j
     integer  :: m
     integer  :: n

! loop index, current QMC sweep count
     integer  :: cstep

! current QMC effective sweep count
     integer  :: nstep

! internal fast cycle numbers, used to control the size of data bins
     integer  :: nfast

! status flag
     integer  :: istat

! accepted QMC flip count
     real(dp) :: accept

! rejected QMC flip count
     real(dp) :: reject

! total QMC flip count
     real(dp) :: tcount

! transition probability between two successive ising configuration
     real(dp) :: markov

! starting time
     real(dp) :: time_begin

! ending time
     real(dp) :: time_end

! time consuming by current iteration
     real(dp) :: time_iter

! time consuming by total iteration
     real(dp) :: time_niter

! real(dp) dummy matrix for pack green's function
     real(dp), allocatable :: msum(:)

! double occupation number, for mpi case
     real(dp), allocatable :: osum(:,:)
     real(dp), allocatable :: oerr(:,:)

! sum of green's function over numprocs processes, for mpi case
     real(dp), allocatable :: gsum(:,:)
     real(dp), allocatable :: gerr(:,:)

! the summed up green's function matrix
     real(dp), allocatable :: tmat(:,:,:)

! allocate memory
     allocate(msum(-ntime+1:ntime-1),  stat=istat)
     if ( istat /= 0 ) then
         call s_print_error('hfqmc_impurity_solver','can not allocate enough memory')
     endif ! back if ( istat /= 0 ) block

! allocate memory
     allocate(osum(norbs,norbs),       stat=istat)
     allocate(oerr(norbs,norbs),       stat=istat)
     if ( istat /= 0 ) then
         call s_print_error('hfqmc_impurity_solver','can not allocate enough memory')
     endif ! back if ( istat /= 0 ) block

! allocate memory
     allocate(gsum(ntime,norbs),       stat=istat)
     allocate(gerr(ntime,norbs),       stat=istat)
     if ( istat /= 0 ) then
         call s_print_error('hfqmc_impurity_solver','can not allocate enough memory')
     endif ! back if ( istat /= 0 ) block

! allocate memory
     allocate(tmat(ntime,ntime,norbs), stat=istat)
     if ( istat /= 0 ) then
         call s_print_error('hfqmc_impurity_solver','can not allocate enough memory')
     endif ! back if ( istat /= 0 ) block

! setup timer
     time_iter = zero
     time_niter = zero

! setup nsweep
! whether it is time to enter QMC data accumulating mode
     if ( iter == 999 ) then
         nsweep = nsweep * 10
     endif ! back if ( iter == 999 ) block

! setup internal iteration parameters
! we need 10 data bins
     nfast = nsweep / 10

! setup statistics counters
     accept = zero
     reject = zero
     tcount = zero

! setup impurity green's function and bath weiss's function matrices
     gmat = zero  ! reset global variables
     wmat = zero  ! reset global variables

     tmat = zero

     msum = zero
     gsum = zero

! setup matrix used by delayed update algorithm
     ktep = 0     ! reset global variables
     atep = zero  ! reset global variables
     btep = zero  ! reset global variables

! setup double occupation number matrix
     osum = zero
     nnmat = zero ! reset global variables

! call hfqmc_solver_init() subroutine to initialize necessary matrices
! note: only the first iteration has the chance to do it. or else the
! imat matrix should be reseted after it is called, and the ising-like
! fields configuration are lost.
     if ( iter == 1 .or. isscf == 1 ) then
         call hfqmc_solver_init() ! init solver-related matrix
!<         ntherm = 100             ! adjust ntherm here for iter = 1
     else
!<         ntherm = 30              ! adjust ntherm here for iter > 1
     endif ! back if ( iter == 1 .or. isscf == 1 ) block

! setup wmat, the bath weiss's function matrix
! the relation between weiss's function and wmat matrix is:
!    wtau($\tau-\tau'$,norbs) ---> wmat($\tau$,$\tau'$,norbs)
     do i=1,norbs
         do j=0,ntime-1
             msum(j) = wtau(j+1,i)
         enddo ! over j={0,ntime-1} loop

         do j=1,ntime-1
             msum(-j) = -msum(ntime-j)
         enddo ! over j={1,ntime-1} loop

         do m=1,ntime
             do n=1,ntime
                 wmat(n,m,i) = msum(n-m)
             enddo ! over n={1,ntime} loop
         enddo ! over m={1,ntime} loop
     enddo ! over i={1,norbs} loop

! setup gmat, the impurity green's function matrix, from wmat matrix using
! cat_clean_update() subroutine. it is a time-consuming subroutine which
! updates green's function from original weiss's function.
! the relation between green's function and gmat matrix is:
!    gtau($\tau-\tau'$,norbs) ---> gmat($\tau$,$\tau'$,norbs)
     do i=1,norbs
         call cat_clean_update(i)
     enddo ! over i={1,norbs} loop

! extract diagonal elements of gmat matrix for delayed update algorithm
     if ( mstep > 1 ) then
         do i=1,norbs
             do j=1,ntime
                 diag(j,i) = gmat(j,j,i)
             enddo ! over j={1,ntime} loop
         enddo ! over i={1,norbs} loop
     endif ! back if ( mstep > 1 ) block

! print main iteration information, only for master node
     if ( myid == master ) then
         write(mystd,'(2X,a)') cname//' >>> HFQMC quantum impurity solver running'
     endif ! back if ( myid == master ) block

!!========================================================================
!!>>> begin main iteration                                             <<<
!!========================================================================

     HFQMC_MAIN_ITERATION: do cstep = 1, nsweep + ntherm

! determine nstep, effective QMC sweep steps
         if ( cstep > ntherm ) then
             nstep = cstep - ntherm
         else
             nstep = 0
         endif ! back if ( cstep > ntherm ) block

! record start time
         if ( mod(nstep, nfast) == 1 ) then
             call cpu_time(time_begin)
         endif ! back if ( mod(nstep, nfast) == 1 ) block

!!========================================================================
!!>>> begin fast loop                                                  <<<
!!========================================================================

         HFQMC_ISING_ITERATION: do m=1,nsing
             HFQMC_SLICE_ITERATION: do n=1,ntime

! evaluate the transition ratio between two ising-like fields configurations
                 call hfqmc_make_detrat(n, m, markov)

! sample the auxiliary ising-like fields
! reject spin-flip update
                 if ( markov < spring_sfmt_stream() ) then

! increase reject 1
                     if ( nstep > 0 ) then
                         reject = reject + one
                         tcount = tcount + one
                     endif ! back if ( nstep > 0 ) block

! accept spin-flip update
                 else

! increase accept 1
                     if ( nstep > 0 ) then
                         accept = accept + one
                         tcount = tcount + one
                     endif ! back if ( nstep > 0 ) block

! update the auxiliary ising-like fields and green's function matrix
                     call hfqmc_make_accept(n, m, cstep)

                 endif ! back if ( markov < spring_sfmt_stream() ) block

             enddo HFQMC_SLICE_ITERATION ! over n={1,ntime} loop
         enddo HFQMC_ISING_ITERATION ! over m={1,nsing} loop

!!========================================================================
!!>>> end fast loop                                                    <<<
!!========================================================================

! reporting quantum impurity solver
!-------------------------------------------------------------------------
! print out QMC trace information, only for master node
         if ( mod(nstep, nfast) == 0 .and. nstep > 0 .and. myid == master ) then
             call hfqmc_print_runtime(iter, nstep, accept, reject, tcount)
         endif ! back if ( mod(nstep, nfast) == 0 .and. nstep > 0 .and. myid == master ) block

! sampling the physical observables
!-------------------------------------------------------------------------
! measure the green's function
! note: here for simpility we only add up to the tmat matrix. in order to
! save the computational time and elimintate the correlation effects between
! two sussceive measurements, we make one measurement per ncarlo cycle.
! noted by li huang, 2007/03/15
         if ( mod(nstep, ncarlo) == 0 .and. nstep > 0 ) then
! if delayed update algorithm is used, it is important to update the status
! of cyclic variables at first
             if ( mstep > 1 ) then
                 call cat_clear_update()
             endif ! back if ( mstep > 1 ) block

             do i=1,norbs
                 do m=1,ntime
                     do n=1,ntime
                         tmat(n,m,i) = tmat(n,m,i) + gmat(n,m,i)
                     enddo ! over n={1,ntime} loop
                 enddo  ! over m={1,ntime} loop
             enddo ! over i={1,norbs} loop
         endif ! back if ( mod(nstep, ncarlo) == 0 .and. nstep > 0 ) block

! measure double occupation number
         if ( mod(nstep, ncarlo) == 0 .and. nstep > 0 ) then
! if delayed update algorithm is used, it is important to update the status
! of cyclic variables at first
             if ( mstep > 1 ) then
                 call cat_clear_update()
             endif ! back if ( mstep > 1 ) block

             do i=1,norbs-1
                 do j=i+1,norbs
                     do m=1,ntime
                         nnmat(i,j) = nnmat(i,j) + ( one - gmat(m,m,i) ) * ( one - gmat(m,m,j) )
                         nnmat(j,i) = nnmat(j,i) + ( one - gmat(m,m,j) ) * ( one - gmat(m,m,i) )
                     enddo ! over m={1,ntime} loop
                 enddo ! over j={i+1,norbs} loop
             enddo ! over i={1,norbs-1} loop
         endif ! back if ( mod(nstep, ncarlo) == 0 .and. nstep > 0 ) block

! reducing immediate results
!-------------------------------------------------------------------------
         if ( mod(nstep, nfast) == 0 .and. nstep > 0 ) then

! build impurity green's function: gtau
             gtau = zero
             do m=1,norbs
                 call hfqmc_make_vertex(m, tmat, msum)
                 do n=1,ntime
                     gtau(n,m) = gtau(n,m) + msum(n-1)
                 enddo ! over n={1,ntime} loop
             enddo ! over m={1,norbs} loop

! collect the impurity green's function data from gtau to gsum
             call hfqmc_reduce_gtau(gsum, gerr)

         endif ! back if ( mod(nstep, nfast) == 0  .and. nstep > 0 ) block

! symmetrizing immediate results
!-------------------------------------------------------------------------
         if ( mod(nstep, nfast) == 0 .and. nstep > 0 ) then

! to deal with and finalize the impurity green's function
             call hfqmc_make_symm(symm, gsum)
             call hfqmc_make_symm(symm, gerr)

! gsum need to be scaled properly before written
             do m=1,norbs
                 do n=1,ntime
                     gsum(n,m) = gsum(n,m) * real(ncarlo) / real(nstep)
                     gerr(n,m) = gerr(n,m) * real(ncarlo) / real(nstep)
                 enddo ! over n={1,ntime} loop
             enddo ! over m={1,norbs} loop

         endif ! back if ( mod(nstep, nfast) == 0  .and. nstep > 0 ) block

! writing immediate results
!-------------------------------------------------------------------------
         if ( mod(nstep, nfast) == 0 .and. nstep > 0 ) then

! write out the impurity green's function, gsum
             if ( myid == master ) then ! only master node can do it
                 call hfqmc_dump_gtau(tmesh, gsum, gerr)
             endif ! back if ( myid == master ) block

         endif ! back if ( mod(nstep, nfast) == 0  .and. nstep > 0 ) block

! timing quantum impurity solver
!-------------------------------------------------------------------------
         if ( mod(nstep, nfast) == 0 .and. nstep > 0 ) then

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

         endif ! back if ( mod(nstep, nfast) == 0  .and. nstep > 0 ) block

     enddo HFQMC_MAIN_ITERATION ! over main loop

!!========================================================================
!!>>> end main iteration                                               <<<
!!========================================================================

! collect the (double) occupation matrix data from nnmat to osum
     call hfqmc_reduce_nmat(osum, oerr)

! collect the impurity green's function data from gtau to gsum
     call hfqmc_reduce_gtau(gsum, gerr)

! symmetrize the impurity green's function (gsum) over spin or over bands
     call hfqmc_make_symm(symm, gsum)
     call hfqmc_make_symm(symm, gerr)

! update original data and calculate the averages simultaneously
     do m=1,norbs
         do n=1,ntime
             gtau(n,m) = gsum(n,m) * real(ncarlo) / real(nsweep)
             gerr(n,m) = gerr(n,m) * real(ncarlo) / real(nsweep)
         enddo ! over n={1,ntime} loop
     enddo ! over m={1,norbs} loop

     do m=1,norbs
         do n=1,norbs
             nnmat(n,m) = osum(n,m) * real(ncarlo) / real( nsweep * ntime )
              oerr(n,m) = oerr(n,m) * real(ncarlo) / real( nsweep * ntime )
         enddo ! over n={1,norbs} loop
     enddo ! over m={1,norbs} loop

! calculate final impurity green's function, bath weiss's function, and
! self energy function at matsubara space by self-consistent equation
     call hfqmc_make_freq()

! calculate impurity occupation number
     call hfqmc_make_nmat()

! write out the final green's function to file, only for master node
     if ( myid == master ) then
         call hfqmc_dump_gtau(tmesh, gtau, gerr)
     endif ! back if ( myid == master ) block

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

! write out occupation number to disk file
     if ( myid == master ) then ! only master node can do it
         call hfqmc_dump_nmat(nmat, nnmat, gerr(1,:), oerr)
     endif ! back if ( myid == master ) block

! print the footer of Hirsch-Fye quantum Monte Carlo quantum impurity solver
     if ( myid == master ) then ! only master node can do it
         write(mystd,'(2X,a)') cname//' >>> HFQMC quantum impurity solver shutdown'
         write(mystd,*)
     endif ! back if ( myid == master ) block

! deallocate memory
     deallocate(msum)
     deallocate(osum)
     deallocate(oerr)
     deallocate(gsum)
     deallocate(gerr)

     deallocate(tmat)

     return
  end subroutine hfqmc_impurity_solver
