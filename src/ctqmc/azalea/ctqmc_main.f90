!=========+=========+=========+=========+=========+=========+=========+>>>
! A test program for dynamical mean field theory (DMFT) self-consistent  !
! engine plus hybridization expansion version continuous time quantum    !
! Monte Carlo (CTQMC) quantum impurity solver                            !
! author  : li huang                                                     !
! version : v2012.08.20T                                                 !
! status  : WARNING: IN TESTING STAGE, USE IT IN YOUR RISK               !
! comment : this impurity solver is based on segment picture formalism   !
!           any question, please contact with huangli712@yahoo.com.cn    !
!=========+=========+=========+=========+=========+=========+=========+>>>

  program ctqmc_main
     use constants
     use control

     use mmpi

     implicit none

! local variables
! loop index
     integer :: iter

! convergence flag
     logical :: convergence

! initialize mpi envirnoment
# if defined (MPI)

! initialize the mpi execution environment
     call mp_init()

! determines the rank of the calling process in the communicator
     call mp_comm_rank(myid)

! determines the size of the group associated with a communicator
     call mp_comm_size(nprocs)

# endif  /* MPI */

! print the running header for continuous time quantum Monte Carlo quantum
! impurity solver and dynamical mean field theory self-consistent engine
     if ( myid == master ) then ! only master node can do it
         call ctqmc_print_header()
     endif

! setup the important parameters for continuous time quantum Monte Carlo
! quantum impurity solver and dynamical mean field theory self-consistent
! engine
     call ctqmc_config()

! print out runtime parameters in summary, only for check
     if ( myid == master ) then ! only master node can do it
         call ctqmc_print_summary()
     endif

! allocate memory and initialize
     call ctqmc_setup_array()

! prepare initial hybridization function, init self-consistent iteration
     call ctqmc_selfer_init()

!-------------------------------------------------------------------------
! note: running mode                                                     !
!-------------------------------------------------------------------------
!    if isscf == 1 .and. isbin == 1                                      !
!        call ctqmc_impurity_solver only, normal mode                    !
!                                                                        !
!    if isscf == 1 .and. isbin == 2                                      !
!        call ctqmc_impurity_solver only, binner mode                    !
!                                                                        !
!    if isscf == 2 .and. isbin == 1                                      !
!        call ctqmc_impurity_solver, normal mode                         !
!        plus                                                            !
!        call ctqmc_dmft_selfer                                          !
!        until convergence                                               !
!                                                                        !
!    if isscf == 2 .and. isbin == 2                                      !
!        call ctqmc_impurity_solver, normal mode                         !
!        plus                                                            !
!        call ctqmc_dmft_selfer                                          !
!        until convergence                                               !
!        plus                                                            !
!        call ctqmc_impurity_solver, binner mode                         !
!-------------------------------------------------------------------------

!=========================================================================
!>>> DMFT ITERATION BEGIN                                              <<<
!=========================================================================

! case A: one-shot non-self-consistent mode
!-------------------------------------------------------------------------
! it is suitable for local density approximation plus dynamical mean field
! theory calculation
     if ( isscf == 1 .and. isbin == 1 ) then

! set the iter number
         iter = niter

! write the iter to screen
         if ( myid == master ) then ! only master node can do it
             write(mystd,'(2X,a,i3,a)') 'AZALEA >>> DMFT iter:', iter, ' <<< SELFING'
         endif

! call the continuous time quantum Monte Carlo quantum impurity solver, to
! build the impurity green's function and self-energy function
         call ctqmc_impurity_solver(iter)

     endif ! back if ( isscf == 1 .and. isbin == 1 ) block

! case B: self-consistent mode
!-------------------------------------------------------------------------
! it is suitable for lattice model hamiltonian plus dynamical mean field
! theory calculation
     DMFT_CTQMC_ITERATION: do iter=1,niter

! check the running mode
         if ( isscf == 1 ) then
             EXIT DMFT_CTQMC_ITERATION ! jump out the iteration
         endif

! write the iter to screen
         if ( myid == master ) then ! only master node can do it
             write(mystd,'(2X,a,i3,a)') 'AZALEA >>> DMFT iter:', iter, ' <<< SELFING'
         endif

! call the continuous time quantum Monte Carlo quantum impurity solver, to
! build the impurity green's function and self-energy function
         call ctqmc_impurity_solver(iter)

! call the self-consistent engine for dynamical mean field theory, to build
! the bath weiss's function and hybridization function
         call ctqmc_dmft_selfer()

! check convergence for dynamical mean field theory iteration
         convergence = .false.
         call ctqmc_dmft_conver(iter, convergence)

! now convergence is achieved
         if ( convergence .eqv. .true. ) then
             EXIT DMFT_CTQMC_ITERATION ! jump out the iteration
         endif

     enddo DMFT_CTQMC_ITERATION ! over iter={1,niter} loop

! case C: binner mode
!-------------------------------------------------------------------------
! perform quantum Monte Carlo data binning
     if ( isbin == 2 ) then

! set the iter number
         iter = 999

! write the iter to screen
         if ( myid == master ) then ! only master node can do it
             write(mystd,'(2X,a,i3,a)') 'AZALEA >>> DMFT iter:', iter, ' <<< BINNING'
         endif

! accumulate the quantum Monte Carlo data
         call ctqmc_impurity_solver(iter)

     endif ! back if ( isbin == 2 ) block

!=========================================================================
!>>> DMFT ITERATION END                                                <<<
!=========================================================================

! deallocate memory and finalize
     call ctqmc_final_array()

! print the footer for continuous time quantum Monte Carlo quantum impurity
! solver and dynamical mean field theory self-consistent engine
     if ( myid == master ) then ! only master node can do it
         call ctqmc_print_footer()
     endif

! finalize mpi envirnoment
# if defined (MPI)

! blocks until all processes have reached this routine
     call mp_barrier()

! terminates mpi execution environment
     call mp_finalize()

# endif  /* MPI */

  end program ctqmc_main
