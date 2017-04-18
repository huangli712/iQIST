!!!=========+=========+=========+=========+=========+=========+=========+!
!!! NARCISSUS @ iQIST                                                    !
!!!                                                                      !
!!! A test program for dynamical mean field theory (DMFT) self-consistent!
!!! engine plus hybridization expansion version continuous time quantum  !
!!! Monte Carlo (CTQMC) quantum impurity solver                          !
!!! author  : Li Huang (at IOP/CAS & SPCLab/CAEP & UNIFR)                !
!!! status  : WARNING: IN TESTING STAGE, USE IT IN YOUR RISK             !
!!! comment : this impurity solver is based on segment picture formalism !
!!!           the dynamical screening effect is included                 !
!!!           any question, please contact with lihuang.dmft@gmail.com   !
!!!=========+=========+=========+=========+=========+=========+=========+!

!!========================================================================
  PROGRAM CTQMC_MAIN !                                                 <<<
!!========================================================================

     use mmpi, only : mp_init       ! init mpi environment
     use mmpi, only : mp_finalize   ! finalize mpi environment
     use mmpi, only : mp_barrier    ! barrier to synchronize the data
     use mmpi, only : mp_comm_rank  ! get index of current process
     use mmpi, only : mp_comm_size  ! get number of processes

     use control, only : isscf      ! self-consistent calculation mode
     use control, only : niter      ! number of self-consistent iteration
     use control, only : nprocs     ! number of processes
     use control, only : myid       ! index of current process
     use control, only : master     ! index of master process

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
     endif ! back if ( myid == master ) block

! setup the parameters for continuous time quantum Monte Carlo quantum
! impurity solver and dynamical mean field theory self-consistent engine
     call ctqmc_config()

! allocate memory and initialize
     call ctqmc_setup_array()

! prepare initial hybridization function, init self-consistent iteration
     call ctqmc_selfer_init()

! print out runtime parameters in summary, only for check
     if ( myid == master ) then ! only master node can do it
         call ctqmc_print_summary()
     endif ! back if ( myid == master ) block

!!========================================================================
!!>>> DMFT ITERATION BEGIN                                             <<<
!!========================================================================

!-------------------------------------------------------------------------
! case A: one-shot non-self-consistent mode
!-------------------------------------------------------------------------
! it is suitable for local density approximation plus dynamical mean field
! theory calculation
!-------------------------------------------------------------------------
     if ( isscf == 1 ) then

! set the iter number
         iter = niter

! write the iter to screen
         if ( myid == master ) then ! only master node can do it
             call ctqmc_print_it_info(iter)
         endif ! back if ( myid == master ) block

! call the continuous time quantum Monte Carlo quantum impurity solver, to
! build the impurity green's function and self-energy function
         call ctqmc_impurity_solver(iter)

     endif ! back if ( isscf == 1 ) block

!-------------------------------------------------------------------------
! case B: self-consistent mode
!-------------------------------------------------------------------------
! it is suitable for lattice model hamiltonian plus dynamical mean field
! theory calculation
!-------------------------------------------------------------------------
     if ( isscf == 2 ) then

     DMFT_CTQMC_ITERATION: &
     do iter=1,niter

! write the iter to screen
         if ( myid == master ) then ! only master node can do it
             call ctqmc_print_it_info(iter)
         endif ! back if ( myid == master ) block

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
         endif ! back if ( convergence .eqv. .true. ) block

     enddo DMFT_CTQMC_ITERATION ! over iter={1,niter} loop

     endif ! back if ( isscf == 2 ) block

!!========================================================================
!!>>> DMFT ITERATION END                                               <<<
!!========================================================================

! deallocate memory and finalize
     call ctqmc_final_array()

! print the footer for continuous time quantum Monte Carlo quantum impurity
! solver and dynamical mean field theory self-consistent engine
     if ( myid == master ) then ! only master node can do it
         call ctqmc_print_footer()
     endif ! back if ( myid == master ) block

! finalize mpi envirnoment
# if defined (MPI)

! blocks until all processes have reached this routine
     call mp_barrier()

! terminates mpi execution environment
     call mp_finalize()

# endif  /* MPI */

!!========================================================================
  END PROGRAM CTQMC_MAIN !                                             <<<
!!========================================================================
