!!!=========+=========+=========+=========+=========+=========+=========+!
!!! iQIST @ NARCISSUS                                                    !
!!!                                                                      !
!!! A highly optimized hybridization expansion version continuous time   !
!!! quantum Monte Carlo (CTQMC) quantum impurity solver plus a classic   !
!!! dynamical mean field theory (DMFT) self-consistent engine            !
!!!                                                                      !
!!! author  : Li Huang (China Academy of Engineering Physics)            !
!!! status  : (WARNING) IN TESTING STAGE, USE IT IN YOUR RISK            !
!!! comment : this impurity solver is based on segment picture formalism !
!!!           any question, please contact with huangli@caep.cn          !
!!!=========+=========+=========+=========+=========+=========+=========+!

!!========================================================================
  PROGRAM CTQMC_MAIN !                                                 <<<
!!========================================================================

     use mmpi, only : mp_init      ! init mpi environment
     use mmpi, only : mp_finalize  ! finalize mpi environment
     use mmpi, only : mp_barrier   ! barrier to synchronize the data
     use mmpi, only : mp_comm_rank ! get index of current process
     use mmpi, only : mp_comm_size ! get number of processes
     !                             !
     use control, only : isscf     ! self-consistent calculation mode
     use control, only : niter     ! number of self-consistent iterations
     use control, only : nprocs    ! number of processes
     use control, only : myid      ! index of current process
     use control, only : master    ! index of master process

     implicit none

!! local variables
     ! loop index
     integer :: iter = 1

     ! convergence flag
     logical :: conv = .false.

!! [body

! initialize mpi envirnoment
# if defined (MPI)

     ! initialize the mpi execution environment
     call mp_init()

     ! determines the rank of the calling process in the communicator
     call mp_comm_rank(myid)

     ! determines the size of the group associated with a communicator
     call mp_comm_size(nprocs)

# endif  /* MPI */

     CTQMC_START: BLOCK

         ! print the welcome messages
         if ( myid == master ) then ! only master node can do it
             call ctqmc_print_header()
         endif ! back if ( myid == master ) block

         ! setup the parameters
         call ctqmc_setup_param()

         ! allocate memory spaces
         call ctqmc_alloc_array()

         ! setup the quantum impurity model
         call ctqmc_setup_model()

         ! print the parameters
         if ( myid == master ) then ! only master node can do it
             call ctqmc_print_summary()
         endif ! back if ( myid == master ) block

     END BLOCK CTQMC_START

!!========================================================================
!!>>> DMFT ITERATION BEGIN                                             <<<
!!========================================================================

     DMFT_CYCLE: do iter=1,niter

         ! write the iter information to screen
         if ( myid == master ) then ! only master node can do it
             call ctqmc_print_it_info(iter)
         endif ! back if ( myid == master ) block

         ! call the quantum impurity solver
         call ctqmc_impurity_solver(iter)

         ! check the self-consistent mode
         if ( isscf == 1 ) then
             EXIT DMFT_CYCLE ! jump out the iteration
         endif ! back if ( isscf == 1 ) block

         ! call the built-in self-consistent engine
         call ctqmc_dmft_selfer()

         ! check whether the convergence is reached
         call ctqmc_dmft_conver(iter, conv)

         ! now the convergence is achieved
         if ( conv .eqv. .true. ) then
             EXIT DMFT_CYCLE ! jump out the iteration
         endif ! back if ( conv .eqv. .true. ) block

     enddo DMFT_CYCLE ! over iter={1,niter} loop

!!========================================================================
!!>>> DMFT ITERATION END                                               <<<
!!========================================================================

     CTQMC_SLEEP: BLOCK

         ! deallocate memory spaces
         call ctqmc_final_array()

         ! print the ending messages
         if ( myid == master ) then ! only master node can do it
             call ctqmc_print_footer()
         endif ! back if ( myid == master ) block

     END BLOCK CTQMC_SLEEP

! finalize mpi envirnoment
# if defined (MPI)

     ! blocks until all processes have reached this routine
     call mp_barrier()

     ! terminates mpi execution environment
     call mp_finalize()

# endif  /* MPI */

!! body]

!!========================================================================
  END PROGRAM CTQMC_MAIN !                                             <<<
!!========================================================================
