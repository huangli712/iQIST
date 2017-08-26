!!!=========+=========+=========+=========+=========+=========+=========+!
!!! AZALEA @ iQIST                                                       !
!!!                                                                      !
!!! A highly optimized diagrammatic framework for dynamical mean field   ! 
!!! theory which can be used to treat non-local correlations in strongly !
!!! correlated systems
!!!                                                                      !
!!! author  : Li Huang (at IOP/CAS & SPCLab/CAEP & UNIFR)                !
!!! status  : (WARNING) IN TESTING STAGE, USE IT IN YOUR RISK            !
!!! comment : now only the dual fermion approach is implemented          !
!!!           any question, please contact with lihuang.dmft@gmail.com   !
!!!=========+=========+=========+=========+=========+=========+=========+!

!!========================================================================
  PROGRAM DT_MAIN !                                                    <<<
!!========================================================================

     use mmpi, only : mp_init      ! init mpi environment
     use mmpi, only : mp_finalize  ! finalize mpi environment
     use mmpi, only : mp_barrier   ! barrier to synchronize the data
     use mmpi, only : mp_comm_rank ! get index of current process
     use mmpi, only : mp_comm_size ! get number of processes

     use control, only : nprocs    ! number of processes
     use control, only : myid      ! index of current process
     use control, only : master    ! index of master process

     implicit none

! initialize mpi envirnoment
# if defined (MPI)

! initialize the mpi execution environment
     call mp_init()

! determines the rank of the calling process in the communicator
     call mp_comm_rank(myid)

! determines the size of the group associated with a communicator
     call mp_comm_size(nprocs)

# endif  /* MPI */

     DMFT_START: BLOCK

! print the welcome messages
         if ( myid == master ) then ! only master node can do it
             call dt_print_header()
         endif ! back if ( myid == master ) block

! setup the parameters
         call dt_setup_param()

! allocate memory spaces
         call dt_alloc_array()

! setup the quantum lattice model
         call dt_setup_model()

! print the runtime parameters
         if ( myid == master ) then ! only master node can do it
             call dt_print_summary()
         endif ! back if ( myid == master ) block

     END BLOCK DMFT_START

!!========================================================================
     call dt_run()
!!========================================================================

     DMFT_SLEEP: BLOCK

! deallocate memory spaces
         call dt_final_array()

! print the ending messages
         if ( myid == master ) then ! only master node can do it
             call dt_print_footer()
         endif ! back if ( myid == master ) block

     END BLOCK DMFT_SLEEP

! finalize mpi envirnoment
# if defined (MPI)

! blocks until all processes have reached this routine
     call mp_barrier()

! terminates mpi execution environment
     call mp_finalize()

# endif  /* MPI */

!!========================================================================
  END PROGRAM DT_MAIN !                                                <<<
!!========================================================================
