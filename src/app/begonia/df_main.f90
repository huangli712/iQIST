
  program df_main
     use mmpi, only : mp_init, mp_finalize
     use mmpi, only : mp_comm_rank, mp_comm_size
     use mmpi, only : mp_barrier

     use df_control, only : nprocs, myid, master

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

! finalize mpi envirnoment
# if defined (MPI)

! blocks until all processes have reached this routine
     call mp_barrier()

! terminates mpi execution environment
     call mp_finalize()

# endif  /* MPI */

     print *, 'Hello World!'
  end program df_main
