
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

     if ( myid == master ) then
         call df_print_header()
     endif

     call df_config()

     if ( myid == master ) then
         call df_print_summary()
     endif

     call df_setup_array()

     call df_mesh_init()
     call df_dmft_init()
     call df_dual_init()
     call df_latt_init()
     call df_vert_init()

     call df_run()

     call df_final_array()

     if ( myid == master ) then
         call df_print_footer()
     endif

! finalize mpi envirnoment
# if defined (MPI)

! blocks until all processes have reached this routine
     call mp_barrier()

! terminates mpi execution environment
     call mp_finalize()

# endif  /* MPI */

  end program df_main
