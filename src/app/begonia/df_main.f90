
  program df_main
     use mmpi, only : mp_init, mp_finalize
     use mmpi, only : mp_comm_rank, mp_comm_size
     use mmpi, only : mp_barrier

     implicit none

     print *, 'Hello World!'
  end program df_main
