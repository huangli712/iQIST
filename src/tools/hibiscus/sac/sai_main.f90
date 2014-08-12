!=========+=========+=========+=========+=========+=========+=========+>>>
! A test program for stochastic analytic continuation method             !
! author  : li huang                                                     !
! version : v2011.08.18T                                                 !
! status  : WARNING: IN TESTING STAGE, USE IT IN YOUR RISK               !
! comment : this code is based originally on Q. S. Wu's code             !
!           any question, please contact with huangli712@yahoo.com.cn    !
!=========+=========+=========+=========+=========+=========+=========+>>>

  program sai_main
     use constants
     use control
     use context

     use mmpi

     implicit none

! local variables
! loop index 
     integer  :: iter

! loop index
     integer  :: dump

! current monte carlo sweep number
     real(dp) :: step

! starting time
     real(dp) :: time_start

! ending time
     real(dp) :: time_end

! initialize mpi envirnoment
# if defined (MPI)

! initialize the mpi execution environment
     call mp_init()

! determines the rank of the calling process in the communicator
     call mp_comm_rank(myid)

! determines the size of the group associated with a communicator
     call mp_comm_size(nprocs)

# endif  /* MPI */

! print the running header for stochastic analytic continuation code
     if ( myid == master ) then ! only master node can do it
         call sai_print_header()
     endif

! setup the important parameters for stochastic analytic continuation code
     call sai_config()

! print out runtime parameters in summary, only for check
     if ( myid == master ) then
         call sai_print_summary()
     endif

! allocate memory and initialize
     call sai_allocate_memory()

! input imaginary time data and related mesh
     call sai_make_init1()

! prepare initial data for stochastic analytic continuation code
     call sai_make_init2()

! warmup the stochastic analytic continuation code, in order to achieve
! equilibrium state quickly
     if ( myid == master ) then ! only master node can do it
         write(mystd,'(4X,a)') 'stochastic analytic continuation warmming'
     endif

     call cpu_time(time_start) ! record starting time
     call sai_warmming()
     call cpu_time(time_end)   ! record ending   time

! print the time information
     if ( myid == master ) then ! only master node can do it
         write(mystd,'(4X,a,f10.3,a)') 'time:', time_end - time_start, 's'
         write(mystd,*)
     endif

! print the monte carlo sampling header
     if ( myid == master ) then ! only master node can do it
         write(mystd,'(2X,a)') 'HIBISCUS >>> SAI stochastic analytic continuation running'
         write(mystd,*)
     endif

! main loop for stochastic analytic continuation code
     step = zero
     SAI_MAIN_LOOP: do iter=1,nstep,ndump

! record start time
         call cpu_time(time_start)

         SAI_DUMP_LOOP: do dump=1,ndump

! increase step by 1
             step = step + one

! perform monte carlo sampling
             call sai_sampling()

! record alpha-resolved image function
             call sai_recording()

         enddo SAI_DUMP_LOOP ! over dump={1,ndump} loop

! record ending time for this iteration
         call cpu_time(time_end)

! reduce alpha-resolved image function from all children nodes
         call sai_reducing()

! dump the statistics data: accept/reject ratio
         if ( myid == master ) then ! only master node can do it
             call sai_dump_aprob(step)
         endif

! dump the alpha-resolved image function
         if ( myid == master ) then ! only master node can do it
             call sai_dump_image(step)
         endif

! it is time to write out the statistics results
         if ( myid == master ) then ! only master node can do it
             call sai_print_runtime(step, time_start, time_end)
         endif

     enddo SAI_MAIN_LOOP ! over iter={1,nstep} loop

! print the footer for stochastic analytic continuation code
     if ( myid == master ) then ! only master node can do it
         write(mystd,'(2X,a)') 'HIBISCUS >>> SAI stochastic analytic continuation shutdown'
         write(mystd,*)
     endif

! deallocate memory and finalize
     call sai_deallocate_memory()

! print the footer for stochastic analytic continuation code
     if ( myid == master ) then ! only master node can do it
         call sai_print_footer()
     endif

! finalize mpi envirnoment
# if defined (MPI)

! blocks until all processes have reached this routine
     call mp_barrier()

! terminates mpi execution environment
     call mp_finalize()

# endif  /* MPI */

  end program sai_main
