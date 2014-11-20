!!!=========+=========+=========+=========+=========+=========+=========+!
!!! HIBISCUS/stoch @ iQIST                                               !
!!!                                                                      !
!!! This tool implements the stochastic analytic continuation method to  !
!!! perform analytical continuation for imaginary time green's function  !
!!! outputed by the hybridization expansion version continuous time      !
!!! quantum Monte Carlo (CT-QMC) or Hirsch-Fye quantum Monte Carlo       !
!!! (HF-QMC) quantum impurity solver                                     !
!!! author  : Li Huang (at IOP/CAS & SPCLab/CAEP & UNIFR)                !
!!! version : v2014.10.11T                                               !
!!! status  : WARNING: IN TESTING STAGE, USE IT IN YOUR RISK             !
!!! comment : this code is based originally on Dr. Q. S. Wu's code       !
!!!           any question, please contact with huangli712@gmail.com     !
!!!=========+=========+=========+=========+=========+=========+=========+!

!!
!!
!! Introduction
!! ============
!!
!! The hibiscus/stoch code is often used to perform the analytical
!! continuation to build spectral function from imaginary-time green's
!! function using the modern stochastic analytic continuation method. In
!! principle, it solves the laplace transformation
!!     G(\tau) = \int kernel A(\omega) d\omega
!! where
!!     kernel = \frac{ \exp{-\tau\omega} }{1.0+\exp{-\beta\omega}}
!! for details of the stochastic analytic continuation method, please
!! refer to:
!!     arXiv:cond-mat/0403055
!!
!! Usage
!! =====
!!
!! # ./sac or bin/sac.x
!!
!! Input
!! =====
!!
!! tau.grn.dat (necessary)
!! sac.in (necessary)
!!
!! Output
!! ======
!!
!! sac.image.dat
!! sac.imsum.dat
!! sac.move.dat
!! sac.swap.dat
!!
!! Documents
!! =========
!!
!! For more details, please go to iqist/doc/manual directory.
!!
!!

  program sac_main
     use constants, only : dp, zero, one, mystd
     use mmpi, only : mp_init, mp_finalize
     use mmpi, only : mp_comm_rank, mp_comm_size
     use mmpi, only : mp_barrier

     use control, only : nprocs, myid, master
     use control, only : nstep, ndump
     use context, only : sac_allocate_memory, sac_deallocate_memory

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
         call sac_print_header()
     endif ! back if ( myid == master ) block

! setup the important parameters for stochastic analytic continuation code
     call sac_config()

! print out runtime parameters in summary, only for check
     if ( myid == master ) then
         call sac_print_summary()
     endif ! back if ( myid == master ) block

! allocate memory and initialize
     call sac_allocate_memory()

! input imaginary time data and related mesh
     call sac_make_init1()

! prepare initial data for stochastic analytic continuation code
     call sac_make_init2()

! warmup the stochastic analytic continuation code, in order to achieve
! equilibrium state quickly
     if ( myid == master ) then ! only master node can do it
         write(mystd,'(4X,a)') 'stochastic analytic continuation warmming'
     endif ! back if ( myid == master ) block

     call cpu_time(time_start) ! record starting time
     call sac_warmming()
     call cpu_time(time_end)   ! record ending   time

! print the time information
     if ( myid == master ) then ! only master node can do it
         write(mystd,'(4X,a,f10.3,a)') 'time:', time_end - time_start, 's'
         write(mystd,*)
     endif ! back if ( myid == master ) block

! print the monte carlo sampling header
     if ( myid == master ) then ! only master node can do it
         write(mystd,'(2X,a)') 'HIBISCUS/stoch >>> stochastic analytic continuation running'
         write(mystd,*)
     endif ! back if ( myid == master ) block

! main loop for stochastic analytic continuation code
     step = zero
     SAC_MAIN_LOOP: do iter=1,nstep,ndump

! record start time
         call cpu_time(time_start)

         SAC_DUMP_LOOP: do dump=1,ndump

! increase step by 1
             step = step + one

! perform monte carlo sampling
             call sac_sampling()

! record alpha-resolved image function
             call sac_recording()

         enddo SAC_DUMP_LOOP ! over dump={1,ndump} loop

! record ending time for this iteration
         call cpu_time(time_end)

! reduce alpha-resolved image function from all children nodes
         call sac_reducing()

! dump the statistics data: accept/reject ratio
         if ( myid == master ) then ! only master node can do it
             call sac_dump_aprob(step)
         endif ! back if ( myid == master ) block

! dump the alpha-resolved image function
         if ( myid == master ) then ! only master node can do it
             call sac_dump_image(step)
         endif ! back if ( myid == master ) block

! it is time to write out the statistics results
         if ( myid == master ) then ! only master node can do it
             call sac_print_runtime(step, time_start, time_end)
         endif ! back if ( myid == master ) block

     enddo SAC_MAIN_LOOP ! over iter={1,nstep} loop

! print the footer for stochastic analytic continuation code
     if ( myid == master ) then ! only master node can do it
         write(mystd,'(2X,a)') 'HIBISCUS/stoch >>> stochastic analytic continuation shutdown'
         write(mystd,*)
     endif ! back if ( myid == master ) block

! deallocate memory and finalize
     call sac_deallocate_memory()

! print the footer for stochastic analytic continuation code
     if ( myid == master ) then ! only master node can do it
         call sac_print_footer()
     endif ! back if ( myid == master ) block

! finalize mpi envirnoment
# if defined (MPI)

! blocks until all processes have reached this routine
     call mp_barrier()

! terminates mpi execution environment
     call mp_finalize()

# endif  /* MPI */

  end program sac_main
