  program test
     use api
     use mmpi
     use constants, only : dp

     implicit none

     type (T_mpi) :: I_mpi
     type (T_segment_azalea) :: I_solver

     integer :: size_t, i
     integer :: mfreq, norbs
     complex(dp), allocatable :: hybf(:)
     complex(dp), allocatable :: grnf(:)
     complex(dp), allocatable :: grnf_s(:)

! allocate memory
     mfreq = 8193
     norbs = 2
     size_t = mfreq * norbs * norbs
     allocate(hybf(size_t))
     allocate(grnf(size_t))
     allocate(grnf_s(size_t))

! initialize the mpi execution environment
     call mp_init()

! determines the rank of the calling process in the communicator
     call mp_comm_rank(I_mpi%myid)

! determines the size of the group associated with a communicator
     call mp_comm_size(I_mpi%nprocs)

! setup I_solver
     I_solver%isscf  = 1
     I_solver%issun  = 1
     I_solver%isspn  = 2
     I_solver%isbin  = 1
     I_solver%nband  = 1
     I_solver%nspin  = 2
     I_solver%norbs  = 2
     I_solver%ncfgs  = 4
     I_solver%niter  = 20
     I_solver%mkink  = 1024
     I_solver%mfreq  = 8193
     I_solver%nfreq  = 128
     I_solver%ntime  = 1024
     I_solver%nflip  = 10000
     I_solver%ntherm = 20000
     I_solver%nsweep = 20000000
     I_solver%nwrite = 2000000
     I_solver%nclean = 20000
     I_solver%nmonte = 100
     I_solver%ncarlo = 100

     I_solver%U     = 4.0
     I_solver%Uc    = 4.0
     I_solver%Uv    = 4.0
     I_solver%Jz    = 0.0
     I_solver%Js    = 0.0
     I_solver%Jp    = 0.0
     I_solver%mune  = 2.0
     I_solver%beta  = 10.0
     I_solver%part  = 0.50
     I_solver%alpha = 0.50

     call init_ctqmc(I_mpi, I_solver)

     do i = 1,20
         call exec_ctqmc(i)
         call get_grnf(size_t, grnf)
         hybf = 0.25 * grnf
         call set_hybf(size_t, hybf)
         print *, 'MAX_ERROR:', maxval(abs(grnf - grnf_s))
         grnf_s = (grnf + grnf_s)/2.0
     enddo

     call stop_ctqmc()

! blocks until all processes have reached this routine
     call mp_barrier()

! terminates mpi execution environment
     call mp_finalize()

! deallocate memory
     deallocate(hybf)
     deallocate(grnf)
     deallocate(grnf_s)

  end program test
