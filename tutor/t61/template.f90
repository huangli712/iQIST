  program test
     use capi
     use mmpi
     use constants, only : dp

     implicit none

     type (T_mpi) :: I_mpi
     type (T_segment_narcissus) :: I_solver

     integer :: i
     integer :: mfreq
     integer :: norbs
     integer :: niter
     integer :: size_t

! hybridization function and green's function
     complex(dp), allocatable :: hybf(:)
     complex(dp), allocatable :: grnf(:)
     complex(dp), allocatable :: grnf_s(:)

! setup parameters
     mfreq = 8193 ! number of matsubara frequency points
     norbs = 2    ! number of orbitals
     niter = 20   ! number of iterations
     size_t = mfreq * norbs * norbs

! allocate memory
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
     I_solver%issun  = 2
     I_solver%isspn  = 1
     I_solver%isbin  = 1
     I_solver%isort  = 1
     I_solver%issus  = 1
     I_solver%isvrt  = 1
     I_solver%isscr  = 1
     I_solver%nband  = 1
     I_solver%nspin  = 2
     I_solver%norbs  = 2
     I_solver%ncfgs  = 4
     I_solver%niter  = 1
     I_solver%lemax  = 32
     I_solver%legrd  = 20001
     I_solver%chmax  = 32
     I_solver%chgrd  = 20001
     I_solver%mkink  = 1024
     I_solver%mfreq  = 8193
     I_solver%nffrq  = 32
     I_solver%nbfrq  = 8
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

! init ctqmc impurity solver
     call cat_init_ctqmc(I_mpi, I_solver)

! try to implement the DMFT self-consistent loop
     do i=1,niter
         call cat_exec_ctqmc(i)
         call cat_get_grnf(size_t, grnf)
         hybf = 0.25 * grnf
         call cat_set_hybf(size_t, hybf)
         print *, 'MAX_ERROR:', maxval(abs(grnf - grnf_s))
         grnf_s = (grnf + grnf_s)/2.0
     enddo ! over i={1,niter} loop

! stop ctqmc impurity solver
     call cat_stop_ctqmc()

! blocks until all processes have reached this routine
     call mp_barrier()

! terminates mpi execution environment
     call mp_finalize()

! deallocate memory
     deallocate(hybf)
     deallocate(grnf)
     deallocate(grnf_s)

  end program test
