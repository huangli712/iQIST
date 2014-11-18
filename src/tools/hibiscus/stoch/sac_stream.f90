!-------------------------------------------------------------------------
! project : hibiscus
! program : sai_config
!           sai_make_init1
!           sai_make_init2
! source  : sai_stream.f90
! type    : subroutine
! author  : li huang (email:huangli712@yahoo.com.cn)
! history : 01/08/2011 by li huang
!           01/09/2011 by li huang
!           01/10/2011 by li huang
! purpose : initialize the stochastic analytic continuation code
! input   :
! output  :
! status  : unstable
! comment :
!-------------------------------------------------------------------------

!>>> setup key parameters for stochastic analytic continuation code
  subroutine sai_config()
     use constants
     use control

     use mmpi

     implicit none

! local variables
! used to check whether the input file (sai.in) exists
     logical :: exists

! setup default control parameters
!-------------------------------------------------------------------------
     ntime = 1024    ! number of imaginary time slice
     nwmax = 128     ! number of frequency point on half axis
     ngrid = 10001   ! number of slice of x in [0,1]
     ngamm = 1024    ! number of r_{\gamma} and a_{\gamma}
     nalph = 10      ! number of alpha parameter used in parallel tempering
     nwarm = 4000    ! maximum number of thermalization steps
     nstep = 2000000 ! maximum number of quantum Monte Carlo sampling steps
     ndump = 20000   ! output period for stochastic analytic continuation code
     ltype = 1       ! measurement scheme, 1 = normal; 2 = legendre scheme
     lemax = 64      ! maximum order for legendre polynomial
     legrd = 20001   ! number of mesh points for legendre polynomial in [-1,1] range
!-------------------------------------------------------------------------
     ainit = 1.00_dp ! initial alpha parameter
     ratio = 2.00_dp ! \alpha_(p+1) / \alpha_p = R
     beta  = 10.0_dp ! inversion of real temperature
     eta1  = 0.02_dp ! lorentz broadening parameter \eta_1
     eta2  = 4E-4_dp ! lorentz broadening parameter \eta_2
     sigma = 1.00_dp ! gauss broadening parameter
     wstep = 0.05_dp ! frequency step, used to build the real frequency mesh
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

! read in input file if possible, only master node can do it
     if ( myid == master ) then
         exists = .false.

! inquire file status: sai.in
         inquire(file = 'sai.in', exist = exists)

! read in parameters, default setting should be overrided
         if ( exists .eqv. .true. ) then
             open(mytmp, file='sai.in', form='formatted', status='unknown')

             read(mytmp,*) ! skip comment lines
             read(mytmp,*)
             read(mytmp,*)
             read(mytmp,*) ntime
             read(mytmp,*) nwmax
             read(mytmp,*) ngrid
             read(mytmp,*) ngamm
             read(mytmp,*) nalph
             read(mytmp,*) nwarm
             read(mytmp,*) nstep
             read(mytmp,*) ndump
             read(mytmp,*) ltype
             read(mytmp,*) lemax
             read(mytmp,*) legrd

             read(mytmp,*) ! skip comment lines
             read(mytmp,*) ainit
             read(mytmp,*) ratio
             read(mytmp,*) beta
             read(mytmp,*) eta1
             read(mytmp,*) sigma
             read(mytmp,*) wstep

             close(mytmp)
         endif ! back if ( exists .eqv. .true. ) block
     endif ! back if ( myid == master ) block

! since config parameters may be updated in master node, it is important
! to broadcast config parameters from root to all children processes
# if defined (MPI)

     call mp_bcast( ntime, master )
     call mp_bcast( nwmax, master )
     call mp_bcast( ngrid, master )
     call mp_bcast( ngamm, master )
     call mp_bcast( nalph, master )
     call mp_bcast( nwarm, master )
     call mp_bcast( nstep, master )
     call mp_bcast( ndump, master )
     call mp_bcast( ltype, master )
     call mp_bcast( lemax, master )
     call mp_bcast( legrd, master )
     call mp_barrier()

     call mp_bcast( ainit, master )
     call mp_bcast( ratio, master )
     call mp_bcast( beta , master )
     call mp_bcast( eta1 , master )
     call mp_bcast( sigma, master )
     call mp_bcast( wstep, master )
     call mp_barrier()

# endif  /* MPI */

! calculate \eta_2
     eta2 = eta1 * eta1

     return
  end subroutine sai_config

!>>> initialize the stochastic analytic continuation code, input original
! imaginary time data and related mesh
  subroutine sai_make_init1()
     use constants
     use control
     use context

     use mmpi

     implicit none

! local variables
! loop index
     integer :: i

! used to check whether the input file (tau.grn.dat) exists
! note: which is ouput from mtau.x
     logical :: exists 

! read in original imaginary time data if available
     if ( myid == master ) then ! only master node can do it
         exists = .false.

! inquire about file's existence
         inquire(file = 'tau.grn.dat', exist = exists)

! find input file: tau.grn.dat, read it
         if ( exists .eqv. .true. ) then

! read in imaginary time function from tau.grn.dat
             open(mytmp, file = 'tau.grn.dat', status = 'unknown')

             do i=1,ntime
                 read(mytmp,*) tmesh(i), G_qmc(i), G_dev(i)
                 if ( abs( G_dev(i) ) < eps6 ) then
                     G_dev(i) = eps6
                 endif
                 G_tau(i) = abs( G_qmc(i) ) / G_dev(i)
             enddo ! over i={1,ntime} loop

             close(mytmp)

         else
             call s_print_error('sac_make_init1','file tau.grn.dat does not exist')
         endif ! back if ( exists .eqv. .true. ) block
     endif ! back if ( myid == master ) block

! since the imaginary time function may be updated in master node, it is
! important to broadcast it from root to all children processes
# if defined (MPI)

! broadcast data
     call mp_bcast(tmesh, master)

! block until all processes have reached here
     call mp_barrier()

! broadcast data
     call mp_bcast(G_qmc, master)
     call mp_bcast(G_tau, master)
     call mp_bcast(G_dev, master)

! block until all processes have reached here
     call mp_barrier()

# endif  /* MPI */

     return
  end subroutine sai_make_init1

!>>> initialize the stochastic analytic continuation code, setup important
! array and variables
  subroutine sai_make_init2()
     use constants
     use control
     use context

     use spring

     implicit none

! local variables
! system time since 1970, Jan 1, used to generate the random number seed
     integer :: system_time

! random number seed for twist generator
     integer :: stream_seed

! init random number generator
     call system_clock(system_time)
     stream_seed = abs( system_time - ( myid * 1981 + 2008 ) * 951049 )
     call spring_sfmt_init(stream_seed)

! generate alpha parameters list
     call sai_make_alpha(alpha)

! generate initial configurations randomly
     call sai_make_rgamm(igamm, rgamm)

! setup frequency mesh
     call sai_make_mesh(wmesh)

! setup default model: flat or gaussian type
     if ( sigma <= zero ) then
         call sai_make_const(model)
     else
         call sai_make_gauss(model)
     endif ! back if ( sigma <= zero ) block

! calculate \phi(\omega)
     call sai_make_fphi(model, F_phi)

! generate a dense grid of x = \phi(\omega)
     call sai_make_grid(wgrid, xgrid)

! generate delta function
     call sai_make_delta(xgrid, F_phi, delta)

! generate legendre polynomial
     call sai_make_ppleg(ppleg, pmesh)

! generate kernel function
     call sai_make_kernel(fkern)

! write out the seed for random number stream, it is useful to reproduce
! the calculation process once fatal error occurs.
     if ( myid == master ) then ! only master node can do it
         write(mystd,'(4X,a)') 'stochastic analytic continuation preparing'
         write(mystd,'(4X,a,i11)') 'seed:', stream_seed
         write(mystd,*)
     endif

     return
  end subroutine sai_make_init2
