!!!-----------------------------------------------------------------------
!!! project : hibiscus/stoch
!!! program : sac_config
!!!           sac_make_init1
!!!           sac_make_init2
!!! source  : sac_stream.f90
!!! type    : subroutines
!!! author  : li huang (email:huangli712@gmail.com)
!!! history : 01/08/2011 by li huang
!!!           01/09/2011 by li huang
!!!           11/19/2014 by li huang
!!! purpose : initialize the stochastic analytic continuation code
!!! status  : unstable
!!! comment :
!!!-----------------------------------------------------------------------

!!>>> sac_config: setup key parameters for stochastic analytic
!!>>> continuation code
  subroutine sac_config()
     use parser, only : p_create, p_parse, p_get, p_destroy
     use mmpi, only : mp_bcast, mp_barrier

     use control ! ALL

     implicit none

! local variables
! used to check whether the input file (sac.in) exists
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

! inquire file status: sac.in
         inquire(file = 'sac.in', exist = exists)

! read in parameters, default setting should be overrided
         if ( exists .eqv. .true. ) then
! create the file parser
             call p_create()
! parse the config file
             call p_parse('sac.in')

! extract parameters
             call p_get('ntime' , ntime )
             call p_get('nwmax' , nwmax )
             call p_get('ngrid' , ngrid )
             call p_get('ngamm' , ngamm )
             call p_get('nalph' , nalph )
             call p_get('nwarm' , nwarm )
             call p_get('nstep' , nstep )
             call p_get('ndump' , ndump )
             call p_get('ltype' , ltype )
             call p_get('lemax' , lemax )
             call p_get('legrd' , legrd )

             call p_get('ainit' , ainit )
             call p_get('ratio' , ratio )
             call p_get('beta'  , beta  )
             call p_get('eta1'  , eta1  )
             call p_get('sigma' , sigma )
             call p_get('wstep' , wstep )

! destroy the parser
             call p_destroy()
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
  end subroutine sac_config

!!>>> sac_make_init1: initialize the stochastic analytic continuation
!!>>> code, input original imaginary time data and related mesh
  subroutine sac_make_init1()
     use constants, only : eps6, mytmp
     use mmpi, only : mp_bcast, mp_barrier

     use control, only : ntime
     use control, only : myid, master
     use context, only : tmesh
     use context, only : G_qmc, G_tau, G_dev

     implicit none

! local variables
! loop index
     integer :: i

! used to check whether the input file (tau.grn.dat) exists
     logical :: exists

! read in original imaginary time data if available
     if ( myid == master ) then ! only master node can do it
         exists = .false.

! inquire about file's existence
         inquire(file = 'tau.grn.dat', exist = exists)
         if ( exists .eqv. .false. ) then
             call s_print_error('sac_make_init1','file tau.grn.dat does not exist')
         endif ! back if ( exists .eqv. .false. ) block

! find input file: tau.grn.dat, read it
! read in imaginary time function from tau.grn.dat
         open(mytmp, file = 'tau.grn.dat', form='formatted', status = 'unknown')

         do i=1,ntime
             read(mytmp,*) tmesh(i), G_qmc(i), G_dev(i)
             if ( abs( G_dev(i) ) < eps6 ) then
                 G_dev(i) = eps6
             endif ! back if ( abs( G_dev(i) ) < eps6 ) block
             G_tau(i) = abs( G_qmc(i) ) / G_dev(i)
         enddo ! over i={1,ntime} loop

         close(mytmp)
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
  end subroutine sac_make_init1

!!>>> sac_make_init2: initialize the stochastic analytic continuation
!!>>> code, setup important arrays and variables
  subroutine sac_make_init2()
     use constants, only : zero, one, mystd
     use spring, only : spring_sfmt_init

     use control, only : nwmax
     use control, only : lemax, legrd
     use control, only : sigma, wstep
     use control, only : myid, master
     use context, only : igamm, rgamm
     use context, only : fkern, ppleg, delta, F_phi, model
     use context, only : wmesh, pmesh, xgrid, wgrid
     use context, only : alpha

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

! setup frequency mesh
     call s_linspace_d(-nwmax * wstep, nwmax * wstep, 2 * nwmax + 1, wmesh)

! build mesh for legendre polynomial in [-1,1]
     call s_linspace_d(-one, one, legrd, pmesh)

! setup legendre polynomial
     call s_legendre(lemax, legrd, pmesh, ppleg)

! setup default model: flat or gaussian type
     if ( sigma <= zero ) then
         call sac_make_const(model)
     else
         call sac_make_gauss(model)
     endif ! back if ( sigma <= zero ) block

! generate alpha parameters list
     call sac_make_alpha(alpha)

! generate initial configurations randomly
     call sac_make_rgamm(igamm, rgamm)

! generate \phi(\omega)
     call sac_make_fphi(model, F_phi)

! generate a dense grid of x = \phi(\omega)
     call sac_make_grid(wgrid, xgrid)

! generate delta function
     call sac_make_delta(xgrid, F_phi, delta)

! generate kernel function
     call sac_make_kernel(fkern)

! write out the seed for random number stream, it is useful to reproduce
! the calculation process once fatal error occurs.
     if ( myid == master ) then ! only master node can do it
         write(mystd,'(4X,a)') 'stochastic analytic continuation preparing'
         write(mystd,'(4X,a,i11)') 'seed:', stream_seed
         write(mystd,*)
     endif ! back if ( myid == master ) block

     return
  end subroutine sac_make_init2
