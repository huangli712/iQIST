!-------------------------------------------------------------------------
! project : azalea
! program : ctqmc_config
!           ctqmc_setup_array
!           ctqmc_selfer_init
!           ctqmc_solver_init
!           ctqmc_final_array
! source  : ctqmc_stream.f90
! type    : subroutine
! author  : li huang (email:huangli712@gmail.com)
! history : 09/16/2009 by li huang
!           09/20/2009 by li huang
!           09/24/2009 by li huang
!           09/27/2009 by li huang
!           10/24/2009 by li huang
!           10/29/2009 by li huang
!           11/01/2009 by li huang
!           11/10/2009 by li huang
!           11/18/2009 by li huang
!           12/01/2009 by li huang
!           12/05/2009 by li huang
!           02/27/2010 by li huang
!           06/08/2010 by li huang
! purpose : initialize and finalize the hybridization expansion version
!           continuous time quantum Monte Carlo (CTQMC) quantum impurity
!           solver and dynamical mean field theory (DMFT) self-consistent
!           engine
! input   :
! output  :
! status  : unstable
! comment :
!-------------------------------------------------------------------------

!>>> setup key parameters for continuous time quantum Monte Carlo quantum
! impurity solver and dynamical mean field theory kernel
  subroutine ctqmc_config()
     use constants
     use control

     use mmpi

     implicit none

! local variables
! used to check whether the input file (solver.ctqmc.in) exists
     logical :: exists

!=========================================================================
! setup dynamical mean field theory self-consistent engine related common variables
!=========================================================================
     isscf  = 1            ! non-self-consistent (1) or self-consistent mode (2)
     issun  = 2            ! without symmetry    (1) or with symmetry   mode (2)
     isspn  = 1            ! spin projection, PM (1) or AFM             mode (2)
     isbin  = 1            ! without binning     (1) or with binning    mode (2)
!-------------------------------------------------------------------------
     nband  = 1            ! number of correlated bands
     nspin  = 2            ! number of spin projection
     norbs  = nspin*nband  ! number of correlated orbitals (= nband * nspin)
     ncfgs  = 2**norbs     ! number of atomic states
     niter  = 20           ! maximum number of DMFT + CTQMC self-consistent iterations
!-------------------------------------------------------------------------
     U      = 4.00_dp      ! U : average Coulomb interaction
     Uc     = 4.00_dp      ! Uc: intraorbital Coulomb interaction
     Uv     = 4.00_dp      ! Uv: interorbital Coulomb interaction, Uv = Uc - 2 * Jz for t2g system
     Jz     = 0.00_dp      ! Jz: Hund's exchange interaction in z axis (Jz = Js = Jp = J)
     Js     = 0.00_dp      ! Js: spin-flip term
     Jp     = 0.00_dp      ! Jp: pair-hopping term
!-------------------------------------------------------------------------
     mune   = 2.00_dp      ! chemical potential or fermi level
     beta   = 8.00_dp      ! inversion of temperature
     part   = 0.50_dp      ! coupling parameter t for Hubbard model
     alpha  = 0.70_dp      ! mixing parameter for self-consistent engine
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

!=========================================================================
! setup continuous time quantum Monte Carlo quantum impurity solver related common variables
!=========================================================================
     mkink  = 1024         ! maximum perturbation expansions order
     mfreq  = 8193         ! maximum number of matsubara frequency
!-------------------------------------------------------------------------
     nfreq  = 128          ! maximum number of matsubara frequency sampling by quantum impurity solver
     ntime  = 1024         ! number of time slice
     nflip  = 20000        ! flip period for spin up and spin down states
     ntherm = 200000       ! maximum number of thermalization steps
     nsweep = 20000000     ! maximum number of quantum Monte Carlo sampling steps
     nwrite = 2000000      ! output period
     nclean = 100000       ! clean update period
     nmonte = 10           ! how often to sampling the gmat and nmat
     ncarlo = 10           ! how often to sampling the gtau and prob
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

! read in input file if possible, only master node can do it
     if ( myid == master ) then
         exists = .false.

! inquire file status: solver.ctqmc.in
         inquire (file = 'solver.ctqmc.in', exist = exists)

! read in parameters, default setting should be overrided
         if ( exists .eqv. .true. ) then
             open(mytmp, file='solver.ctqmc.in', form='formatted', status='unknown')

             read(mytmp,*)
             read(mytmp,*)
             read(mytmp,*)
!------------------------------------------------------------------------+
             read(mytmp,*) isscf                                         !
             read(mytmp,*) issun                                         !
             read(mytmp,*) isspn                                         !
             read(mytmp,*) isbin                                         !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^+

             read(mytmp,*)
!------------------------------------------------------------------------+
             read(mytmp,*) nband                                         !
             read(mytmp,*) nspin                                         !
             read(mytmp,*) norbs                                         !
             read(mytmp,*) ncfgs                                         !
             read(mytmp,*) niter                                         !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^+

             read(mytmp,*)
!------------------------------------------------------------------------+
             read(mytmp,*) U                                             !
             read(mytmp,*) Uc                                            !
             read(mytmp,*) Uv                                            !
             read(mytmp,*) Jz                                            !
             read(mytmp,*) Js                                            !
             read(mytmp,*) Jp                                            !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^+

             read(mytmp,*)
!------------------------------------------------------------------------+
             read(mytmp,*) mune                                          !
             read(mytmp,*) beta                                          !
             read(mytmp,*) part                                          !
             read(mytmp,*) alpha                                         !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^+

             read(mytmp,*)
!------------------------------------------------------------------------+
             read(mytmp,*) mkink                                         !
             read(mytmp,*) mfreq                                         !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^+

             read(mytmp,*)
!------------------------------------------------------------------------+
             read(mytmp,*) nfreq                                         !
             read(mytmp,*) ntime                                         !
             read(mytmp,*) nflip                                         !
             read(mytmp,*) ntherm                                        !
             read(mytmp,*) nsweep                                        !
             read(mytmp,*) nwrite                                        !
             read(mytmp,*) nclean                                        !
             read(mytmp,*) nmonte                                        !
             read(mytmp,*) ncarlo                                        !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^+

             close(mytmp)
         endif ! back if ( exists .eqv. .true. ) block
     endif ! back if ( myid == master ) block

! since config parameters may be updated in master node, it is important
! to broadcast config parameters from root to all children processes
# if defined (MPI)

!------------------------------------------------------------------------+
     call mp_bcast( isscf , master )                                     !
     call mp_bcast( issun , master )                                     !
     call mp_bcast( isspn , master )                                     !
     call mp_bcast( isbin , master )                                     !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^+
     call mp_barrier()

!------------------------------------------------------------------------+
     call mp_bcast( nband , master )                                     !
     call mp_bcast( nspin , master )                                     !
     call mp_bcast( norbs , master )                                     !
     call mp_bcast( ncfgs , master )                                     !
     call mp_bcast( niter , master )                                     !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^+
     call mp_barrier()

!------------------------------------------------------------------------+
     call mp_bcast( U     , master )                                     !
     call mp_bcast( Uc    , master )                                     !
     call mp_bcast( Uv    , master )                                     !
     call mp_bcast( Jz    , master )                                     !
     call mp_bcast( Js    , master )                                     !
     call mp_bcast( Jp    , master )                                     !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^+
     call mp_barrier()

!------------------------------------------------------------------------+
     call mp_bcast( mune  , master )                                     !
     call mp_bcast( beta  , master )                                     !
     call mp_bcast( part  , master )                                     !
     call mp_bcast( alpha , master )                                     !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^+
     call mp_barrier()

!------------------------------------------------------------------------+
     call mp_bcast( mkink , master )                                     !
     call mp_bcast( mfreq , master )                                     !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^+
     call mp_barrier()

!------------------------------------------------------------------------+
     call mp_bcast( nfreq , master )                                     !
     call mp_bcast( ntime , master )                                     !
     call mp_bcast( nflip , master )                                     !
     call mp_bcast( ntherm, master )                                     !
     call mp_bcast( nsweep, master )                                     !
     call mp_bcast( nwrite, master )                                     !
     call mp_bcast( nclean, master )                                     !
     call mp_bcast( nmonte, master )                                     !
     call mp_bcast( ncarlo, master )                                     !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^+
     call mp_barrier()

# endif  /* MPI */

     return
  end subroutine ctqmc_config

!>>> allocate memory for global variables and then initialize them
  subroutine ctqmc_setup_array()
     use context

     implicit none

! allocate memory for context module
     call ctqmc_allocate_memory_clur()

     call ctqmc_allocate_memory_mesh()
     call ctqmc_allocate_memory_meat()
     call ctqmc_allocate_memory_umat()
     call ctqmc_allocate_memory_mmat()

     call ctqmc_allocate_memory_gmat()
     call ctqmc_allocate_memory_wmat()
     call ctqmc_allocate_memory_smat()

     return
  end subroutine ctqmc_setup_array

!>>> initialize the continuous time quantum Monte Carlo quantum impurity
! solver plus dynamical mean field theory self-consistent engine
  subroutine ctqmc_selfer_init()
     use constants
     use control
     use context

     use mmpi

     implicit none

! local variables
! loop index
     integer  :: i
     integer  :: j
     integer  :: k

! used to check whether the input file (solver.hyb.in or solver.eimp.in) exists
     logical  :: exists

! dummy real variables
     real(dp) :: rtmp
     real(dp) :: r1, r2
     real(dp) :: i1, i2

! build identity: unity
     !!unity = czero
     !!do i=1,norbs
     !!    unity(i,i) = cone
     !!enddo ! over i={1,norbs} loop

! build imaginary time tau mesh: tmesh
     do i=1,ntime
         tmesh(i) = zero + ( beta - zero ) / real(ntime - 1) * real(i - 1)
     enddo ! over i={1,ntime} loop

! build matsubara frequency mesh: rmesh
     do j=1,mfreq
         rmesh(j) = ( two * real(j - 1) + one ) * ( pi / beta )
     enddo ! over j={1,mfreq} loop

! build matsubara frequency mesh: cmesh
     !!do k=1,mfreq
     !!    cmesh(k) = czi * ( two * real(k - 1) + one ) * ( pi / beta )
     !!enddo ! over k={1,mfreq} loop

! build initial green's function: i * 2.0 * ( w - sqrt(w*w + 1) )
! using the analytical equation at non-interaction limit, and then
! build initial hybridization function using self-consistent condition
     do i=1,mfreq
         call s_identity_z( norbs, hybf(i,:,:) )
         hybf(i,:,:) = hybf(i,:,:) * (part**2) * (czi*two) * ( rmesh(i) - sqrt( rmesh(i)**2 + one ) )
     enddo ! over i={1,mfreq} loop

! read in initial hybridization function if available
!-------------------------------------------------------------------------
     if ( myid == master ) then ! only master node can do it
         exists = .false.

! inquire about file's existence
         inquire (file = 'solver.hyb.in', exist = exists)

! find input file: solver.hyb.in, read it
         if ( exists .eqv. .true. ) then

             hybf = czero ! reset it to zero

! read in hybridization function from solver.hyb.in
             open(mytmp, file='solver.hyb.in', form='formatted', status='unknown')
             do i=1,nband
                 do j=1,mfreq
                     read(mytmp,*) k, rtmp, r1, i1, r2, i2
                     hybf(j,i,i) = dcmplx(r1,i1)             ! spin up part
                     hybf(j,i+nband,i+nband) = dcmplx(r2,i2) ! spin dn part
                 enddo ! over j={1,mfreq} loop
                 read(mytmp,*) ! skip two lines
                 read(mytmp,*)
             enddo ! over i={1,nband} loop
             close(mytmp)

         endif ! back if ( exists .eqv. .true. ) block
     endif ! back if ( myid == master ) block

! write out the hybridization function
     if ( myid == master ) then ! only master node can do it
         call ctqmc_dump_hybf(rmesh, hybf)
     endif

! since the hybridization function may be updated in master node, it is
! important to broadcast it from root to all children processes
# if defined (MPI)

! broadcast data
     call mp_bcast(hybf, master)

! block until all processes have reached here
     call mp_barrier()

# endif  /* MPI */

! setup initial symm
     symm = 1

! setup initial eimp
     eimp = zero

! read in impurity level and orbital symmetry if available
!-------------------------------------------------------------------------
     if ( myid == master ) then ! only master node can do it
         exists = .false.

! inquire about file's existence
         inquire (file = 'solver.eimp.in', exist = exists)

! find input file: solver.eimp.in, read it
         if ( exists .eqv. .true. ) then

! read in impurity level from solver.eimp.in
             open(mytmp, file='solver.eimp.in', form='formatted', status='unknown')
             do i=1,norbs
                 read(mytmp,*) k, eimp(i), symm(i)
             enddo ! over i={1,norbs} loop
             close(mytmp)

         endif ! back if ( exists .eqv. .true. ) block
     endif ! back if ( myid == master ) block

! broadcast eimp and symm from master node to all children nodes
# if defined (MPI)

! broadcast data
     call mp_bcast(eimp, master)

! broadcast data
     call mp_bcast(symm, master)

! block until all processes have reached here
     call mp_barrier()

# endif  /* MPI */

     return
  end subroutine ctqmc_selfer_init

!>>> initialize the continuous time quantum Monte Carlo quantum impurity solver
  subroutine ctqmc_solver_init()
     use constants
     use control
     use context

     use stack
     use spring

     implicit none

! local variables
! loop index
     integer :: i
     integer :: j

! system time since 1970, Jan 1, used to generate the random number seed
     integer :: system_time

! random number seed for twist generator
     integer :: stream_seed

! init random number generator
     call system_clock(system_time)
     stream_seed = abs( system_time - ( myid * 1981 + 2008 ) * 951049 )
     call spring_sfmt_init(stream_seed)

! init empty_s and empty_e stack structure
     do i=1,norbs
         call istack_clean( empty_s(i) )
         call istack_clean( empty_e(i) )
     enddo ! over i={1,norbs} loop

     do i=1,norbs
         do j=mkink,1,-1
             call istack_push( empty_s(i), j )
             call istack_push( empty_e(i), j )
         enddo ! over j={mkink,1} loop
     enddo ! over i={1,norbs} loop

! init statistics variables
     insert_tcount = zero
     insert_accept = zero
     insert_reject = zero

     remove_tcount = zero
     remove_accept = zero
     remove_reject = zero

     lshift_tcount = zero
     lshift_accept = zero
     lshift_reject = zero

     rshift_tcount = zero
     rshift_accept = zero
     rshift_reject = zero

     reswap_tcount = zero
     reswap_accept = zero
     reswap_reject = zero

     reflip_tcount = zero
     reflip_accept = zero
     reflip_reject = zero

! init global variables
     ckink   = 0
     cstat   = 0

! init hist  array
     hist    = 0

! init rank  array
     rank    = 0

! init stts  array
! stts = 0 : null occupation case
! stts = 1 : partial occupation case, segment scheme
! stts = 2 : partial occupation case, anti-segment scheme
! stts = 3 : full occupation case
     stts    = 0

! init index array
     index_s = 0
     index_e = 0

! init time  array
     time_s  = zero
     time_e  = zero

! init probability for atomic states
     prob    = zero

! init auxiliary physical observables
     paux    = zero

! init occupation number array
     nmat    = zero
     nnmat   = zero

! init M-matrix related array
     mmat    = zero
     lspace  = zero
     rspace  = zero

! init imaginary time impurity green's function array
     gtau    = zero

! init imaginary time bath weiss's function array
     wtau    = zero

! init exponent array exp_s and exp_e
     exp_s   = czero
     exp_e   = czero

! init G-matrix related array
     gmat    = czero
     lsaves  = czero
     rsaves  = czero

! init impurity green's function array
     grnf    = czero

! init bath weiss's function array
     wssf    = czero

! init self-energy function array
! note: sig1 should not be reinitialized here, since it is used to keep
! the persistency of self-energy function
!<     sig1    = czero
     sig2    = czero

! calculate two-index pair interaction, uumat
     call ctqmc_make_uumat(uumat)

! fourier transformation hybridization function from matsubara frequency
! space to imaginary time space
     call ctqmc_four_hybf(hybf, htau)

! symmetrize the hybridization function on imaginary time axis if needed
     if ( issun == 2 .or. isspn == 1 ) then
         call ctqmc_symm_gtau(symm, htau)
     endif

! calculate the 2nd-derivates of htau, which is used in spline subroutines
     call ctqmc_make_hsed(tmesh, htau, hsed)

! write out the hybridization function on imaginary time axis
     if ( myid == master ) then ! only master node can do it
         call ctqmc_dump_htau(tmesh, htau)
     endif

! write out the seed for random number stream, it is useful to reproduce
! the calculation process once fatal error occurs.
     if ( myid == master ) then ! only master node can do it
         write(mystd,'(4X,a,i11)') 'seed:', stream_seed
     endif

     return
  end subroutine ctqmc_solver_init

!>>> garbage collection for this program, please refer to ctqmc_setup_array
  subroutine ctqmc_final_array()
     use context

     implicit none

! deallocate memory for context module
     call ctqmc_deallocate_memory_clur()

     call ctqmc_deallocate_memory_mesh()
     call ctqmc_deallocate_memory_meat()
     call ctqmc_deallocate_memory_umat()
     call ctqmc_deallocate_memory_mmat()

     call ctqmc_deallocate_memory_gmat()
     call ctqmc_deallocate_memory_wmat()
     call ctqmc_deallocate_memory_smat()

     return
  end subroutine ctqmc_final_array
