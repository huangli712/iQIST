!!!-----------------------------------------------------------------------
!!! project : narcissus
!!! program : ctqmc_setup_param
!!!           ctqmc_setup_array
!!!           ctqmc_setup_model
!!!           ctqmc_solver_init
!!!           ctqmc_final_array
!!! source  : ctqmc_stream.f90
!!! type    : subroutines
!!! author  : li huang (email:lihuang.dmft@gmail.com)
!!! history : 09/16/2009 by li huang (created)
!!!           04/18/2017 by li huang (last modified)
!!! purpose : initialize and finalize the hybridization expansion version
!!!           continuous time quantum Monte Carlo (CTQMC) quantum impurity
!!!           solver and dynamical mean field theory (DMFT) self-consistent
!!!           engine
!!! status  : unstable
!!! comment :
!!!-----------------------------------------------------------------------

!!>>> @ctqmc_setup_param
!!>>> setup key parameters for continuous time quantum Monte Carlo quantum
!!>>> impurity solver and dynamical mean field theory kernel
  subroutine ctqmc_setup_param()
     use parser, only : p_create, p_parse, p_get, p_destroy
     use mmpi, only : mp_bcast, mp_barrier

     use control ! ALL

     implicit none

! local variables
! used to check whether the input file (solver.ctqmc.in) exists
     logical :: exists

!!========================================================================
!!>>> setup general control flags                                      <<<
!!========================================================================
     isscf  = 2            ! non-self-consistent (1) or self-consistent mode (2)
     issun  = 2            ! without symmetry    (1) or with symmetry   mode (2)
     isspn  = 1            ! spin projection, PM (1) or AFM             mode (2)
     isbin  = 2            ! without binning     (1) or with binning    mode (2)
     isort  = 1            ! normal measurement  (1) or legendre polynomial  (2) or chebyshev polynomial (3)
     issus  = 1            ! without suscept.    (1) or with susceptibility  (2)
     isvrt  = 1            ! without vertex      (1) or with vertex function (2)
     isscr  = 1            ! normal (1) or holstein-hubbard (2) or plasmon pole (3) or ohmic model (4)
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

!!========================================================================
!!>>> setup common variables for quantum impurity model                <<<
!!========================================================================
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
     lc     = 1.00_dp      ! lc: screening strength
     wc     = 1.00_dp      ! wc: screening frequency
!-------------------------------------------------------------------------
     mune   = 2.00_dp      ! chemical potential or fermi level
     beta   = 8.00_dp      ! inversion of temperature
     part   = 0.50_dp      ! coupling parameter t for Hubbard model
     alpha  = 0.70_dp      ! mixing parameter for self-consistent engine
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

!!========================================================================
!!>>> setup common variables for quantum impurity solver               <<<
!!========================================================================
     lemax  = 32           ! maximum order for legendre polynomial
     legrd  = 20001        ! number of mesh points for legendre polynomial
     chmax  = 32           ! maximum order for chebyshev polynomial
     chgrd  = 20001        ! number of mesh points for chebyshev polynomial
!-------------------------------------------------------------------------
     mkink  = 1024         ! maximum perturbation expansions order
     mfreq  = 8193         ! maximum number of matsubara frequency
!-------------------------------------------------------------------------
     nffrq  = 32           ! number of matsubara frequency for the two-particle green's function
     nbfrq  = 8            ! number of bosonic frequncy for the two-particle green's function
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
! create the file parser
             call p_create()

! parse the config file
             call p_parse('solver.ctqmc.in')

! extract parameters
             call p_get('isscf' , isscf )
             call p_get('issun' , issun )
             call p_get('isspn' , isspn )
             call p_get('isbin' , isbin )
             call p_get('isort' , isort )
             call p_get('issus' , issus )
             call p_get('isvrt' , isvrt )
             call p_get('isscr' , isscr )

             call p_get('nband' , nband )
             call p_get('nspin' , nspin )
             call p_get('norbs' , norbs )
             call p_get('ncfgs' , ncfgs )
             call p_get('niter' , niter )

             call p_get('U'     , U     )
             call p_get('Uc'    , Uc    )
             call p_get('Uv'    , Uv    )
             call p_get('Jz'    , Jz    )
             call p_get('Js'    , Js    )
             call p_get('Jp'    , Jp    )
             call p_get('lc'    , lc    )
             call p_get('wc'    , wc    )

             call p_get('mune'  , mune  )
             call p_get('beta'  , beta  )
             call p_get('part'  , part  )
             call p_get('alpha' , alpha )

             call p_get('lemax' , lemax )
             call p_get('legrd' , legrd )
             call p_get('chmax' , chmax )
             call p_get('chgrd' , chgrd )

             call p_get('mkink' , mkink )
             call p_get('mfreq' , mfreq )

             call p_get('nffrq' , nffrq )
             call p_get('nbfrq' , nbfrq )
             call p_get('nfreq' , nfreq )
             call p_get('ntime' , ntime )
             call p_get('nflip' , nflip )
             call p_get('ntherm', ntherm)
             call p_get('nsweep', nsweep)
             call p_get('nwrite', nwrite)
             call p_get('nclean', nclean)
             call p_get('nmonte', nmonte)
             call p_get('ncarlo', ncarlo)

! destroy the parser
             call p_destroy()
         endif ! back if ( exists .eqv. .true. ) block
     endif ! back if ( myid == master ) block

! since config parameters may be updated in master node, it is important
! to broadcast config parameters from root to all children processes
# if defined (MPI)

     call mp_bcast( isscf , master )
     call mp_bcast( issun , master )
     call mp_bcast( isspn , master )
     call mp_bcast( isbin , master )
     call mp_bcast( isort , master )
     call mp_bcast( issus , master )
     call mp_bcast( isvrt , master )
     call mp_bcast( isscr , master )
     call mp_barrier()

     call mp_bcast( nband , master )
     call mp_bcast( nspin , master )
     call mp_bcast( norbs , master )
     call mp_bcast( ncfgs , master )
     call mp_bcast( niter , master )
     call mp_barrier()

     call mp_bcast( U     , master )
     call mp_bcast( Uc    , master )
     call mp_bcast( Uv    , master )
     call mp_bcast( Jz    , master )
     call mp_bcast( Js    , master )
     call mp_bcast( Jp    , master )
     call mp_bcast( lc    , master )
     call mp_bcast( wc    , master )
     call mp_barrier()

     call mp_bcast( mune  , master )
     call mp_bcast( beta  , master )
     call mp_bcast( part  , master )
     call mp_bcast( alpha , master )
     call mp_barrier()

     call mp_bcast( lemax , master )
     call mp_bcast( legrd , master )
     call mp_bcast( chmax , master )
     call mp_bcast( chgrd , master )
     call mp_barrier()

     call mp_bcast( mkink , master )
     call mp_bcast( mfreq , master )
     call mp_barrier()

     call mp_bcast( nffrq , master )
     call mp_bcast( nbfrq , master )
     call mp_bcast( nfreq , master )
     call mp_bcast( ntime , master )
     call mp_bcast( nflip , master )
     call mp_bcast( ntherm, master )
     call mp_bcast( nsweep, master )
     call mp_bcast( nwrite, master )
     call mp_bcast( nclean, master )
     call mp_bcast( nmonte, master )
     call mp_bcast( ncarlo, master )
     call mp_barrier()

# endif  /* MPI */

     return
  end subroutine ctqmc_setup_param

!!>>> @ctqmc_setup_array
!!>>> allocate memory for global variables and then initialize them
  subroutine ctqmc_setup_array()
     use context ! ALL

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

!!>>> @ctqmc_setup_model
!!>>> setup impurity model for continuous time quantum Monte Carlo quantum
!!>>> impurity solver and dynamical mean field theory kernel
  subroutine ctqmc_setup_model()
     use constants, only : dp, zero, one, two, pi, czi, czero, mytmp
     use mmpi, only : mp_bcast, mp_barrier

     use control, only : isscr
     use control, only : norbs
     use control, only : lemax, legrd, chmax, chgrd
     use control, only : mfreq
     use control, only : ntime
     use control, only : beta, part
     use control, only : myid, master
     use context, only : tmesh, rmesh, pmesh, qmesh, ppleg, qqche
     use context, only : symm, eimp, ktau, ptau, uumat
     use context, only : hybf

     implicit none

! local variables
! loop index
     integer  :: i
     integer  :: j
     integer  :: k
     integer  :: l

! used to check whether the input file (solver.hyb.in or solver.eimp.in
! or solver.umat.in or solver.ktau.in) exists
     logical  :: exists

! dummy real variables
     real(dp) :: rtmp
     real(dp) :: r1, r2
     real(dp) :: i1, i2

! build imaginary time tau mesh: tmesh
     call s_linspace_d(zero, beta, ntime, tmesh)

! build matsubara frequency mesh: rmesh
     call s_linspace_d(pi / beta, (two * mfreq - one) * (pi / beta), mfreq, rmesh)

! build mesh for legendre polynomial in [-1,1]
     call s_linspace_d(-one, one, legrd, pmesh)

! build mesh for chebyshev polynomial in [-1,1]
     call s_linspace_d(-one, one, chgrd, qmesh)

! build legendre polynomial in [-1,1]
     call s_legendre(lemax, legrd, pmesh, ppleg)

! build chebyshev polynomial in [-1,1]
! note: it is second kind chebyshev polynomial
     call s_chebyshev(chmax, chgrd, qmesh, qqche)

! build initial green's function: i * 2.0 * ( w - sqrt(w*w + 1) )
! using the analytical equation at non-interaction limit, and then
! build initial hybridization function using self-consistent condition
     do i=1,mfreq
         call s_identity_z( norbs, hybf(i,:,:) )
         hybf(i,:,:) = hybf(i,:,:) * (part**2) * (czi*two)
         hybf(i,:,:) = hybf(i,:,:) * ( rmesh(i) - sqrt( rmesh(i)**2 + one ) )
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
             do i=1,norbs
                 do j=1,mfreq
                     read(mytmp,*) k, rtmp, r1, i1, r2, i2
                     hybf(j,i,i) = dcmplx(r1,i1)
                 enddo ! over j={1,mfreq} loop
                 read(mytmp,*) ! skip two lines
                 read(mytmp,*)
             enddo ! over i={1,norbs} loop
             close(mytmp)

         endif ! back if ( exists .eqv. .true. ) block
     endif ! back if ( myid == master ) block

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

! calculate two-index Coulomb interaction, uumat
     call ctqmc_make_uumat(uumat)

! read in two-index Coulomb interaction if available
!-------------------------------------------------------------------------
     if ( myid == master ) then ! only master node can do it
         exists = .false.

! inquire about file's existence
         inquire (file = 'solver.umat.in', exist = exists)

! find input file: solver.umat.in, read it
         if ( exists .eqv. .true. ) then

! read in Coulomb interaction matrix from solver.umat.in
             open(mytmp, file='solver.umat.in', form='formatted', status='unknown')
             do i=1,norbs
                 do j=1,norbs
                     read(mytmp,*) k, l, rtmp
                     uumat(k,l) = rtmp
                 enddo ! over j={1,norbs} loop
             enddo ! over i={1,norbs} loop
             close(mytmp)

         endif ! back if ( exists .eqv. .true. ) block
     endif ! back if ( myid == master ) block

! broadcast uumat from master node to all children nodes
# if defined (MPI)

! broadcast data
     call mp_bcast(uumat, master)

! block until all processes have reached here
     call mp_barrier()

# endif  /* MPI */

! setup initial ktau
     ktau = zero

! setup initial ptau
     ptau = zero

! read in initial screening function and its derivates if available
!-------------------------------------------------------------------------
     if ( myid == master ) then ! only master node can do it
         exists = .false.

! inquire about file's existence
         inquire (file = 'solver.ktau.in', exist = exists)

! find input file: solver.ktau.in, read it
         if ( exists .eqv. .true. ) then

! read in screening function and its derivates from solver.ktau.in
             open(mytmp, file='solver.ktau.in', form='formatted', status='unknown')
             read(mytmp,*) ! skip one line
             do i=1,ntime
                 read(mytmp,*) rtmp, ktau(i), ptau(i)
             enddo ! over i={1,ntime} loop
             close(mytmp)

         else
             if ( isscr == 99 ) then
                 call s_print_error('ctqmc_selfer_init','solver.ktau.in does not exist')
             endif ! back if ( isscr == 99 ) block
         endif ! back if ( exists .eqv. .true. ) block
     endif ! back if ( myid == master ) block

! since the screening function and its derivates may be updated in master
! node, it is important to broadcast it from root to all children processes
# if defined (MPI)

! broadcast data
     call mp_bcast(ktau, master)

! broadcast data
     call mp_bcast(ptau, master)

! block until all processes have reached here
     call mp_barrier()

# endif  /* MPI */

! FINAL STEP
!-------------------------------------------------------------------------
! shift the Coulomb interaction matrix and chemical potential if retarded
! interaction or the so-called dynamical screening effect is considered
     call ctqmc_make_shift(uumat, one)

     return
  end subroutine ctqmc_setup_model

!!>>> @ctqmc_solver_init
!!>>> initialize the key variables for continuous time quantum Monte Carlo
!!>>> quantum impurity solver
  subroutine ctqmc_solver_init()
     use constants, only : zero, czero
     use spring, only : spring_sfmt_init
     use stack, only : istack_clean, istack_push

     use control ! ALL
     use context ! ALL

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

! for stack data structure
!-------------------------------------------------------------------------
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

! for integer variables
!-------------------------------------------------------------------------
! init global variables
     ckink   = 0
     cstat   = 0

! for real variables
!-------------------------------------------------------------------------
! init statistics variables
     insert_tcount = zero; insert_accept = zero; insert_reject = zero
     remove_tcount = zero; remove_accept = zero; remove_reject = zero
     lshift_tcount = zero; lshift_accept = zero; lshift_reject = zero
     rshift_tcount = zero; rshift_accept = zero; rshift_reject = zero
     reflip_tcount = zero; reflip_accept = zero; reflip_reject = zero

! for integer arrays
!-------------------------------------------------------------------------
! init index array
     index_s = 0
     index_e = 0

! init rank  array
     rank    = 0

! init stts  array
! stts = 0 : null occupation case
! stts = 1 : partial occupation case, segment scheme
! stts = 2 : partial occupation case, anti-segment scheme
! stts = 3 : full occupation case
     stts    = 0

! for real arrays
!-------------------------------------------------------------------------
! init time  array
     time_s  = zero
     time_e  = zero

! init hist  array
     hist    = zero

! init auxiliary physical observables
     paux    = zero

! init probability for atomic states
     prob    = zero

! init occupation number array
     nmat    = zero
     nnmat   = zero

! init < k^2 > - < k >^2 array
     kmat    = zero
     kkmat   = zero

! init fidelity susceptibility array
     lmat    = zero
     rmat    = zero
     lrmat   = zero

! init powers of local magnetization array
     szpow   = zero

! init spin-spin correlation function
     schi    = zero
     sschi   = zero
     ssfom   = zero

! init orbital-orbital correlation function
     ochi    = zero
     oochi   = zero
     oofom   = zero

! init two-particle green's function
     g2_re   = zero
     g2_im   = zero
     h2_re   = zero
     h2_im   = zero

! init particle-particle pair susceptibility
     ps_re   = zero
     ps_im   = zero

! init prefactor for improved estimator
     pref    = zero

! init M-matrix related array
     mmat    = zero
     lspace  = zero
     rspace  = zero

! init imaginary time impurity green's function array
     gtau    = zero
     ftau    = zero

! init imaginary time bath weiss's function array
     wtau    = zero

! for complex arrays
!-------------------------------------------------------------------------
! init exponent array exp_s and exp_e
     exp_s   = czero
     exp_e   = czero

! init G-matrix related array
     gmat    = czero
     lsaves  = czero
     rsaves  = czero

! init impurity green's function array
     grnf    = czero
     frnf    = czero

! init bath weiss's function array
     wssf    = czero

! init self-energy function array
! note: sig1 should not be reinitialized here, since it is used to keep
! the persistency of self-energy function
!<     sig1    = czero
     sig2    = czero

! for the other variables/arrays
!-------------------------------------------------------------------------
! fourier transformation hybridization function from matsubara frequency
! space to imaginary time space
     call ctqmc_four_hybf(hybf, htau)

! symmetrize the hybridization function on imaginary time axis if needed
     if ( issun == 2 .or. isspn == 1 ) then
         call ctqmc_symm_gtau(symm, htau)
     endif ! back if ( issun == 2 .or. isspn == 1 ) block

! calculate the 2nd-derivates of htau, which is used in spline subroutines
     call ctqmc_eval_hsed(tmesh, htau, hsed)

! calculate the 2nd-derivates of ktau, which is used in spline subroutines
     call ctqmc_eval_ksed(tmesh, ktau, ksed)

! calculate the 2nd-derivates of ptau, which is used in spline subroutines
     call ctqmc_eval_ksed(tmesh, ptau, psed)

! dump the necessary files
!-------------------------------------------------------------------------
! write out the hybridization function in matsubara frequency axis
     if ( myid == master ) then ! only master node can do it
         call ctqmc_dump_hybf(rmesh, hybf)
     endif ! back if ( myid == master ) block

! write out the hybridization function on imaginary time axis
     if ( myid == master ) then ! only master node can do it
         call ctqmc_dump_htau(tmesh, htau)
     endif ! back if ( myid == master ) block

! write out the screening function and its derivates
     if ( myid == master ) then ! only master node can do it
         call ctqmc_dump_ktau(tmesh, ktau, ptau, ksed, psed)
     endif ! back if ( myid == master ) block

! write out the seed for random number stream, it is useful to reproduce
! the calculation process once fatal error occurs.
     if ( myid == master ) then ! only master node can do it
         write(mystd,'(4X,a,i11)') 'seed:', stream_seed
     endif ! back if ( myid == master ) block

     return
  end subroutine ctqmc_solver_init

!!>>> @ctqmc_final_array
!!>>> garbage collection for this code, please refer to ctqmc_setup_array
  subroutine ctqmc_final_array()
     use context ! ALL

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
