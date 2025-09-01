!!!-----------------------------------------------------------------------
!!! project : iqist @ narcissus
!!! program : ctqmc_setup_param
!!!           ctqmc_setup_model
!!!           ctqmc_input_mesh_
!!!           ctqmc_input_hybf_
!!!           ctqmc_input_eimp_
!!!           ctqmc_input_umat_
!!!           ctqmc_input_ktau_
!!!           ctqmc_alloc_array
!!!           ctqmc_reset_array
!!!           ctqmc_final_array
!!! source  : ctqmc_stream.f90
!!! type    : subroutines
!!! author  : li huang (email:huangli@caep.cn)
!!! history : 09/16/2009 by li huang (created)
!!!           05/06/2025 by li huang (last modified)
!!! purpose : initialize and finalize the hybridization expansion version
!!!           continuous time quantum Monte Carlo (CTQMC) quantum impurity
!!!           solver and dynamical mean field theory (DMFT) self-consistent
!!!           computational engine.
!!! status  : unstable
!!! comment :
!!!-----------------------------------------------------------------------

!!========================================================================
!!>>> config quantum impurity solver                                   <<<
!!========================================================================

!!
!! @sub ctqmc_setup_param
!!
!! setup key parameters for continuous time quantum Monte Carlo quantum
!! impurity solver and dynamical mean field theory kernel
!!
  subroutine ctqmc_setup_param()
     use parser, only : p_create
     use parser, only : p_parse
     use parser, only : p_get
     use parser, only : p_destroy

     use mmpi, only : mp_bcast
     use mmpi, only : mp_barrier

     use control ! ALL

     implicit none

!! local variables
     ! used to check whether the input file (solver.ctqmc.in) exists
     logical :: exists

!! [body

     ! setup general control flags
     !--------------------------------------------------------------------
     isscf  = 1         ! self-consistent scheme: one-shot mode
     isscr  = 1         ! dynamic interaction: static interaction
     isbnd  = 1         ! symmetry (band part): not to be symmetrized
     isspn  = 1         ! symmetry (spin part): not to be symmetrized
     iswor  = 1         ! worm algorithm: disable
     isort  = 1         ! advanced basis: standard representation
     isobs  = 1         ! various physical observables: nothing
     issus  = 1         ! charge/spin susceptibility: nothing
     isvrt  = 1         ! two-particle green's functions: nothing
     !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

     ! setup common variables for dynamical mean field theory
     !--------------------------------------------------------------------
     niter  = 20        ! maximum number of self-consistent iterations
     !--------------------------------------------------------------------
     alpha  = 0.70_dp   ! mixing parameter for self-consistent iterations
     !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

     ! setup common variables for quantum impurity model
     !--------------------------------------------------------------------
     nband  = 1         ! number of correlated bands
     nspin  = 2         ! number of spin projections
     norbs  = 2         ! number of correlated orbitals
     ncfgs  = 4         ! number of atomic eigenstates
     !--------------------------------------------------------------------
     Uc     = 4.00_dp   ! intra-orbital Coulomb interaction
     Jz     = 0.00_dp   ! Hund's exchange interaction in z axis
     lc     = 1.00_dp   ! screening strength
     wc     = 1.00_dp   ! screening frequency
     !--------------------------------------------------------------------
     mune   = 2.00_dp   ! chemical potential or fermi level
     beta   = 8.00_dp   ! inversion of temperature
     part   = 0.50_dp   ! hopping parameter t for Hubbard model
     !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

     ! setup common variables for quantum impurity solver
     !--------------------------------------------------------------------
     lemax  = 32        ! maximum expansion order for legendre polynomial
     legrd  = 20001     ! number of mesh points for legendre polynomial
     svmax  = 32        ! maximum expansion order for svd polynomial
     svgrd  = 2001      ! number of mesh points for svd polynomial
     !--------------------------------------------------------------------
     mkink  = 1024      ! maximum perturbation expansion order
     mfreq  = 8193      ! maximum number of matsubara frequency points
     !--------------------------------------------------------------------
     nffrq  = 32        ! number of fermionic frequency
     nbfrq  = 8         ! number of bosonic frequncy
     nfreq  = 128       ! number of sampled matsubara frequency points
     ntime  = 1024      ! number of time slices
     !--------------------------------------------------------------------
     nflip  = 20000     ! flip period for spin up and spin down states
     ntherm = 200000    ! number of thermalization steps
     nsweep = 20000000  ! number of Monte Carlo sweeping steps
     nwrite = 2000000   ! output period
     nclean = 100000    ! clean update period
     nmonte = 10        ! how often to sample the observables
     ncarlo = 10        ! how often to sample the observables
     !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

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
             call p_get('isscr' , isscr )
             call p_get('isbnd' , isbnd )
             call p_get('isspn' , isspn )
             call p_get('iswor' , iswor )
             call p_get('isort' , isort )
             call p_get('isobs' , isobs )
             call p_get('issus' , issus )
             call p_get('isvrt' , isvrt )

             call p_get('niter' , niter )

             call p_get('alpha' , alpha )

             call p_get('nband' , nband )
             call p_get('nspin' , nspin )
             call p_get('norbs' , norbs )
             call p_get('ncfgs' , ncfgs )

             call p_get('Uc'    , Uc    )
             call p_get('Jz'    , Jz    )
             call p_get('lc'    , lc    )
             call p_get('wc'    , wc    )

             call p_get('mune'  , mune  )
             call p_get('beta'  , beta  )
             call p_get('part'  , part  )

             call p_get('lemax' , lemax )
             call p_get('legrd' , legrd )
             call p_get('svmax' , svmax )
             call p_get('svgrd' , svgrd )

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

! since config parameters may be updated in master node, it is crucial
! to broadcast config parameters from root to all children processes
# if defined (MPI)

     call mp_bcast( isscf , master )
     call mp_bcast( isscr , master )
     call mp_bcast( isbnd , master )
     call mp_bcast( isspn , master )
     call mp_bcast( iswor , master )
     call mp_bcast( isort , master )
     call mp_bcast( isobs , master )
     call mp_bcast( issus , master )
     call mp_bcast( isvrt , master )
     call mp_barrier()

     call mp_bcast( niter , master )
     call mp_barrier()

     call mp_bcast( alpha , master )
     call mp_barrier()

     call mp_bcast( nband , master )
     call mp_bcast( nspin , master )
     call mp_bcast( norbs , master )
     call mp_bcast( ncfgs , master )
     call mp_barrier()

     call mp_bcast( Uc    , master )
     call mp_bcast( Jz    , master )
     call mp_bcast( lc    , master )
     call mp_bcast( wc    , master )
     call mp_barrier()

     call mp_bcast( mune  , master )
     call mp_bcast( beta  , master )
     call mp_bcast( part  , master )
     call mp_barrier()

     call mp_bcast( lemax , master )
     call mp_bcast( legrd , master )
     call mp_bcast( svmax , master )
     call mp_bcast( svgrd , master )
     call mp_barrier()

     call mp_bcast( mkink , master )
     call mp_bcast( mfreq , master )
     call mp_barrier()

     call mp_bcast( nffrq , master )
     call mp_bcast( nbfrq , master )
     call mp_bcast( nfreq , master )
     call mp_bcast( ntime , master )
     call mp_barrier()

     call mp_bcast( nflip , master )
     call mp_bcast( ntherm, master )
     call mp_bcast( nsweep, master )
     call mp_bcast( nwrite, master )
     call mp_bcast( nclean, master )
     call mp_bcast( nmonte, master )
     call mp_bcast( ncarlo, master )
     call mp_barrier()

# endif  /* MPI */

!! body]

     return
  end subroutine ctqmc_setup_param

!!
!! @sub ctqmc_setup_model
!!
!! setup quantum impurity model for continuous time quantum Monte Carlo
!! quantum impurity solver and dynamical mean field theory kernel
!!
  subroutine ctqmc_setup_model()
     implicit none

!! [body

     ! build various meshes (tmesh, rmesh, lmesh, and rep_l, etc)
     call ctqmc_input_mesh_()

     ! build initial hybridization function (hybf)
     call ctqmc_input_hybf_()

     ! build symmetry vector and impurity level (symm and eimp)
     call ctqmc_input_eimp_()

     ! build Coulomb interaction matrix (umat)
     call ctqmc_input_umat_()

     ! build dynamic interaction if necessary (ktau and ptau)
     call ctqmc_input_ktau_()

!! body]

     return
  end subroutine ctqmc_setup_model

!!========================================================================
!!>>> config quantum impurity model                                    <<<
!!========================================================================

!!
!! @sub ctqmc_input_mesh_
!!
!! try to create various meshes, including time mesh, frequency mesh etc
!!
  subroutine ctqmc_input_mesh_()
     use constants, only : pi, zero, one, two

     use control, only : isort
     use control, only : lemax, legrd
     use control, only : svmax, svgrd
     use control, only : mfreq
     use control, only : ntime
     use control, only : beta

     use context, only : tmesh, rmesh
     use context, only : lmesh, rep_l
     use context, only : smesh, rep_s

     implicit none

!! [body

     ! build imaginary time mesh: tmesh
     call s_linspace_d(zero, beta, ntime, tmesh)

     ! build matsubara frequency mesh: rmesh
     call s_linspace_d(pi / beta, (two * mfreq - one) * (pi / beta), mfreq, rmesh)

     ! build mesh for legendre orthogonal polynomial in [-1,1]
     if ( isort == 2 ) then
         call s_linspace_d(-one, one, legrd, lmesh)
     endif ! back if ( isort == 2 ) block

     ! build mesh for svd orthogonal polynomial in [-1,1]
     if ( isort == 3 ) then
         call s_linspace_d(-one, one, svgrd, smesh)
     endif ! back if ( isort == 3 ) block

     ! build legendre orthogonal polynomial in [-1,1]
     if ( isort == 2 ) then
         call s_leg_basis(lemax, legrd, lmesh, rep_l)
     endif ! back if ( isort == 2 ) block

     ! build svd orthogonal polynomial in [-1,1]
     ! .false. means fermionic kernel, and .true. means bosonic kernel
     if ( isort == 3 ) then
         call s_svd_basis(svmax, svgrd, smesh, rep_s, .false., beta)
     endif ! back if ( isort == 3 ) block

!! body]

     return
  end subroutine ctqmc_input_mesh_

!!
!! @sub ctqmc_input_hybf_
!!
!! try to build initial hybridization function from solver.hyb.in
!!
  subroutine ctqmc_input_hybf_()
     use constants, only : dp
     use constants, only : one, two
     use constants, only : czi, czero
     use constants, only : mytmp

     use mmpi, only : mp_bcast
     use mmpi, only : mp_barrier

     use control, only : norbs
     use control, only : mfreq
     use control, only : part
     use control, only : myid, master

     use context, only : rmesh
     use context, only : hybf

     implicit none

!! local variables
     ! loop index
     integer  :: i
     integer  :: j
     integer  :: k

     ! used to check whether the input file (solver.hyb.in) exists
     logical  :: exists

     ! dummy real variables
     real(dp) :: rtmp
     real(dp) :: r1, r2
     real(dp) :: i1, i2

!! [body

     ! build initial green's function using the analytical expression at
     ! non-interaction limit:
     !
     !     G = i * 2.0 * ( w - sqrt(w*w + 1) ),
     !
     ! and then build initial hybridization function using self-consistent
     ! condition for bethe lattice:
     !
     !     \Delta = t^2 * G,
     !
     ! where t is the hopping parameter
     do i=1,mfreq
         call s_identity_z( norbs, hybf(i,:,:) )
         hybf(i,:,:) = hybf(i,:,:) * (part**2) * (czi*two)
         hybf(i,:,:) = hybf(i,:,:) * ( rmesh(i) - sqrt( rmesh(i)**2 + one ) )
     enddo ! over i={1,mfreq} loop

     ! read in initial hybridization function if available
     !--------------------------------------------------------------------
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
     !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

! since the hybridization function may be updated in master node, it is
! important to broadcast it from root to all children processes
# if defined (MPI)

     ! broadcast data
     call mp_bcast(hybf, master)

     ! block until all processes have reached here
     call mp_barrier()

# endif  /* MPI */

!! body]

     return
  end subroutine ctqmc_input_hybf_

!!
!! @sub ctqmc_input_eimp_
!!
!! try to build orbital symmetry and impurity level from solver.eimp.in
!!
  subroutine ctqmc_input_eimp_()
     use constants, only : zero
     use constants, only : mytmp

     use mmpi, only : mp_bcast
     use mmpi, only : mp_barrier

     use control, only : norbs
     use control, only : myid, master

     use context, only : symm, eimp

     implicit none

!! local variables
     ! loop index
     integer  :: i
     integer  :: j

     ! used to check whether the input file (solver.eimp.in) exists
     logical  :: exists

!! [body

     ! setup initial symm
     symm = 1

     ! setup initial eimp
     eimp = zero

     ! read in impurity level and orbital symmetry if available
     !--------------------------------------------------------------------
     if ( myid == master ) then ! only master node can do it
         exists = .false.

         ! inquire about file's existence
         inquire (file = 'solver.eimp.in', exist = exists)

         ! find input file: solver.eimp.in, read it
         if ( exists .eqv. .true. ) then

             ! read in impurity level from solver.eimp.in
             open(mytmp, file='solver.eimp.in', form='formatted', status='unknown')
             do i=1,norbs
                 read(mytmp,*) j, eimp(i), symm(i)
             enddo ! over i={1,norbs} loop
             close(mytmp)

         endif ! back if ( exists .eqv. .true. ) block
     endif ! back if ( myid == master ) block
     !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

! broadcast eimp and symm from master node to all children nodes
# if defined (MPI)

     ! broadcast data
     call mp_bcast(eimp, master)

     ! broadcast data
     call mp_bcast(symm, master)

     ! block until all processes have reached here
     call mp_barrier()

# endif  /* MPI */

!! body]

     return
  end subroutine ctqmc_input_eimp_

!!
!! @sub ctqmc_input_umat_
!!
!! try to build the Coulomb interaction matrix from solver.umat.in
!!
  subroutine ctqmc_input_umat_()
     use constants, only : dp
     use constants, only : mytmp

     use mmpi, only : mp_bcast
     use mmpi, only : mp_barrier

     use control, only : norbs
     use control, only : myid, master

     use context, only : umat

     implicit none

!! local variables
     ! loop index
     integer  :: i
     integer  :: j
     integer  :: k
     integer  :: l

     ! used to check whether the input file (solver.umat.in) exists
     logical  :: exists

     ! dummy real variables
     real(dp) :: rtmp

!! [body

     ! calculate two-index Coulomb interaction, umat
     call ctqmc_make_umat(umat)

     ! read in two-index Coulomb interaction if available
     !--------------------------------------------------------------------
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
                     umat(k,l) = rtmp
                 enddo ! over j={1,norbs} loop
             enddo ! over i={1,norbs} loop
             close(mytmp)

         endif ! back if ( exists .eqv. .true. ) block
     endif ! back if ( myid == master ) block
     !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

! broadcast umat from master node to all children nodes
# if defined (MPI)

     ! broadcast data
     call mp_bcast(umat, master)

     ! block until all processes have reached here
     call mp_barrier()

# endif  /* MPI */

!! body]

     return
  end subroutine ctqmc_input_umat_

!!
!! @sub ctqmc_input_ktau_
!!
!! try to build the dynamic interaction model from solver.ktau.in
!!
  subroutine ctqmc_input_ktau_()
     use constants, only : dp
     use constants, only : zero, one
     use constants, only : mytmp

     use mmpi, only : mp_bcast
     use mmpi, only : mp_barrier

     use control, only : isscr
     use control, only : ntime
     use control, only : myid, master

     use context, only : ktau, ptau
     use context, only : umat

     implicit none

!! local variables
     ! loop index
     integer  :: i

     ! used to check whether the input file (solver.ktau.in) exists
     logical  :: exists

     ! dummy real variables
     real(dp) :: rtmp

!! [body

     ! setup initial ktau
     ktau = zero

     ! setup initial ptau
     ptau = zero

     ! check isscr, if the interaction is static, return directly
     if ( isscr == 1 ) RETURN

     ! read in initial screening function and its derivates if available
     !--------------------------------------------------------------------
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
             if ( isscr == 4 ) then
                 call s_print_error('ctqmc_input_ktau_','solver.ktau.in does not exist')
             endif ! back if ( isscr == 4 ) block
         endif ! back if ( exists .eqv. .true. ) block
     endif ! back if ( myid == master ) block
     !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

! since the screening function and its derivates may be updated in master
! node, it is important to broadcast them from root to all children processes
# if defined (MPI)

     ! broadcast data
     call mp_bcast(ktau, master)

     ! broadcast data
     call mp_bcast(ptau, master)

     ! block until all processes have reached here
     call mp_barrier()

# endif  /* MPI */

     ! shift the Coulomb interaction matrix and chemical potential if
     ! retarded interaction or the so-called dynamic screening effect
     ! is considered. this step is very important to get the correct
     ! occupation number
     call ctqmc_make_lift(umat, one)

!! body]

     return
  end subroutine ctqmc_input_ktau_

!!========================================================================
!!>>> manage memory for quantum impurity solver                        <<<
!!========================================================================

!!
!! @sub ctqmc_alloc_array
!!
!! allocate memory for global variables and then initialize them
!!
  subroutine ctqmc_alloc_array()
     use context ! ALL

     implicit none

!! [body

     ! allocate memory for context module
     call cat_alloc_clur()
     !
     call cat_alloc_mesh()
     call cat_alloc_meat()
     call cat_alloc_umat()
     call cat_alloc_mmat()
     !
     call cat_alloc_gmat()
     call cat_alloc_wmat()
     call cat_alloc_smat()

!! body]

     return
  end subroutine ctqmc_alloc_array

!!
!! @sub ctqmc_reset_array
!!
!! reset the key variables for continuous time quantum Monte Carlo quantum
!! impurity solver
!!
  subroutine ctqmc_reset_array()
     use constants, only : zero
     use constants, only : czero
     use constants, only : mystd

     use spring, only : spring_sfmt_init
     use stack, only : istack_clean
     use stack, only : istack_push

     use control, only : myid, master

     use context ! ALL

     implicit none

!! local variables
     ! loop index
     integer :: i
     integer :: j

     ! system time since 1970, Jan 1
     ! it is used to generate the random number seed
     integer :: system_time

     ! random number seed for twist generator
     integer :: stream_seed

!! [body

     ! init random number generator
     call system_clock(system_time)
     stream_seed = abs( system_time - ( myid * 1981 + 2008 ) * 951049 )
     call spring_sfmt_init(stream_seed)

     !>>> ctqmc_core module
     !--------------------------------------------------------------------

     ! init global variables
     ckink = 0
     cstat = 0

     ! init statistics variables
     ins_t = zero; ins_a = zero; ins_r = zero
     rmv_t = zero; rmv_a = zero; rmv_r = zero
     lsh_t = zero; lsh_a = zero; lsh_r = zero
     rsh_t = zero; rsh_a = zero; rsh_r = zero
     rfl_t = zero; rfl_a = zero; rfl_r = zero

     !>>> ctqmc_clur module
     !--------------------------------------------------------------------

     ! init index
     index_s = 0
     index_e = 0

     ! init time
     time_s  = zero
     time_e  = zero

     ! init exponent
     exp_s   = czero
     exp_e   = czero

     ! init stack
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

     !>>> ctqmc_mesh module
     !--------------------------------------------------------------------

     ! the variables have been initialized at ctqmc_setup_model()

     !>>> ctqmc_meat module
     !--------------------------------------------------------------------

     ! init autocorrelation function
     ac_v = zero
     ac_f = zero

     ! init histogram
     hist = zero

     ! init probability for atomic eigenstates
     prob = zero

     ! init auxiliary physical observables
     paux = zero

     ! init occupation number
     nimp = zero
     nmat = zero

     ! init kinetic energy fluctuation
     knop = zero
     kmat = zero

     ! init fidelity susceptibility
     lnop = zero
     rnop = zero
     lrmm = zero

     ! init powers of local magnetization
     szpw = zero

     ! init spin-spin correlation function
     schi = zero
     sp_t = zero
     sp_w = zero

     ! init charge-charge correlation function
     cchi = zero
     ch_t = zero
     ch_w = zero

     ! init two-particle green's function
     g2ph = czero
     g2pp = czero

     ! init two-particle vertex function
     h2ph = czero
     h2pp = czero

     !>>> ctqmc_umat module
     !--------------------------------------------------------------------

     ! some variables have been initialized at ctqmc_setup_model()

     ! init rank
     rank = 0

     ! init stts
     stts = 0

     ! init prefactor for improved estimator
     pref = zero

     !>>> ctqmc_mmat module
     !--------------------------------------------------------------------

     ! init M-matrix related array
     mmat   = zero
     lspace = zero
     rspace = zero

     ! init G-matrix related array
     gmat   = czero
     lsaves = czero
     rsaves = czero

     !>>> ctqmc_gmat module
     !--------------------------------------------------------------------

     ! init imaginary time impurity green's function
     gtau = zero
     ftau = zero

     ! init matsubara impurity green's function
     grnf = czero
     frnf = czero

     !>>> ctqmc_wmat module
     !--------------------------------------------------------------------

     ! some variables have been initialized at ctqmc_setup_model()

     ! init imaginary time bath weiss's function
     wtau = zero

     ! init matsubara bath weiss's function
     wssf = czero

     !>>> ctqmc_smat module
     !--------------------------------------------------------------------

     ! sig1 should not be reinitialized here
     ! it is used to keep the persistency of self-energy function

     ! init self-energy function
!<   sig1 = czero
     sig2 = czero

     !>>> postprocess hybridization function
     !--------------------------------------------------------------------

     ! fourier hybridization function from frequency space to time space
     call ctqmc_four_hybf(hybf, htau)

     ! symmetrize the hybridization function on imaginary time axis
     call ctqmc_symm_gtau(symm, htau)

     ! calculate the 2nd-derivates of htau
     ! it is used in spline interpolation
     call ctqmc_eval_hsed(htau, hsed)

     !>>> postprocess dynamic interaction
     !--------------------------------------------------------------------

     ! calculate the 2nd-derivates of ktau
     ! it is used in spline interpolation
     call ctqmc_eval_ksed(ktau, ksed)

     ! calculate the 2nd-derivates of ptau
     ! it is used in spline interpolation
     call ctqmc_eval_ksed(ptau, psed)

     !>>> dump the necessary files
     !--------------------------------------------------------------------

     ! write out the hybridization function
     if ( myid == master ) then ! only master node can do it
         call ctqmc_dump_hybf(hybf)
         call ctqmc_dump_htau(htau)
     endif ! back if ( myid == master ) block

     ! write out the screening function and its derivates
     if ( myid == master ) then ! only master node can do it
         call ctqmc_dump_ktau(ktau, ptau, ksed, psed)
     endif ! back if ( myid == master ) block

     ! write out the seed for random number stream, it is useful to
     ! reproduce the calculation process once fatal error occurs.
     if ( myid == master ) then ! only master node can do it
         write(mystd,'(4X,a,i11)') 'seed:', stream_seed
     endif ! back if ( myid == master ) block

!! body]

     return
  end subroutine ctqmc_reset_array

!!
!! @sub ctqmc_final_array
!!
!! garbage collection for this code, please refer to ctqmc_alloc_array()
!!
  subroutine ctqmc_final_array()
     use context ! ALL

     implicit none

!! [body

     ! deallocate memory for context module
     call cat_free_clur()
     !
     call cat_free_mesh()
     call cat_free_meat()
     call cat_free_umat()
     call cat_free_mmat()
     !
     call cat_free_gmat()
     call cat_free_wmat()
     call cat_free_smat()

!! body]

     return
  end subroutine ctqmc_final_array
