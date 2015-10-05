!!!-----------------------------------------------------------------------
!!! project : camellia
!!! program : ctqmc_config
!!!           ctqmc_setup_array
!!!           ctqmc_selfer_init
!!!           ctqmc_solver_init
!!!           ctqmc_final_array
!!! source  : ctqmc_stream.f90
!!! type    : subroutines
!!! author  : li huang (email:lihuang.dmft@gmail.com)
!!! history : 09/16/2009 by li huang (created)
!!!           08/17/2015 by li huang (last modified)
!!! purpose : initialize and finalize the hybridization expansion version
!!!           continuous time quantum Monte Carlo (CTQMC) quantum impurity
!!!           solver and dynamical mean field theory (DMFT) self-consistent
!!!           engine
!!! status  : unstable
!!! comment :
!!!-----------------------------------------------------------------------

!!>>> ctqmc_config: setup key parameters for continuous time quantum Monte
!!>>> Carlo quantum impurity solver and dynamical mean field theory kernel
  subroutine ctqmc_config()
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
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

!!========================================================================
!!>>> setup common variables for quantum impurity model                <<<
!!========================================================================
     nband  = 1            ! number of correlated bands
     nspin  = 2            ! number of spin projection
     norbs  = nspin*nband  ! number of correlated orbitals (= nband * nspin)
     ncfgs  = 2**norbs     ! number of atomic states
     nzero  = 128          ! maximum number of non-zero elements in sparse matrix style
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
     nvect  = 4            ! number of selected eigenvectors 
     nleja  = 64           ! maximum number of real leja points
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

             call p_get('nband' , nband )
             call p_get('nspin' , nspin )
             call p_get('norbs' , norbs )
             call p_get('ncfgs' , ncfgs )
             call p_get('nzero' , nzero )
             call p_get('niter' , niter )

             call p_get('U'     , U     )
             call p_get('Uc'    , Uc    )
             call p_get('Uv'    , Uv    )
             call p_get('Jz'    , Jz    )
             call p_get('Js'    , Js    )
             call p_get('Jp'    , Jp    )

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
             call p_get('nvect' , nvect )
             call p_get('nleja' , nleja )
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

!------------------------------------------------------------------------+
     call mp_bcast( isscf , master )                                     !
     call mp_bcast( issun , master )                                     !
     call mp_bcast( isspn , master )                                     !
     call mp_bcast( isbin , master )                                     !
     call mp_bcast( isort , master )                                     !
     call mp_bcast( isvrt , master )                                     !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^+
     call mp_barrier()

!------------------------------------------------------------------------+
     call mp_bcast( nband , master )                                     !
     call mp_bcast( nspin , master )                                     !
     call mp_bcast( norbs , master )                                     !
     call mp_bcast( ncfgs , master )                                     !
     call mp_bcast( nvect , master )                                     !
     call mp_bcast( nhmat , master )                                     !
     call mp_bcast( nfmat , master )                                     !
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
     call mp_bcast( lemax , master )                                     !
     call mp_bcast( legrd , master )                                     !
     call mp_bcast( chmax , master )                                     !
     call mp_bcast( chgrd , master )                                     !
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
     call mp_bcast( nleja , master )                                     !
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

! init the real leja points algorithm
     call leja_setup_param(ncfgs, nhmat, nleja)

     return
  end subroutine ctqmc_config

!>>> allocate memory for global variables and then initialize them
  subroutine ctqmc_setup_array()
     use context

     use leja

     implicit none

! allocate memory for leja module
     call leja_setup_array()

! allocate memory for context module
     call ctqmc_allocate_memory_clur()
     call ctqmc_allocate_memory_flvr()

     call ctqmc_allocate_memory_umat()
     call ctqmc_allocate_memory_fmat()
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

     use leja
     use sparse

     use mmpi

     implicit none

! local variables
! loop index
     integer  :: i
     integer  :: j
     integer  :: k

! dummy integer variables
     integer  :: j1, j2, j3

! used to check whether the input file (solver.hyb.in or solver.eimp.in) exists
     logical  :: exists

! dummy real variables
     real(dp) :: rtmp
     real(dp) :: r1, r2
     real(dp) :: i1, i2

! build identity: unity
     unity = czero
     do i=1,norbs
         unity(i,i) = cone
     enddo ! over i={1,norbs} loop

! build mesh for legendre polynomial in [-1,1]
     do i=1,legrd
         pmesh(i) = real(i - 1) * two / real(legrd - 1) - one
     enddo ! over i={1,legrd} loop

! build mesh for chebyshev polynomial in [-1,1]
     do i=1,chgrd
         qmesh(i) = real(i - 1) * two / real(chgrd - 1) - one
     enddo ! over i={1,chgrd} loop

! build imaginary time tau mesh: tmesh
     do i=1,ntime
         tmesh(i) = zero + ( beta - zero ) / real(ntime - 1) * real(i - 1)
     enddo ! over i={1,ntime} loop

! build matsubara frequency mesh: rmesh
     do j=1,mfreq
         rmesh(j) = ( two * real(j - 1) + one ) * ( pi / beta )
     enddo ! over j={1,mfreq} loop

! build matsubara frequency mesh: cmesh
     do k=1,mfreq
         cmesh(k) = czi * ( two * real(k - 1) + one ) * ( pi / beta )
     enddo ! over k={1,mfreq} loop

! build legendre polynomial in [-1,1]
     if ( lemax <= 2 ) then
         call ctqmc_print_error('ctqmc_selfer_init','lemax must be larger than 2')
     endif

     do i=1,legrd
         ppleg(i,1) = one
         ppleg(i,2) = pmesh(i)
         do j=3,lemax
             k = j - 1
             ppleg(i,j) = ( real(2*k-1) * pmesh(i) * ppleg(i,j-1) - real(k-1) * ppleg(i,j-2) ) / real(k)
         enddo ! over j={3,lemax} loop
     enddo ! over i={1,legrd} loop

! build chebyshev polynomial in [-1,1]
! note: it is second kind chebyshev polynomial
     if ( chmax <= 2 ) then
         call ctqmc_print_error('ctqmc_selfer_init','chmax must be larger than 2')
     endif

     do i=1,chgrd
         qqche(i,1) = one
         qqche(i,2) = two * qmesh(i)
         do j=3,chmax
             qqche(i,j) = two * qmesh(i) * qqche(i,j-1) - qqche(i,j-2)
         enddo ! over j={3,chmax} loop
     enddo ! over i={1,chgrd} loop

! build initial green's function: i * 2.0 * ( w - sqrt(w*w + 1) )
! using the analytical equation at non-interaction limit, and then
! build initial hybridization function using self-consistent condition
     do i=1,mfreq
         hybf(i,:,:) = unity * (part**2) * (czi*two) * ( rmesh(i) - sqrt( rmesh(i)**2 + one ) )
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

! setup initial eigs, naux, and saux
     eigs = zero
     naux = zero
     saux = zero

! setup initial op_c and op_d matrix
     op_c = zero
     op_d = zero

! setup initial hmat, vmat, and wmat matrix
     hmat = zero
     vmat = zero
     wmat = zero

! read in initial F matrix if available
!-------------------------------------------------------------------------
     if ( myid == master ) then ! only master node can do it
         exists = .false.

! inquire about file's existence
         inquire (file = 'atom.cix', exist = exists)

! find input file: atom.cix, read it
! file atom.cix is necessary, the code can not run without it
         if ( exists .eqv. .true. ) then

! open data file
             open(mytmp, file='atom.cix', form='formatted', status='unknown')

! read in eigenvalues for local hamiltonian matrix from atom.cix
             read(mytmp,*) ! skip one line
             do i=1,ncfgs
                 read(mytmp,*) k, eigs(i), naux(i), saux(i)
             enddo ! over i={1,ncfgs} loop

! read in F matrix from atom.cix
             read(mytmp,*) ! skip one line
             do i=1,norbs
                 do j=1,ncfgs
                     do k=1,ncfgs
                         read(mytmp,*) j1, j2, j3, op_d(k,j,i)
                     enddo ! over k={1,ncfgs} loop
                 enddo ! over j={1,ncfgs} loop
             enddo ! over i={1,norbs} loop

! read in local hamiltonian matrix from atom.cix
             read(mytmp,*) ! skip one line
             do i=1,ncfgs
                 do j=1,ncfgs
                     read(mytmp,*) j1, j2, hmat(j,i)
                 enddo ! over j={1,ncfgs} loop
             enddo ! over i={1,ncfgs} loop

! read in eigenvectors for local hamiltonian matrix from atom.cix
             read(mytmp,*) ! skip one line
             do i=1,ncfgs
                 do j=1,ncfgs
                     read(mytmp,*) j1, j2, vmat(j,i)
                 enddo ! over j={1,ncfgs} loop
             enddo ! over i={1,ncfgs} loop

! close data file
             close(mytmp)

! add the contribution from chemical potential to eigenvalues
!<             do i=1,ncfgs
!<                 eigs(i) = eigs(i) - mune * naux(i)
!<             enddo ! over i={1,ncfgs} loop

! substract the eigenvalues zero point, here we store the eigen energy
! zero point in U
             r1 = minval(eigs)
             r2 = maxval(eigs)
             U  = r1 + one ! here we choose the minimum as zero point
             do i=1,ncfgs
                 eigs(i) = eigs(i) - U
             enddo ! over i={1,ncfgs} loop

! check eigs
! note: \infity - \infity is undefined, which return NaN
             do i=1,ncfgs
                 if ( isnan( exp( - beta * eigs(i) ) - exp( - beta * eigs(i) ) ) ) then
                     call ctqmc_print_error('ctqmc_selfer_init','NaN error, please adjust the zero base of eigs')
                 endif
             enddo ! over i={1,ncfgs} loop

! check whether hmat is a real symmetric matrix
             do i=1,ncfgs-1
                 do j=i,ncfgs
                     if ( hmat(i,j) /= hmat(j,i) ) then
                         call ctqmc_print_error('ctqmc_selfer_init','hmat is not a real symmetric matrix')
                     endif
                 enddo ! over j={i,ncfgs} loop
             enddo ! over i={1,ncfgs-1} loop

! calculate op_c from op_d
             do i=1,norbs
                 op_c(:,:,i) = transpose( op_d(:,:,i) )
             enddo ! over i={1,norbs} loop

! calculate wmat from vmat, vmat = A, wmat = A^{T}
             wmat = transpose( vmat )

         else
             call ctqmc_print_error('ctqmc_selfer_init','file atom.cix does not exist')
         endif ! back if ( exists .eqv. .true. ) block
     endif ! back if ( myid == master ) block

! broadcast U, hmat, vmat, wmat, op_c, op_d, et al from master node to all children nodes
# if defined (MPI)

! broadcast data
     call mp_bcast(U,    master)

! block until all processes have reached here
     call mp_barrier()

! broadcast data
     call mp_bcast(hmat, master)

! block until all processes have reached here
     call mp_barrier()

! broadcast data
     call mp_bcast(vmat, master)
     call mp_bcast(wmat, master)

! block until all processes have reached here
     call mp_barrier()

! broadcast data
     call mp_bcast(op_c, master)
     call mp_bcast(op_d, master)

! block until all processes have reached here
     call mp_barrier()

! broadcast data
     call mp_bcast(eigs, master)

! block until all processes have reached here
     call mp_barrier()

! broadcast data
     call mp_bcast(naux, master)
     call mp_bcast(saux, master)

! block until all processes have reached here
     call mp_barrier()

# endif  /* MPI */

! now all the processes have one copies of hmat
! convert hmat from dense-stored matrix form to row-stored sparse matrix
     call sparse_dns_to_csr( ncfgs, ncfgs, nhmat, hmat, sop_h, sop_jh, sop_ih )

! now all the processes have one copies of op_c and op_d
! convert op_c from dense-stored matrix form to row-stored sparse matrix
     do i=1,norbs
         call sparse_dns_to_csr( ncfgs, ncfgs, nfmat, op_c(:,:,i), sop_c(:,i), sop_jc(:,i), sop_ic(:,i) )
     enddo ! over i={1,norbs} loop

! convert op_d from dense-stored matrix form to row-stored sparse matrix
     do i=1,norbs
         call sparse_dns_to_csr( ncfgs, ncfgs, nfmat, op_d(:,:,i), sop_d(:,i), sop_jd(:,i), sop_id(:,i) )
     enddo ! over i={1,norbs} loop

! now we rotate op_c and op_d matrix from occupation number basis to eigen basis
! they are necessary in the calculations of (double) occupation numbers
! note: since the data contained in hmat are transfered into op_h already,
! it is not useful any more, then we utilize it as a dummy matrix hereafter
     hmat = zero
     do i=1,norbs
         call dgemm('N', 'N', ncfgs, ncfgs, ncfgs, one, wmat, ncfgs, op_c(:,:,i), ncfgs, zero, hmat, ncfgs)
         call dgemm('N', 'N', ncfgs, ncfgs, ncfgs, one, hmat, ncfgs, vmat, ncfgs, zero, op_c(:,:,i), ncfgs)

         call dgemm('N', 'N', ncfgs, ncfgs, ncfgs, one, wmat, ncfgs, op_d(:,:,i), ncfgs, zero, hmat, ncfgs)
         call dgemm('N', 'N', ncfgs, ncfgs, ncfgs, one, hmat, ncfgs, vmat, ncfgs, zero, op_d(:,:,i), ncfgs)
     enddo ! over i={1,norbs} loop

! prepare necessary data for real leja points algorithm
     call leja_build_spmat(ncfgs, nhmat, -sop_h, sop_jh, sop_ih)

! note: we can not deallocate op_c and op_d to release the memory at here,
! since op_d is still used at ctqmc_make_hub1() subroutine

     return
  end subroutine ctqmc_selfer_init

!>>> initialize the continuous time quantum Monte Carlo quantum impurity solver
  subroutine ctqmc_solver_init()
     use constants
     use control
     use context

     use stack
     use sparse
     use spring

     implicit none

! local variables
! loop index
     integer :: i
     integer :: j

! status flag
     integer :: istat

! system time since 1970, Jan 1, used to generate the random number seed
     integer :: system_time

! random number seed for twist generator
     integer :: stream_seed

! dummy arrays
     real(dp), allocatable :: raux(:,:)
     real(dp), allocatable :: taux(:,:)

! allocate memory
     allocate(raux(ncfgs,ncfgs), stat=istat)
     if ( istat /= 0 ) then
         call ctqmc_print_error('ctqmc_solver_init','can not allocate enough memory')
     endif

     allocate(taux(ncfgs,ncfgs), stat=istat)
     if ( istat /= 0 ) then
         call ctqmc_print_error('ctqmc_solver_init','can not allocate enough memory')
     endif

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

! init empty_v stack structure
     call istack_clean( empty_v )
     do j=mkink,1,-1
         call istack_push( empty_v, j )
     enddo ! over j={mkink,1} loop

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

     reflip_tcount = zero
     reflip_accept = zero
     reflip_reject = zero

! init global variables
     ckink   = 0
     csign   = 1
     cnegs   = 0
     caves   = 0

! init hist  array
     hist    = 0

! init rank  array
     rank    = 0

! init index array
     index_s = 0
     index_e = 0

     index_t = 0
     index_v = 0

! init type  array
     type_v  = 1

! init flvr  array
     flvr_v  = 1

! init time  array
     time_s  = zero
     time_e  = zero

     time_v  = zero

! init probability for atomic states
     prob    = zero
     ddmat   = zero

! init auxiliary physical observables
     paux    = zero

! init spin-spin correlation function
     schi    = zero
     sschi   = zero

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

! init exponent array expt_v
     expt_v  = zero

! init exponent array expt_t
     expt_t  = zero

! init exponent array expt_z
! expt_z is used to store e^{ -\beta \cdot H } persistently
     do i=1,ncfgs
         expt_z(i) = exp( - eigs(i) * beta )
     enddo ! over i={1,ncfgs} loop

! init matrix_ntrace and matrix_ptrace, they need to be scaled to avoid
! numerical instability
     matrix_ntrace = sum( expt_z )
     matrix_ptrace = sum( expt_z )

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

! init op_n, < c^{\dag} c >,
! which are used to calculate occupation number
     do i=1,norbs
         call dgemm('N', 'N', ncfgs, ncfgs, ncfgs, one, op_c(:,:,i), ncfgs, op_d(:,:,i), ncfgs, zero, hmat, ncfgs)
         call sparse_dns_to_csr( ncfgs, ncfgs, nfmat, hmat, sop_n(:,i), sop_jn(:,i), sop_in(:,i) )
     enddo ! over i={1,norbs} loop
     hmat    = zero ! do not forget to reset it

! init op_m, < c^{\dag} c c^{\dag} c >,
! which are used to calculate double occupation number
! note: here we use op_a and op_b as dummy matrix temporarily
     do i=1,norbs-1
         do j=i+1,norbs
             call dgemm('N', 'N', ncfgs, ncfgs, ncfgs, one, op_c(:,:,i), ncfgs, op_d(:,:,i), ncfgs, zero, raux, ncfgs)
             call dgemm('N', 'N', ncfgs, ncfgs, ncfgs, one, op_c(:,:,j), ncfgs, op_d(:,:,j), ncfgs, zero, taux, ncfgs)
             call dgemm('N', 'N', ncfgs, ncfgs, ncfgs, one, raux, ncfgs, taux, ncfgs, zero, hmat, ncfgs)
             call sparse_dns_to_csr( ncfgs, ncfgs, nfmat, hmat, sop_m(:,i,j), sop_jm(:,i,j), sop_im(:,i,j) )

             call dgemm('N', 'N', ncfgs, ncfgs, ncfgs, one, op_c(:,:,j), ncfgs, op_d(:,:,j), ncfgs, zero, raux, ncfgs)
             call dgemm('N', 'N', ncfgs, ncfgs, ncfgs, one, op_c(:,:,i), ncfgs, op_d(:,:,i), ncfgs, zero, taux, ncfgs)
             call dgemm('N', 'N', ncfgs, ncfgs, ncfgs, one, raux, ncfgs, taux, ncfgs, zero, hmat, ncfgs)
             call sparse_dns_to_csr( ncfgs, ncfgs, nfmat, hmat, sop_m(:,j,i), sop_jm(:,j,i), sop_im(:,j,i) )
         enddo ! over j={i+1,norbs} loop
     enddo ! over i={1,norbs-1} loop
     hmat    = zero ! do not forget to reset it

! fourier transformation hybridization function from matsubara frequency
! space to imaginary time space
     call ctqmc_fourier_hybf(hybf, htau)

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

! deallocate memory
     deallocate(raux)
     deallocate(taux)

     return
  end subroutine ctqmc_solver_init

!>>> garbage collection for this program, please refer to ctqmc_setup_array
  subroutine ctqmc_final_array()
     use context

     use leja

     implicit none

! deallocate memory for leja module
     call leja_final_array()

! deallocate memory for context module
     call ctqmc_deallocate_memory_clur()
     call ctqmc_deallocate_memory_flvr()

     call ctqmc_deallocate_memory_umat()
     call ctqmc_deallocate_memory_fmat()
     call ctqmc_deallocate_memory_mmat()

     call ctqmc_deallocate_memory_gmat()
     call ctqmc_deallocate_memory_wmat()
     call ctqmc_deallocate_memory_smat()

     return
  end subroutine ctqmc_final_array
