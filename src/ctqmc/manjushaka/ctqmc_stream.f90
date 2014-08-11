!-------------------------------------------------------------------------
! project : lavender
! program : ctqmc_config
!           ctqmc_setup_array
!           ctqmc_selfer_init
!           ctqmc_solver_init
!           ctqmc_final_array
! source  : ctqmc_stream.f90
! type    : subroutine
! author  : li huang (email:huangli712@yahoo.com.cn)
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
     isscf  = 2            ! non-self-consistent (1) or self-consistent mode (2)
     issun  = 2            ! without symmetry    (1) or with symmetry   mode (2)
     isspn  = 1            ! spin projection, PM (1) or AFM             mode (2)
     isbin  = 2            ! without binning     (1) or with binning    mode (2)
     isort  = 1            ! normal measurement  (1) or legendre polynomial  (2) or chebyshev polynomial (3)
     isvrt  = 1            ! without vertex      (1) or with vertex function (2)
     itrun  = 1            ! how to truncate the Hilbert space
     idoub  = 1            ! whether to measure the double occupancy number
!-------------------------------------------------------------------------
     nband  = 1            ! number of correlated bands
     nspin  = 2            ! number of spin projection
     norbs  = nspin*nband  ! number of correlated orbitals (= nband * nspin)
     ncfgs  = 2**norbs     ! number of atomic states
     niter  = 20           ! maximum number of DMFT + CTQMC self-consistent iterations
     nmini  = 0            ! minimum of occupancy number 
     nmaxi  = norbs        ! maximum of occupancy number 
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
     npart  = 16           ! number of parts that the imaginary time axis is split
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
             read(mytmp,*) isort                                         !
             read(mytmp,*) isvrt                                         !
             read(mytmp,*) itrun                                         !
             read(mytmp,*) idoub                                         !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^+

             read(mytmp,*)
!------------------------------------------------------------------------+
             read(mytmp,*) nband                                         !
             read(mytmp,*) nspin                                         !
             read(mytmp,*) norbs                                         !
             read(mytmp,*) ncfgs                                         !
             read(mytmp,*) niter                                         !
             read(mytmp,*) nmini                                         !
             read(mytmp,*) nmaxi                                         !
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
             read(mytmp,*) lemax                                         !
             read(mytmp,*) legrd                                         !
             read(mytmp,*) chmax                                         !
             read(mytmp,*) chgrd                                         !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^+

             read(mytmp,*)
!------------------------------------------------------------------------+
             read(mytmp,*) mkink                                         !
             read(mytmp,*) mfreq                                         !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^+

             read(mytmp,*)
!------------------------------------------------------------------------+
             read(mytmp,*) nffrq                                         !
             read(mytmp,*) nbfrq                                         !
             read(mytmp,*) nfreq                                         !
             read(mytmp,*) ntime                                         !
             read(mytmp,*) npart                                         !
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
     call mp_bcast( isort , master )                                     !
     call mp_bcast( isvrt , master )                                     !
     call mp_bcast( itrun , master )                                     !
     call mp_bcast( idoub , master )                                     !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^+
     call mp_barrier()

!------------------------------------------------------------------------+
     call mp_bcast( nband , master )                                     !
     call mp_bcast( nspin , master )                                     !
     call mp_bcast( norbs , master )                                     !
     call mp_bcast( ncfgs , master )                                     !
     call mp_bcast( niter , master )                                     !
     call mp_bcast( nmini , master )                                     !
     call mp_bcast( nmaxi , master )                                     !
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
     call mp_bcast( nffrq , master )                                     !
     call mp_bcast( nbfrq , master )                                     !
     call mp_bcast( nfreq , master )                                     !
     call mp_bcast( ntime , master )                                     !
     call mp_bcast( npart , master )                                     !
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
     call ctqmc_allocate_memory_flvr()

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

     use m_sector
     use m_npart

     implicit none

! local variables
! loop index
     integer  :: i,j,k,ii

! dummy integer variables
     integer  :: j1

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

!-------------------------------------------------------------------------
! first, read the information of sectors
     eigs = zero
     naux = zero

     if (myid == master) then ! only master node can do it
         exists = .false.
! inquire about file 'atom.cix', this file is necessary, the code can not run without it
         inquire (file = 'atom.cix', exist = exists)

! find 'atom.cix', read it 
         if ( exists .eqv. .true. ) then
             open(mytmp, file='atom.cix', form='formatted', status='unknown')
! skip three header lines
             read(mytmp,*) 
             read(mytmp,*) 
             read(mytmp,*) 
! read the total number of sectors, maximum dimension of sectors, and average dimension of sectors
             read(mytmp,*) nsectors, max_dim_sect, ave_dim_sect

! after we know the total number of sectors, we can allocate memory for array sectors and parts
             call ctqmc_allocate_memory_sect()

! read each sector's information
             do i=1, nsectors
                 read(mytmp,*) ! skip the header

! read the dimension, total number of electrons, number of fermion operators, and start index of this sector
                 read(mytmp,*) j1, sectors(i)%ndim, sectors(i)%nelectron, sectors(i)%nops, sectors(i)%istart

! allocate the memory for sectors(i)
                 call alloc_one_sector(sectors(i))

! read the next_sector index
                 read(mytmp,*) ! skip the header
                 do j=1, sectors(i)%nops
                     read(mytmp,*) j1, sectors(i)%next_sector(j,0), sectors(i)%next_sector(j,1)  
                 enddo

! read the eigenvalue of this sector
                 read(mytmp,*) ! skip the header
                 do j=1, sectors(i)%ndim
                     read(mytmp,*) j1, sectors(i)%myeigval(j)
                 enddo
             enddo
             close(mytmp) 
         else
             call ctqmc_print_error('ctqmc_selfer_init','file atom.cix does not exist')
         endif ! back if ( exists .eqv. .true. ) block

!-------------------------------------------------------------------------
! read the fmat
         exists = .false.
! file 'atom.fmat' is necessary, the code can not run without it
         inquire (file = 'atom.fmat', exist = exists)
! find file 'atom.fmat', read it
         if ( exists .eqv. .true. ) then
             open(mytmp, file='atom.fmat', form='unformatted', status='unknown')
             do i=1, nsectors
                 do j=1, sectors(i)%nops
                     do k=0,1
                         ii = sectors(i)%next_sector(j,k)
                         if (ii == -1) cycle
                         sectors(i)%myfmat(j,k)%n = sectors(ii)%ndim
                         sectors(i)%myfmat(j,k)%m = sectors(i)%ndim
                         call alloc_one_fmat(sectors(i)%myfmat(j,k))
                         read(mytmp) sectors(i)%myfmat(j,k)%item
                     enddo 
                 enddo 
             enddo 
             close(mytmp)
         else
             call ctqmc_print_error('ctqmc_selfer_init','file atom.fmat does not exist')
         endif
     endif ! back if ( myid == master ) block

# if defined (MPI)
! block until all processes have reached here
     call mp_barrier()

     call mp_bcast(nsectors,     master) 
     call mp_bcast(max_dim_sect, master) 
     call mp_bcast(ave_dim_sect, master) 

     if (myid /= master ) then
         call ctqmc_allocate_memory_sect()
     endif

     do i=1, nsectors
         call mp_barrier()
         call mp_bcast(sectors(i)%ndim,        master) 
         call mp_bcast(sectors(i)%nelectron,   master) 
         call mp_bcast(sectors(i)%nops,        master) 
         call mp_bcast(sectors(i)%istart,      master) 
         if ( myid /= master) then
             call alloc_one_sector(sectors(i))
         endif
         call mp_bcast(sectors(i)%next_sector, master) 
         call mp_bcast(sectors(i)%myeigval,    master) 
     enddo
     call mp_barrier()

     do i=1, nsectors
         do j=1, sectors(i)%nops
             do k=0,1
                 ii = sectors(i)%next_sector(j,k)
                 if (ii == -1) cycle
                 if ( myid /= master) then
                     sectors(i)%myfmat(j,k)%n = sectors(ii)%ndim
                     sectors(i)%myfmat(j,k)%m = sectors(i)%ndim
                     call alloc_one_fmat(sectors(i)%myfmat(j,k))
                 endif
                 call mp_barrier()
                 call mp_bcast(sectors(i)%myfmat(j,k)%item, master) 
             enddo 
         enddo 
     enddo 
     call mp_barrier()

# endif  /* MPI */

!-------------------------------------------------------------------------
! add the contribution from chemical potential to eigenvalues
     j1 = 0
     do i=1,nsectors
         do j=1, sectors(i)%ndim
             j1 = j1 + 1
             sectors(i)%myeigval(j) = sectors(i)%myeigval(j) - mune * sectors(i)%nelectron
             eigs(j1) = sectors(i)%myeigval(j)  
             naux(j1) = sectors(i)%nelectron
         enddo
     enddo 

! dump eigs for reference
     if (myid == master) then
         open(mytmp, file='solver.eigs.dat')
         do i=1, nsectors
             write(mytmp, '(2i5,F20.10)') i, sectors(i)%nelectron, minval(sectors(i)%myeigval)
         enddo
         close(mytmp)
     endif

! substract the eigenvalues zero point, here we store the eigen energy zero point in U
     r1 = minval(eigs)
     r2 = maxval(eigs)
! here we choose the minimum as zero point
     U  = r1 + one              
     do i=1,ncfgs
         eigs(i) = eigs(i) - U
     enddo 

! check eigs
! note: \infity - \infity is undefined, which return NaN
     do i=1,ncfgs
         if ( isnan( exp( - beta * eigs(i) ) - exp( - beta * eigs(i) ) ) ) then
             call ctqmc_print_error('ctqmc_selfer_init','NaN error, please &
                                               adjust the zero base of eigs')
         endif
     enddo 

     return
  end subroutine ctqmc_selfer_init

!>>> initialize the continuous time quantum Monte Carlo quantum impurity solver
  subroutine ctqmc_solver_init()
     use constants
     use control
     use context

     use stack
     use spring

     use m_sector
     use m_npart

     implicit none

! local variables
! loop index
     integer  :: i
     integer  :: j

! system time since 1970, Jan 1, used to generate the random number seed
     integer  :: system_time

! random number seed for twist generator
     integer  :: stream_seed

! init random number generator
     call system_clock(system_time)
     stream_seed = abs( system_time - ( myid * 1981 + 2008 ) * 951049 )
     !stream_seed = 123456
     call spring_sfmt_init(stream_seed)
     call random_seed()

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

! init orbital-orbital correlation function
     ochi    = zero
     oochi   = zero

! init two-particle green's function
     g2_re   = zero
     g2_im   = zero

! init vertex function
     h2_re   = zero
     h2_im   = zero

! init occupation number array
     nmat    = zero
     nnmat   = zero

! init M-matrix related array
     mmat    = zero
     lspace  = zero
     rspace  = zero

! init imaginary time impurity green's function array
     gtau    = zero
     ftau    = zero

! init imaginary time bath weiss's function array
     wtau    = zero

! init exponent array expt_v
     expt_v  = zero

! init exponent array expt_t
! expt_t(:,1) : used to store trial  e^{ -(\beta - \tau_n) \cdot H }
! expt_t(:,2) : used to store normal e^{ -(\beta - \tau_n) \cdot H }
! expt_t(:,3) : used to store e^{ -\beta \cdot H } persistently
! expt_t(:,4) : conserved, not used so far
     do i=1,ncfgs
         expt_t(i, 1) = exp( - eigs(i) * beta )
         expt_t(i, 2) = exp( - eigs(i) * beta )
         expt_t(i, 3) = exp( - eigs(i) * beta )
         expt_t(i, 4) = exp( - eigs(i) * beta )
     enddo ! over i={1,ncfgs} loop

! init matrix_ntrace and matrix_ptrace
     matrix_ntrace = sum( expt_t(:, 1) )
     matrix_ptrace = sum( expt_t(:, 2) )

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

! truncate the Hilbert space here
     call ctqmc_make_trunc()

! allocate final_product, occu, double_occu for un-truncated sectors 
     call ctqmc_allocate_memory_occu()
   
! build occu, double_occu for un-truncated sectors
     call ctqmc_make_occu()

! allocate memory for npart
     call ctqmc_allocate_memory_part()

     num_prod = zero
     is_save = 1
     is_copy = .false.
     col_copy = 0
     ops = 0
     ope = 0

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

     return
  end subroutine ctqmc_solver_init

!>>> garbage collection for this program, please refer to ctqmc_setup_array
  subroutine ctqmc_final_array()
     use context

     use m_sector
     use m_npart

     implicit none

! deallocate memory for context module
     call ctqmc_deallocate_memory_clur()
     call ctqmc_deallocate_memory_flvr()

     call ctqmc_deallocate_memory_umat()
     call ctqmc_deallocate_memory_mmat()

     call ctqmc_deallocate_memory_gmat()
     call ctqmc_deallocate_memory_wmat()
     call ctqmc_deallocate_memory_smat()

     call ctqmc_deallocate_memory_sect()

     return
  end subroutine ctqmc_final_array
