!-------------------------------------------------------------------------
! project : begonia
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
!-------------------------------------------------------------------------
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

!=========================================================================
! setup continuous time quantum Monte Carlo quantum impurity solver related common variables
!=========================================================================
     mkink  = 1024         ! maximum perturbation expansions order
     mfreq  = 8193         ! maximum number of matsubara frequency
!-------------------------------------------------------------------------
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
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^+

             read(mytmp,*)
!------------------------------------------------------------------------+
             read(mytmp,*) nband                                         !
             read(mytmp,*) nspin                                         !
             read(mytmp,*) norbs                                         !
             read(mytmp,*) ncfgs                                         !
             read(mytmp,*) nzero                                         !
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
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^+
     call mp_barrier()

!------------------------------------------------------------------------+
     call mp_bcast( nband , master )                                     !
     call mp_bcast( nspin , master )                                     !
     call mp_bcast( norbs , master )                                     !
     call mp_bcast( ncfgs , master )                                     !
     call mp_bcast( nzero , master )                                     !
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

! read in initial F matrix if available
!-------------------------------------------------------------------------
!     if ( myid == master ) then ! only master node can do it
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
! read the total number of sectors
             read(mytmp,*) nsect
! after we know the total number of sectors, we can allocate memory for array sect
             call ctqmc_allocate_memory_sect()

! read the data for each sector
             do i=1, nsect
                 sect(i)%indx = i
                 sect(i)%nops = norbs
                 read(mytmp,*) ! skip the header
! read the dimension, total number of electrons of this sector
                 read(mytmp,*) j1, sect(i)%ndim, sect(i)%nelectron 
! here, we allocate the memory for sect(i)
                 call alloc_one_sector(sect(i))
! read the next_sector 
                 read(mytmp,*) ! skip the header
                 do j=1, sect(i)%nops
                     read(mytmp,*) j1, sect(i)%next_sector(j,0), sect(i)%next_sector(j,1)  
                 enddo
! read the eigenvalue of this sector
                 read(mytmp,*) ! skip the header
                 do j=1, sect(i)%ndim
                     read(mytmp,*) j1, sect(i)%eigenvalue(j)
                 enddo
             enddo

! read the fmat for each sector
             read(mytmp, *) ! skip the header
             do i=1, nsect
                 read(mytmp,*) ! skip the header
! for each sector, loop over all the operators
                 do j=1, sect(i)%nops 
                     read(mytmp,*) ! skip the header  
! for each operator, read the fmat
                     do k=0, 1
! first, for the annihilation operator
                         ii = sect(i)%next_sector(j,k)
                         if ( ii /= 0 ) then
                             sect(i)%myfmat(j,k)%n = sect(ii)%ndim
                             sect(i)%myfmat(j,k)%m = sect(i)%ndim
! here, we should first allocate the fmat, and nullify it
                             call alloc_one_fmat(myfmat(j,k))
                             call nullify_one_fmat(myfmat(j,k))
! now, we read the fmat item
                             read(mytmp,*) ! skip the header
                             do col=1, sect(i)%ndim
                                 do row=1, sect(ii)%ndim
                                     read(mytmp, *) j1, j2, sect(i)%myfmat(j,0)%item(row, col)
                                 enddo
                             enddo
                         endif ! back to if ( ii/= 0 ) then block
                     enddo ! back to k=0,1 block
                 enddo  ! back to j=1, sect(i)%nops block
             enddo ! back to i=1, nsect block

! close data file
             close(mytmp)

! add the contribution from chemical potential to eigenvalues
             do i=1,nsect
                 sect(i)%eigenvalue = sect(i)%eigenvalue - mune * sect(i)%nelectron
             enddo ! over i={1,nsect} loop

! substract the eigenvalues zero point, here we store the eigen energy
! zero point in U
             r1 = minval(sect(i)%eigenvalue)
             r2 = maxval(sect(i)%eigenvalue)
             do i=2, nsect
                 real_tmp1 = minval(sect(i)%eigenvalue)
                 real_tmp2 = maxval(sect(i)%eigenvalue)
                 if (real_tmp1 < r1) then 
                     r1 = real_tmp1
                 endif
                 if (real_tmp2 > r2) then
                     r2 = real_tmp2
                 endif
             enddo
             U  = r1 + one ! here we choose the minimum as zero point
             do i=1,nsect
                 sect(i)%eigenvalue = sect(i)%eigenvalue - U
             enddo ! over i={1,nsect} loop

! check eigs
! note: \infity - \infity is undefined, which return NaN
             do i=1,nsect
                 do j=1, sect(i)%ndim
                     if ( isnan( exp( - beta * sect(i)%eigenvalue(j) ) - exp( - beta * sect(i)%eigenvalue(j) ) ) ) then
                         call ctqmc_print_error('ctqmc_selfer_init','NaN error, please adjust the zero base of eigs')
                     endif
                 enddo
             enddo ! over i={1,nsect} loop

         else
             call ctqmc_print_error('ctqmc_selfer_init','file atom.cix does not exist')
         endif ! back if ( exists .eqv. .true. ) block
!>>>     endif ! back if ( myid == master ) block

# if defined (MPI)
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
     use sparse
     use spring

     implicit none

! local variables
! loop index
     integer  :: i
     integer  :: j

! system time since 1970, Jan 1, used to generate the random number seed
     integer  :: system_time

! random number seed for twist generator
     integer  :: stream_seed

! dummy sparse matrix in CSR format
     integer  :: sop_it(ncfgs+1)
     integer  :: sop_jt(nzero)
     real(dp) :: sop_t(nzero)

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
         call sparse_csr_mm_csr(   ncfgs, ncfgs, ncfgs, nzero,           &
                                   sop_c(:,i), sop_jc(:,i), sop_ic(:,i), &
                                   sop_d(:,i), sop_jd(:,i), sop_id(:,i), &
                                   sop_t     , sop_jt     , sop_it      )

         call sparse_csr_cp_csr( ncfgs, nzero, sop_t, sop_jt, sop_it, sop_n(:,i), sop_jn(:,i), sop_in(:,i) )
     enddo ! over i={1,norbs} loop

! init op_m, < c^{\dag} c c^{\dag} c >,
! which are used to calculate double occupation number
! note: here we use op_a and op_b as dummy matrix temporarily
     do i=1,norbs-1
         do j=i+1,norbs
             call sparse_csr_mm_csr(   ncfgs, ncfgs, ncfgs, nzero,       &
                                   sop_c(:,i), sop_jc(:,i), sop_ic(:,i), &
                                   sop_d(:,i), sop_jd(:,i), sop_id(:,i), &
                                   sop_a(:,1), sop_ja(:,1), sop_ia(:,1) )
             call sparse_csr_mm_csr(   ncfgs, ncfgs, ncfgs, nzero,       &
                                   sop_c(:,j), sop_jc(:,j), sop_ic(:,j), &
                                   sop_d(:,j), sop_jd(:,j), sop_id(:,j), &
                                   sop_b(:,1), sop_jb(:,1), sop_ib(:,1) )

             call sparse_csr_mm_csr(   ncfgs, ncfgs, ncfgs, nzero,       &
                                   sop_a(:,1), sop_ja(:,1), sop_ia(:,1), &
                                   sop_b(:,1), sop_jb(:,1), sop_ib(:,1), &
                                   sop_t     , sop_jt     , sop_it      )
             call sparse_csr_cp_csr( ncfgs, nzero, sop_t, sop_jt, sop_it, sop_m(:,i,j), sop_jm(:,i,j), sop_im(:,i,j) )

             call sparse_csr_mm_csr(   ncfgs, ncfgs, ncfgs, nzero,       &
                                   sop_b(:,1), sop_jb(:,1), sop_ib(:,1), &
                                   sop_a(:,1), sop_ja(:,1), sop_ia(:,1), &
                                   sop_t     , sop_jt     , sop_it      )
             call sparse_csr_cp_csr( ncfgs, nzero, sop_t, sop_jt, sop_it, sop_m(:,j,i), sop_jm(:,j,i), sop_im(:,j,i) )
         enddo ! over j={i+1,norbs} loop
     enddo ! over i={1,norbs-1} loop

! reinit sparse matrix op_a (sop_a, sop_ja, sop_ia)
! reinit sparse matrix op_b (sop_b, sop_jb, sop_ib)
! the related dense matrix should be an identity matrix
     do i=1,npart
         call sparse_uni_to_csr( ncfgs, nzero, sop_a(:,i), sop_ja(:,i), sop_ia(:,i) )
         call sparse_uni_to_csr( ncfgs, nzero, sop_b(:,i), sop_jb(:,i), sop_ib(:,i) )
     enddo ! over i={1,norbs} loop

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

     implicit none

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
