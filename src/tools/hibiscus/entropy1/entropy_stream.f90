!-------------------------------------------------------------------------
! project : hibiscus
! program : entropy_config
!           entropy_make_init1
!           entropy_make_init2
! source  : entropy_stream.f90
! type    : subroutine
! author  : li huang (email:huangli712@yahoo.com.cn)
! history : 01/08/2011 by li huang
!           01/09/2011 by li huang
!           01/10/2011 by li huang
!           01/20/2011 by li huang
!           01/26/2011 by li huang
! purpose : initialize the classic maximum entropy method code
! input   :
! output  :
! status  : unstable
! comment :
!-------------------------------------------------------------------------

!>>> setup key parameters for classic maximum entropy method code
  subroutine entropy_config()
     use constants
     use control

     use mmpi

     implicit none

! local variables
! used to check whether the input file (entropy.in) exists
     logical :: exists

! setup default control parameters
!-------------------------------------------------------------------------
     ntime = 129       ! number of imaginary time slice
     nwmax = 200       ! number of frequency point on half axis
     niter = 20        ! number of cycles for maximum entropy method
     ntune = 20        ! number of smooth runs for maximum entropy method
     nstep = 4000      ! number of annealing steps per maximum entropy method cycle
     nband = 1         ! number of bands
     norbs = 2         ! number of orbitals
     ntype = 1         ! type of default model, if ntype = 0, gaussion type, if ntype = 1, flat type
!-------------------------------------------------------------------------
     ainit = 1200._dp  ! initial alpha parameter
     devia = 0.001_dp  ! it is the deviation from the green's function
     beta  = 10.00_dp  ! inversion of real temperature
     sigma = 1.600_dp  ! gauss broadening parameter
     wstep = 0.025_dp  ! frequency step, used to build the frequency mesh
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

! read in input file if possible, only master node can do it
     if ( myid == master ) then
         exists = .false.

! inquire file status: entropy.in
         inquire(file = 'entropy.in', exist = exists)

! read in parameters, default setting should be overrided
         if ( exists .eqv. .true. ) then
             open(mytmp, file='entropy.in', form='formatted', status='unknown')

             read(mytmp,*) ! skip comment lines
             read(mytmp,*)
             read(mytmp,*)
             read(mytmp,*) ntime
             read(mytmp,*) nwmax
             read(mytmp,*) niter
             read(mytmp,*) ntune
             read(mytmp,*) nstep
             read(mytmp,*) nband
             read(mytmp,*) norbs
             read(mytmp,*) ntype

             read(mytmp,*) ! skip comment lines
             read(mytmp,*) ainit
             read(mytmp,*) devia
             read(mytmp,*) beta
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
     call mp_bcast( niter, master )
     call mp_bcast( ntune, master )
     call mp_bcast( nstep, master )
     call mp_bcast( nband, master )
     call mp_bcast( norbs, master )
     call mp_bcast( ntype, master )
     call mp_barrier()

     call mp_bcast( ainit, master )
     call mp_bcast( devia, master )
     call mp_bcast( beta , master )
     call mp_bcast( sigma, master )
     call mp_bcast( wstep, master )
     call mp_barrier()

# endif  /* MPI */

     return
  end subroutine entropy_config

!>>> initialize the classic maximum entropy method code, input original
! imaginary time data and related mesh
  subroutine entropy_make_init1(tmesh, G_qmc, G_dev)
     use constants
     use control

     use mmpi

     implicit none

! external arguments
! time slice data
     real(dp), intent(out) :: tmesh(ntime)

! imaginary green's function data
     real(dp), intent(out) :: G_qmc(ntime,norbs)

! error bar data
     real(dp), intent(out) :: G_dev(ntime,norbs)

! local variables
! loop index for bands
     integer :: i

! loop index for time slice
     integer :: j

! used to check whether the input file (tau.grn.dat) exists
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

             do i=1,nband
                 do j=1,ntime

! read in data
                     read(mytmp,*) tmesh(j), G_qmc(j,i), G_dev(j,i), G_qmc(j,i+nband), G_dev(j,i+nband)

! deal with error bar data
! it is based gaussian model
                     if( ntype == 0 ) then
                         G_dev(j,i)       = G_dev(j,i)       * devia
                         G_dev(j,i+nband) = G_dev(j,i+nband) * devia
! it is based flat model
                     else
                         G_dev(j,i)       = devia
                         G_dev(j,i+nband) = devia
                     endif ! back if ( ntype == 0 ) block

                     if ( abs(G_dev(j,i)) < 0.00001_dp ) then
                         G_dev(j,i)       = 0.00001_dp
                         G_dev(j,i+nband) = 0.00001_dp
                     endif
                     G_dev(j,i)       = one / G_dev(j,i)**2
                     G_dev(j,i+nband) = one / G_dev(j,i+nband)**2

                 enddo ! over j={1,ntime} loop
                 read(mytmp,*) ! skip two lines
                 read(mytmp,*)
             enddo ! over i={1,nband} loop

! close data file
             close(mytmp)

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
     call mp_bcast(G_dev, master)

! block until all processes have reached here
     call mp_barrier()

# endif  /* MPI */

     return
  end subroutine entropy_make_init1

!>>> initialize the classic maximum entropy method code, setup important
! array and variables
  subroutine entropy_make_init2(tmesh, wmesh, model, fnorm, fkern)
     use constants
     use control

     use spring

     implicit none

! external arguments
! imaginary time mesh
     real(dp), intent(in)  :: tmesh(ntime)

! real frequency grid
     real(dp), intent(out) :: wmesh(-nwmax:nwmax)

! default model D(\omega)
     real(dp), intent(out) :: model(-nwmax:nwmax)

! normalization function
     real(dp), intent(out) :: fnorm(-nwmax:nwmax,3)

! fermion kernel function
     real(dp), intent(out) :: fkern(-nwmax:nwmax,ntime)

! local variables
! system time since 1970, Jan 1, used to generate the random number seed
     integer :: system_time

! random number seed for twist generator
     integer :: stream_seed

! init random number generator
     call system_clock(system_time)
     stream_seed = abs( system_time - ( myid * 1981 + 2008 ) * 951049 )
     call spring_sfmt_init(stream_seed)

! build real frequency grid
     call s_linspace_d(-nwmax * wstep, nwmax * wstep, 2 * nwmax + 1, wmesh)

! calculate fermion kernel function
     call entropy_make_fkern(tmesh, wmesh, fkern)

! calculate normalization function f0, f1, and f2
     call entropy_make_fnorm(wmesh, fnorm)

! build default model function
     call entropy_make_model(wmesh, model)

! normalize the model
     call entropy_make_normal(one, fnorm(:,1), model)

     return
  end subroutine entropy_make_init2
