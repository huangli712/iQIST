!!!-----------------------------------------------------------------------
!!! project : CSML (Common Service Modules Library)
!!! program : mmpi
!!! source  : m_mpi.f90
!!! type    : module
!!! author  : li huang (email:lihuang.dmft@gmail.com)
!!! history : 08/09/2006 by li huang (created)
!!!           04/19/2017 by li huang (last modified)
!!! purpose : define my own mpi calls, inspired by famous quantum espresso
!!!           code. we note that the original mpi interfaces/subroutines
!!!           are rather complicated for newbies, thus we try to wrap the
!!!           most important and useful (not all) mpi subroutines (in our
!!!           opinion) in this module to facilite the usage of mpi.
!!! status  : unstable
!!! comment : this module has been tested under the following environment:
!!!               mpich1    1.2.7p1
!!!               mpich2    1.2.1p1
!!!               mpich     3.0.3, 3.0.4, 3.1.4
!!!               mvapich2  1.2.0p1
!!!               openmpi   1.6.4, 1.7.1, 1.8.3
!!!               intel mpi 3.2.0
!!!-----------------------------------------------------------------------

!!
!!
!! Introduction
!! ============
!!
!! This module wraps the most useful mpi calls by using generic programming
!! techniques. It supports most of the collective operations (such as BCAST,
!! GATHER, REDUCE, etc.). These subroutines are
!!     MPI_INIT(),
!!     MPI_FINALIZE(),
!!     MPI_WTIME(),
!!     MPI_WTICK(),
!!     MPI_BARRIER(),
!!     MPI_DIMS_CREATE(),
!!     MPI_CART_CREATE(),
!!     MPI_CART_COORDS(),
!!     MPI_COMM_SPLIT(),
!!     MPI_COMM_RANK(),
!!     MPI_COMM_SIZE(),
!!     MPI_GET_PROCESSOR_NAME(),
!!     MPI_BCAST(),
!!     MPI_GATHER(),
!!     MPI_GATHERV(),
!!     MPI_ALLGATHER(),
!!     MPI_ALLGATHERV(),
!!     MPI_REDUCE(),
!!     MPI_ALLREDUCE(),
!! etc. However, none of the point-to-point operations is supported. in the
!! module, we also try to implement a light-weight error handler. enjoy it!
!!
!! Usage
!! =====
!!
!! 1. include mpi support
!! ----------------------
!!
!! use mmpi
!!
!! pay attention to the module name. it is mmpi, instead of mpi.
!!
!! 2. init mpi environment
!! -----------------------
!!
!! call mp_init() ! init mpi environment
!! call mp_comm_rank(myid) ! get current process it
!! call mp_comm_size(nprocs) ! get number of processes
!!
!! 3. broadcast data
!! -----------------
!!
!! real(dp), allocatable :: real_data(:,:,:)
!! integer, allocatable :: int_data(:)
!! complex(dp), allocatable :: cmplx_data(:,:,:,:)
!!
!! call mp_bcast(real_data, master)
!! call mp_bcast(int_data, master)
!! call mp_bcast(cmplx_data, master)
!!
!! here master == 0 which means the master node/root process.
!!
!! 4. all-reduce data
!! ------------------
!!
!! real(dp), allocatable :: real_data(:)
!! real(dp), allocatable :: real_data_mpi(:)
!!
!! call mp_allreduce(real_data, real_data_mpi) ! all-readuce data
!! real_data = real_data_mpi / number_of_processes ! calculate the average
!!
!! 5. setup barrier
!! ----------------
!!
!! call mp_barrier()
!!
!! 6. finialize mpi environment
!! ----------------------------
!!
!! call mp_finalize()
!!
!!

!!>>> whether the compiler support mpi environment, i.e, mpif90
# if defined (MPI)

  module mmpi
     use mpi

     implicit none

!!========================================================================
!!>>> declare global constants                                         <<<
!!========================================================================

! dp: number precision, double precision for reals
     integer, private, parameter :: dp    = kind(1.0d0)

! mystd: device descriptor, console output
     integer, private, parameter :: mystd = 6

!!========================================================================
!!>>> declare mpi constants (datatypes)                                <<<
!!========================================================================

! m_log: datatype, boolean
     integer, private, parameter :: m_log = MPI_LOGICAL

! m_int: datatype, integer
     integer, private, parameter :: m_int = MPI_INTEGER

! m_rdp: datatype, double precision float
     integer, private, parameter :: m_rdp = MPI_DOUBLE_PRECISION

! m_cdp: datatype, double precision complex
     integer, private, parameter :: m_cdp = MPI_DOUBLE_COMPLEX

!!========================================================================
!!>>> declare common constants                                         <<<
!!========================================================================

! ndims: number of cartesian dimensions, and we need a 2D grid
     integer, private, parameter :: ndims = 2

! reorder: ranking may be reordered (true) or not (false) in a new grid
     logical, private, parameter :: reorder = .false.

! periods: specifying whether the grid is periodic or not in each dimens
     logical, private, parameter :: periods(ndims) = .false.

!!========================================================================
!!>>> declare common variables                                         <<<
!!========================================================================

! ierror: error code for mpi subroutines
     integer, private :: ierror

! istat: status code for mpi subroutines
     integer, private :: istat

! opera: default operator
     integer, private :: opera

! group: default communicator
     integer, private :: group

! isize: size of elements
     integer, private :: isize

!-------------------------------------------------------------------------

! mpi_comm_cart: communicator for cartesian topology
     integer, public, save :: mpi_comm_cart

! mpi_comm_row: rowwise-striped communicator for cartesian topology
     integer, public, save :: mpi_comm_row

! mpi_comm_col: columnwise-striped communicator for cartesian topology
     integer, public, save :: mpi_comm_col

!!========================================================================
!!>>> declare accessibility for module routines                        <<<
!!========================================================================

!!>>> mpi information operation
     public :: mp_info

!!>>> mpi environment operation

! initialize mpi environment
     public :: mp_init

! finalize mpi environment
     public :: mp_finalize

! obtain id of current process
     public :: mp_comm_rank

! obtain the size of process group
     public :: mp_comm_size

! query machine name
     public :: mp_processor

!!>>> mpi cartesian topology operation

! creates a division of processors in a cartesian grid
     public :: mp_dims_create

! makes a new communicator to which topology is cartesian
     public :: mp_cart_create

! determines process coords in cartesian topology given rank in group
     public :: mp_cart_coords

! creates new communicators based on colors and keys for rowwise grid
     public :: mp_comm_split_row

! creates new communicators based on colors and keys for columnwise grid
     public :: mp_comm_split_col

!!>>> synchronics operations

! manually block until all processes reach
     public :: mp_barrier

!!>>> time handler

! obtain time consuming by parallelized program
     public :: mp_wtime

! obtain time precision of mp_wtime
     public :: mp_wtick

!!>>> broadcasting operations

! broadcasting bool
     private :: mp_bcast_log0

! broadcasting bool(:)
     private :: mp_bcast_log1

! broadcasting bool(:,:)
     private :: mp_bcast_log2

! broadcasting bool(:,:,:)
     private :: mp_bcast_log3

! broadcasting bool(:,:,:,:)
     private :: mp_bcast_log4

! broadcasting bool(:,:,:,:,:)
     private :: mp_bcast_log5

! broadcasting int
     private :: mp_bcast_int0

! broadcasting int(:)
     private :: mp_bcast_int1

! broadcasting int(:,:)
     private :: mp_bcast_int2

! broadcasting int(:,:,:)
     private :: mp_bcast_int3

! broadcasting int(:,:,:,:)
     private :: mp_bcast_int4

! broadcasting int(:,:,:,:,:)
     private :: mp_bcast_int5

! broadcasting real
     private :: mp_bcast_rdp0

! broadcasting real(:)
     private :: mp_bcast_rdp1

! broadcasting real(:,:)
     private :: mp_bcast_rdp2

! broadcasting real(:,:,:)
     private :: mp_bcast_rdp3

! broadcasting real(:,:,:,:)
     private :: mp_bcast_rdp4

! broadcasting real(:,:,:,:,:)
     private :: mp_bcast_rdp5

! broadcasting complex
     private :: mp_bcast_cdp0

! broadcasting complex(:)
     private :: mp_bcast_cdp1

! broadcasting complex(:,:)
     private :: mp_bcast_cdp2

! broadcasting complex(:,:,:)
     private :: mp_bcast_cdp3

! broadcasting complex(:,:,:,:)
     private :: mp_bcast_cdp4

! broadcasting complex(:,:,:,:,:)
     private :: mp_bcast_cdp5

!!>>> gathering operations

! gathering int(:)
     private :: mp_gather_int1

! gathering int(:,:)
     private :: mp_gather_int2

! gathering int(:,:,:)
     private :: mp_gather_int3

! gathering int(:,:,:,:)
     private :: mp_gather_int4

! gathering int(:,:,:,:,:)
     private :: mp_gather_int5

! gathering real(:)
     private :: mp_gather_rdp1

! gathering real(:,:)
     private :: mp_gather_rdp2

! gathering real(:,:,:)
     private :: mp_gather_rdp3

! gathering real(:,:,:,:)
     private :: mp_gather_rdp4

! gathering real(:,:,:,:,:)
     private :: mp_gather_rdp5

! gathering complex(:)
     private :: mp_gather_cdp1

! gathering complex(:,:)
     private :: mp_gather_cdp2

! gathering complex(:,:,:)
     private :: mp_gather_cdp3

! gathering complex(:,:,:,:)
     private :: mp_gather_cdp4

! gathering complex(:,:,:,:,:)
     private :: mp_gather_cdp5

!!>>> gatherving operations

! gatherving int(:)
     private :: mp_gatherv_int1

! gatherving int(:,:)
     private :: mp_gatherv_int2

! gatherving int(:,:,:)
     private :: mp_gatherv_int3

! gatherving int(:,:,:,:)
     private :: mp_gatherv_int4

! gatherving int(:,:,:,:,:)
     private :: mp_gatherv_int5

! gatherving real(:)
     private :: mp_gatherv_rdp1

! gatherving real(:,:)
     private :: mp_gatherv_rdp2

! gatherving real(:,:,:)
     private :: mp_gatherv_rdp3

! gatherving real(:,:,:,:)
     private :: mp_gatherv_rdp4

! gatherving real(:,:,:,:,:)
     private :: mp_gatherv_rdp5

! gatherving complex(:)
     private :: mp_gatherv_cdp1

! gatherving complex(:,:)
     private :: mp_gatherv_cdp2

! gatherving complex(:,:,:)
     private :: mp_gatherv_cdp3

! gatherving complex(:,:,:,:)
     private :: mp_gatherv_cdp4

! gatherving complex(:,:,:,:,:)
     private :: mp_gatherv_cdp5

!!>>> allgathering operations

! allgathering int(:)
     private :: mp_allgather_int1

! allgathering int(:,:)
     private :: mp_allgather_int2

! allgathering int(:,:,:)
     private :: mp_allgather_int3

! allgathering int(:,:,:,:)
     private :: mp_allgather_int4

! allgathering int(:,:,:,:,:)
     private :: mp_allgather_int5

! allgathering real(:)
     private :: mp_allgather_rdp1

! allgathering real(:,:)
     private :: mp_allgather_rdp2

! allgathering real(:,:,:)
     private :: mp_allgather_rdp3

! allgathering real(:,:,:,:)
     private :: mp_allgather_rdp4

! allgathering real(:,:,:,:,:)
     private :: mp_allgather_rdp5

! allgathering complex(:)
     private :: mp_allgather_cdp1

! allgathering complex(:,:)
     private :: mp_allgather_cdp2

! allgathering complex(:,:,:)
     private :: mp_allgather_cdp3

! allgathering complex(:,:,:,:)
     private :: mp_allgather_cdp4

! allgathering complex(:,:,:,:,:)
     private :: mp_allgather_cdp5

!!>>> allgatherving operations

! allgatherving int(:)
     private :: mp_allgatherv_int1

! allgatherving int(:,:)
     private :: mp_allgatherv_int2

! allgatherving int(:,:,:)
     private :: mp_allgatherv_int3

! allgatherving int(:,:,:,:)
     private :: mp_allgatherv_int4

! allgatherving int(:,:,:,:,:)
     private :: mp_allgatherv_int5

! allgatherving real(:)
     private :: mp_allgatherv_rdp1

! allgatherving real(:,:)
     private :: mp_allgatherv_rdp2

! allgatherving real(:,:,:)
     private :: mp_allgatherv_rdp3

! allgatherving real(:,:,:,:)
     private :: mp_allgatherv_rdp4

! allgatherving real(:,:,:,:,:)
     private :: mp_allgatherv_rdp5

! allgatherving complex(:)
     private :: mp_allgatherv_cdp1

! allgatherving complex(:,:)
     private :: mp_allgatherv_cdp2

! allgatherving complex(:,:,:)
     private :: mp_allgatherv_cdp3

! allgatherving complex(:,:,:,:)
     private :: mp_allgatherv_cdp4

! allgatherving complex(:,:,:,:,:)
     private :: mp_allgatherv_cdp5

!!>>> reducing operations

! readucing int
     private :: mp_reduce_int0

! reducing int(:)
     private :: mp_reduce_int1

! reducing int(:,:)
     private :: mp_reduce_int2

! reducing int(:,:,:)
     private :: mp_reduce_int3

! reducing int(:,:,:,:)
     private :: mp_reduce_int4

! reducing int(:,:,:,:,:)
     private :: mp_reduce_int5

! reducing real
     private :: mp_reduce_rdp0

! reducing real(:)
     private :: mp_reduce_rdp1

! reducing real(:,:)
     private :: mp_reduce_rdp2

! reducing real(:,:,:)
     private :: mp_reduce_rdp3

! reducing real(:,:,:,:)
     private :: mp_reduce_rdp4

! reducing real(:,:,:,:,:)
     private :: mp_reduce_rdp5

! reducing complex
     private :: mp_reduce_cdp0

! reducing complex(:)
     private :: mp_reduce_cdp1

! reducing complex(:,:)
     private :: mp_reduce_cdp2

! reducing complex(:,:,:)
     private :: mp_reduce_cdp3

! reducing complex(:,:,:,:)
     private :: mp_reduce_cdp4

! reducing complex(:,:,:,:,:)
     private :: mp_reduce_cdp5

!!>>> allreducing operations

! allreducing int
     private :: mp_allreduce_int0

! allreducing int(:)
     private :: mp_allreduce_int1

! allreducing int(:,:)
     private :: mp_allreduce_int2

! allreducing int(:,:,:)
     private :: mp_allreduce_int3

! allreducing int(:,:,:,:)
     private :: mp_allreduce_int4

! allreducing int(:,:,:,:,:)
     private :: mp_allreduce_int5

! allreducing real
     private :: mp_allreduce_rdp0

! allreducing real(:)
     private :: mp_allreduce_rdp1

! allreducing real(:,:)
     private :: mp_allreduce_rdp2

! allreducing real(:,:,:)
     private :: mp_allreduce_rdp3

! allreducing real(:,:,:,:)
     private :: mp_allreduce_rdp4

! allreducing real(:,:,:,:,:)
     private :: mp_allreduce_rdp5

! allreducing complex
     private :: mp_allreduce_cdp0

! allreducing complex(:)
     private :: mp_allreduce_cdp1

! allreducing complex(:,:)
     private :: mp_allreduce_cdp2

! allreducing complex(:,:,:)
     private :: mp_allreduce_cdp3

! allreducing complex(:,:,:,:)
     private :: mp_allreduce_cdp4

! allreducing complex(:,:,:,:,:)
     private :: mp_allreduce_cdp5

!!>>> error handler

! echo the error information
     private :: mp_error

!!========================================================================
!!>>> declare interface and module procedure                           <<<
!!========================================================================

! we define these interfaces to implement the so called "generic" software
! engineering technique

!!>>> mpi_bcast subroutines
     public :: mp_bcast
     interface mp_bcast
         module procedure mp_bcast_log0
         module procedure mp_bcast_log1
         module procedure mp_bcast_log2
         module procedure mp_bcast_log3
         module procedure mp_bcast_log4
         module procedure mp_bcast_log5

         module procedure mp_bcast_int0
         module procedure mp_bcast_int1
         module procedure mp_bcast_int2
         module procedure mp_bcast_int3
         module procedure mp_bcast_int4
         module procedure mp_bcast_int5

         module procedure mp_bcast_rdp0
         module procedure mp_bcast_rdp1
         module procedure mp_bcast_rdp2
         module procedure mp_bcast_rdp3
         module procedure mp_bcast_rdp4
         module procedure mp_bcast_rdp5

         module procedure mp_bcast_cdp0
         module procedure mp_bcast_cdp1
         module procedure mp_bcast_cdp2
         module procedure mp_bcast_cdp3
         module procedure mp_bcast_cdp4
         module procedure mp_bcast_cdp5
     end interface mp_bcast

!!>>> mpi_gather subroutines
     public :: mp_gather
     interface mp_gather
         module procedure mp_gather_int1
         module procedure mp_gather_int2
         module procedure mp_gather_int3
         module procedure mp_gather_int4
         module procedure mp_gather_int5

         module procedure mp_gather_rdp1
         module procedure mp_gather_rdp2
         module procedure mp_gather_rdp3
         module procedure mp_gather_rdp4
         module procedure mp_gather_rdp5

         module procedure mp_gather_cdp1
         module procedure mp_gather_cdp2
         module procedure mp_gather_cdp3
         module procedure mp_gather_cdp4
         module procedure mp_gather_cdp5
     end interface mp_gather

!!>>> mpi_gatherv subroutines
     public :: mp_gatherv
     interface mp_gatherv
         module procedure mp_gatherv_int1
         module procedure mp_gatherv_int2
         module procedure mp_gatherv_int3
         module procedure mp_gatherv_int4
         module procedure mp_gatherv_int5

         module procedure mp_gatherv_rdp1
         module procedure mp_gatherv_rdp2
         module procedure mp_gatherv_rdp3
         module procedure mp_gatherv_rdp4
         module procedure mp_gatherv_rdp5

         module procedure mp_gatherv_cdp1
         module procedure mp_gatherv_cdp2
         module procedure mp_gatherv_cdp3
         module procedure mp_gatherv_cdp4
         module procedure mp_gatherv_cdp5
     end interface mp_gatherv

!!>>> mpi_allgather subroutines
     public :: mp_allgather
     interface mp_allgather
         module procedure mp_allgather_int1
         module procedure mp_allgather_int2
         module procedure mp_allgather_int3
         module procedure mp_allgather_int4
         module procedure mp_allgather_int5

         module procedure mp_allgather_rdp1
         module procedure mp_allgather_rdp2
         module procedure mp_allgather_rdp3
         module procedure mp_allgather_rdp4
         module procedure mp_allgather_rdp5

         module procedure mp_allgather_cdp1
         module procedure mp_allgather_cdp2
         module procedure mp_allgather_cdp3
         module procedure mp_allgather_cdp4
         module procedure mp_allgather_cdp5
     end interface mp_allgather

!!>>> mpi_allgatherv subroutines
     public :: mp_allgatherv
     interface mp_allgatherv
         module procedure mp_allgatherv_int1
         module procedure mp_allgatherv_int2
         module procedure mp_allgatherv_int3
         module procedure mp_allgatherv_int4
         module procedure mp_allgatherv_int5

         module procedure mp_allgatherv_rdp1
         module procedure mp_allgatherv_rdp2
         module procedure mp_allgatherv_rdp3
         module procedure mp_allgatherv_rdp4
         module procedure mp_allgatherv_rdp5

         module procedure mp_allgatherv_cdp1
         module procedure mp_allgatherv_cdp2
         module procedure mp_allgatherv_cdp3
         module procedure mp_allgatherv_cdp4
         module procedure mp_allgatherv_cdp5
     end interface mp_allgatherv

!!>>> mpi_reduce subroutines
     public :: mp_reduce
     interface mp_reduce
         module procedure mp_reduce_int0
         module procedure mp_reduce_int1
         module procedure mp_reduce_int2
         module procedure mp_reduce_int3
         module procedure mp_reduce_int4
         module procedure mp_reduce_int5

         module procedure mp_reduce_rdp0
         module procedure mp_reduce_rdp1
         module procedure mp_reduce_rdp2
         module procedure mp_reduce_rdp3
         module procedure mp_reduce_rdp4
         module procedure mp_reduce_rdp5

         module procedure mp_reduce_cdp0
         module procedure mp_reduce_cdp1
         module procedure mp_reduce_cdp2
         module procedure mp_reduce_cdp3
         module procedure mp_reduce_cdp4
         module procedure mp_reduce_cdp5
     end interface mp_reduce

!!>>> mpi_allreduce subroutines
     public :: mp_allreduce
     interface mp_allreduce
         module procedure mp_allreduce_int0
         module procedure mp_allreduce_int1
         module procedure mp_allreduce_int2
         module procedure mp_allreduce_int3
         module procedure mp_allreduce_int4
         module procedure mp_allreduce_int5

         module procedure mp_allreduce_rdp0
         module procedure mp_allreduce_rdp1
         module procedure mp_allreduce_rdp2
         module procedure mp_allreduce_rdp3
         module procedure mp_allreduce_rdp4
         module procedure mp_allreduce_rdp5

         module procedure mp_allreduce_cdp0
         module procedure mp_allreduce_cdp1
         module procedure mp_allreduce_cdp2
         module procedure mp_allreduce_cdp3
         module procedure mp_allreduce_cdp4
         module procedure mp_allreduce_cdp5
     end interface mp_allreduce

  contains

!!========================================================================
!!>>> MPI information operations                                       <<<
!!========================================================================

!!
!! @sub mp_info
!!
!! return the current information about mpi environment
!!
     subroutine mp_info()
         implicit none

# if defined (OMPI)

# define STR_MPI 'Good news, your compiler is mpi-compatible (openmpi).'
         write(mystd,'(a)') STR_MPI

# else   /* OMPI */

# define STR_MPI 'Good news, your compiler is mpi-compatible (mpich).'
         write(mystd,'(a)') STR_MPI

# endif  /* OMPI */

         return
     end subroutine mp_info

!!========================================================================
!!>>> MPI initialize and finalize operations                           <<<
!!========================================================================

!!
!! @sub mp_init
!!
!! initialize mpi environment
!!
     subroutine mp_init()
         implicit none

! invoke related mpi subroutines
         call MPI_INIT(ierror)

! handler for return code
         call mp_error('mp_init', ierror)

         return
     end subroutine mp_init

!!
!! @sub mp_finalize
!!
!! finalize mpi environment
!!
     subroutine mp_finalize()
         implicit none

! invoke related mpi subroutines
         call MPI_FINALIZE(ierror)

! handler for return code
         call mp_error('mp_finalize', ierror)

         return
     end subroutine mp_finalize

!!========================================================================
!!>>> MPI setup operations                                             <<<
!!========================================================================

!!
!! @sub mp_comm_rank
!!
!! determine the rank of the current process
!!
     subroutine mp_comm_rank(myid, gid)
         implicit none

! external arguments
         integer, intent(out) :: myid
         integer, optional, intent(in) :: gid

! set current communicator
         if ( present(gid) .eqv. .true. ) then
             group = gid
         else
             group = MPI_COMM_WORLD
         endif ! back if ( present(gid) .eqv. .true. ) block

! invoke related mpi subroutines
         call MPI_COMM_RANK(group, myid, ierror)

! handler for return code
         call mp_error('mp_comm_rank', ierror)

         return
     end subroutine mp_comm_rank

!!
!! @sub mp_comm_size
!!
!! evaluate the number of processes in current communicator
!!
     subroutine mp_comm_size(nprocs, gid)
         implicit none

! external arguments
         integer, intent(out) :: nprocs
         integer, optional, intent(in) :: gid

! set current communicator
         if ( present(gid) .eqv. .true. ) then
             group = gid
         else
             group = MPI_COMM_WORLD
         endif ! back if ( present(gid) .eqv. .true. ) block

! invoke related mpi subroutines
         call MPI_COMM_SIZE(group, nprocs, ierror)

! handler for return code
         call mp_error('mp_comm_size', ierror)

         return
     end subroutine mp_comm_size

!!
!! @sub mp_processor
!!
!! determine the current workstation's name
!!
     subroutine mp_processor(workstation)
         implicit none

! external arguments
         character(len=MPI_MAX_PROCESSOR_NAME), intent(out) :: workstation

! invoke related mpi subroutines
         call MPI_GET_PROCESSOR_NAME(workstation, istat, ierror)

! handler for return code
         call mp_error('mp_processor', ierror)

         return
     end subroutine mp_processor

!!========================================================================
!!>>> MPI cartesian topology operations                                <<<
!!========================================================================

!!
!! @sub mp_dims_create
!!
!! creates a division of processors in a cartesian grid
!!
     subroutine mp_dims_create(nprocs, dims)
         implicit none

! external arguments
         integer, intent(in) :: nprocs
         integer, intent(inout) :: dims(ndims)

! invoke related mpi subroutines
         call MPI_DIMS_CREATE(nprocs, ndims, dims, ierror)

! handler for return code
         call mp_error('mp_dims_create', ierror)

         return
     end subroutine mp_dims_create

!!
!! @sub mp_cart_create
!!
!! makes a new communicator to which topology is cartesian
!!
     subroutine mp_cart_create(dims)
         implicit none

! external arguments
         integer, intent(in) :: dims(ndims)

! invoke related mpi subroutines
! note: mpi_comm_cart should be overwriten in output
         call MPI_CART_CREATE(MPI_COMM_WORLD, ndims, dims, periods, reorder, mpi_comm_cart, ierror)

! handler for return code
         call mp_error('mp_cart_create', ierror)

         return
     end subroutine mp_cart_create

!!
!! @sub mp_cart_coords
!!
!! determines process coords in cartesian topology
!!
     subroutine mp_cart_coords(myid, cx, cy)
         implicit none

! external arguments
         integer, intent(in) :: myid
         integer, intent(out) :: cx
         integer, intent(out) :: cy

! local variables
         integer :: coords(ndims)

! invoke related mpi subroutines
         call MPI_CART_COORDS(mpi_comm_cart, myid, ndims, coords, ierror)

! copy coordinates to cx and cy, in principle, the dimension of coords
! is ndims (now ndims is equal to 2)
         cx = coords(1)
         cy = coords(2)

! handler for return code
         call mp_error('mp_cart_coords', ierror)

         return
     end subroutine mp_cart_coords

!!
!! @sub mp_comm_split_row
!!
!! creates new communicators based on colors and keys
!!
     subroutine mp_comm_split_row(color, key)
         implicit none

! external arguments
         integer :: color
         integer :: key

! invoke related mpi subroutines
! note: mpi_comm_row should be overwritten in output
         call MPI_COMM_SPLIT(mpi_comm_cart, color, key, mpi_comm_row, ierror)

! handler for return code
         call mp_error('mp_comm_split_row', ierror)

         return
     end subroutine mp_comm_split_row

!!
!! @sub mp_comm_split_col
!!
!! creates new communicators based on colors and keys
!!
     subroutine mp_comm_split_col(color, key)
         implicit none

! external arguments
         integer :: color
         integer :: key

! invoke related mpi subroutines
! note: mpi_comm_col should be overwritten in output
         call MPI_COMM_SPLIT(mpi_comm_cart, color, key, mpi_comm_col, ierror)

! handler for return code
         call mp_error('mp_comm_split_col', ierror)

         return
     end subroutine mp_comm_split_col

!!========================================================================
!!>>> MPI barrier operations                                           <<<
!!========================================================================

!!
!! @sub mp_barrier
!!
!! blocks until all process have reached this routine
!!
     subroutine mp_barrier(gid)
         implicit none

! external arguments
         integer, optional, intent(in) :: gid

! set current communicator
         if ( present(gid) .eqv. .true. ) then
             group = gid
         else
             group = MPI_COMM_WORLD
         endif ! back if ( present(gid) .eqv. .true. ) block

! invoke related mpi subroutines
         call MPI_BARRIER(group, ierror)

! handler for return code
         call mp_error('mp_barrier', ierror)

         return
     end subroutine mp_barrier

!!========================================================================
!!>>> MPI time operations                                              <<<
!!========================================================================

!!
!! @sub mp_wtime
!!
!! returns an elapsed time on the calling processor
!!
     subroutine mp_wtime(time)
         implicit none

! external arguments
         real(dp), intent(out) :: time

! invoke related mpi subroutines
         time = MPI_WTIME()

         return
     end subroutine mp_wtime

!!
!! @sub mp_wtick
!!
!! returns the resolution of MPI_Wtime
!!
     subroutine mp_wtick(tick)
         implicit none

! external arguments
         real(dp), intent(out) :: tick

! invoke related mpi subroutines
         tick = MPI_WTICK()

         return
     end subroutine mp_wtick

!!========================================================================
!!>>> MPI collective operations: broadcasting                          <<<
!!========================================================================

!!
!! @sub mp_bcast_log0
!!
!! broadcasts bool from the process with rank "root"
!!
     subroutine mp_bcast_log0(data, root, gid)
         implicit none

! external arguments
         logical, intent(in) :: data
         integer, intent(in) :: root
         integer, optional, intent(in) :: gid

! set current communicator
         if ( present(gid) .eqv. .true. ) then
             group = gid
         else
             group = MPI_COMM_WORLD
         endif ! back if ( present(gid) .eqv. .true. ) block

! barrier until all processes reach here
         call mp_barrier(group)

! invoke realted MPI subroutines
         call MPI_BCAST(data, 1, m_log, root, group, ierror)

! handler for return code
         call mp_error('mp_bcast_log0', ierror)

         return
     end subroutine mp_bcast_log0

!!
!! @sub mp_bcast_log1
!!
!! broadcasts bool(:) from the process with rank "root"
!!
     subroutine mp_bcast_log1(data, root, gid)
         implicit none

! external arguments
         logical, intent(in) :: data(:)
         integer, intent(in) :: root
         integer, optional, intent(in) :: gid

! set current communicator
         if ( present(gid) .eqv. .true. ) then
             group = gid
         else
             group = MPI_COMM_WORLD
         endif ! back if ( present(gid) .eqv. .true. ) block

! barrier until all processes reach here
         call mp_barrier(group)

! setup element count
         isize = size(data)

! invoke realted MPI subroutines
         call MPI_BCAST(data, isize, m_log, root, group, ierror)

! handler for return code
         call mp_error('mp_bcast_log1', ierror)

         return
     end subroutine mp_bcast_log1

!!
!! @sub mp_bcast_log2
!!
!! broadcasts bool2 from the process with rank "root"
!!
     subroutine mp_bcast_log2(data, root, gid)
         implicit none

! external arguments
         logical, intent(in) :: data(:,:)
         integer, intent(in) :: root
         integer, optional, intent(in) :: gid

! set current communicator
         if ( present(gid) .eqv. .true. ) then
             group = gid
         else
             group = MPI_COMM_WORLD
         endif ! back if ( present(gid) .eqv. .true. ) block

! barrier until all processes reach here
         call mp_barrier(group)

! setup element count
         isize = size(data)

! invoke realted MPI subroutines
         call MPI_BCAST(data, isize, m_log, root, group, ierror)

! handler for return code
         call mp_error('mp_bcast_log2', ierror)

         return
     end subroutine mp_bcast_log2

!!
!! @sub mp_bcast_log3
!!
!! broadcasts bool3 from the process with rank "root"
!!
     subroutine mp_bcast_log3(data, root, gid)
         implicit none

! external arguments
         logical, intent(in) :: data(:,:,:)
         integer, intent(in) :: root
         integer, optional, intent(in) :: gid

! set current communicator
         if ( present(gid) .eqv. .true. ) then
             group = gid
         else
             group = MPI_COMM_WORLD
         endif ! back if ( present(gid) .eqv. .true. ) block

! barrier until all processes reach here
         call mp_barrier(group)

! setup element count
         isize = size(data)

! invoke realted MPI subroutines
         call MPI_BCAST(data, isize, m_log, root, group, ierror)

! handler for return code
         call mp_error('mp_bcast_log3', ierror)

         return
     end subroutine mp_bcast_log3

!!
!! @sub mp_bcast_log4
!!
!! broadcasts bool4 from the process with rank "root"
!!
     subroutine mp_bcast_log4(data, root, gid)
         implicit none

! external arguments
         logical, intent(in) :: data(:,:,:,:)
         integer, intent(in) :: root
         integer, optional, intent(in) :: gid

! set current communicator
         if ( present(gid) .eqv. .true. ) then
             group = gid
         else
             group = MPI_COMM_WORLD
         endif ! back if ( present(gid) .eqv. .true. ) block

! barrier until all processes reach here
         call mp_barrier(group)

! setup element count
         isize = size(data)

! invoke realted MPI subroutines
         call MPI_BCAST(data, isize, m_log, root, group, ierror)

! handler for return code
         call mp_error('mp_bcast_log4', ierror)

         return
     end subroutine mp_bcast_log4

!!
!! @sub mp_bcast_log5
!!
!! broadcasts bool5 from the process with rank "root"
!!
     subroutine mp_bcast_log5(data, root, gid)
         implicit none

! external arguments
         logical, intent(in) :: data(:,:,:,:,:)
         integer, intent(in) :: root
         integer, optional, intent(in) :: gid

! set current communicator
         if ( present(gid) .eqv. .true. ) then
             group = gid
         else
             group = MPI_COMM_WORLD
         endif ! back if ( present(gid) .eqv. .true. ) block

! barrier until all processes reach here
         call mp_barrier(group)

! setup element count
         isize = size(data)

! invoke realted MPI subroutines
         call MPI_BCAST(data, isize, m_log, root, group, ierror)

! handler for return code
         call mp_error('mp_bcast_log5', ierror)

         return
     end subroutine mp_bcast_log5

!!
!! @sub mp_bcast_int0
!!
!! broadcasts int from the process with rank "root"
!!
     subroutine mp_bcast_int0(data, root, gid)
         implicit none

! external arguments
         integer, intent(in) :: data
         integer, intent(in) :: root
         integer, optional, intent(in) :: gid

! set current communicator
         if ( present(gid) .eqv. .true. ) then
             group = gid
         else
             group = MPI_COMM_WORLD
         endif ! back if ( present(gid) .eqv. .true. ) block

! barrier until all processes reach here
         call mp_barrier(group)

! invoke realted MPI subroutines
         call MPI_BCAST(data, 1, m_int, root, group, ierror)

! handler for return code
         call mp_error('mp_bcast_int0', ierror)

         return
     end subroutine mp_bcast_int0

!!
!! @sub mp_bcast_int1
!!
!! broadcasts int(:) from the process with rank "root"
!!
     subroutine mp_bcast_int1(data, root, gid)
         implicit none

! external arguments
         integer, intent(in) :: data(:)
         integer, intent(in) :: root
         integer, optional, intent(in) :: gid

! set current communicator
         if ( present(gid) .eqv. .true. ) then
             group = gid
         else
             group = MPI_COMM_WORLD
         endif ! back if ( present(gid) .eqv. .true. ) block

! barrier until all processes reach here
         call mp_barrier(group)

! setup element count
         isize = size(data)

! invoke realted MPI subroutines
         call MPI_BCAST(data, isize, m_int, root, group, ierror)

! handler for return code
         call mp_error('mp_bcast_int1', ierror)

         return
     end subroutine mp_bcast_int1

!!
!! @sub mp_bcast_int2
!!
!! broadcasts int2 from the process with rank "root"
!!
     subroutine mp_bcast_int2(data, root, gid)
         implicit none

! external arguments
         integer, intent(in) :: data(:,:)
         integer, intent(in) :: root
         integer, optional, intent(in) :: gid

! set current communicator
         if ( present(gid) .eqv. .true. ) then
             group = gid
         else
             group = MPI_COMM_WORLD
         endif ! back if ( present(gid) .eqv. .true. ) block

! barrier until all processes reach here
         call mp_barrier(group)

! setup element count
         isize = size(data)

! invoke realted MPI subroutines
         call MPI_BCAST(data, isize, m_int, root, group, ierror)

! handler for return code
         call mp_error('mp_bcast_int2', ierror)

         return
     end subroutine mp_bcast_int2

!!
!! @sub mp_bcast_int3
!!
!! broadcasts int3 from the process with rank "root"
!!
     subroutine mp_bcast_int3(data, root, gid)
         implicit none

! external arguments
         integer, intent(in) :: data(:,:,:)
         integer, intent(in) :: root
         integer, optional, intent(in) :: gid

! set current communicator
         if ( present(gid) .eqv. .true. ) then
             group = gid
         else
             group = MPI_COMM_WORLD
         endif ! back if ( present(gid) .eqv. .true. ) block

! barrier until all processes reach here
         call mp_barrier(group)

! setup element count
         isize = size(data)

! invoke realted MPI subroutines
         call MPI_BCAST(data, isize, m_int, root, group, ierror)

! handler for return code
         call mp_error('mp_bcast_int3', ierror)

         return
     end subroutine mp_bcast_int3

!!
!! @sub mp_bcast_int4
!!
!! broadcasts int4 from the process with rank "root"
!!
     subroutine mp_bcast_int4(data, root, gid)
         implicit none

! external arguments
         integer, intent(in) :: data(:,:,:,:)
         integer, intent(in) :: root
         integer, optional, intent(in) :: gid

! set current communicator
         if ( present(gid) .eqv. .true. ) then
             group = gid
         else
             group = MPI_COMM_WORLD
         endif ! back if ( present(gid) .eqv. .true. ) block

! barrier until all processes reach here
         call mp_barrier(group)

! setup element count
         isize = size(data)

! invoke realted MPI subroutines
         call MPI_BCAST(data, isize, m_int, root, group, ierror)

! handler for return code
         call mp_error('mp_bcast_int4', ierror)

         return
     end subroutine mp_bcast_int4

!!
!! @sub mp_bcast_int5
!!
!! broadcasts int5 from the process with rank "root"
!!
     subroutine mp_bcast_int5(data, root, gid)
         implicit none

! external arguments
         integer, intent(in) :: data(:,:,:,:,:)
         integer, intent(in) :: root
         integer, optional, intent(in) :: gid

! set current communicator
         if ( present(gid) .eqv. .true. ) then
             group = gid
         else
             group = MPI_COMM_WORLD
         endif ! back if ( present(gid) .eqv. .true. ) block

! barrier until all processes reach here
         call mp_barrier(group)

! setup element count
         isize = size(data)

! invoke realted MPI subroutines
         call MPI_BCAST(data, isize, m_int, root, group, ierror)

! handler for return code
         call mp_error('mp_bcast_int5', ierror)

         return
     end subroutine mp_bcast_int5

!!
!! @sub mp_bcast_rdp0
!!
!! broadcasts real from the process with rank "root"
!!
     subroutine mp_bcast_rdp0(data, root, gid)
         implicit none

! external arguments
         real(dp), intent(in) :: data
         integer, intent(in) :: root
         integer, optional, intent(in) :: gid

! set current communicator
         if ( present(gid) .eqv. .true. ) then
             group = gid
         else
             group = MPI_COMM_WORLD
         endif ! back if ( present(gid) .eqv. .true. ) block

! barrier until all processes reach here
         call mp_barrier(group)

! invoke related mpi subroutines
         call MPI_BCAST(data, 1, m_rdp, root, group, ierror)

! handler for return code
         call mp_error('mp_bcast_rdp0', ierror)

         return
     end subroutine mp_bcast_rdp0

!!
!! @sub mp_bcast_rdp1
!!
!! broadcasts real(:) from the process with rank "root"
!!
     subroutine mp_bcast_rdp1(data, root, gid)
         implicit none

! external arguments
         real(dp), intent(in) :: data(:)
         integer, intent(in) :: root
         integer, optional, intent(in) :: gid

! set current communicator
         if ( present(gid) .eqv. .true. ) then
             group = gid
         else
             group = MPI_COMM_WORLD
         endif ! back if ( present(gid) .eqv. .true. ) block

! barrier until all processes reach here
         call mp_barrier(group)

! setup element count
         isize = size(data)

! invoke related mpi subroutines
         call MPI_BCAST(data, isize, m_rdp, root, group, ierror)

! handler for return code
         call mp_error('mp_bcast_rdp1', ierror)

         return
     end subroutine mp_bcast_rdp1

!!
!! @sub mp_bcast_rdp2
!!
!! broadcasts real2 from the process with rank "root"
!!
     subroutine mp_bcast_rdp2(data, root, gid)
         implicit none

! external arguments
         real(dp), intent(in) :: data(:,:)
         integer, intent(in) :: root
         integer, optional, intent(in) :: gid

! set current communicator
         if ( present(gid) .eqv. .true. ) then
             group = gid
         else
             group = MPI_COMM_WORLD
         endif ! back if ( present(gid) .eqv. .true. ) block

! barrier until all processes reach here
         call mp_barrier(group)

! setup element count
         isize = size(data)

! invoke related mpi subroutines
         call MPI_BCAST(data, isize, m_rdp, root, group, ierror)

! handler for return code
         call mp_error('mp_bcast_rdp2', ierror)

         return
     end subroutine mp_bcast_rdp2

!!
!! @sub mp_bcast_rdp3
!!
!! broadcasts real3 from the process with rank "root"
!!
     subroutine mp_bcast_rdp3(data, root, gid)
         implicit none

! external arguments
         real(dp), intent(in) :: data(:,:,:)
         integer, intent(in) :: root
         integer, optional, intent(in) :: gid

! set current communicator
         if ( present(gid) .eqv. .true. ) then
             group = gid
         else
             group = MPI_COMM_WORLD
         endif ! back if ( present(gid) .eqv. .true. ) block

! barrier until all processes reach here
         call mp_barrier(group)

! setup element count
         isize = size(data)

! invoke related mpi subroutines
         call MPI_BCAST(data, isize, m_rdp, root, group, ierror)

! handler for return code
         call mp_error('mp_bcast_rdp3', ierror)

         return
     end subroutine mp_bcast_rdp3

!!
!! @sub mp_bcast_rdp4
!!
!! broadcasts real4 from the process with rank "root"
!!
     subroutine mp_bcast_rdp4(data, root, gid)
         implicit none

! external arguments
         real(dp), intent(in) :: data(:,:,:,:)
         integer, intent(in) :: root
         integer, optional, intent(in) :: gid

! set current communicator
         if ( present(gid) .eqv. .true. ) then
             group = gid
         else
             group = MPI_COMM_WORLD
         endif ! back if ( present(gid) .eqv. .true. ) block

! barrier until all processes reach here
         call mp_barrier(group)

! setup element count
         isize = size(data)

! invoke related mpi subroutines
         call MPI_BCAST(data, isize, m_rdp, root, group, ierror)

! handler for return code
         call mp_error('mp_bcast_rdp4', ierror)

         return
     end subroutine mp_bcast_rdp4

!!
!! @sub mp_bcast_rdp5
!!
!! broadcasts real5 from the process with rank "root"
!!
     subroutine mp_bcast_rdp5(data, root, gid)
         implicit none

! external arguments
         real(dp), intent(in) :: data(:,:,:,:,:)
         integer, intent(in) :: root
         integer, optional, intent(in) :: gid

! set current communicator
         if ( present(gid) .eqv. .true. ) then
             group = gid
         else
             group = MPI_COMM_WORLD
         endif ! back if ( present(gid) .eqv. .true. ) block

! barrier until all processes reach here
         call mp_barrier(group)

! setup element count
         isize = size(data)

! invoke related mpi subroutines
         call MPI_BCAST(data, isize, m_rdp, root, group, ierror)

! handler for return code
         call mp_error('mp_bcast_rdp5', ierror)

         return
     end subroutine mp_bcast_rdp5

!!
!! @sub mp_bcast_cdp0
!!
!! broadcasts complex from the process with rank "root"
!!
     subroutine mp_bcast_cdp0(data, root, gid)
         implicit none

! external arguments
         complex(dp), intent(in) :: data
         integer, intent(in) :: root
         integer, optional, intent(in) :: gid

! set current communicator
         if ( present(gid) .eqv. .true. ) then
             group = gid
         else
             group = MPI_COMM_WORLD
         endif ! back if ( present(gid) .eqv. .true. ) block

! barrier until all processes reach here
         call mp_barrier(group)

! invoke related mpi subroutines
         call MPI_BCAST(data, 1, m_cdp, root, group, ierror)

! handler for return code
         call mp_error('mp_bcast_cdp0', ierror)

         return
     end subroutine mp_bcast_cdp0

!!
!! @sub mp_bcast_cdp1
!!
!! broadcasts complex(:) from the process with rank "root"
!!
     subroutine mp_bcast_cdp1(data, root, gid)
         implicit none

! external arguments
         complex(dp), intent(in) :: data(:)
         integer, intent(in) :: root
         integer, optional, intent(in) :: gid

! set current communicator
         if ( present(gid) .eqv. .true. ) then
             group = gid
         else
             group = MPI_COMM_WORLD
         endif ! back if ( present(gid) .eqv. .true. ) block

! barrier until all processes reach here
         call mp_barrier(group)

! setup element count
         isize = size(data)

! invoke related mpi subroutines
         call MPI_BCAST(data, isize, m_cdp, root, group, ierror)

! handler for return code
         call mp_error('mp_bcast_cdp1', ierror)

         return
     end subroutine mp_bcast_cdp1

!!
!! @sub mp_bcast_cdp2
!!
!! broadcasts complex2 from the process with rank "root"
!!
     subroutine mp_bcast_cdp2(data, root, gid)
         implicit none

! external arguments
         complex(dp), intent(in) :: data(:,:)
         integer, intent(in) :: root
         integer, optional, intent(in) :: gid

! set current communicator
         if ( present(gid) .eqv. .true. ) then
             group = gid
         else
             group = MPI_COMM_WORLD
         endif ! back if ( present(gid) .eqv. .true. ) block

! barrier until all processes reach here
         call mp_barrier(group)

! setup element count
         isize = size(data)

! invoke related mpi subroutines
         call MPI_BCAST(data, isize, m_cdp, root, group, ierror)

! handler for return code
         call mp_error('mp_bcast_cdp2', ierror)

         return
     end subroutine mp_bcast_cdp2

!!
!! @sub mp_bcast_cdp3
!!
!! broadcasts complex3 from the process with rank "root"
!!
     subroutine mp_bcast_cdp3(data, root, gid)
         implicit none

! external arguments
         complex(dp), intent(in) :: data(:,:,:)
         integer, intent(in) :: root
         integer, optional, intent(in) :: gid

! set current communicator
         if ( present(gid) .eqv. .true. ) then
             group = gid
         else
             group = MPI_COMM_WORLD
         endif ! back if ( present(gid) .eqv. .true. ) block

! barrier until all processes reach here
         call mp_barrier(group)

! setup element count
         isize = size(data)

! invoke related mpi subroutines
         call MPI_BCAST(data, isize, m_cdp, root, group, ierror)

! handler for return code
         call mp_error('mp_bcast_cdp3', ierror)

         return
     end subroutine mp_bcast_cdp3

!!
!! @sub mp_bcast_cdp4
!!
!! broadcasts complex4 from the process with rank "root"
!!
     subroutine mp_bcast_cdp4(data, root, gid)
         implicit none

! external arguments
         complex(dp), intent(in) :: data(:,:,:,:)
         integer, intent(in) :: root
         integer, optional, intent(in) :: gid

! set current communicator
         if ( present(gid) .eqv. .true. ) then
             group = gid
         else
             group = MPI_COMM_WORLD
         endif ! back if ( present(gid) .eqv. .true. ) block

! barrier until all processes reach here
         call mp_barrier(group)

! setup element count
         isize = size(data)

! invoke related mpi subroutines
         call MPI_BCAST(data, isize, m_cdp, root, group, ierror)

! handler for return code
         call mp_error('mp_bcast_cdp4', ierror)

         return
     end subroutine mp_bcast_cdp4

!!
!! @sub mp_bcast_cdp5
!!
!! broadcasts complex5 from the process with rank "root"
!!
     subroutine mp_bcast_cdp5(data, root, gid)
         implicit none

! external arguments
         complex(dp), intent(in) :: data(:,:,:,:,:)
         integer, intent(in) :: root
         integer, optional, intent(in) :: gid

! set current communicator
         if ( present(gid) .eqv. .true. ) then
             group = gid
         else
             group = MPI_COMM_WORLD
         endif ! back if ( present(gid) .eqv. .true. ) block

! barrier until all processes reach here
         call mp_barrier(group)

! setup element count
         isize = size(data)

! invoke related mpi subroutines
         call MPI_BCAST(data, isize, m_cdp, root, group, ierror)

! handler for return code
         call mp_error('mp_bcast_cdp5', ierror)

         return
     end subroutine mp_bcast_cdp5

!!========================================================================
!!>>> MPI collective operations : gathering                            <<<
!!========================================================================

!!
!! @sub mp_gather_int1
!!
!! gather integer data from every processes to rank 0
!!
     subroutine mp_gather_int1(send, data, root, gid)
         implicit none

! external arguments
         integer, intent(in) :: send(:)
         integer, intent(inout) :: data(:)

         integer, intent(in) :: root
         integer, optional, intent(in) :: gid

! set current communicator
         if ( present(gid) .eqv. .true. ) then
             group = gid
         else
             group = MPI_COMM_WORLD
         endif ! back if ( present(gid) .eqv. .true. ) block

! barrier until all processes reach here
         call mp_barrier(group)

! setup element count
         isize = size(send)

! invoke related mpi subroutines
         call MPI_GATHER(send, isize, m_int, data, isize, m_int, root, group, ierror)

! handler for return code
         call mp_error('mp_gather_int1', ierror)

         return
     end subroutine mp_gather_int1

!!
!! @sub mp_gather_int2
!!
!! gather integer data from every processes to rank 0
!!
     subroutine mp_gather_int2(send, data, root, gid)
         implicit none

! external arguments
         integer, intent(in) :: send(:,:)
         integer, intent(inout) :: data(:,:)

         integer, intent(in) :: root
         integer, optional, intent(in) :: gid

! set current communicator
         if ( present(gid) .eqv. .true. ) then
             group = gid
         else
             group = MPI_COMM_WORLD
         endif ! back if ( present(gid) .eqv. .true. ) block

! barrier until all processes reach here
         call mp_barrier(group)

! setup element count
         isize = size(send)

! invoke related mpi subroutines
         call MPI_GATHER(send, isize, m_int, data, isize, m_int, root, group, ierror)

! handler for return code
         call mp_error('mp_gather_int2', ierror)

         return
     end subroutine mp_gather_int2

!!
!! @sub mp_gather_int3
!!
!! gather integer data from every processes to rank 0
!!
     subroutine mp_gather_int3(send, data, root, gid)
         implicit none

! external arguments
         integer, intent(in) :: send(:,:,:)
         integer, intent(inout) :: data(:,:,:)

         integer, intent(in) :: root
         integer, optional, intent(in) :: gid

! set current communicator
         if ( present(gid) .eqv. .true. ) then
             group = gid
         else
             group = MPI_COMM_WORLD
         endif ! back if ( present(gid) .eqv. .true. ) block

! barrier until all processes reach here
         call mp_barrier(group)

! setup element count
         isize = size(send)

! invoke related mpi subroutines
         call MPI_GATHER(send, isize, m_int, data, isize, m_int, root, group, ierror)

! handler for return code
         call mp_error('mp_gather_int3', ierror)

         return
     end subroutine mp_gather_int3

!!
!! @sub mp_gather_int4
!!
!! gather integer data from every processes to rank 0
!!
     subroutine mp_gather_int4(send, data, root, gid)
         implicit none

! external arguments
         integer, intent(in) :: send(:,:,:,:)
         integer, intent(inout) :: data(:,:,:,:)

         integer, intent(in) :: root
         integer, optional, intent(in) :: gid

! set current communicator
         if ( present(gid) .eqv. .true. ) then
             group = gid
         else
             group = MPI_COMM_WORLD
         endif ! back if ( present(gid) .eqv. .true. ) block

! barrier until all processes reach here
         call mp_barrier(group)

! setup element count
         isize = size(send)

! invoke related mpi subroutines
         call MPI_GATHER(send, isize, m_int, data, isize, m_int, root, group, ierror)

! handler for return code
         call mp_error('mp_gather_int4', ierror)

         return
     end subroutine mp_gather_int4

!!
!! @sub mp_gather_int5
!!
!! gather integer data from every processes to rank 0
!!
     subroutine mp_gather_int5(send, data, root, gid)
         implicit none

! external arguments
         integer, intent(in) :: send(:,:,:,:,:)
         integer, intent(inout) :: data(:,:,:,:,:)

         integer, intent(in) :: root
         integer, optional, intent(in) :: gid

! set current communicator
         if ( present(gid) .eqv. .true. ) then
             group = gid
         else
             group = MPI_COMM_WORLD
         endif ! back if ( present(gid) .eqv. .true. ) block

! barrier until all processes reach here
         call mp_barrier(group)

! setup element count
         isize = size(send)

! invoke related mpi subroutines
         call MPI_GATHER(send, isize, m_int, data, isize, m_int, root, group, ierror)

! handler for return code
         call mp_error('mp_gather_int5', ierror)

         return
     end subroutine mp_gather_int5

!!
!! @sub mp_gather_rdp1
!!
!! gather real(dp) data from every processes to rank 0
!!
     subroutine mp_gather_rdp1(send, data, root, gid)
         implicit none

! external arguments
         real(dp), intent(in) :: send(:)
         real(dp), intent(inout) :: data(:)

         integer, intent(in) :: root
         integer, optional, intent(in) :: gid

! set current communicator
         if ( present(gid) .eqv. .true. ) then
             group = gid
         else
             group = MPI_COMM_WORLD
         endif ! back if ( present(gid) .eqv. .true. ) block

! barrier until all processes reach here
         call mp_barrier(group)

! setup element count
         isize = size(send)

! invoke related mpi subroutines
         call MPI_GATHER(send, isize, m_rdp, data, isize, m_rdp, root, group, ierror)

! handler for return code
         call mp_error('mp_gather_rdp1', ierror)

         return
     end subroutine mp_gather_rdp1

!!
!! @sub mp_gather_rdp2
!!
!! gather real(dp) data from every processes to rank 0
!!
     subroutine mp_gather_rdp2(send, data, root, gid)
         implicit none

! external arguments
         real(dp), intent(in) :: send(:,:)
         real(dp), intent(inout) :: data(:,:)

         integer, intent(in) :: root
         integer, optional, intent(in) :: gid

! set current communicator
         if ( present(gid) .eqv. .true. ) then
             group = gid
         else
             group = MPI_COMM_WORLD
         endif ! back if ( present(gid) .eqv. .true. ) block

! barrier until all processes reach here
         call mp_barrier(group)

! setup element count
         isize = size(send)

! invoke related mpi subroutines
         call MPI_GATHER(send, isize, m_rdp, data, isize, m_rdp, root, group, ierror)

! handler for return code
         call mp_error('mp_gather_rdp2', ierror)

         return
     end subroutine mp_gather_rdp2

!!
!! @sub mp_gather_rdp3
!!
!! gather real(dp) data from every processes to rank 0
!!
     subroutine mp_gather_rdp3(send, data, root, gid)
         implicit none

! external arguments
         real(dp), intent(in) :: send(:,:,:)
         real(dp), intent(inout) :: data(:,:,:)

         integer, intent(in) :: root
         integer, optional, intent(in) :: gid

! set current communicator
         if ( present(gid) .eqv. .true. ) then
             group = gid
         else
             group = MPI_COMM_WORLD
         endif ! back if ( present(gid) .eqv. .true. ) block

! barrier until all processes reach here
         call mp_barrier(group)

! setup element count
         isize = size(send)

! invoke related mpi subroutines
         call MPI_GATHER(send, isize, m_rdp, data, isize, m_rdp, root, group, ierror)

! handler for return code
         call mp_error('mp_gather_rdp3', ierror)

         return
     end subroutine mp_gather_rdp3

!!
!! @sub mp_gather_rdp4
!!
!! gather real(dp) data from every processes to rank 0
!!
     subroutine mp_gather_rdp4(send, data, root, gid)
         implicit none

! external arguments
         real(dp), intent(in) :: send(:,:,:,:)
         real(dp), intent(inout) :: data(:,:,:,:)

         integer, intent(in) :: root
         integer, optional, intent(in) :: gid

! set current communicator
         if ( present(gid) .eqv. .true. ) then
             group = gid
         else
             group = MPI_COMM_WORLD
         endif ! back if ( present(gid) .eqv. .true. ) block

! barrier until all processes reach here
         call mp_barrier(group)

! setup element count
         isize = size(send)

! invoke related mpi subroutines
         call MPI_GATHER(send, isize, m_rdp, data, isize, m_rdp, root, group, ierror)

! handler for return code
         call mp_error('mp_gather_rdp4', ierror)

         return
     end subroutine mp_gather_rdp4

!!
!! @sub mp_gather_rdp5
!!
!! gather real(dp) data from every processes to rank 0
!!
     subroutine mp_gather_rdp5(send, data, root, gid)
         implicit none

! external arguments
         real(dp), intent(in) :: send(:,:,:,:,:)
         real(dp), intent(inout) :: data(:,:,:,:,:)

         integer, intent(in) :: root
         integer, optional, intent(in) :: gid

! set current communicator
         if ( present(gid) .eqv. .true. ) then
             group = gid
         else
             group = MPI_COMM_WORLD
         endif ! back if ( present(gid) .eqv. .true. ) block

! barrier until all processes reach here
         call mp_barrier(group)

! setup element count
         isize = size(send)

! invoke related mpi subroutines
         call MPI_GATHER(send, isize, m_rdp, data, isize, m_rdp, root, group, ierror)

! handler for return code
         call mp_error('mp_gather_rdp5', ierror)

         return
     end subroutine mp_gather_rdp5

!!
!! @sub mp_gather_cdp1
!!
!! gather complex(dp) data from every processes to rank 0
!!
     subroutine mp_gather_cdp1(send, data, root, gid)
         implicit none

! external arguments
         complex(dp), intent(in) :: send(:)
         complex(dp), intent(inout) :: data(:)

         integer, intent(in) :: root
         integer, optional, intent(in) :: gid

! set current communicator
         if ( present(gid) .eqv. .true. ) then
             group = gid
         else
             group = MPI_COMM_WORLD
         endif ! back if ( present(gid) .eqv. .true. ) block

! barrier until all processes reach here
         call mp_barrier(group)

! setup element count
         isize = size(send)

! invoke related mpi subroutines
         call MPI_GATHER(send, isize, m_cdp, data, isize, m_cdp, root, group, ierror)

! handler for return code
         call mp_error('mp_gather_cdp1', ierror)

         return
     end subroutine mp_gather_cdp1

!!
!! @sub mp_gather_cdp2
!!
!! gather complex(dp) data from every processes to rank 0
!!
     subroutine mp_gather_cdp2(send, data, root, gid)
         implicit none

! external arguments
         complex(dp), intent(in) :: send(:,:)
         complex(dp), intent(inout) :: data(:,:)

         integer, intent(in) :: root
         integer, optional, intent(in) :: gid

! set current communicator
         if ( present(gid) .eqv. .true. ) then
             group = gid
         else
             group = MPI_COMM_WORLD
         endif ! back if ( present(gid) .eqv. .true. ) block

! barrier until all processes reach here
         call mp_barrier(group)

! setup element count
         isize = size(send)

! invoke related mpi subroutines
         call MPI_GATHER(send, isize, m_cdp, data, isize, m_cdp, root, group, ierror)

! handler for return code
         call mp_error('mp_gather_cdp2', ierror)

         return
     end subroutine mp_gather_cdp2

!!
!! @sub mp_gather_cdp3
!!
!! gather complex(dp) data from every processes to rank 0
!!
     subroutine mp_gather_cdp3(send, data, root, gid)
         implicit none

! external arguments
         complex(dp), intent(in) :: send(:,:,:)
         complex(dp), intent(inout) :: data(:,:,:)

         integer, intent(in) :: root
         integer, optional, intent(in) :: gid

! set current communicator
         if ( present(gid) .eqv. .true. ) then
             group = gid
         else
             group = MPI_COMM_WORLD
         endif ! back if ( present(gid) .eqv. .true. ) block

! barrier until all processes reach here
         call mp_barrier(group)

! setup element count
         isize = size(send)

! invoke related mpi subroutines
         call MPI_GATHER(send, isize, m_cdp, data, isize, m_cdp, root, group, ierror)

! handler for return code
         call mp_error('mp_gather_cdp3', ierror)

         return
     end subroutine mp_gather_cdp3

!!
!! @sub mp_gather_cdp4
!!
!! gather complex(dp) data from every processes to rank 0
!!
     subroutine mp_gather_cdp4(send, data, root, gid)
         implicit none

! external arguments
         complex(dp), intent(in) :: send(:,:,:,:)
         complex(dp), intent(inout) :: data(:,:,:,:)

         integer, intent(in) :: root
         integer, optional, intent(in) :: gid

! set current communicator
         if ( present(gid) .eqv. .true. ) then
             group = gid
         else
             group = MPI_COMM_WORLD
         endif ! back if ( present(gid) .eqv. .true. ) block

! barrier until all processes reach here
         call mp_barrier(group)

! setup element count
         isize = size(send)

! invoke related mpi subroutines
         call MPI_GATHER(send, isize, m_cdp, data, isize, m_cdp, root, group, ierror)

! handler for return code
         call mp_error('mp_gather_cdp4', ierror)

         return
     end subroutine mp_gather_cdp4

!!
!! @sub mp_gather_cdp5
!!
!! gather complex(dp) data from every processes to rank 0
!!
     subroutine mp_gather_cdp5(send, data, root, gid)
         implicit none

! external arguments
         complex(dp), intent(in) :: send(:,:,:,:,:)
         complex(dp), intent(inout) :: data(:,:,:,:,:)

         integer, intent(in) :: root
         integer, optional, intent(in) :: gid

! set current communicator
         if ( present(gid) .eqv. .true. ) then
             group = gid
         else
             group = MPI_COMM_WORLD
         endif ! back if ( present(gid) .eqv. .true. ) block

! barrier until all processes reach here
         call mp_barrier(group)

! setup element count
         isize = size(send)

! invoke related mpi subroutines
         call MPI_GATHER(send, isize, m_cdp, data, isize, m_cdp, root, group, ierror)

! handler for return code
         call mp_error('mp_gather_cdp5', ierror)

         return
     end subroutine mp_gather_cdp5

!!========================================================================
!!>>> MPI collective operations : gatherving                           <<<
!!========================================================================

!!
!! @sub mp_gatherv_int1
!!
!! gather integer data from every processes to rank 0
!!
     subroutine mp_gatherv_int1(send, data, recv, disp, root, gid)
         implicit none

! external arguments
         integer, intent(in) :: send(:)
         integer, intent(inout) :: data(:)

         integer, intent(in) :: recv(:)
         integer, intent(in) :: disp(:)

         integer, intent(in) :: root
         integer, optional, intent(in) :: gid

! set current communicator
         if ( present(gid) .eqv. .true. ) then
             group = gid
         else
             group = MPI_COMM_WORLD
         endif ! back if ( present(gid) .eqv. .true. ) block

! barrier until all processes reach here
         call mp_barrier(group)

! setup element count
         isize = size(send)

! invoke related mpi subroutines
         call MPI_GATHERV(send, isize, m_int, data, recv, disp, m_int, root, group, ierror)

! handler for return code
         call mp_error('mp_gatherv_int1', ierror)

         return
     end subroutine mp_gatherv_int1

!!
!! @sub mp_gatherv_int2
!!
!! gather integer data from every processes to rank 0
!!
     subroutine mp_gatherv_int2(send, data, recv, disp, root, gid)
         implicit none

! external arguments
         integer, intent(in) :: send(:,:)
         integer, intent(inout) :: data(:,:)

         integer, intent(in) :: recv(:)
         integer, intent(in) :: disp(:)

         integer, intent(in) :: root
         integer, optional, intent(in) :: gid

! set current communicator
         if ( present(gid) .eqv. .true. ) then
             group = gid
         else
             group = MPI_COMM_WORLD
         endif ! back if ( present(gid) .eqv. .true. ) block

! barrier until all processes reach here
         call mp_barrier(group)

! setup element count
         isize = size(send)

! invoke related mpi subroutines
         call MPI_GATHERV(send, isize, m_int, data, recv, disp, m_int, root, group, ierror)

! handler for return code
         call mp_error('mp_gatherv_int2', ierror)

         return
     end subroutine mp_gatherv_int2

!!
!! @sub mp_gatherv_int3
!!
!! gather integer data from every processes to rank 0
!!
     subroutine mp_gatherv_int3(send, data, recv, disp, root, gid)
         implicit none

! external arguments
         integer, intent(in) :: send(:,:,:)
         integer, intent(inout) :: data(:,:,:)

         integer, intent(in) :: recv(:)
         integer, intent(in) :: disp(:)

         integer, intent(in) :: root
         integer, optional, intent(in) :: gid

! set current communicator
         if ( present(gid) .eqv. .true. ) then
             group = gid
         else
             group = MPI_COMM_WORLD
         endif ! back if ( present(gid) .eqv. .true. ) block

! barrier until all processes reach here
         call mp_barrier(group)

! setup element count
         isize = size(send)

! invoke related mpi subroutines
         call MPI_GATHERV(send, isize, m_int, data, recv, disp, m_int, root, group, ierror)

! handler for return code
         call mp_error('mp_gatherv_int3', ierror)

         return
     end subroutine mp_gatherv_int3

!!
!! @sub mp_gatherv_int4
!!
!! gather integer data from every processes to rank 0
!!
     subroutine mp_gatherv_int4(send, data, recv, disp, root, gid)
         implicit none

! external arguments
         integer, intent(in) :: send(:,:,:,:)
         integer, intent(inout) :: data(:,:,:,:)

         integer, intent(in) :: recv(:)
         integer, intent(in) :: disp(:)

         integer, intent(in) :: root
         integer, optional, intent(in) :: gid

! set current communicator
         if ( present(gid) .eqv. .true. ) then
             group = gid
         else
             group = MPI_COMM_WORLD
         endif ! back if ( present(gid) .eqv. .true. ) block

! barrier until all processes reach here
         call mp_barrier(group)

! setup element count
         isize = size(send)

! invoke related mpi subroutines
         call MPI_GATHERV(send, isize, m_int, data, recv, disp, m_int, root, group, ierror)

! handler for return code
         call mp_error('mp_gatherv_int4', ierror)

         return
     end subroutine mp_gatherv_int4

!!
!! @sub mp_gatherv_int5
!!
!! gather integer data from every processes to rank 0
!!
     subroutine mp_gatherv_int5(send, data, recv, disp, root, gid)
         implicit none

! external arguments
         integer, intent(in) :: send(:,:,:,:,:)
         integer, intent(inout) :: data(:,:,:,:,:)

         integer, intent(in) :: recv(:)
         integer, intent(in) :: disp(:)

         integer, intent(in) :: root
         integer, optional, intent(in) :: gid

! set current communicator
         if ( present(gid) .eqv. .true. ) then
             group = gid
         else
             group = MPI_COMM_WORLD
         endif ! back if ( present(gid) .eqv. .true. ) block

! barrier until all processes reach here
         call mp_barrier(group)

! setup element count
         isize = size(send)

! invoke related mpi subroutines
         call MPI_GATHERV(send, isize, m_int, data, recv, disp, m_int, root, group, ierror)

! handler for return code
         call mp_error('mp_gatherv_int5', ierror)

         return
     end subroutine mp_gatherv_int5

!!
!! @sub mp_gatherv_rdp1
!!
!! gather real(dp) data from every processes to rank 0
!!
     subroutine mp_gatherv_rdp1(send, data, recv, disp, root, gid)
         implicit none

! external arguments
         real(dp), intent(in) :: send(:)
         real(dp), intent(inout) :: data(:)

         integer, intent(in) :: recv(:)
         integer, intent(in) :: disp(:)

         integer, intent(in) :: root
         integer, optional, intent(in) :: gid

! set current communicator
         if ( present(gid) .eqv. .true. ) then
             group = gid
         else
             group = MPI_COMM_WORLD
         endif ! back if ( present(gid) .eqv. .true. ) block

! barrier until all processes reach here
         call mp_barrier(group)

! setup element count
         isize = size(send)

! invoke related mpi subroutines
         call MPI_GATHERV(send, isize, m_rdp, data, recv, disp, m_rdp, root, group, ierror)

! handler for return code
         call mp_error('mp_gatherv_rdp1', ierror)

         return
     end subroutine mp_gatherv_rdp1

!!
!! @sub mp_gatherv_rdp2
!!
!! gather real(dp) data from every processes to rank 0
!!
     subroutine mp_gatherv_rdp2(send, data, recv, disp, root, gid)
         implicit none

! external arguments
         real(dp), intent(in) :: send(:,:)
         real(dp), intent(inout) :: data(:,:)

         integer, intent(in) :: recv(:)
         integer, intent(in) :: disp(:)

         integer, intent(in) :: root
         integer, optional, intent(in) :: gid

! set current communicator
         if ( present(gid) .eqv. .true. ) then
             group = gid
         else
             group = MPI_COMM_WORLD
         endif ! back if ( present(gid) .eqv. .true. ) block

! barrier until all processes reach here
         call mp_barrier(group)

! setup element count
         isize = size(send)

! invoke related mpi subroutines
         call MPI_GATHERV(send, isize, m_rdp, data, recv, disp, m_rdp, root, group, ierror)

! handler for return code
         call mp_error('mp_gatherv_rdp2', ierror)

         return
     end subroutine mp_gatherv_rdp2

!!
!! @sub mp_gatherv_rdp3
!!
!! gather real(dp) data from every processes to rank 0
!!
     subroutine mp_gatherv_rdp3(send, data, recv, disp, root, gid)
         implicit none

! external arguments
         real(dp), intent(in) :: send(:,:,:)
         real(dp), intent(inout) :: data(:,:,:)

         integer, intent(in) :: recv(:)
         integer, intent(in) :: disp(:)

         integer, intent(in) :: root
         integer, optional, intent(in) :: gid

! set current communicator
         if ( present(gid) .eqv. .true. ) then
             group = gid
         else
             group = MPI_COMM_WORLD
         endif ! back if ( present(gid) .eqv. .true. ) block

! barrier until all processes reach here
         call mp_barrier(group)

! setup element count
         isize = size(send)

! invoke related mpi subroutines
         call MPI_GATHERV(send, isize, m_rdp, data, recv, disp, m_rdp, root, group, ierror)

! handler for return code
         call mp_error('mp_gatherv_rdp3', ierror)

         return
     end subroutine mp_gatherv_rdp3

!!
!! @sub mp_gatherv_rdp4
!!
!! gather real(dp) data from every processes to rank 0
!!
     subroutine mp_gatherv_rdp4(send, data, recv, disp, root, gid)
         implicit none

! external arguments
         real(dp), intent(in) :: send(:,:,:,:)
         real(dp), intent(inout) :: data(:,:,:,:)

         integer, intent(in) :: recv(:)
         integer, intent(in) :: disp(:)

         integer, intent(in) :: root
         integer, optional, intent(in) :: gid

! set current communicator
         if ( present(gid) .eqv. .true. ) then
             group = gid
         else
             group = MPI_COMM_WORLD
         endif ! back if ( present(gid) .eqv. .true. ) block

! barrier until all processes reach here
         call mp_barrier(group)

! setup element count
         isize = size(send)

! invoke related mpi subroutines
         call MPI_GATHERV(send, isize, m_rdp, data, recv, disp, m_rdp, root, group, ierror)

! handler for return code
         call mp_error('mp_gatherv_rdp4', ierror)

         return
     end subroutine mp_gatherv_rdp4

!!
!! @sub mp_gatherv_rdp5
!!
!! gather real(dp) data from every processes to rank 0
!!
     subroutine mp_gatherv_rdp5(send, data, recv, disp, root, gid)
         implicit none

! external arguments
         real(dp), intent(in) :: send(:,:,:,:,:)
         real(dp), intent(inout) :: data(:,:,:,:,:)

         integer, intent(in) :: recv(:)
         integer, intent(in) :: disp(:)

         integer, intent(in) :: root
         integer, optional, intent(in) :: gid

! set current communicator
         if ( present(gid) .eqv. .true. ) then
             group = gid
         else
             group = MPI_COMM_WORLD
         endif ! back if ( present(gid) .eqv. .true. ) block

! barrier until all processes reach here
         call mp_barrier(group)

! setup element count
         isize = size(send)

! invoke related mpi subroutines
         call MPI_GATHERV(send, isize, m_rdp, data, recv, disp, m_rdp, root, group, ierror)

! handler for return code
         call mp_error('mp_gatherv_rdp5', ierror)

         return
     end subroutine mp_gatherv_rdp5

!!
!! @sub mp_gatherv_cdp1
!!
!! gather complex(dp) data from every processes to rank 0
!!
     subroutine mp_gatherv_cdp1(send, data, recv, disp, root, gid)
         implicit none

! external arguments
         complex(dp), intent(in) :: send(:)
         complex(dp), intent(inout) :: data(:)

         integer, intent(in) :: recv(:)
         integer, intent(in) :: disp(:)

         integer, intent(in) :: root
         integer, optional, intent(in) :: gid

! set current communicator
         if ( present(gid) .eqv. .true. ) then
             group = gid
         else
             group = MPI_COMM_WORLD
         endif ! back if ( present(gid) .eqv. .true. ) block

! barrier until all processes reach here
         call mp_barrier(group)

! setup element count
         isize = size(send)

! invoke related mpi subroutines
         call MPI_GATHERV(send, isize, m_cdp, data, recv, disp, m_cdp, root, group, ierror)

! handler for return code
         call mp_error('mp_gatherv_cdp1', ierror)

         return
     end subroutine mp_gatherv_cdp1

!!
!! @sub mp_gatherv_cdp2
!!
!! gather complex(dp) data from every processes to rank 0
!!
     subroutine mp_gatherv_cdp2(send, data, recv, disp, root, gid)
         implicit none

! external arguments
         complex(dp), intent(in) :: send(:,:)
         complex(dp), intent(inout) :: data(:,:)

         integer, intent(in) :: recv(:)
         integer, intent(in) :: disp(:)

         integer, intent(in) :: root
         integer, optional, intent(in) :: gid

! set current communicator
         if ( present(gid) .eqv. .true. ) then
             group = gid
         else
             group = MPI_COMM_WORLD
         endif ! back if ( present(gid) .eqv. .true. ) block

! barrier until all processes reach here
         call mp_barrier(group)

! setup element count
         isize = size(send)

! invoke related mpi subroutines
         call MPI_GATHERV(send, isize, m_cdp, data, recv, disp, m_cdp, root, group, ierror)

! handler for return code
         call mp_error('mp_gatherv_cdp2', ierror)

         return
     end subroutine mp_gatherv_cdp2

!!
!! @sub mp_gatherv_cdp3
!!
!! gather complex(dp) data from every processes to rank 0
!!
     subroutine mp_gatherv_cdp3(send, data, recv, disp, root, gid)
         implicit none

! external arguments
         complex(dp), intent(in) :: send(:,:,:)
         complex(dp), intent(inout) :: data(:,:,:)

         integer, intent(in) :: recv(:)
         integer, intent(in) :: disp(:)

         integer, intent(in) :: root
         integer, optional, intent(in) :: gid

! set current communicator
         if ( present(gid) .eqv. .true. ) then
             group = gid
         else
             group = MPI_COMM_WORLD
         endif ! back if ( present(gid) .eqv. .true. ) block

! barrier until all processes reach here
         call mp_barrier(group)

! setup element count
         isize = size(send)

! invoke related mpi subroutines
         call MPI_GATHERV(send, isize, m_cdp, data, recv, disp, m_cdp, root, group, ierror)

! handler for return code
         call mp_error('mp_gatherv_cdp3', ierror)

         return
     end subroutine mp_gatherv_cdp3

!!
!! @sub mp_gatherv_cdp4
!!
!! gather complex(dp) data from every processes to rank 0
!!
     subroutine mp_gatherv_cdp4(send, data, recv, disp, root, gid)
         implicit none

! external arguments
         complex(dp), intent(in) :: send(:,:,:,:)
         complex(dp), intent(inout) :: data(:,:,:,:)

         integer, intent(in) :: recv(:)
         integer, intent(in) :: disp(:)

         integer, intent(in) :: root
         integer, optional, intent(in) :: gid

! set current communicator
         if ( present(gid) .eqv. .true. ) then
             group = gid
         else
             group = MPI_COMM_WORLD
         endif ! back if ( present(gid) .eqv. .true. ) block

! barrier until all processes reach here
         call mp_barrier(group)

! setup element count
         isize = size(send)

! invoke related mpi subroutines
         call MPI_GATHERV(send, isize, m_cdp, data, recv, disp, m_cdp, root, group, ierror)

! handler for return code
         call mp_error('mp_gatherv_cdp4', ierror)

         return
     end subroutine mp_gatherv_cdp4

!!
!! @sub mp_gatherv_cdp5
!!
!! gather complex(dp) data from every processes to rank 0
!!
     subroutine mp_gatherv_cdp5(send, data, recv, disp, root, gid)
         implicit none

! external arguments
         complex(dp), intent(in) :: send(:,:,:,:,:)
         complex(dp), intent(inout) :: data(:,:,:,:,:)

         integer, intent(in) :: recv(:)
         integer, intent(in) :: disp(:)

         integer, intent(in) :: root
         integer, optional, intent(in) :: gid

! set current communicator
         if ( present(gid) .eqv. .true. ) then
             group = gid
         else
             group = MPI_COMM_WORLD
         endif ! back if ( present(gid) .eqv. .true. ) block

! barrier until all processes reach here
         call mp_barrier(group)

! setup element count
         isize = size(send)

! invoke related mpi subroutines
         call MPI_GATHERV(send, isize, m_cdp, data, recv, disp, m_cdp, root, group, ierror)

! handler for return code
         call mp_error('mp_gatherv_cdp5', ierror)

         return
     end subroutine mp_gatherv_cdp5

!!========================================================================
!!>>> MPI collective operations: allgathering                          <<<
!!========================================================================

!!
!! @sub mp_allgather_int1
!!
!! gather integer data from all processes and then redistribute it to
!! all processes
!!
     subroutine mp_allgather_int1(send, data, gid)
         implicit none

! external arguments
         integer, intent(in) :: send(:)
         integer, intent(inout) :: data(:)

         integer, optional, intent(in) :: gid

! set current communicator
         if ( present(gid) .eqv. .true. ) then
             group = gid
         else
             group = MPI_COMM_WORLD
         endif ! back if ( present(gid) .eqv. .true. ) block

! barrier until all processes reach here
         call mp_barrier(group)

! setup element count
         isize = size(send)

! invoke related mpi subroutines
         call MPI_ALLGATHER(send, isize, m_int, data, isize, m_int, group, ierror)

! handler for return code
         call mp_error('mp_allgather_int1', ierror)

         return
     end subroutine mp_allgather_int1

!!
!! @sub mp_allgather_int2
!!
!! gather integer data from all processes and then redistribute it to
!! all processes
!!
     subroutine mp_allgather_int2(send, data, gid)
         implicit none

! external arguments
         integer, intent(in) :: send(:,:)
         integer, intent(inout) :: data(:,:)

         integer, optional, intent(in) :: gid

! set current communicator
         if ( present(gid) .eqv. .true. ) then
             group = gid
         else
             group = MPI_COMM_WORLD
         endif ! back if ( present(gid) .eqv. .true. ) block

! barrier until all processes reach here
         call mp_barrier(group)

! setup element count
         isize = size(send)

! invoke related mpi subroutines
         call MPI_ALLGATHER(send, isize, m_int, data, isize, m_int, group, ierror)

! handler for return code
         call mp_error('mp_allgather_int2', ierror)

         return
     end subroutine mp_allgather_int2

!!
!! @sub mp_allgather_int3
!!
!! gather integer data from all processes and then redistribute it to
!! all processes
!!
     subroutine mp_allgather_int3(send, data, gid)
         implicit none

! external arguments
         integer, intent(in) :: send(:,:,:)
         integer, intent(inout) :: data(:,:,:)

         integer, optional, intent(in) :: gid

! set current communicator
         if ( present(gid) .eqv. .true. ) then
             group = gid
         else
             group = MPI_COMM_WORLD
         endif ! back if ( present(gid) .eqv. .true. ) block

! barrier until all processes reach here
         call mp_barrier(group)

! setup element count
         isize = size(send)

! invoke related mpi subroutines
         call MPI_ALLGATHER(send, isize, m_int, data, isize, m_int, group, ierror)

! handler for return code
         call mp_error('mp_allgather_int3', ierror)

         return
     end subroutine mp_allgather_int3

!!
!! @sub mp_allgather_int4
!!
!! gather integer data from all processes and then redistribute it to
!! all processes
!!
     subroutine mp_allgather_int4(send, data, gid)
         implicit none

! external arguments
         integer, intent(in) :: send(:,:,:,:)
         integer, intent(inout) :: data(:,:,:,:)

         integer, optional, intent(in) :: gid

! set current communicator
         if ( present(gid) .eqv. .true. ) then
             group = gid
         else
             group = MPI_COMM_WORLD
         endif ! back if ( present(gid) .eqv. .true. ) block

! barrier until all processes reach here
         call mp_barrier(group)

! setup element count
         isize = size(send)

! invoke related mpi subroutines
         call MPI_ALLGATHER(send, isize, m_int, data, isize, m_int, group, ierror)

! handler for return code
         call mp_error('mp_allgather_int4', ierror)

         return
     end subroutine mp_allgather_int4

!!
!! @sub mp_allgather_int5
!!
!! gather integer data from all processes and then redistribute it to
!! all processes
!!
     subroutine mp_allgather_int5(send, data, gid)
         implicit none

! external arguments
         integer, intent(in) :: send(:,:,:,:,:)
         integer, intent(inout) :: data(:,:,:,:,:)

         integer, optional, intent(in) :: gid

! set current communicator
         if ( present(gid) .eqv. .true. ) then
             group = gid
         else
             group = MPI_COMM_WORLD
         endif ! back if ( present(gid) .eqv. .true. ) block

! barrier until all processes reach here
         call mp_barrier(group)

! setup element count
         isize = size(send)

! invoke related mpi subroutines
         call MPI_ALLGATHER(send, isize, m_int, data, isize, m_int, group, ierror)

! handler for return code
         call mp_error('mp_allgather_int5', ierror)

         return
     end subroutine mp_allgather_int5

!!
!! @sub mp_allgather_rdp1
!!
!! gather real(dp) data from all processes and then redistribute it to
!! all processes
!!
     subroutine mp_allgather_rdp1(send, data, gid)
         implicit none

! external arguments
         real(dp), intent(in) :: send(:)
         real(dp), intent(inout) :: data(:)

         integer, optional, intent(in) :: gid

! set current communicator
         if ( present(gid) .eqv. .true. ) then
             group = gid
         else
             group = MPI_COMM_WORLD
         endif ! back if ( present(gid) .eqv. .true. ) block

! barrier until all processes reach here
         call mp_barrier(group)

! setup element count
         isize = size(send)

! invoke related mpi subroutines
         call MPI_ALLGATHER(send, isize, m_rdp, data, isize, m_rdp, group, ierror)

! handler for return code
         call mp_error('mp_allgather_rdp1', ierror)

         return
     end subroutine mp_allgather_rdp1

!!
!! @sub mp_allgather_rdp2
!!
!! gather real(dp) data from all processes and then redistribute it to
!! all processes
!!
     subroutine mp_allgather_rdp2(send, data, gid)
         implicit none

! external arguments
         real(dp), intent(in) :: send(:,:)
         real(dp), intent(inout) :: data(:,:)

         integer, optional, intent(in) :: gid

! set current communicator
         if ( present(gid) .eqv. .true. ) then
             group = gid
         else
             group = MPI_COMM_WORLD
         endif ! back if ( present(gid) .eqv. .true. ) block

! barrier until all processes reach here
         call mp_barrier(group)

! setup element count
         isize = size(send)

! invoke related mpi subroutines
         call MPI_ALLGATHER(send, isize, m_rdp, data, isize, m_rdp, group, ierror)

! handler for return code
         call mp_error('mp_allgather_rdp2', ierror)

         return
     end subroutine mp_allgather_rdp2

!!
!! @sub mp_allgather_rdp3
!!
!! gather real(dp) data from all processes and then redistribute it to
!! all processes
!!
     subroutine mp_allgather_rdp3(send, data, gid)
         implicit none

! external arguments
         real(dp), intent(in) :: send(:,:,:)
         real(dp), intent(inout) :: data(:,:,:)

         integer, optional, intent(in) :: gid

! set current communicator
         if ( present(gid) .eqv. .true. ) then
             group = gid
         else
             group = MPI_COMM_WORLD
         endif ! back if ( present(gid) .eqv. .true. ) block

! barrier until all processes reach here
         call mp_barrier(group)

! setup element count
         isize = size(send)

! invoke related mpi subroutines
         call MPI_ALLGATHER(send, isize, m_rdp, data, isize, m_rdp, group, ierror)

! handler for return code
         call mp_error('mp_allgather_rdp3', ierror)

         return
     end subroutine mp_allgather_rdp3

!!
!! @sub mp_allgather_rdp4
!!
!! gather real(dp) data from all processes and then redistribute it to
!! all processes
!!
     subroutine mp_allgather_rdp4(send, data, gid)
         implicit none

! external arguments
         real(dp), intent(in) :: send(:,:,:,:)
         real(dp), intent(inout) :: data(:,:,:,:)

         integer, optional, intent(in) :: gid

! set current communicator
         if ( present(gid) .eqv. .true. ) then
             group = gid
         else
             group = MPI_COMM_WORLD
         endif ! back if ( present(gid) .eqv. .true. ) block

! barrier until all processes reach here
         call mp_barrier(group)

! setup element count
         isize = size(send)

! invoke related mpi subroutines
         call MPI_ALLGATHER(send, isize, m_rdp, data, isize, m_rdp, group, ierror)

! handler for return code
         call mp_error('mp_allgather_rdp4', ierror)

         return
     end subroutine mp_allgather_rdp4

!!
!! @sub mp_allgather_rdp5
!!
!! gather real(dp) data from all processes and then redistribute it to
!! all processes
!!
     subroutine mp_allgather_rdp5(send, data, gid)
         implicit none

! external arguments
         real(dp), intent(in) :: send(:,:,:,:,:)
         real(dp), intent(inout) :: data(:,:,:,:,:)

         integer, optional, intent(in) :: gid

! set current communicator
         if ( present(gid) .eqv. .true. ) then
             group = gid
         else
             group = MPI_COMM_WORLD
         endif ! back if ( present(gid) .eqv. .true. ) block

! barrier until all processes reach here
         call mp_barrier(group)

! setup element count
         isize = size(send)

! invoke related mpi subroutines
         call MPI_ALLGATHER(send, isize, m_rdp, data, isize, m_rdp, group, ierror)

! handler for return code
         call mp_error('mp_allgather_rdp5', ierror)

         return
     end subroutine mp_allgather_rdp5

!!
!! @sub mp_allgather_cdp1
!!
!! gather complex(dp) data from all processes and then redistribute it to
!! all processes
!!
     subroutine mp_allgather_cdp1(send, data, gid)
         implicit none

! external arguments
         complex(dp), intent(in) :: send(:)
         complex(dp), intent(inout) :: data(:)

         integer, optional, intent(in) :: gid

! set current communicator
         if ( present(gid) .eqv. .true. ) then
             group = gid
         else
             group = MPI_COMM_WORLD
         endif ! back if ( present(gid) .eqv. .true. ) block

! barrier until all processes reach here
         call mp_barrier(group)

! setup element count
         isize = size(send)

! invoke related mpi subroutines
         call MPI_ALLGATHER(send, isize, m_cdp, data, isize, m_cdp, group, ierror)

! handler for return code
         call mp_error('mp_allgather_cdp1', ierror)

         return
     end subroutine mp_allgather_cdp1

!!
!! @sub mp_allgather_cdp2
!!
!! gather complex(dp) data from all processes and then redistribute it to
!! all processes
!!
     subroutine mp_allgather_cdp2(send, data, gid)
         implicit none

! external arguments
         complex(dp), intent(in) :: send(:,:)
         complex(dp), intent(inout) :: data(:,:)

         integer, optional, intent(in) :: gid

! set current communicator
         if ( present(gid) .eqv. .true. ) then
             group = gid
         else
             group = MPI_COMM_WORLD
         endif ! back if ( present(gid) .eqv. .true. ) block

! barrier until all processes reach here
         call mp_barrier(group)

! setup element count
         isize = size(send)

! invoke related mpi subroutines
         call MPI_ALLGATHER(send, isize, m_cdp, data, isize, m_cdp, group, ierror)

! handler for return code
         call mp_error('mp_allgather_cdp2', ierror)

         return
     end subroutine mp_allgather_cdp2

!!
!! @sub mp_allgather_cdp3
!!
!! gather complex(dp) data from all processes and then redistribute it to
!! all processes
!!
     subroutine mp_allgather_cdp3(send, data, gid)
         implicit none

! external arguments
         complex(dp), intent(in) :: send(:,:,:)
         complex(dp), intent(inout) :: data(:,:,:)

         integer, optional, intent(in) :: gid

! set current communicator
         if ( present(gid) .eqv. .true. ) then
             group = gid
         else
             group = MPI_COMM_WORLD
         endif ! back if ( present(gid) .eqv. .true. ) block

! barrier until all processes reach here
         call mp_barrier(group)

! setup element count
         isize = size(send)

! invoke related mpi subroutines
         call MPI_ALLGATHER(send, isize, m_cdp, data, isize, m_cdp, group, ierror)

! handler for return code
         call mp_error('mp_allgather_cdp3', ierror)

         return
     end subroutine mp_allgather_cdp3

!!
!! @sub mp_allgather_cdp4
!!
!! gather complex(dp) data from all processes and then redistribute it to
!! all processes
!!
     subroutine mp_allgather_cdp4(send, data, gid)
         implicit none

! external arguments
         complex(dp), intent(in) :: send(:,:,:,:)
         complex(dp), intent(inout) :: data(:,:,:,:)

         integer, optional, intent(in) :: gid

! set current communicator
         if ( present(gid) .eqv. .true. ) then
             group = gid
         else
             group = MPI_COMM_WORLD
         endif ! back if ( present(gid) .eqv. .true. ) block

! barrier until all processes reach here
         call mp_barrier(group)

! setup element count
         isize = size(send)

! invoke related mpi subroutines
         call MPI_ALLGATHER(send, isize, m_cdp, data, isize, m_cdp, group, ierror)

! handler for return code
         call mp_error('mp_allgather_cdp4', ierror)

         return
     end subroutine mp_allgather_cdp4

!!
!! @sub mp_allgather_cdp5
!!
!! gather complex(dp) data from all processes and then redistribute it to
!! all processes
!!
     subroutine mp_allgather_cdp5(send, data, gid)
         implicit none

! external arguments
         complex(dp), intent(in) :: send(:,:,:,:,:)
         complex(dp), intent(inout) :: data(:,:,:,:,:)

         integer, optional, intent(in) :: gid

! set current communicator
         if ( present(gid) .eqv. .true. ) then
             group = gid
         else
             group = MPI_COMM_WORLD
         endif ! back if ( present(gid) .eqv. .true. ) block

! barrier until all processes reach here
         call mp_barrier(group)

! setup element count
         isize = size(send)

! invoke related mpi subroutines
         call MPI_ALLGATHER(send, isize, m_cdp, data, isize, m_cdp, group, ierror)

! handler for return code
         call mp_error('mp_allgather_cdp5', ierror)

         return
     end subroutine mp_allgather_cdp5

!!========================================================================
!!>>> MPI collective operations: allgatherving                         <<<
!!========================================================================

!!
!! @sub mp_allgatherv_int1
!!
!! gather integer data from all processes and then redistribute it to
!! all processes
!!
     subroutine mp_allgatherv_int1(send, data, recv, disp, gid)
         implicit none

! external arguments
         integer, intent(in) :: send(:)
         integer, intent(inout) :: data(:)

         integer, intent(in) :: recv(:)
         integer, intent(in) :: disp(:)

         integer, optional, intent(in) :: gid

! set current communicator
         if ( present(gid) .eqv. .true. ) then
             group = gid
         else
             group = MPI_COMM_WORLD
         endif ! back if ( present(gid) .eqv. .true. ) block

! barrier until all processes reach here
         call mp_barrier(group)

! setup element count
         isize = size(send)

! invoke related mpi subroutines
         call MPI_ALLGATHERV(send, isize, m_int, data, recv, disp, m_int, group, ierror)

! handler for return code
         call mp_error('mp_allgatherv_int1', ierror)

         return
     end subroutine mp_allgatherv_int1

!!
!! @sub mp_allgatherv_int2
!!
!! gather integer data from all processes and then redistribute it to
!! all processes
!!
     subroutine mp_allgatherv_int2(send, data, recv, disp, gid)
         implicit none

! external arguments
         integer, intent(in) :: send(:,:)
         integer, intent(inout) :: data(:,:)

         integer, intent(in) :: recv(:)
         integer, intent(in) :: disp(:)

         integer, optional, intent(in) :: gid

! set current communicator
         if ( present(gid) .eqv. .true. ) then
             group = gid
         else
             group = MPI_COMM_WORLD
         endif ! back if ( present(gid) .eqv. .true. ) block

! barrier until all processes reach here
         call mp_barrier(group)

! setup element count
         isize = size(send)

! invoke related mpi subroutines
         call MPI_ALLGATHERV(send, isize, m_int, data, recv, disp, m_int, group, ierror)

! handler for return code
         call mp_error('mp_allgatherv_int2', ierror)

         return
     end subroutine mp_allgatherv_int2

!!
!! @sub mp_allgatherv_int3
!!
!! gather integer data from all processes and then redistribute it to
!! all processes
!!
     subroutine mp_allgatherv_int3(send, data, recv, disp, gid)
         implicit none

! external arguments
         integer, intent(in) :: send(:,:,:)
         integer, intent(inout) :: data(:,:,:)

         integer, intent(in) :: recv(:)
         integer, intent(in) :: disp(:)

         integer, optional, intent(in) :: gid

! set current communicator
         if ( present(gid) .eqv. .true. ) then
             group = gid
         else
             group = MPI_COMM_WORLD
         endif ! back if ( present(gid) .eqv. .true. ) block

! barrier until all processes reach here
         call mp_barrier(group)

! setup element count
         isize = size(send)

! invoke related mpi subroutines
         call MPI_ALLGATHERV(send, isize, m_int, data, recv, disp, m_int, group, ierror)

! handler for return code
         call mp_error('mp_allgatherv_int3', ierror)

         return
     end subroutine mp_allgatherv_int3

!!
!! @sub mp_allgatherv_int4
!!
!! gather integer data from all processes and then redistribute it to
!! all processes
!!
     subroutine mp_allgatherv_int4(send, data, recv, disp, gid)
         implicit none

! external arguments
         integer, intent(in) :: send(:,:,:,:)
         integer, intent(inout) :: data(:,:,:,:)

         integer, intent(in) :: recv(:)
         integer, intent(in) :: disp(:)

         integer, optional, intent(in) :: gid

! set current communicator
         if ( present(gid) .eqv. .true. ) then
             group = gid
         else
             group = MPI_COMM_WORLD
         endif ! back if ( present(gid) .eqv. .true. ) block

! barrier until all processes reach here
         call mp_barrier(group)

! setup element count
         isize = size(send)

! invoke related mpi subroutines
         call MPI_ALLGATHERV(send, isize, m_int, data, recv, disp, m_int, group, ierror)

! handler for return code
         call mp_error('mp_allgatherv_int4', ierror)

         return
     end subroutine mp_allgatherv_int4

!!
!! @sub mp_allgatherv_int5
!!
!! gather integer data from all processes and then redistribute it to
!! all processes
!!
     subroutine mp_allgatherv_int5(send, data, recv, disp, gid)
         implicit none

! external arguments
         integer, intent(in) :: send(:,:,:,:,:)
         integer, intent(inout) :: data(:,:,:,:,:)

         integer, intent(in) :: recv(:)
         integer, intent(in) :: disp(:)

         integer, optional, intent(in) :: gid

! set current communicator
         if ( present(gid) .eqv. .true. ) then
             group = gid
         else
             group = MPI_COMM_WORLD
         endif ! back if ( present(gid) .eqv. .true. ) block

! barrier until all processes reach here
         call mp_barrier(group)

! setup element count
         isize = size(send)

! invoke related mpi subroutines
         call MPI_ALLGATHERV(send, isize, m_int, data, recv, disp, m_int, group, ierror)

! handler for return code
         call mp_error('mp_allgatherv_int5', ierror)

         return
     end subroutine mp_allgatherv_int5

!!
!! @sub mp_allgatherv_rdp1
!!
!! gather real(dp) data from all processes and then redistribute it to
!! all processes
!!
     subroutine mp_allgatherv_rdp1(send, data, recv, disp, gid)
         implicit none

! external arguments
         real(dp), intent(in) :: send(:)
         real(dp), intent(inout) :: data(:)

         integer, intent(in) :: recv(:)
         integer, intent(in) :: disp(:)

         integer, optional, intent(in) :: gid

! set current communicator
         if ( present(gid) .eqv. .true. ) then
             group = gid
         else
             group = MPI_COMM_WORLD
         endif ! back if ( present(gid) .eqv. .true. ) block

! barrier until all processes reach here
         call mp_barrier(group)

! setup element count
         isize = size(send)

! invoke related mpi subroutines
         call MPI_ALLGATHERV(send, isize, m_rdp, data, recv, disp, m_rdp, group, ierror)

! handler for return code
         call mp_error('mp_allgatherv_rdp1', ierror)

         return
     end subroutine mp_allgatherv_rdp1

!!
!! @sub mp_allgatherv_rdp2
!!
!! gather real(dp) data from all processes and then redistribute it to
!! all processes
!!
     subroutine mp_allgatherv_rdp2(send, data, recv, disp, gid)
         implicit none

! external arguments
         real(dp), intent(in) :: send(:,:)
         real(dp), intent(inout) :: data(:,:)

         integer, intent(in) :: recv(:)
         integer, intent(in) :: disp(:)

         integer, optional, intent(in) :: gid

! set current communicator
         if ( present(gid) .eqv. .true. ) then
             group = gid
         else
             group = MPI_COMM_WORLD
         endif ! back if ( present(gid) .eqv. .true. ) block

! barrier until all processes reach here
         call mp_barrier(group)

! setup element count
         isize = size(send)

! invoke related mpi subroutines
         call MPI_ALLGATHERV(send, isize, m_rdp, data, recv, disp, m_rdp, group, ierror)

! handler for return code
         call mp_error('mp_allgatherv_rdp2', ierror)

         return
     end subroutine mp_allgatherv_rdp2

!!
!! @sub mp_allgatherv_rdp3
!!
!! gather real(dp) data from all processes and then redistribute it to
!! all processes
!!
     subroutine mp_allgatherv_rdp3(send, data, recv, disp, gid)
         implicit none

! external arguments
         real(dp), intent(in) :: send(:,:,:)
         real(dp), intent(inout) :: data(:,:,:)

         integer, intent(in) :: recv(:)
         integer, intent(in) :: disp(:)

         integer, optional, intent(in) :: gid

! set current communicator
         if ( present(gid) .eqv. .true. ) then
             group = gid
         else
             group = MPI_COMM_WORLD
         endif ! back if ( present(gid) .eqv. .true. ) block

! barrier until all processes reach here
         call mp_barrier(group)

! setup element count
         isize = size(send)

! invoke related mpi subroutines
         call MPI_ALLGATHERV(send, isize, m_rdp, data, recv, disp, m_rdp, group, ierror)

! handler for return code
         call mp_error('mp_allgatherv_rdp3', ierror)

         return
     end subroutine mp_allgatherv_rdp3

!!
!! @sub mp_allgatherv_rdp4
!!
!! gather real(dp) data from all processes and then redistribute it to
!! all processes
!!
     subroutine mp_allgatherv_rdp4(send, data, recv, disp, gid)
         implicit none

! external arguments
         real(dp), intent(in) :: send(:,:,:,:)
         real(dp), intent(inout) :: data(:,:,:,:)

         integer, intent(in) :: recv(:)
         integer, intent(in) :: disp(:)

         integer, optional, intent(in) :: gid

! set current communicator
         if ( present(gid) .eqv. .true. ) then
             group = gid
         else
             group = MPI_COMM_WORLD
         endif ! back if ( present(gid) .eqv. .true. ) block

! barrier until all processes reach here
         call mp_barrier(group)

! setup element count
         isize = size(send)

! invoke related mpi subroutines
         call MPI_ALLGATHERV(send, isize, m_rdp, data, recv, disp, m_rdp, group, ierror)

! handler for return code
         call mp_error('mp_allgatherv_rdp4', ierror)

         return
     end subroutine mp_allgatherv_rdp4

!!
!! @sub mp_allgatherv_rdp5
!!
!! gather real(dp) data from all processes and then redistribute it to
!! all processes
!!
     subroutine mp_allgatherv_rdp5(send, data, recv, disp, gid)
         implicit none

! external arguments
         real(dp), intent(in) :: send(:,:,:,:,:)
         real(dp), intent(inout) :: data(:,:,:,:,:)

         integer, intent(in) :: recv(:)
         integer, intent(in) :: disp(:)

         integer, optional, intent(in) :: gid

! set current communicator
         if ( present(gid) .eqv. .true. ) then
             group = gid
         else
             group = MPI_COMM_WORLD
         endif ! back if ( present(gid) .eqv. .true. ) block

! barrier until all processes reach here
         call mp_barrier(group)

! setup element count
         isize = size(send)

! invoke related mpi subroutines
         call MPI_ALLGATHERV(send, isize, m_rdp, data, recv, disp, m_rdp, group, ierror)

! handler for return code
         call mp_error('mp_allgatherv_rdp5', ierror)

         return
     end subroutine mp_allgatherv_rdp5

!!
!! @sub mp_allgatherv_cdp1
!!
!! gather complex(dp) data from all processes and then redistribute it to
!! all processes
!!
     subroutine mp_allgatherv_cdp1(send, data, recv, disp, gid)
         implicit none

! external arguments
         complex(dp), intent(in) :: send(:)
         complex(dp), intent(inout) :: data(:)

         integer, intent(in) :: recv(:)
         integer, intent(in) :: disp(:)

         integer, optional, intent(in) :: gid

! set current communicator
         if ( present(gid) .eqv. .true. ) then
             group = gid
         else
             group = MPI_COMM_WORLD
         endif ! back if ( present(gid) .eqv. .true. ) block

! barrier until all processes reach here
         call mp_barrier(group)

! setup element count
         isize = size(send)

! invoke related mpi subroutines
         call MPI_ALLGATHERV(send, isize, m_cdp, data, recv, disp, m_cdp, group, ierror)

! handler for return code
         call mp_error('mp_allgatherv_cdp1', ierror)

         return
     end subroutine mp_allgatherv_cdp1

!!
!! @sub mp_allgatherv_cdp2
!!
!! gather complex(dp) data from all processes and then redistribute it to
!! all processes
!!
     subroutine mp_allgatherv_cdp2(send, data, recv, disp, gid)
         implicit none

! external arguments
         complex(dp), intent(in) :: send(:,:)
         complex(dp), intent(inout) :: data(:,:)

         integer, intent(in) :: recv(:)
         integer, intent(in) :: disp(:)

         integer, optional, intent(in) :: gid

! set current communicator
         if ( present(gid) .eqv. .true. ) then
             group = gid
         else
             group = MPI_COMM_WORLD
         endif ! back if ( present(gid) .eqv. .true. ) block

! barrier until all processes reach here
         call mp_barrier(group)

! setup element count
         isize = size(send)

! invoke related mpi subroutines
         call MPI_ALLGATHERV(send, isize, m_cdp, data, recv, disp, m_cdp, group, ierror)

! handler for return code
         call mp_error('mp_allgatherv_cdp2', ierror)

         return
     end subroutine mp_allgatherv_cdp2

!!
!! @sub mp_allgatherv_cdp3
!!
!! gather complex(dp) data from all processes and then redistribute it to
!! all processes
!!
     subroutine mp_allgatherv_cdp3(send, data, recv, disp, gid)
         implicit none

! external arguments
         complex(dp), intent(in) :: send(:,:,:)
         complex(dp), intent(inout) :: data(:,:,:)

         integer, intent(in) :: recv(:)
         integer, intent(in) :: disp(:)

         integer, optional, intent(in) :: gid

! set current communicator
         if ( present(gid) .eqv. .true. ) then
             group = gid
         else
             group = MPI_COMM_WORLD
         endif ! back if ( present(gid) .eqv. .true. ) block

! barrier until all processes reach here
         call mp_barrier(group)

! setup element count
         isize = size(send)

! invoke related mpi subroutines
         call MPI_ALLGATHERV(send, isize, m_cdp, data, recv, disp, m_cdp, group, ierror)

! handler for return code
         call mp_error('mp_allgatherv_cdp3', ierror)

         return
     end subroutine mp_allgatherv_cdp3

!!
!! @sub mp_allgatherv_cdp4
!!
!! gather complex(dp) data from all processes and then redistribute it to
!! all processes
!!
     subroutine mp_allgatherv_cdp4(send, data, recv, disp, gid)
         implicit none

! external arguments
         complex(dp), intent(in) :: send(:,:,:,:)
         complex(dp), intent(inout) :: data(:,:,:,:)

         integer, intent(in) :: recv(:)
         integer, intent(in) :: disp(:)

         integer, optional, intent(in) :: gid

! set current communicator
         if ( present(gid) .eqv. .true. ) then
             group = gid
         else
             group = MPI_COMM_WORLD
         endif ! back if ( present(gid) .eqv. .true. ) block

! barrier until all processes reach here
         call mp_barrier(group)

! setup element count
         isize = size(send)

! invoke related mpi subroutines
         call MPI_ALLGATHERV(send, isize, m_cdp, data, recv, disp, m_cdp, group, ierror)

! handler for return code
         call mp_error('mp_allgatherv_cdp4', ierror)

         return
     end subroutine mp_allgatherv_cdp4

!!
!! @sub mp_allgatherv_cdp5
!!
!! gather complex(dp) data from all processes and then redistribute it to
!! all processes
!!
     subroutine mp_allgatherv_cdp5(send, data, recv, disp, gid)
         implicit none

! external arguments
         complex(dp), intent(in) :: send(:,:,:,:,:)
         complex(dp), intent(inout) :: data(:,:,:,:,:)

         integer, intent(in) :: recv(:)
         integer, intent(in) :: disp(:)

         integer, optional, intent(in) :: gid

! set current communicator
         if ( present(gid) .eqv. .true. ) then
             group = gid
         else
             group = MPI_COMM_WORLD
         endif ! back if ( present(gid) .eqv. .true. ) block

! barrier until all processes reach here
         call mp_barrier(group)

! setup element count
         isize = size(send)

! invoke related mpi subroutines
         call MPI_ALLGATHERV(send, isize, m_cdp, data, recv, disp, m_cdp, group, ierror)

! handler for return code
         call mp_error('mp_allgatherv_cdp5', ierror)

         return
     end subroutine mp_allgatherv_cdp5

!!========================================================================
!!>>> MPI collective operations: reducing                              <<<
!!========================================================================

!!
!! @sub mp_reduce_int0
!!
!! reduce 1 integer from all processes
!!
     subroutine mp_reduce_int0(source, data, root, mop, gid)
         implicit none

! external arguments
         integer, intent(in) :: source
         integer, intent(inout) :: data
         integer, intent(in) :: root
         integer, optional, intent(in) :: mop
         integer, optional, intent(in) :: gid

! set current operator
         if ( present(mop) .eqv. .true. ) then
             opera = mop
         else
             opera = MPI_SUM
         endif ! back if ( present(mop) .eqv. .true. ) block

! set current communicator
         if ( present(gid) .eqv. .true. ) then
             group = gid
         else
             group = MPI_COMM_WORLD
         endif ! back if ( present(gid) .eqv. .true. ) block

! barrier until all processes reach here
         call mp_barrier(group)

! invoke related mpi subroutines
         call MPI_REDUCE(source, data, 1, m_int, opera, root, group, ierror)

! handler for return code
         call mp_error('mp_reduce_int0', ierror)

         return
     end subroutine mp_reduce_int0

!!
!! @sub mp_reduce_int1
!!
!! reduce integer vector from all processes
!!
     subroutine mp_reduce_int1(source, data, root, mop, gid)
         implicit none

! external arguments
         integer, intent(in) :: source(:)
         integer, intent(inout) :: data(:)
         integer, intent(in) :: root
         integer, optional, intent(in) :: mop
         integer, optional, intent(in) :: gid

! set current operator
         if ( present(mop) .eqv. .true. ) then
             opera = mop
         else
             opera = MPI_SUM
         endif ! back if ( present(mop) .eqv. .true. ) block

! set current communicator
         if ( present(gid) .eqv. .true. ) then
             group = gid
         else
             group = MPI_COMM_WORLD
         endif ! back if ( present(gid) .eqv. .true. ) block

! barrier until all processes reach here
         call mp_barrier(group)

! setup element count
         isize = size(source)

! invoke related mpi subroutines
         call MPI_REDUCE(source, data, isize, m_int, opera, root, group, ierror)

! handler for return code
         call mp_error('mp_reduce_int1', ierror)

         return
     end subroutine mp_reduce_int1

!!
!! @sub mp_reduce_int2
!!
!! reduce integer matrix from all processes
!!
     subroutine mp_reduce_int2(source, data, root, mop, gid)
         implicit none

! external arguments
         integer, intent(in) :: source(:,:)
         integer, intent(inout) :: data(:,:)
         integer, intent(in) :: root
         integer, optional, intent(in) :: mop
         integer, optional, intent(in) :: gid

! set current operator
         if ( present(mop) .eqv. .true. ) then
             opera = mop
         else
             opera = MPI_SUM
         endif ! back if ( present(mop) .eqv. .true. ) block

! set current communicator
         if ( present(gid) .eqv. .true. ) then
             group = gid
         else
             group = MPI_COMM_WORLD
         endif ! back if ( present(gid) .eqv. .true. ) block

! barrier until all processes reach here
         call mp_barrier(group)

! setup element count
         isize = size(source)

! invoke related mpi subroutines
         call MPI_REDUCE(source, data, isize, m_int, opera, root, group, ierror)

! handler for return code
         call mp_error('mp_reduce_int2', ierror)

         return
     end subroutine mp_reduce_int2

!!
!! @sub mp_reduce_int3
!!
!! reduce integer matrix from all processes
!!
     subroutine mp_reduce_int3(source, data, root, mop, gid)
         implicit none

! external arguments
         integer, intent(in) :: source(:,:,:)
         integer, intent(inout) :: data(:,:,:)
         integer, intent(in) :: root
         integer, optional, intent(in) :: mop
         integer, optional, intent(in) :: gid

! set current operator
         if ( present(mop) .eqv. .true. ) then
             opera = mop
         else
             opera = MPI_SUM
         endif ! back if ( present(mop) .eqv. .true. ) block

! set current communicator
         if ( present(gid) .eqv. .true. ) then
             group = gid
         else
             group = MPI_COMM_WORLD
         endif ! back if ( present(gid) .eqv. .true. ) block

! barrier until all processes reach here
         call mp_barrier(group)

! setup element count
         isize = size(source)

! invoke related mpi subroutines
         call MPI_REDUCE(source, data, isize, m_int, opera, root, group, ierror)

! handler for return code
         call mp_error('mp_reduce_int3', ierror)

         return
     end subroutine mp_reduce_int3

!!
!! @sub mp_reduce_int4
!!
!! reduce integer matrix from all processes
!!
     subroutine mp_reduce_int4(source, data, root, mop, gid)
         implicit none

! external arguments
         integer, intent(in) :: source(:,:,:,:)
         integer, intent(inout) :: data(:,:,:,:)
         integer, intent(in) :: root
         integer, optional, intent(in) :: mop
         integer, optional, intent(in) :: gid

! set current operator
         if ( present(mop) .eqv. .true. ) then
             opera = mop
         else
             opera = MPI_SUM
         endif ! back if ( present(mop) .eqv. .true. ) block

! set current communicator
         if ( present(gid) .eqv. .true. ) then
             group = gid
         else
             group = MPI_COMM_WORLD
         endif ! back if ( present(gid) .eqv. .true. ) block

! barrier until all processes reach here
         call mp_barrier(group)

! setup element count
         isize = size(source)

! invoke related mpi subroutines
         call MPI_REDUCE(source, data, isize, m_int, opera, root, group, ierror)

! handler for return code
         call mp_error('mp_reduce_int4', ierror)

         return
     end subroutine mp_reduce_int4

!!
!! @sub mp_reduce_int5
!!
!! reduce integer matrix from all processes
!!
     subroutine mp_reduce_int5(source, data, root, mop, gid)
         implicit none

! external arguments
         integer, intent(in) :: source(:,:,:,:,:)
         integer, intent(inout) :: data(:,:,:,:,:)
         integer, intent(in) :: root
         integer, optional, intent(in) :: mop
         integer, optional, intent(in) :: gid

! set current operator
         if ( present(mop) .eqv. .true. ) then
             opera = mop
         else
             opera = MPI_SUM
         endif ! back if ( present(mop) .eqv. .true. ) block

! set current communicator
         if ( present(gid) .eqv. .true. ) then
             group = gid
         else
             group = MPI_COMM_WORLD
         endif ! back if ( present(gid) .eqv. .true. ) block

! barrier until all processes reach here
         call mp_barrier(group)

! setup element count
         isize = size(source)

! invoke related mpi subroutines
         call MPI_REDUCE(source, data, isize, m_int, opera, root, group, ierror)

! handler for return code
         call mp_error('mp_reduce_int5', ierror)

         return
     end subroutine mp_reduce_int5

!!
!! @sub mp_reduce_rdp0
!!
!! reduce 1 real(dp) from all processes
!!
     subroutine mp_reduce_rdp0(source, data, root, mop, gid)
         implicit none

! external arguments
         real(dp), intent(in) :: source
         real(dp), intent(inout) :: data
         integer, intent(in) :: root
         integer, optional, intent(in) :: mop
         integer, optional, intent(in) :: gid

! set current operator
         if ( present(mop) .eqv. .true. ) then
             opera = mop
         else
             opera = MPI_SUM
         endif ! back if ( present(mop) .eqv. .true. ) block

! set current communicator
         if ( present(gid) .eqv. .true. ) then
             group = gid
         else
             group = MPI_COMM_WORLD
         endif ! back if ( present(gid) .eqv. .true. ) block

! barrier until all processes reach here
         call mp_barrier(group)

! invoke related mpi subroutines
         call MPI_REDUCE(source, data, 1, m_rdp, opera, root, group, ierror)

! handler for return code
         call mp_error('mp_reduce_rdp0', ierror)

         return
     end subroutine mp_reduce_rdp0

!!
!! @sub mp_reduce_rdp1
!!
!! reduce real(dp) vector from all processes
!!
     subroutine mp_reduce_rdp1(source, data, root, mop, gid)
         implicit none

! external arguments
         real(dp), intent(in) :: source(:)
         real(dp), intent(inout) :: data(:)
         integer, intent(in) :: root
         integer, optional, intent(in) :: mop
         integer, optional, intent(in) :: gid

! set current operator
         if ( present(mop) .eqv. .true. ) then
             opera = mop
         else
             opera = MPI_SUM
         endif ! back if ( present(mop) .eqv. .true. ) block

! set current communicator
         if ( present(gid) .eqv. .true. ) then
             group = gid
         else
             group = MPI_COMM_WORLD
         endif ! back if ( present(gid) .eqv. .true. ) block

! barrier until all processes reach here
         call mp_barrier(group)

! setup element count
         isize = size(source)

! invoke related mpi subroutines
         call MPI_REDUCE(source, data, isize, m_rdp, opera, root, group, ierror)

! handler for return code
         call mp_error('mp_reduce_rdp1', ierror)

         return
     end subroutine mp_reduce_rdp1

!!
!! @sub mp_reduce_rdp2
!!
!! reduce real(dp) matrix from all processes
!!
     subroutine mp_reduce_rdp2(source, data, root, mop, gid)
         implicit none

! external arguments
         real(dp), intent(in) :: source(:,:)
         real(dp), intent(inout) :: data(:,:)
         integer, intent(in) :: root
         integer, optional, intent(in) :: mop
         integer, optional, intent(in) :: gid

! set current operator
         if ( present(mop) .eqv. .true. ) then
             opera = mop
         else
             opera = MPI_SUM
         endif ! back if ( present(mop) .eqv. .true. ) block

! set current communicator
         if ( present(gid) .eqv. .true. ) then
             group = gid
         else
             group = MPI_COMM_WORLD
         endif ! back if ( present(gid) .eqv. .true. ) block

! barrier until all processes reach here
         call mp_barrier(group)

! setup element count
         isize = size(source)

! invoke related mpi subroutines
         call MPI_REDUCE(source, data, isize, m_rdp, opera, root, group, ierror)

! handler for return code
         call mp_error('mp_reduce_rdp2', ierror)

         return
     end subroutine mp_reduce_rdp2

!!
!! @sub mp_reduce_rdp3
!!
!! reduce real(dp) matrix from all processes
!!
     subroutine mp_reduce_rdp3(source, data, root, mop, gid)
         implicit none

! external arguments
         real(dp), intent(in) :: source(:,:,:)
         real(dp), intent(inout) :: data(:,:,:)
         integer, intent(in) :: root
         integer, optional, intent(in) :: mop
         integer, optional, intent(in) :: gid

! set current operator
         if ( present(mop) .eqv. .true. ) then
             opera = mop
         else
             opera = MPI_SUM
         endif ! back if ( present(mop) .eqv. .true. ) block

! set current communicator
         if ( present(gid) .eqv. .true. ) then
             group = gid
         else
             group = MPI_COMM_WORLD
         endif ! back if ( present(gid) .eqv. .true. ) block

! barrier until all processes reach here
         call mp_barrier(group)

! setup element count
         isize = size(source)

! invoke related mpi subroutines
         call MPI_REDUCE(source, data, isize, m_rdp, opera, root, group, ierror)

! handler for return code
         call mp_error('mp_reduce_rdp3', ierror)

         return
     end subroutine mp_reduce_rdp3

!!
!! @sub mp_reduce_rdp4
!!
!! reduce real(dp) matrix from all processes
!!
     subroutine mp_reduce_rdp4(source, data, root, mop, gid)
         implicit none

! external arguments
         real(dp), intent(in) :: source(:,:,:,:)
         real(dp), intent(inout) :: data(:,:,:,:)
         integer, intent(in) :: root
         integer, optional, intent(in) :: mop
         integer, optional, intent(in) :: gid

! set current operator
         if ( present(mop) .eqv. .true. ) then
             opera = mop
         else
             opera = MPI_SUM
         endif ! back if ( present(mop) .eqv. .true. ) block

! set current communicator
         if ( present(gid) .eqv. .true. ) then
             group = gid
         else
             group = MPI_COMM_WORLD
         endif ! back if ( present(gid) .eqv. .true. ) block

! barrier until all processes reach here
         call mp_barrier(group)

! setup element count
         isize = size(source)

! invoke related mpi subroutines
         call MPI_REDUCE(source, data, isize, m_rdp, opera, root, group, ierror)

! handler for return code
         call mp_error('mp_reduce_rdp4', ierror)

         return
     end subroutine mp_reduce_rdp4

!!
!! @sub mp_reduce_rdp5
!!
!! reduce real(dp) matrix from all processes
!!
     subroutine mp_reduce_rdp5(source, data, root, mop, gid)
         implicit none

! external arguments
         real(dp), intent(in) :: source(:,:,:,:,:)
         real(dp), intent(inout) :: data(:,:,:,:,:)
         integer, intent(in) :: root
         integer, optional, intent(in) :: mop
         integer, optional, intent(in) :: gid

! set current operator
         if ( present(mop) .eqv. .true. ) then
             opera = mop
         else
             opera = MPI_SUM
         endif ! back if ( present(mop) .eqv. .true. ) block

! set current communicator
         if ( present(gid) .eqv. .true. ) then
             group = gid
         else
             group = MPI_COMM_WORLD
         endif ! back if ( present(gid) .eqv. .true. ) block

! barrier until all processes reach here
         call mp_barrier(group)

! setup element count
         isize = size(source)

! invoke related mpi subroutines
         call MPI_REDUCE(source, data, isize, m_rdp, opera, root, group, ierror)

! handler for return code
         call mp_error('mp_reduce_rdp5', ierror)

         return
     end subroutine mp_reduce_rdp5

!!
!! @sub mp_reduce_cdp0
!!
!! reduce 1 complex(dp) from all processes
!!
     subroutine mp_reduce_cdp0(source, data, root, mop, gid)
         implicit none

! external arguments
         complex(dp), intent(in) :: source
         complex(dp), intent(inout) :: data
         integer, intent(in) :: root
         integer, optional, intent(in) :: mop
         integer, optional, intent(in) :: gid

! set current operator
         if ( present(mop) .eqv. .true. ) then
             opera = mop
         else
             opera = MPI_SUM
         endif ! back if ( present(mop) .eqv. .true. ) block

! set current communicator
         if ( present(gid) .eqv. .true. ) then
             group = gid
         else
             group = MPI_COMM_WORLD
         endif ! back if ( present(gid) .eqv. .true. ) block

! barrier until all processes reach here
         call mp_barrier(group)

! invoke related mpi subroutines
         call MPI_REDUCE(source, data, 1, m_cdp, opera, root, group, ierror)

! handler for return code
         call mp_error('mp_reduce_cdp0', ierror)

         return
     end subroutine mp_reduce_cdp0

!!
!! @sub mp_reduce_cdp1
!!
!! reduce complex(dp) vector from all processes
!!
     subroutine mp_reduce_cdp1(source, data, root, mop, gid)
         implicit none

! external arguments
         complex(dp), intent(in) :: source(:)
         complex(dp), intent(inout) :: data(:)
         integer, intent(in) :: root
         integer, optional, intent(in) :: mop
         integer, optional, intent(in) :: gid

! set current operator
         if ( present(mop) .eqv. .true. ) then
             opera = mop
         else
             opera = MPI_SUM
         endif ! back if ( present(mop) .eqv. .true. ) block

! set current communicator
         if ( present(gid) .eqv. .true. ) then
             group = gid
         else
             group = MPI_COMM_WORLD
         endif ! back if ( present(gid) .eqv. .true. ) block

! barrier until all processes reach here
         call mp_barrier(group)

! setup element count
         isize = size(source)

! invoke related mpi subroutines
         call MPI_REDUCE(source, data, isize, m_cdp, opera, root, group, ierror)

! handler for return code
         call mp_error('mp_reduce_cdp1', ierror)

         return
     end subroutine mp_reduce_cdp1

!!
!! @sub mp_reduce_cdp2
!!
!! reduce complex(dp) matrix from all processes
!!
     subroutine mp_reduce_cdp2(source, data, root, mop, gid)
         implicit none

! external arguments
         complex(dp), intent(in) :: source(:,:)
         complex(dp), intent(inout) :: data(:,:)
         integer, intent(in) :: root
         integer, optional, intent(in) :: mop
         integer, optional, intent(in) :: gid

! set current operator
         if ( present(mop) .eqv. .true. ) then
             opera = mop
         else
             opera = MPI_SUM
         endif ! back if ( present(mop) .eqv. .true. ) block

! set current communicator
         if ( present(gid) .eqv. .true. ) then
             group = gid
         else
             group = MPI_COMM_WORLD
         endif ! back if ( present(gid) .eqv. .true. ) block

! barrier until all processes reach here
         call mp_barrier(group)

! setup element count
         isize = size(source)

! invoke related mpi subroutines
         call MPI_REDUCE(source, data, isize, m_cdp, opera, root, group, ierror)

! handler for return code
         call mp_error('mp_reduce_cdp2', ierror)

         return
     end subroutine mp_reduce_cdp2

!!
!! @sub mp_reduce_cdp3
!!
!! reduce complex(dp) matrix from all processes
!!
     subroutine mp_reduce_cdp3(source, data, root, mop, gid)
         implicit none

! external arguments
         complex(dp), intent(in) :: source(:,:,:)
         complex(dp), intent(inout) :: data(:,:,:)
         integer, intent(in) :: root
         integer, optional, intent(in) :: mop
         integer, optional, intent(in) :: gid

! set current operator
         if ( present(mop) .eqv. .true. ) then
             opera = mop
         else
             opera = MPI_SUM
         endif ! back if ( present(mop) .eqv. .true. ) block

! set current communicator
         if ( present(gid) .eqv. .true. ) then
             group = gid
         else
             group = MPI_COMM_WORLD
         endif ! back if ( present(gid) .eqv. .true. ) block

! barrier until all processes reach here
         call mp_barrier(group)

! setup element count
         isize = size(source)

! invoke related mpi subroutines
         call MPI_REDUCE(source, data, isize, m_cdp, opera, root, group, ierror)

! handler for return code
         call mp_error('mp_reduce_cdp3', ierror)

         return
     end subroutine mp_reduce_cdp3

!!
!! @sub mp_reduce_cdp4
!!
!! reduce complex(dp) matrix from all processes
!!
     subroutine mp_reduce_cdp4(source, data, root, mop, gid)
         implicit none

! external arguments
         complex(dp), intent(in) :: source(:,:,:,:)
         complex(dp), intent(inout) :: data(:,:,:,:)
         integer, intent(in) :: root
         integer, optional, intent(in) :: mop
         integer, optional, intent(in) :: gid

! set current operator
         if ( present(mop) .eqv. .true. ) then
             opera = mop
         else
             opera = MPI_SUM
         endif ! back if ( present(mop) .eqv. .true. ) block

! set current communicator
         if ( present(gid) .eqv. .true. ) then
             group = gid
         else
             group = MPI_COMM_WORLD
         endif ! back if ( present(gid) .eqv. .true. ) block

! barrier until all processes reach here
         call mp_barrier(group)

! setup element count
         isize = size(source)

! invoke related mpi subroutines
         call MPI_REDUCE(source, data, isize, m_cdp, opera, root, group, ierror)

! handler for return code
         call mp_error('mp_reduce_cdp4', ierror)

         return
     end subroutine mp_reduce_cdp4

!!
!! @sub mp_reduce_cdp5
!!
!! reduce complex(dp) matrix from all processes
!!
     subroutine mp_reduce_cdp5(source, data, root, mop, gid)
         implicit none

! external arguments
         complex(dp), intent(in) :: source(:,:,:,:,:)
         complex(dp), intent(inout) :: data(:,:,:,:,:)
         integer, intent(in) :: root
         integer, optional, intent(in) :: mop
         integer, optional, intent(in) :: gid

! set current operator
         if ( present(mop) .eqv. .true. ) then
             opera = mop
         else
             opera = MPI_SUM
         endif ! back if ( present(mop) .eqv. .true. ) block

! set current communicator
         if ( present(gid) .eqv. .true. ) then
             group = gid
         else
             group = MPI_COMM_WORLD
         endif ! back if ( present(gid) .eqv. .true. ) block

! barrier until all processes reach here
         call mp_barrier(group)

! setup element count
         isize = size(source)

! invoke related mpi subroutines
         call MPI_REDUCE(source, data, isize, m_cdp, opera, root, group, ierror)

! handler for return code
         call mp_error('mp_reduce_cdp5', ierror)

         return
     end subroutine mp_reduce_cdp5

!!========================================================================
!!>>> MPI collective operations: allreducing                           <<<
!!========================================================================

!!
!! @sub mp_allreduce_int0
!!
!! reduce 1 integer from all processes
!!
     subroutine mp_allreduce_int0(source, data, mop, gid)
         implicit none

! external arguments
         integer, intent(in) :: source
         integer, intent(inout) :: data
         integer, optional, intent(in) :: mop
         integer, optional, intent(in) :: gid

! set current operator
         if ( present(mop) .eqv. .true. ) then
             opera = mop
         else
             opera = MPI_SUM
         endif ! back if ( present(mop) .eqv. .true. ) block

! set current communicator
         if ( present(gid) .eqv. .true. ) then
             group = gid
         else
             group = MPI_COMM_WORLD
         endif ! back if ( present(gid) .eqv. .true. ) block

! barrier until all processes reach here
         call mp_barrier(group)

! invoke related mpi subroutines
         call MPI_ALLREDUCE(source, data, 1, m_int, opera, group, ierror)

! handler for return code
         call mp_error('mp_allreduce_int0', ierror)

         return
     end subroutine mp_allreduce_int0

!!
!! @sub mp_allreduce_int1
!!
!! reduce integer vector from all processes
!!
     subroutine mp_allreduce_int1(source, data, mop, gid)
         implicit none

! external arguments
         integer, intent(in) :: source(:)
         integer, intent(inout) :: data(:)
         integer, optional, intent(in) :: mop
         integer, optional, intent(in) :: gid

! set current operator
         if ( present(mop) .eqv. .true. ) then
             opera = mop
         else
             opera = MPI_SUM
         endif ! back if ( present(mop) .eqv. .true. ) block

! set current communicator
         if ( present(gid) .eqv. .true. ) then
             group = gid
         else
             group = MPI_COMM_WORLD
         endif ! back if ( present(gid) .eqv. .true. ) block

! barrier until all processes reach here
         call mp_barrier(group)

! setup element count
         isize = size(source)

! invoke related mpi subroutines
         call MPI_ALLREDUCE(source, data, isize, m_int, opera, group, ierror)

! handler for return code
         call mp_error('mp_allreduce_int1', ierror)

         return
     end subroutine mp_allreduce_int1

!!
!! @sub mp_allreduce_int2
!!
!! reduce integer matrix from all processes
!!
     subroutine mp_allreduce_int2(source, data, mop, gid)
         implicit none

! external arguments
         integer, intent(in) :: source(:,:)
         integer, intent(inout) :: data(:,:)
         integer, optional, intent(in) :: mop
         integer, optional, intent(in) :: gid

! set current operator
         if ( present(mop) .eqv. .true. ) then
             opera = mop
         else
             opera = MPI_SUM
         endif ! back if ( present(mop) .eqv. .true. ) block

! set current communicator
         if ( present(gid) .eqv. .true. ) then
             group = gid
         else
             group = MPI_COMM_WORLD
         endif ! back if ( present(gid) .eqv. .true. ) block

! barrier until all processes reach here
         call mp_barrier(group)

! setup element count
         isize = size(source)

! invoke related mpi subroutines
         call MPI_ALLREDUCE(source, data, isize, m_int, opera, group, ierror)

! handler for return code
         call mp_error('mp_allreduce_int2', ierror)

         return
     end subroutine mp_allreduce_int2

!!
!! @sub mp_allreduce_int3
!!
!! reduce integer matrix from all processes
!!
     subroutine mp_allreduce_int3(source, data, mop, gid)
         implicit none

! external arguments
         integer, intent(in) :: source(:,:,:)
         integer, intent(inout) :: data(:,:,:)
         integer, optional, intent(in) :: mop
         integer, optional, intent(in) :: gid

! set current operator
         if ( present(mop) .eqv. .true. ) then
             opera = mop
         else
             opera = MPI_SUM
         endif ! back if ( present(mop) .eqv. .true. ) block

! set current communicator
         if ( present(gid) .eqv. .true. ) then
             group = gid
         else
             group = MPI_COMM_WORLD
         endif ! back if ( present(gid) .eqv. .true. ) block

! barrier until all processes reach here
         call mp_barrier(group)

! setup element count
         isize = size(source)

! invoke related mpi subroutines
         call MPI_ALLREDUCE(source, data, isize, m_int, opera, group, ierror)

! handler for return code
         call mp_error('mp_allreduce_int3', ierror)

         return
     end subroutine mp_allreduce_int3

!!
!! @sub mp_allreduce_int4
!!
!! reduce integer matrix from all processes
!!
     subroutine mp_allreduce_int4(source, data, mop, gid)
         implicit none

! external arguments
         integer, intent(in) :: source(:,:,:,:)
         integer, intent(inout) :: data(:,:,:,:)
         integer, optional, intent(in) :: mop
         integer, optional, intent(in) :: gid

! set current operator
         if ( present(mop) .eqv. .true. ) then
             opera = mop
         else
             opera = MPI_SUM
         endif ! back if ( present(mop) .eqv. .true. ) block

! set current communicator
         if ( present(gid) .eqv. .true. ) then
             group = gid
         else
             group = MPI_COMM_WORLD
         endif ! back if ( present(gid) .eqv. .true. ) block

! barrier until all processes reach here
         call mp_barrier(group)

! setup element count
         isize = size(source)

! invoke related mpi subroutines
         call MPI_ALLREDUCE(source, data, isize, m_int, opera, group, ierror)

! handler for return code
         call mp_error('mp_allreduce_int4', ierror)

         return
     end subroutine mp_allreduce_int4

!!
!! @sub mp_allreduce_int5
!!
!! reduce integer matrix from all processes
!!
     subroutine mp_allreduce_int5(source, data, mop, gid)
         implicit none

! external arguments
         integer, intent(in) :: source(:,:,:,:,:)
         integer, intent(inout) :: data(:,:,:,:,:)
         integer, optional, intent(in) :: mop
         integer, optional, intent(in) :: gid

! set current operator
         if ( present(mop) .eqv. .true. ) then
             opera = mop
         else
             opera = MPI_SUM
         endif ! back if ( present(mop) .eqv. .true. ) block

! set current communicator
         if ( present(gid) .eqv. .true. ) then
             group = gid
         else
             group = MPI_COMM_WORLD
         endif ! back if ( present(gid) .eqv. .true. ) block

! barrier until all processes reach here
         call mp_barrier(group)

! setup element count
         isize = size(source)

! invoke related mpi subroutines
         call MPI_ALLREDUCE(source, data, isize, m_int, opera, group, ierror)

! handler for return code
         call mp_error('mp_allreduce_int5', ierror)

         return
     end subroutine mp_allreduce_int5

!!
!! @sub mp_allreduce_rdp0
!!
!! reduce 1 real(dp) from all processes
!!
     subroutine mp_allreduce_rdp0(source, data, mop, gid)
         implicit none

! external arguments
         real(dp), intent(in) :: source
         real(dp), intent(inout) :: data
         integer, optional, intent(in) :: mop
         integer, optional, intent(in) :: gid

! set current operator
         if ( present(mop) .eqv. .true. ) then
             opera = mop
         else
             opera = MPI_SUM
         endif ! back if ( present(mop) .eqv. .true. ) block

! set current communicator
         if ( present(gid) .eqv. .true. ) then
             group = gid
         else
             group = MPI_COMM_WORLD
         endif ! back if ( present(gid) .eqv. .true. ) block

! barrier until all processes reach here
         call mp_barrier(group)

! invoke related mpi subroutines
         call MPI_ALLREDUCE(source, data, 1, m_rdp, opera, group, ierror)

! handler for return code
         call mp_error('mp_allreduce_rdp0', ierror)

         return
     end subroutine mp_allreduce_rdp0

!!
!! @sub mp_allreduce_rdp1
!!
!! reduce real(dp) vector from all processes
!!
     subroutine mp_allreduce_rdp1(source, data, mop, gid)
         implicit none

! external arguments
         real(dp), intent(in) :: source(:)
         real(dp), intent(inout) :: data(:)
         integer, optional, intent(in) :: mop
         integer, optional, intent(in) :: gid

! set current operator
         if ( present(mop) .eqv. .true. ) then
             opera = mop
         else
             opera = MPI_SUM
         endif ! back if ( present(mop) .eqv. .true. ) block

! set current communicator
         if ( present(gid) .eqv. .true. ) then
             group = gid
         else
             group = MPI_COMM_WORLD
         endif ! back if ( present(gid) .eqv. .true. ) block

! barrier until all processes reach here
         call mp_barrier(group)

! setup element count
         isize = size(source)

! invoke related mpi subroutines
         call MPI_ALLREDUCE(source, data, isize, m_rdp, opera, group, ierror)

! handler for return code
         call mp_error('mp_allreduce_rdp1', ierror)

         return
     end subroutine mp_allreduce_rdp1

!!
!! @sub mp_allreduce_rdp2
!!
!! reduce real(dp) matrix from all processes
!!
     subroutine mp_allreduce_rdp2(source, data, mop, gid)
         implicit none

! external arguments
         real(dp), intent(in) :: source(:,:)
         real(dp), intent(inout) :: data(:,:)
         integer, optional, intent(in) :: mop
         integer, optional, intent(in) :: gid

! set current operator
         if ( present(mop) .eqv. .true. ) then
             opera = mop
         else
             opera = MPI_SUM
         endif ! back if ( present(mop) .eqv. .true. ) block

! set current communicator
         if ( present(gid) .eqv. .true. ) then
             group = gid
         else
             group = MPI_COMM_WORLD
         endif ! back if ( present(gid) .eqv. .true. ) block

! barrier until all processes reach here
         call mp_barrier(group)

! setup element count
         isize = size(source)

! invoke related mpi subroutines
         call MPI_ALLREDUCE(source, data, isize, m_rdp, opera, group, ierror)

! handler for return code
         call mp_error('mp_allreduce_rdp2', ierror)

         return
     end subroutine mp_allreduce_rdp2

!!
!! @sub mp_allreduce_rdp3
!!
!! reduce real(dp) matrix from all processes
!!
     subroutine mp_allreduce_rdp3(source, data, mop, gid)
         implicit none

! external arguments
         real(dp), intent(in) :: source(:,:,:)
         real(dp), intent(inout) :: data(:,:,:)
         integer, optional, intent(in) :: mop
         integer, optional, intent(in) :: gid

! set current operator
         if ( present(mop) .eqv. .true. ) then
             opera = mop
         else
             opera = MPI_SUM
         endif ! back if ( present(mop) .eqv. .true. ) block

! set current communicator
         if ( present(gid) .eqv. .true. ) then
             group = gid
         else
             group = MPI_COMM_WORLD
         endif ! back if ( present(gid) .eqv. .true. ) block

! barrier until all processes reach here
         call mp_barrier(group)

! setup element count
         isize = size(source)

! invoke related mpi subroutines
         call MPI_ALLREDUCE(source, data, isize, m_rdp, opera, group, ierror)

! handler for return code
         call mp_error('mp_allreduce_rdp3', ierror)

         return
     end subroutine mp_allreduce_rdp3

!!
!! @sub mp_allreduce_rdp4
!!
!! reduce real(dp) matrix from all processes
!!
     subroutine mp_allreduce_rdp4(source, data, mop, gid)
         implicit none

! external arguments
         real(dp), intent(in) :: source(:,:,:,:)
         real(dp), intent(inout) :: data(:,:,:,:)
         integer, optional, intent(in) :: mop
         integer, optional, intent(in) :: gid

! set current operator
         if ( present(mop) .eqv. .true. ) then
             opera = mop
         else
             opera = MPI_SUM
         endif ! back if ( present(mop) .eqv. .true. ) block

! set current communicator
         if ( present(gid) .eqv. .true. ) then
             group = gid
         else
             group = MPI_COMM_WORLD
         endif ! back if ( present(gid) .eqv. .true. ) block

! barrier until all processes reach here
         call mp_barrier(group)

! setup element count
         isize = size(source)

! invoke related mpi subroutines
         call MPI_ALLREDUCE(source, data, isize, m_rdp, opera, group, ierror)

! handler for return code
         call mp_error('mp_allreduce_rdp4', ierror)

         return
     end subroutine mp_allreduce_rdp4

!!
!! @sub mp_allreduce_rdp5
!!
!! reduce real(dp) matrix from all processes
!!
     subroutine mp_allreduce_rdp5(source, data, mop, gid)
         implicit none

! external arguments
         real(dp), intent(in) :: source(:,:,:,:,:)
         real(dp), intent(inout) :: data(:,:,:,:,:)
         integer, optional, intent(in) :: mop
         integer, optional, intent(in) :: gid

! set current operator
         if ( present(mop) .eqv. .true. ) then
             opera = mop
         else
             opera = MPI_SUM
         endif ! back if ( present(mop) .eqv. .true. ) block

! set current communicator
         if ( present(gid) .eqv. .true. ) then
             group = gid
         else
             group = MPI_COMM_WORLD
         endif ! back if ( present(gid) .eqv. .true. ) block

! barrier until all processes reach here
         call mp_barrier(group)

! setup element count
         isize = size(source)

! invoke related mpi subroutines
         call MPI_ALLREDUCE(source, data, isize, m_rdp, opera, group, ierror)

! handler for return code
         call mp_error('mp_allreduce_rdp5', ierror)

         return
     end subroutine mp_allreduce_rdp5

!!
!! @sub mp_allreduce_cdp0
!!
!! reduce 1 complex(dp) from all processes
!!
     subroutine mp_allreduce_cdp0(source, data, mop, gid)
         implicit none

! external arguments
         complex(dp), intent(in) :: source
         complex(dp), intent(inout) :: data
         integer, optional, intent(in) :: mop
         integer, optional, intent(in) :: gid

! set current operator
         if ( present(mop) .eqv. .true. ) then
             opera = mop
         else
             opera = MPI_SUM
         endif ! back if ( present(mop) .eqv. .true. ) block

! set current communicator
         if ( present(gid) .eqv. .true. ) then
             group = gid
         else
             group = MPI_COMM_WORLD
         endif ! back if ( present(gid) .eqv. .true. ) block

! barrier until all processes reach here
         call mp_barrier(group)

! invoke related mpi subroutines
         call MPI_ALLREDUCE(source, data, 1, m_cdp, opera, group, ierror)

! handler for return code
         call mp_error('mp_allreduce_cdp0', ierror)

         return
     end subroutine mp_allreduce_cdp0

!!
!! @sub mp_allreduce_cdp1
!!
!! reduce complex(dp) vector from all processes
!!
     subroutine mp_allreduce_cdp1(source, data, mop, gid)
         implicit none

! external arguments
         complex(dp), intent(in) :: source(:)
         complex(dp), intent(inout) :: data(:)
         integer, optional, intent(in) :: mop
         integer, optional, intent(in) :: gid

! set current operator
         if ( present(mop) .eqv. .true. ) then
             opera = mop
         else
             opera = MPI_SUM
         endif ! back if ( present(mop) .eqv. .true. ) block

! set current communicator
         if ( present(gid) .eqv. .true. ) then
             group = gid
         else
             group = MPI_COMM_WORLD
         endif ! back if ( present(gid) .eqv. .true. ) block

! barrier until all processes reach here
         call mp_barrier(group)

! setup element count
         isize = size(source)

! invoke related mpi subroutines
         call MPI_ALLREDUCE(source, data, isize, m_cdp, opera, group, ierror)

! handler for return code
         call mp_error('mp_allreduce_cdp1', ierror)

         return
     end subroutine mp_allreduce_cdp1

!!
!! @sub mp_allreduce_cdp2
!!
!! reduce complex(dp) matrix from all processes
!!
     subroutine mp_allreduce_cdp2(source, data, mop, gid)
         implicit none

! external arguments
         complex(dp), intent(in) :: source(:,:)
         complex(dp), intent(inout) :: data(:,:)
         integer, optional, intent(in) :: mop
         integer, optional, intent(in) :: gid

! set current operator
         if ( present(mop) .eqv. .true. ) then
             opera = mop
         else
             opera = MPI_SUM
         endif ! back if ( present(mop) .eqv. .true. ) block

! set current communicator
         if ( present(gid) .eqv. .true. ) then
             group = gid
         else
             group = MPI_COMM_WORLD
         endif ! back if ( present(gid) .eqv. .true. ) block

! barrier until all processes reach here
         call mp_barrier(group)

! setup element count
         isize = size(source)

! invoke related mpi subroutines
         call MPI_ALLREDUCE(source, data, isize, m_cdp, opera, group, ierror)

! handler for return code
         call mp_error('mp_allreduce_cdp2', ierror)

         return
     end subroutine mp_allreduce_cdp2

!!
!! @sub mp_allreduce_cdp3
!!
!! reduce complex(dp) matrix from all processes
!!
     subroutine mp_allreduce_cdp3(source, data, mop, gid)
         implicit none

! external arguments
         complex(dp), intent(in) :: source(:,:,:)
         complex(dp), intent(inout) :: data(:,:,:)
         integer, optional, intent(in) :: mop
         integer, optional, intent(in) :: gid

! set current operator
         if ( present(mop) .eqv. .true. ) then
             opera = mop
         else
             opera = MPI_SUM
         endif ! back if ( present(mop) .eqv. .true. ) block

! set current communicator
         if ( present(gid) .eqv. .true. ) then
             group = gid
         else
             group = MPI_COMM_WORLD
         endif ! back if ( present(gid) .eqv. .true. ) block

! barrier until all processes reach here
         call mp_barrier(group)

! setup element count
         isize = size(source)

! invoke related mpi subroutines
         call MPI_ALLREDUCE(source, data, isize, m_cdp, opera, group, ierror)

! handler for return code
         call mp_error('mp_allreduce_cdp3', ierror)

         return
     end subroutine mp_allreduce_cdp3

!!
!! @sub mp_allreduce_cdp4
!!
!! reduce complex(dp) matrix from all processes
!!
     subroutine mp_allreduce_cdp4(source, data, mop, gid)
         implicit none

! external arguments
         complex(dp), intent(in) :: source(:,:,:,:)
         complex(dp), intent(inout) :: data(:,:,:,:)
         integer, optional, intent(in) :: mop
         integer, optional, intent(in) :: gid

! set current operator
         if ( present(mop) .eqv. .true. ) then
             opera = mop
         else
             opera = MPI_SUM
         endif ! back if ( present(mop) .eqv. .true. ) block

! set current communicator
         if ( present(gid) .eqv. .true. ) then
             group = gid
         else
             group = MPI_COMM_WORLD
         endif ! back if ( present(gid) .eqv. .true. ) block

! barrier until all processes reach here
         call mp_barrier(group)

! setup element count
         isize = size(source)

! invoke related mpi subroutines
         call MPI_ALLREDUCE(source, data, isize, m_cdp, opera, group, ierror)

! handler for return code
         call mp_error('mp_allreduce_cdp4', ierror)

         return
     end subroutine mp_allreduce_cdp4

!!
!! @sub mp_allreduce_cdp5
!!
!! reduce complex(dp) matrix from all processes
!!
     subroutine mp_allreduce_cdp5(source, data, mop, gid)
         implicit none

! external arguments
         complex(dp), intent(in) :: source(:,:,:,:,:)
         complex(dp), intent(inout) :: data(:,:,:,:,:)
         integer, optional, intent(in) :: mop
         integer, optional, intent(in) :: gid

! set current operator
         if ( present(mop) .eqv. .true. ) then
             opera = mop
         else
             opera = MPI_SUM
         endif ! back if ( present(mop) .eqv. .true. ) block

! set current communicator
         if ( present(gid) .eqv. .true. ) then
             group = gid
         else
             group = MPI_COMM_WORLD
         endif ! back if ( present(gid) .eqv. .true. ) block

! barrier until all processes reach here
         call mp_barrier(group)

! setup element count
         isize = size(source)

! invoke related mpi subroutines
         call MPI_ALLREDUCE(source, data, isize, m_cdp, opera, group, ierror)

! handler for return code
         call mp_error('mp_allreduce_cdp5', ierror)

         return
     end subroutine mp_allreduce_cdp5

!!========================================================================
!!>>> MPI handler for return code                                      <<<
!!========================================================================

!!
!! @sub mp_error
!!
!! deal with the return code of MPI subroutine
!!
# define STR_ERR_COMM      'invalid communicator in mpi call.'
# define STR_ERR_COUNT     'invalid count in mpi call.'
# define STR_ERR_TYPE      'invalid datatype in mpi call.'
# define STR_ERR_BUFFER    'invalid buffer in mpi call.'
# define STR_ERR_ROOT      'invalid root in mpi call.'
# define STR_ERR_ARG       'invalid argument in mpi call.'
# define STR_ERR_TAG       'invalid tag in mpi call.'
# define STR_ERR_RANK      'invalid rank in mpi call.'
# define STR_ERR_GROUP     'null group passed to mpi call.'
# define STR_ERR_OP        'invalid operation in mpi call.'
# define STR_ERR_TOPOLOGY  'invalid topology in mpi call.'
# define STR_ERR_DIMS      'illegal dimension argument in mpi call.'
# define STR_ERR_UNKNOWN   'unknown error in mpi call.'
# define STR_ERR_TRUNCATE  'message truncated on receive in mpi call.'
# define STR_ERR_OTHER     'other error in mpi call.'
# define STR_ERR_INTERN    'internal error code in mpi call.'
# define STR_ERR_IN_STATUS 'look in status for error value.'
# define STR_ERR_PENDING   'pending request in mpi call.'
# define STR_ERR_REQUEST   'illegal mpi_request handle in mpi call.'
# define STR_ERR_LASTCODE  'last error code in mpi call.'
     subroutine mp_error(sub, err)
         implicit none

! external arguments
! subroutine name
         character(len=*), intent(in) :: sub

! error no.
         integer, intent(in) :: err

         select case (err)

             case (MPI_SUCCESS)
                 return

             case (MPI_ERR_COMM)
                 write(mystd,'(2a)') sub, STR_ERR_COMM

             case (MPI_ERR_COUNT)
                 write(mystd,'(2a)') sub, STR_ERR_COUNT

             case (MPI_ERR_TYPE)
                 write(mystd,'(2a)') sub, STR_ERR_TYPE

             case (MPI_ERR_BUFFER)
                 write(mystd,'(2a)') sub, STR_ERR_BUFFER

             case (MPI_ERR_ROOT)
                 write(mystd,'(2a)') sub, STR_ERR_ROOT

             case (MPI_ERR_ARG)
                 write(mystd,'(2a)') sub, STR_ERR_ARG

             case (MPI_ERR_TAG)
                 write(mystd,'(2a)') sub, STR_ERR_TAG

             case (MPI_ERR_RANK)
                 write(mystd,'(2a)') sub, STR_ERR_RANK

             case (MPI_ERR_GROUP)
                 write(mystd,'(2a)') sub, STR_ERR_GROUP

             case (MPI_ERR_OP)
                 write(mystd,'(2a)') sub, STR_ERR_OP

             case (MPI_ERR_TOPOLOGY)
                 write(mystd,'(2a)') sub, STR_ERR_TOPOLOGY

             case (MPI_ERR_DIMS)
                 write(mystd,'(2a)') sub, STR_ERR_DIMS

             case (MPI_ERR_UNKNOWN)
                 write(mystd,'(2a)') sub, STR_ERR_UNKNOWN

             case (MPI_ERR_TRUNCATE)
                 write(mystd,'(2a)') sub, STR_ERR_TRUNCATE

             case (MPI_ERR_OTHER)
                 write(mystd,'(2a)') sub, STR_ERR_OTHER

             case (MPI_ERR_INTERN)
                 write(mystd,'(2a)') sub, STR_ERR_INTERN

             case (MPI_ERR_IN_STATUS)
                 write(mystd,'(2a)') sub, STR_ERR_IN_STATUS

             case (MPI_ERR_PENDING)
                 write(mystd,'(2a)') sub, STR_ERR_PENDING

             case (MPI_ERR_REQUEST)
                 write(mystd,'(2a)') sub, STR_ERR_REQUEST

             case (MPI_ERR_LASTCODE)
                 write(mystd,'(2a)') sub, STR_ERR_LASTCODE

             case default
                 return

         end select

         return
     end subroutine mp_error

  end module mmpi

!!>>> current used compiler is not mpif90
# else   /* MPI */

  module mmpi
     implicit none

!!========================================================================
!!>>> declare global constants                                         <<<
!!========================================================================

! dp: number precision, double precision for reals
     integer, private, parameter :: dp = kind(1.0d0)

! mystd: device descriptor, console output
     integer, private, parameter :: mystd = 6

!!========================================================================
!!>>> declare accessibility for module routines                        <<<
!!========================================================================

!!>>> mpi information operation
     public :: mp_info

!!>>> mpi environment operation
     public :: mp_init
     public :: mp_finalize
     public :: mp_comm_rank
     public :: mp_comm_size

!!>>> synchronics operations
     public :: mp_barrier

!!>>> broadcasting operations
     public :: mp_bcast

!!>>> gathering operations
     public :: mp_gather

!!>>> gatherving operations
     public :: mp_gatherv

!!>>> allgathering operations
     public :: mp_allgather

!!>>> allgatherving operations
     public :: mp_allgatherv

!!>>> reducing operations
     public :: mp_reduce

!!>>> allreducing operations
     public :: mp_allreduce

  contains

!!========================================================================
!!>>> MPI information operations                                       <<<
!!========================================================================

!!
!! @sub mp_info
!!
!! return the current information about mpi environment
!!
     subroutine mp_info()
         implicit none

         write(mystd,'(a)') 'Bad news, your compiler is mpi-incompatible.'

         return
     end subroutine mp_info

!!========================================================================
!!>>> MPI initialize and finalize operations                           <<<
!!========================================================================

!!
!! @sub mp_init
!!
!! initialize mpi environment
!!
     subroutine mp_init()
         implicit none

         return
     end subroutine mp_init

!!
!! @sub mp_finalize
!!
!! finalize mpi environment
!!
     subroutine mp_finalize()
         implicit none

         return
     end subroutine mp_finalize

!!========================================================================
!!>>> MPI setup operations                                             <<<
!!========================================================================

!!
!! @sub mp_comm_rank
!!
!! determine the rank of the current process
!!
     subroutine mp_comm_rank()
         implicit none

         return
     end subroutine mp_comm_rank

!!
!! @sub mp_comm_size
!!
!! evaluate the number of processes in current communicator
!!
     subroutine mp_comm_size()
         implicit none

         return
     end subroutine mp_comm_size

!!========================================================================
!!>>> MPI barrier operations                                           <<<
!!========================================================================

!!
!! @sub mp_barrier
!!
!! blocks until all process have reached this routine
!!
     subroutine mp_barrier()
         implicit none

         return
     end subroutine mp_barrier

!!========================================================================
!!>>> MPI collective operations: broadcasting                          <<<
!!========================================================================

!!
!! @sub mp_bcast
!!
!! broadcasts sth. from the process with rank "root"
!!
     subroutine mp_bcast()
         implicit none

         return
     end subroutine mp_bcast

!!========================================================================
!!>>> MPI collective operations : gathering                            <<<
!!========================================================================

!!
!! @sub mp_gather_int1
!!
!! gather sth. from every processes to rank 0
!!
     subroutine mp_gather()
         implicit none

         return
     end subroutine mp_gather

!!========================================================================
!!>>> MPI collective operations : gatherving                           <<<
!!========================================================================

!!
!! @sub mp_gatherv
!!
!! gather sth. from every processes to rank 0
!!
     subroutine mp_gatherv()
         implicit none

         return
     end subroutine mp_gatherv

!!========================================================================
!!>>> MPI collective operations: allgathering                          <<<
!!========================================================================

!!
!! @sub mp_allgather
!!
!! gather sth. from all processes and then redistribute it to all processes
!!
     subroutine mp_allgather()
         implicit none

         return
     end subroutine mp_allgather

!!========================================================================
!!>>> MPI collective operations: allgatherving                         <<<
!!========================================================================

!!
!! @sub mp_allgatherv
!!
!! gather sth. from all processes and then redistribute it to all processes
!!
     subroutine mp_allgatherv()
         implicit none

         return
     end subroutine mp_allgatherv

!!========================================================================
!!>>> MPI collective operations: reducing                              <<<
!!========================================================================

!!
!! @sub mp_reduce
!!
!! reduce sth. from all processes
!!
     subroutine mp_reduce()
         implicit none

         return
     end subroutine mp_reduce

!!========================================================================
!!>>> MPI collective operations: allreducing                           <<<
!!========================================================================

!!
!! @sub mp_allreduce
!!
!! reduce sth. from all processes
!!
     subroutine mp_allreduce()
         implicit none

         return
     end subroutine mp_allreduce

  end module mmpi

# endif  /* MPI */
