!-------------------------------------------------------------------------
! project : azalea@fantasy
! program : mmpi
! source  : mod_mpi.f90
! type    : module
! author  : li huang (email:huangli712yahoo.com.cn)
! history : 08/09/2006 by li huang
!           08/14/2006 by li huang
!           08/18/2006 by li huang
!           08/21/2006 by li huang
!           08/25/2006 by li huang
!           08/29/2006 by li huang
!           09/11/2006 by li huang
!           09/13/2006 by li huang
!           09/26/2006 by li huang
!           10/31/2006 by li huang
!           11/15/2006 by li huang
!           11/26/2006 by li huang
!           12/10/2006 by li huang
!           02/20/2007 by li huang
!           02/28/2007 by li huang
!           03/02/2007 by li huang
!           03/14/2007 by li huang
!           04/17/2007 by li huang
!           04/19/2007 by li huang
!           01/11/2008 by li huang
!           10/24/2008 by li huang
!           12/15/2008 by li huang
!           04/27/2009 by li huang
!           07/29/2009 by li huang
!           08/22/2009 by li huang
!           12/18/2009 by li huang
!           02/27/2010 by li huang
! purpose : to define my own mpi calls, inspired by quantum-espresso code.
!           note that the original mpi subroutines are rather complicated
!           for newbies. thus we try to wrap the most important and useful
!           (not all) mpi subroutines, including
!               MPI_INIT(),
!               MPI_FINALIZE(),
!               MPI_WTIME(),
!               MPI_WTICK(),
!               MPI_BARRIER(),
!               MPI_DIMS_CREATE(),
!               MPI_CART_CREATE(),
!               MPI_CART_COORDS(),
!               MPI_COMM_SPLIT(),
!               MPI_COMM_RANK(),
!               MPI_COMM_SIZE(),
!               MPI_GET_PROCESSOR_NAME(),
!               MPI_BCAST(),
!               MPI_GATHER(),
!               MPI_GATHERV(),
!               MPI_ALLGATHER(),
!               MPI_ALLGATHERV(),
!               MPI_REDUCE(),
!               MPI_ALLREDUCE(),
!           etc, to facilite the usage of mpi. in the module, we also
!           implement a new light-weight error handler. enjoy it!
! status  : unstable
! comment : this module has been tested under the following environment:
!               mpich1    1.2.7p1
!               mpich2    1.2.1p1
!               mvapich2  1.2.0p1
!               intel mpi 3.2.0
!-------------------------------------------------------------------------

! whether the compiler support mpi environment, i.e, mpif90
# if defined (MPI)

  module mmpi
     use mpi

     implicit none

!-------------------------------------------------------------------------
!::: declare global constants                                          :::
!-------------------------------------------------------------------------

! dp: number precision, double precision for reals
     integer, private, parameter :: dp = kind(1.0d0)

! mystd: device descriptor, console output
     integer, private, parameter :: mystd = 6

!-------------------------------------------------------------------------
!::: declare mpi constants (datatypes)                                 :::
!-------------------------------------------------------------------------

! mpi_log: datatype, boolean
     integer, private, parameter :: mpi_log    = MPI_LOGICAL

! mpi_mint: datatype, integer
     integer, private, parameter :: mpi_mint   = MPI_INTEGER

! mpi_dreal: datatype, double precision float
     integer, private, parameter :: mpi_dreal  = MPI_DOUBLE_PRECISION

! mpi_dcmplx: datatype, double precision complex
     integer, private, parameter :: mpi_dcmplx = MPI_DOUBLE_COMPLEX

!-------------------------------------------------------------------------
!::: declare common constants                                          :::
!-------------------------------------------------------------------------

! ndims: number of cartesian dimensions, and we need a 2D grid
     integer, private, parameter :: ndims = 2

! reorder: ranking may be reordered (true) or not (false) in a new grid
     logical, private, parameter :: reorder = .false.

! periods: specifying whether the grid is periodic or not in each dimens
     logical, private, parameter :: periods(ndims) = .false.

!-------------------------------------------------------------------------
!::: declare common variables                                          :::
!-------------------------------------------------------------------------

! ierror: error code for mpi subroutines
     integer, private :: ierror

! istat: status code for mpi subroutines
     integer, private :: istat

! group: default communicator
     integer, private :: group

! isize: size of elements
     integer, private :: isize

!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

! mpi_comm_cart: communicator for cartesian topology
     integer, public, save :: mpi_comm_cart

! mpi_comm_row: rowwise-striped communicator for cartesian topology
     integer, public, save :: mpi_comm_row

! mpi_comm_col: columnwise-striped communicator for cartesian topology
     integer, public, save :: mpi_comm_col

!-------------------------------------------------------------------------
!::: declare accessibility for module routines                         :::
!-------------------------------------------------------------------------

!>>> mpi information operation
     public :: mp_info

!>>> mpi environment operation

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

!>>> mpi cartesian topology operation

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

!>>> synchronics operations

! manually block until all processes reach
     public :: mp_barrier

!>>> time handler

! obtain time consuming by parallelized program
     public :: mp_wtime

! obtain time precision of mp_wtime
     public :: mp_wtick

!>>> broadcasting operations

! broadcasting bool
     private :: mp_bcast_bool

! broadcasting bool(:)
     private :: mp_bcast_bool1

! broadcasting bool(:,:)
     private :: mp_bcast_bool2

! broadcasting int
     private :: mp_bcast_int

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
     private :: mp_bcast_rdp

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
     private :: mp_bcast_cdp

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

!>>> gathering operations

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

!>>> gatherving operations

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

!>>> allgathering operations

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

!>>> allgatherving operations

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

!>>> reducing operations

! readucing int
     private :: mp_reduce_int

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
     private :: mp_reduce_rdp

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
     private :: mp_reduce_cdp

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

!>>> allreducing operations

! allreducing int
     private :: mp_allreduce_int

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
     private :: mp_allreduce_rdp

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
     private :: mp_allreduce_cdp

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

!>>> error handler

! echo the error information
     private :: mp_error

!-------------------------------------------------------------------------
!::: declare interface and module procedure                            :::
!-------------------------------------------------------------------------

! we define these interfaces to implement the so called "generic" software
! engineering technique

! mpi_bcast subroutines
     public :: mp_bcast
     interface mp_bcast
         module procedure mp_bcast_bool
         module procedure mp_bcast_bool1
         module procedure mp_bcast_bool2

         module procedure mp_bcast_int
         module procedure mp_bcast_int1
         module procedure mp_bcast_int2
         module procedure mp_bcast_int3
         module procedure mp_bcast_int4
         module procedure mp_bcast_int5

         module procedure mp_bcast_rdp
         module procedure mp_bcast_rdp1
         module procedure mp_bcast_rdp2
         module procedure mp_bcast_rdp3
         module procedure mp_bcast_rdp4
         module procedure mp_bcast_rdp5

         module procedure mp_bcast_cdp
         module procedure mp_bcast_cdp1
         module procedure mp_bcast_cdp2
         module procedure mp_bcast_cdp3
         module procedure mp_bcast_cdp4
         module procedure mp_bcast_cdp5
     end interface mp_bcast

! mpi_gather subroutines
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

! mpi_gatherv subroutines
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

! mpi_allgather subroutines
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

! mpi_allgatherv subroutines
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

! mpi_reduce subroutines
     public :: mp_reduce
     interface mp_reduce
         module procedure mp_reduce_int
         module procedure mp_reduce_int1
         module procedure mp_reduce_int2
         module procedure mp_reduce_int3
         module procedure mp_reduce_int4
         module procedure mp_reduce_int5

         module procedure mp_reduce_rdp
         module procedure mp_reduce_rdp1
         module procedure mp_reduce_rdp2
         module procedure mp_reduce_rdp3
         module procedure mp_reduce_rdp4
         module procedure mp_reduce_rdp5

         module procedure mp_reduce_cdp
         module procedure mp_reduce_cdp1
         module procedure mp_reduce_cdp2
         module procedure mp_reduce_cdp3
         module procedure mp_reduce_cdp4
         module procedure mp_reduce_cdp5
     end interface mp_reduce

! mpi_allreduce subroutines
     public :: mp_allreduce
     interface mp_allreduce
         module procedure mp_allreduce_int
         module procedure mp_allreduce_int1
         module procedure mp_allreduce_int2
         module procedure mp_allreduce_int3
         module procedure mp_allreduce_int4
         module procedure mp_allreduce_int5

         module procedure mp_allreduce_rdp
         module procedure mp_allreduce_rdp1
         module procedure mp_allreduce_rdp2
         module procedure mp_allreduce_rdp3
         module procedure mp_allreduce_rdp4
         module procedure mp_allreduce_rdp5

         module procedure mp_allreduce_cdp
         module procedure mp_allreduce_cdp1
         module procedure mp_allreduce_cdp2
         module procedure mp_allreduce_cdp3
         module procedure mp_allreduce_cdp4
         module procedure mp_allreduce_cdp5
     end interface mp_allreduce

  contains

!-------------------------------------------------------------------------
!::: MPI information operations                                        :::
!-------------------------------------------------------------------------

! mp_info: return the current information about mpi environment
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

!-------------------------------------------------------------------------
!::: MPI initialize and finalize operations                            :::
!-------------------------------------------------------------------------

! mp_init: initialize mpi environment
     subroutine mp_init()
         implicit none

! invoke related mpi subroutines
         call MPI_INIT(ierror)

! handler for return code
         call mp_error('mp_init', ierror)

         return
     end subroutine mp_init

! mp_finalize: finalize mpi environment
     subroutine mp_finalize()
         implicit none

! invoke related mpi subroutines
         call MPI_FINALIZE(ierror)

! handler for return code
         call mp_error('mp_finalize', ierror)

         return
     end subroutine mp_finalize

!-------------------------------------------------------------------------
!::: MPI setup operations                                              :::
!-------------------------------------------------------------------------

! mp_comm_rank: determine the rank of the current process
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
         endif

! invoke related mpi subroutines
         call MPI_COMM_RANK(group, myid, ierror)

! handler for return code
         call mp_error('mp_comm_rank', ierror)

         return
     end subroutine mp_comm_rank

! mp_comm_size: evaluate the number of processes in current communicator
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
         endif

! invoke related mpi subroutines
         call MPI_COMM_SIZE(group, nprocs, ierror)

! handler for return code
         call mp_error('mp_comm_size', ierror)

         return
     end subroutine mp_comm_size

! mp_processor: determine the current workstation's name
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

!-------------------------------------------------------------------------
!::: MPI cartesian topology operations                                 :::
!-------------------------------------------------------------------------

! mp_dims_create: creates a division of processors in a cartesian grid
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

! mp_cart_create: makes a new communicator to which topology is cartesian
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

! mp_cart_coords: determines process coords in cartesian topology
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

! mp_comm_split_row: creates new communicators based on colors and keys
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

! mp_comm_split_col: creates new communicators based on colors and keys
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

!-------------------------------------------------------------------------
!::: MPI barrier operations                                            :::
!-------------------------------------------------------------------------

! mp_barrier: blocks until all process have reached this routine
     subroutine mp_barrier(gid)
         implicit none

! external arguments
         integer, optional, intent(in) :: gid

! set current communicator
         if ( present(gid) .eqv. .true. ) then
             group = gid
         else
             group = MPI_COMM_WORLD
         endif

! invoke related mpi subroutines
         call MPI_BARRIER(group, ierror)

! handler for return code
         call mp_error('mp_barrier', ierror)

         return
     end subroutine mp_barrier

!-------------------------------------------------------------------------
!::: MPI time operations                                               :::
!-------------------------------------------------------------------------

! mp_wtime: returns an elapsed time on the calling processor
     subroutine mp_wtime(time)
         implicit none

! external arguments
         real(dp), intent(out) :: time

! invoke related mpi subroutines
         time = MPI_WTIME()

         return
     end subroutine mp_wtime

! mp_wtick: returns the resolution of MPI_Wtime
     subroutine mp_wtick(tick)
         implicit none

! external arguments
         real(dp), intent(out) :: tick

! invoke related mpi subroutines
         tick = MPI_WTICK()

         return
     end subroutine mp_wtick

!-------------------------------------------------------------------------
!::: MPI collective operations: broadcasting                           :::
!-------------------------------------------------------------------------

! mp_bcast_bool: broadcasts bool from the process with rank "root"
     subroutine mp_bcast_bool(data, root, gid)
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
         endif

! barrier until all processes reach here
         call mp_barrier(group)

! invoke realted MPI subroutines
         call MPI_BCAST(data, 1, mpi_log, root, group, ierror)

! handler for return code
         call mp_error('mp_bcast_bool', ierror)

         return
     end subroutine mp_bcast_bool

! mp_bcast_bool1: broadcasts bool(:) from the process with rank "root"
     subroutine mp_bcast_bool1(data, root, gid)
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
         endif

! barrier until all processes reach here
         call mp_barrier(group)

! setup element count
         isize = size(data)

! invoke realted MPI subroutines
         call MPI_BCAST(data, isize, mpi_log, root, group, ierror)

! handler for return code
         call mp_error('mp_bcast_bool1', ierror)

         return
     end subroutine mp_bcast_bool1

! mp_bcast_bool2: broadcasts bool(:,:) from the process with rank "root"
     subroutine mp_bcast_bool2(data, root, gid)
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
         endif

! barrier until all processes reach here
         call mp_barrier(group)

! setup element count
         isize = size(data)

! invoke realted MPI subroutines
         call MPI_BCAST(data, isize, mpi_log, root, group, ierror)

! handler for return code
         call mp_error('mp_bcast_bool2', ierror)

         return
     end subroutine mp_bcast_bool2

! mp_bcast_int: broadcasts int from the process with rank "root"
     subroutine mp_bcast_int(data, root, gid)
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
         endif

! barrier until all processes reach here
         call mp_barrier(group)

! invoke realted MPI subroutines
         call MPI_BCAST(data, 1, mpi_mint, root, group, ierror)

! handler for return code
         call mp_error('mp_bcast_int', ierror)

         return
     end subroutine mp_bcast_int

! mp_bcast_int1: broadcasts int(:) from the process with rank "root"
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
         endif

! barrier until all processes reach here
         call mp_barrier(group)

! setup element count
         isize = size(data)

! invoke realted MPI subroutines
         call MPI_BCAST(data, isize, mpi_mint, root, group, ierror)

! handler for return code
         call mp_error('mp_bcast_int1', ierror)

         return
     end subroutine mp_bcast_int1

! mp_bcast_int2: broadcasts int(:,:) from the process with rank "root"
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
         endif

! barrier until all processes reach here
         call mp_barrier(group)

! setup element count
         isize = size(data)

! invoke realted MPI subroutines
         call MPI_BCAST(data, isize, mpi_mint, root, group, ierror)

! handler for return code
         call mp_error('mp_bcast_int2', ierror)

         return
     end subroutine mp_bcast_int2

! mp_bcast_int3: broadcasts int(:,:,:) from the process with rank "root"
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
         endif

! barrier until all processes reach here
         call mp_barrier(group)

! setup element count
         isize = size(data)

! invoke realted MPI subroutines
         call MPI_BCAST(data, isize, mpi_mint, root, group, ierror)

! handler for return code
         call mp_error('mp_bcast_int3', ierror)

         return
     end subroutine mp_bcast_int3

! mp_bcast_int4: broadcasts int4 from the process with rank "root"
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
         endif

! barrier until all processes reach here
         call mp_barrier(group)

! setup element count
         isize = size(data)

! invoke realted MPI subroutines
         call MPI_BCAST(data, isize, mpi_mint, root, group, ierror)

! handler for return code
         call mp_error('mp_bcast_int4', ierror)

         return
     end subroutine mp_bcast_int4

! mp_bcast_int5: broadcasts int5 from the process with rank "root"
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
         endif

! barrier until all processes reach here
         call mp_barrier(group)

! setup element count
         isize = size(data)

! invoke realted MPI subroutines
         call MPI_BCAST(data, isize, mpi_mint, root, group, ierror)

! handler for return code
         call mp_error('mp_bcast_int5', ierror)

         return
     end subroutine mp_bcast_int5

! mp_bcast_rdp: broadcasts real from the process with rank "root"
     subroutine mp_bcast_rdp(data, root, gid)
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
         endif

! barrier until all processes reach here
         call mp_barrier(group)

! invoke related mpi subroutines
         call MPI_BCAST(data, 1, mpi_dreal, root, group, ierror)

! handler for return code
         call mp_error('mp_bcast_rdp', ierror)

         return
     end subroutine mp_bcast_rdp

! mp_bcast_rdp1: broadcasts real(:) from the process with rank "root"
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
         endif

! barrier until all processes reach here
         call mp_barrier(group)

! setup element count
         isize = size(data)

! invoke related mpi subroutines
         call MPI_BCAST(data, isize, mpi_dreal, root, group, ierror)

! handler for return code
         call mp_error('mp_bcast_rdp1', ierror)

         return
     end subroutine mp_bcast_rdp1

! mp_bcast_rdp2: broadcasts real(:,:) from the process with rank "root"
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
         endif

! barrier until all processes reach here
         call mp_barrier(group)

! setup element count
         isize = size(data)

! invoke related mpi subroutines
         call MPI_BCAST(data, isize, mpi_dreal, root, group, ierror)

! handler for return code
         call mp_error('mp_bcast_rdp2', ierror)

         return
     end subroutine mp_bcast_rdp2

! mp_bcast_rdp3: broadcasts real(:,:,:) from the process with rank "root"
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
         endif

! barrier until all processes reach here
         call mp_barrier(group)

! setup element count
         isize = size(data)

! invoke related mpi subroutines
         call MPI_BCAST(data, isize, mpi_dreal, root, group, ierror)

! handler for return code
         call mp_error('mp_bcast_rdp3', ierror)

         return
     end subroutine mp_bcast_rdp3

! mp_bcast_rdp4: broadcasts real4 from the process with rank "root"
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
         endif

! barrier until all processes reach here
         call mp_barrier(group)

! setup element count
         isize = size(data)

! invoke related mpi subroutines
         call MPI_BCAST(data, isize, mpi_dreal, root, group, ierror)

! handler for return code
         call mp_error('mp_bcast_rdp4', ierror)

         return
     end subroutine mp_bcast_rdp4

! mp_bcast_rdp5: broadcasts real5 from the process with rank "root"
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
         endif

! barrier until all processes reach here
         call mp_barrier(group)

! setup element count
         isize = size(data)

! invoke related mpi subroutines
         call MPI_BCAST(data, isize, mpi_dreal, root, group, ierror)

! handler for return code
         call mp_error('mp_bcast_rdp5', ierror)

         return
     end subroutine mp_bcast_rdp5

! mp_bcast_cdp: broadcasts complex from the process with rank "root"
     subroutine mp_bcast_cdp(data, root, gid)
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
         endif

! barrier until all processes reach here
         call mp_barrier(group)

! invoke related mpi subroutines
         call MPI_BCAST(data, 1, mpi_dcmplx, root, group, ierror)

! handler for return code
         call mp_error('mp_bcast_cdp', ierror)

         return
     end subroutine mp_bcast_cdp

! mp_bcast_cdp1: broadcasts complex(:) from the process with rank "root"
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
         endif

! barrier until all processes reach here
         call mp_barrier(group)

! setup element count
         isize = size(data)

! invoke related mpi subroutines
         call MPI_BCAST(data, isize, mpi_dcmplx, root, group, ierror)

! handler for return code
         call mp_error('mp_bcast_cdp1', ierror)

         return
     end subroutine mp_bcast_cdp1

! mp_bcast_cdp2: broadcasts complex2 from the process with rank "root"
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
         endif

! barrier until all processes reach here
         call mp_barrier(group)

! setup element count
         isize = size(data)

! invoke related mpi subroutines
         call MPI_BCAST(data, isize, mpi_dcmplx, root, group, ierror)

! handler for return code
         call mp_error('mp_bcast_cdp2', ierror)

         return
     end subroutine mp_bcast_cdp2

! mp_bcast_cdp3: broadcasts complex3 from the process with rank "root"
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
         endif

! barrier until all processes reach here
         call mp_barrier(group)

! setup element count
         isize = size(data)

! invoke related mpi subroutines
         call MPI_BCAST(data, isize, mpi_dcmplx, root, group, ierror)

! handler for return code
         call mp_error('mp_bcast_cdp3', ierror)

         return
     end subroutine mp_bcast_cdp3

! mp_bcast_cdp4: broadcasts complex4 from the process with rank "root"
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
         endif

! barrier until all processes reach here
         call mp_barrier(group)

! setup element count
         isize = size(data)

! invoke related mpi subroutines
         call MPI_BCAST(data, isize, mpi_dcmplx, root, group, ierror)

! handler for return code
         call mp_error('mp_bcast_cdp4', ierror)

         return
     end subroutine mp_bcast_cdp4

! mp_bcast_cdp5: broadcasts complex5 from the process with rank "root"
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
         endif

! barrier until all processes reach here
         call mp_barrier(group)

! setup element count
         isize = size(data)

! invoke related mpi subroutines
         call MPI_BCAST(data, isize, mpi_dcmplx, root, group, ierror)

! handler for return code
         call mp_error('mp_bcast_cdp5', ierror)

         return
     end subroutine mp_bcast_cdp5

!-------------------------------------------------------------------------
!::: MPI collective operations : gathering                             :::
!-------------------------------------------------------------------------

! mp_gather_int1: gather integer data from every processes to rank 0
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
         endif

! barrier until all processes reach here
         call mp_barrier(group)

! setup element count
         isize = size(send)

! invoke related mpi subroutines
         call MPI_GATHER(send, isize, mpi_mint, data, isize, mpi_mint, root, group, ierror)

! handler for return code
         call mp_error('mp_gather_int1', ierror)

         return
     end subroutine mp_gather_int1

! mp_gather_int2: gather integer data from every processes to rank 0
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
         endif

! barrier until all processes reach here
         call mp_barrier(group)

! setup element count
         isize = size(send)

! invoke related mpi subroutines
         call MPI_GATHER(send, isize, mpi_mint, data, isize, mpi_mint, root, group, ierror)

! handler for return code
         call mp_error('mp_gather_int2', ierror)

         return
     end subroutine mp_gather_int2

! mp_gather_int3: gather integer data from every processes to rank 0
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
         endif

! barrier until all processes reach here
         call mp_barrier(group)

! setup element count
         isize = size(send)

! invoke related mpi subroutines
         call MPI_GATHER(send, isize, mpi_mint, data, isize, mpi_mint, root, group, ierror)

! handler for return code
         call mp_error('mp_gather_int3', ierror)

         return
     end subroutine mp_gather_int3

! mp_gather_int4: gather integer data from every processes to rank 0
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
         endif

! barrier until all processes reach here
         call mp_barrier(group)

! setup element count
         isize = size(send)

! invoke related mpi subroutines
         call MPI_GATHER(send, isize, mpi_mint, data, isize, mpi_mint, root, group, ierror)

! handler for return code
         call mp_error('mp_gather_int4', ierror)

         return
     end subroutine mp_gather_int4

! mp_gather_int5: gather integer data from every processes to rank 0
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
         endif

! barrier until all processes reach here
         call mp_barrier(group)

! setup element count
         isize = size(send)

! invoke related mpi subroutines
         call MPI_GATHER(send, isize, mpi_mint, data, isize, mpi_mint, root, group, ierror)

! handler for return code
         call mp_error('mp_gather_int5', ierror)

         return
     end subroutine mp_gather_int5

! mp_gather_rdp1: gather real(dp) data from every processes to rank 0
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
         endif

! barrier until all processes reach here
         call mp_barrier(group)

! setup element count
         isize = size(send)

! invoke related mpi subroutines
         call MPI_GATHER(send, isize, mpi_dreal, data, isize, mpi_dreal, root, group, ierror)

! handler for return code
         call mp_error('mp_gather_rdp1', ierror)

         return
     end subroutine mp_gather_rdp1

! mp_gather_rdp2: gather real(dp) data from every processes to rank 0
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
         endif

! barrier until all processes reach here
         call mp_barrier(group)

! setup element count
         isize = size(send)

! invoke related mpi subroutines
         call MPI_GATHER(send, isize, mpi_dreal, data, isize, mpi_dreal, root, group, ierror)

! handler for return code
         call mp_error('mp_gather_rdp2', ierror)

         return
     end subroutine mp_gather_rdp2

! mp_gather_rdp3: gather real(dp) data from every processes to rank 0
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
         endif

! barrier until all processes reach here
         call mp_barrier(group)

! setup element count
         isize = size(send)

! invoke related mpi subroutines
         call MPI_GATHER(send, isize, mpi_dreal, data, isize, mpi_dreal, root, group, ierror)

! handler for return code
         call mp_error('mp_gather_rdp3', ierror)

         return
     end subroutine mp_gather_rdp3

! mp_gather_rdp4: gather real(dp) data from every processes to rank 0
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
         endif

! barrier until all processes reach here
         call mp_barrier(group)

! setup element count
         isize = size(send)

! invoke related mpi subroutines
         call MPI_GATHER(send, isize, mpi_dreal, data, isize, mpi_dreal, root, group, ierror)

! handler for return code
         call mp_error('mp_gather_rdp4', ierror)

         return
     end subroutine mp_gather_rdp4

! mp_gather_rdp5: gather real(dp) data from every processes to rank 0
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
         endif

! barrier until all processes reach here
         call mp_barrier(group)

! setup element count
         isize = size(send)

! invoke related mpi subroutines
         call MPI_GATHER(send, isize, mpi_dreal, data, isize, mpi_dreal, root, group, ierror)

! handler for return code
         call mp_error('mp_gather_rdp5', ierror)

         return
     end subroutine mp_gather_rdp5

! mp_gather_cdp1: gather complex(dp) data from every processes to rank 0
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
         endif

! barrier until all processes reach here
         call mp_barrier(group)

! setup element count
         isize = size(send)

! invoke related mpi subroutines
         call MPI_GATHER(send, isize, mpi_dcmplx, data, isize, mpi_dcmplx, root, group, ierror)

! handler for return code
         call mp_error('mp_gather_cdp1', ierror)

         return
     end subroutine mp_gather_cdp1

! mp_gather_cdp2: gather complex(dp) data from every processes to rank 0
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
         endif

! barrier until all processes reach here
         call mp_barrier(group)

! setup element count
         isize = size(send)

! invoke related mpi subroutines
         call MPI_GATHER(send, isize, mpi_dcmplx, data, isize, mpi_dcmplx, root, group, ierror)

! handler for return code
         call mp_error('mp_gather_cdp2', ierror)

         return
     end subroutine mp_gather_cdp2

! mp_gather_cdp3: gather complex(dp) data from every processes to rank 0
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
         endif

! barrier until all processes reach here
         call mp_barrier(group)

! setup element count
         isize = size(send)

! invoke related mpi subroutines
         call MPI_GATHER(send, isize, mpi_dcmplx, data, isize, mpi_dcmplx, root, group, ierror)

! handler for return code
         call mp_error('mp_gather_cdp3', ierror)

         return
     end subroutine mp_gather_cdp3

! mp_gather_cdp4: gather complex(dp) data from every processes to rank 0
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
         endif

! barrier until all processes reach here
         call mp_barrier(group)

! setup element count
         isize = size(send)

! invoke related mpi subroutines
         call MPI_GATHER(send, isize, mpi_dcmplx, data, isize, mpi_dcmplx, root, group, ierror)

! handler for return code
         call mp_error('mp_gather_cdp4', ierror)

         return
     end subroutine mp_gather_cdp4

! mp_gather_cdp5: gather complex(dp) data from every processes to rank 0
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
         endif

! barrier until all processes reach here
         call mp_barrier(group)

! setup element count
         isize = size(send)

! invoke related mpi subroutines
         call MPI_GATHER(send, isize, mpi_dcmplx, data, isize, mpi_dcmplx, root, group, ierror)

! handler for return code
         call mp_error('mp_gather_cdp5', ierror)

         return
     end subroutine mp_gather_cdp5

!-------------------------------------------------------------------------
!::: MPI collective operations : gatherving                            :::
!-------------------------------------------------------------------------

! mp_gatherv_int1: gather integer data from every processes to rank 0
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
         endif

! barrier until all processes reach here
         call mp_barrier(group)

! setup element count
         isize = size(send)

! invoke related mpi subroutines
         call MPI_GATHERV(send, isize, mpi_mint, data, recv, disp, mpi_mint, root, group, ierror)

! handler for return code
         call mp_error('mp_gatherv_int1', ierror)

         return
     end subroutine mp_gatherv_int1

! mp_gatherv_int2: gather integer data from every processes to rank 0
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
         endif

! barrier until all processes reach here
         call mp_barrier(group)

! setup element count
         isize = size(send)

! invoke related mpi subroutines
         call MPI_GATHERV(send, isize, mpi_mint, data, recv, disp, mpi_mint, root, group, ierror)

! handler for return code
         call mp_error('mp_gatherv_int2', ierror)

         return
     end subroutine mp_gatherv_int2

! mp_gatherv_int3: gather integer data from every processes to rank 0
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
         endif

! barrier until all processes reach here
         call mp_barrier(group)

! setup element count
         isize = size(send)

! invoke related mpi subroutines
         call MPI_GATHERV(send, isize, mpi_mint, data, recv, disp, mpi_mint, root, group, ierror)

! handler for return code
         call mp_error('mp_gatherv_int3', ierror)

         return
     end subroutine mp_gatherv_int3

! mp_gatherv_int4: gather integer data from every processes to rank 0
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
         endif

! barrier until all processes reach here
         call mp_barrier(group)

! setup element count
         isize = size(send)

! invoke related mpi subroutines
         call MPI_GATHERV(send, isize, mpi_mint, data, recv, disp, mpi_mint, root, group, ierror)

! handler for return code
         call mp_error('mp_gatherv_int4', ierror)

         return
     end subroutine mp_gatherv_int4

! mp_gatherv_int5: gather integer data from every processes to rank 0
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
         endif

! barrier until all processes reach here
         call mp_barrier(group)

! setup element count
         isize = size(send)

! invoke related mpi subroutines
         call MPI_GATHERV(send, isize, mpi_mint, data, recv, disp, mpi_mint, root, group, ierror)

! handler for return code
         call mp_error('mp_gatherv_int5', ierror)

         return
     end subroutine mp_gatherv_int5

! mp_gatherv_rdp1: gather real data from every processes to rank 0
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
         endif

! barrier until all processes reach here
         call mp_barrier(group)

! setup element count
         isize = size(send)

! invoke related mpi subroutines
         call MPI_GATHERV(send, isize, mpi_dreal, data, recv, disp, mpi_dreal, root, group, ierror)

! handler for return code
         call mp_error('mp_gatherv_rdp1', ierror)

         return
     end subroutine mp_gatherv_rdp1

! mp_gatherv_rdp2: gather real data from every processes to rank 0
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
         endif

! barrier until all processes reach here
         call mp_barrier(group)

! setup element count
         isize = size(send)

! invoke related mpi subroutines
         call MPI_GATHERV(send, isize, mpi_dreal, data, recv, disp, mpi_dreal, root, group, ierror)

! handler for return code
         call mp_error('mp_gatherv_rdp2', ierror)

         return
     end subroutine mp_gatherv_rdp2

! mp_gatherv_rdp3: gather real data from every processes to rank 0
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
         endif

! barrier until all processes reach here
         call mp_barrier(group)

! setup element count
         isize = size(send)

! invoke related mpi subroutines
         call MPI_GATHERV(send, isize, mpi_dreal, data, recv, disp, mpi_dreal, root, group, ierror)

! handler for return code
         call mp_error('mp_gatherv_rdp3', ierror)

         return
     end subroutine mp_gatherv_rdp3

! mp_gatherv_rdp4: gather real data from every processes to rank 0
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
         endif

! barrier until all processes reach here
         call mp_barrier(group)

! setup element count
         isize = size(send)

! invoke related mpi subroutines
         call MPI_GATHERV(send, isize, mpi_dreal, data, recv, disp, mpi_dreal, root, group, ierror)

! handler for return code
         call mp_error('mp_gatherv_rdp4', ierror)

         return
     end subroutine mp_gatherv_rdp4

! mp_gatherv_rdp5: gather real data from every processes to rank 0
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
         endif

! barrier until all processes reach here
         call mp_barrier(group)

! setup element count
         isize = size(send)

! invoke related mpi subroutines
         call MPI_GATHERV(send, isize, mpi_dreal, data, recv, disp, mpi_dreal, root, group, ierror)

! handler for return code
         call mp_error('mp_gatherv_rdp5', ierror)

         return
     end subroutine mp_gatherv_rdp5

! mp_gatherv_cdp1: gather complex data from every processes to rank 0
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
         endif

! barrier until all processes reach here
         call mp_barrier(group)

! setup element count
         isize = size(send)

! invoke related mpi subroutines
         call MPI_GATHERV(send, isize, mpi_dcmplx, data, recv, disp, mpi_dcmplx, root, group, ierror)

! handler for return code
         call mp_error('mp_gatherv_cdp1', ierror)

         return
     end subroutine mp_gatherv_cdp1

! mp_gatherv_cdp2: gather complex data from every processes to rank 0
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
         endif

! barrier until all processes reach here
         call mp_barrier(group)

! setup element count
         isize = size(send)

! invoke related mpi subroutines
         call MPI_GATHERV(send, isize, mpi_dcmplx, data, recv, disp, mpi_dcmplx, root, group, ierror)

! handler for return code
         call mp_error('mp_gatherv_cdp2', ierror)

         return
     end subroutine mp_gatherv_cdp2

! mp_gatherv_cdp3: gather complex data from every processes to rank 0
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
         endif

! barrier until all processes reach here
         call mp_barrier(group)

! setup element count
         isize = size(send)

! invoke related mpi subroutines
         call MPI_GATHERV(send, isize, mpi_dcmplx, data, recv, disp, mpi_dcmplx, root, group, ierror)

! handler for return code
         call mp_error('mp_gatherv_cdp3', ierror)

         return
     end subroutine mp_gatherv_cdp3

! mp_gatherv_cdp4: gather complex data from every processes to rank 0
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
         endif

! barrier until all processes reach here
         call mp_barrier(group)

! setup element count
         isize = size(send)

! invoke related mpi subroutines
         call MPI_GATHERV(send, isize, mpi_dcmplx, data, recv, disp, mpi_dcmplx, root, group, ierror)

! handler for return code
         call mp_error('mp_gatherv_cdp4', ierror)

         return
     end subroutine mp_gatherv_cdp4

! mp_gatherv_cdp5: gather complex data from every processes to rank 0
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
         endif

! barrier until all processes reach here
         call mp_barrier(group)

! setup element count
         isize = size(send)

! invoke related mpi subroutines
         call MPI_GATHERV(send, isize, mpi_dcmplx, data, recv, disp, mpi_dcmplx, root, group, ierror)

! handler for return code
         call mp_error('mp_gatherv_cdp5', ierror)

         return
     end subroutine mp_gatherv_cdp5

!-------------------------------------------------------------------------
!::: MPI collective operations: allgathering                           :::
!-------------------------------------------------------------------------

! mp_allgather_int1: gather integer data from all processes and then
! redistribute it to all processes
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
         endif

! barrier until all processes reach here
         call mp_barrier(group)

! setup element count
         isize = size(send)

! invoke related mpi subroutines
         call MPI_ALLGATHER(send, isize, mpi_mint, data, isize, mpi_mint, group, ierror)

! handler for return code
         call mp_error('mp_allgather_int1', ierror)

         return
     end subroutine mp_allgather_int1

! mp_allgather_int2: gather integer data from all processes and then
! redistribute it to all processes
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
         endif

! barrier until all processes reach here
         call mp_barrier(group)

! setup element count
         isize = size(send)

! invoke related mpi subroutines
         call MPI_ALLGATHER(send, isize, mpi_mint, data, isize, mpi_mint, group, ierror)

! handler for return code
         call mp_error('mp_allgather_int2', ierror)

         return
     end subroutine mp_allgather_int2

! mp_allgather_int3: gather integer data from all processes and then
! redistribute it to all processes
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
         endif

! barrier until all processes reach here
         call mp_barrier(group)

! setup element count
         isize = size(send)

! invoke related mpi subroutines
         call MPI_ALLGATHER(send, isize, mpi_mint, data, isize, mpi_mint, group, ierror)

! handler for return code
         call mp_error('mp_allgather_int3', ierror)

         return
     end subroutine mp_allgather_int3

! mp_allgather_int4: gather integer data from all processes and then
! redistribute it to all processes
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
         endif

! barrier until all processes reach here
         call mp_barrier(group)

! setup element count
         isize = size(send)

! invoke related mpi subroutines
         call MPI_ALLGATHER(send, isize, mpi_mint, data, isize, mpi_mint, group, ierror)

! handler for return code
         call mp_error('mp_allgather_int4', ierror)

         return
     end subroutine mp_allgather_int4

! mp_allgather_int5: gather integer data from all processes and then
! redistribute it to all processes
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
         endif

! barrier until all processes reach here
         call mp_barrier(group)

! setup element count
         isize = size(send)

! invoke related mpi subroutines
         call MPI_ALLGATHER(send, isize, mpi_mint, data, isize, mpi_mint, group, ierror)

! handler for return code
         call mp_error('mp_allgather_int5', ierror)

         return
     end subroutine mp_allgather_int5

! mp_allgather_rdp1: gather real data from all processes and then
! redistribute it to all processes
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
         endif

! barrier until all processes reach here
         call mp_barrier(group)

! setup element count
         isize = size(send)

! invoke related mpi subroutines
         call MPI_ALLGATHER(send, isize, mpi_dreal, data, isize, mpi_dreal, group, ierror)

! handler for return code
         call mp_error('mp_allgather_rdp1', ierror)

         return
     end subroutine mp_allgather_rdp1

! mp_allgather_rdp2: gather real data from all processes and then
! redistribute it to all processes
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
         endif

! barrier until all processes reach here
         call mp_barrier(group)

! setup element count
         isize = size(send)

! invoke related mpi subroutines
         call MPI_ALLGATHER(send, isize, mpi_dreal, data, isize, mpi_dreal, group, ierror)

! handler for return code
         call mp_error('mp_allgather_rdp2', ierror)

         return
     end subroutine mp_allgather_rdp2

! mp_allgather_rdp3: gather real data from all processes and then
! redistribute it to all processes
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
         endif

! barrier until all processes reach here
         call mp_barrier(group)

! setup element count
         isize = size(send)

! invoke related mpi subroutines
         call MPI_ALLGATHER(send, isize, mpi_dreal, data, isize, mpi_dreal, group, ierror)

! handler for return code
         call mp_error('mp_allgather_rdp3', ierror)

         return
     end subroutine mp_allgather_rdp3

! mp_allgather_rdp4: gather real data from all processes and then
! redistribute it to all processes
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
         endif

! barrier until all processes reach here
         call mp_barrier(group)

! setup element count
         isize = size(send)

! invoke related mpi subroutines
         call MPI_ALLGATHER(send, isize, mpi_dreal, data, isize, mpi_dreal, group, ierror)

! handler for return code
         call mp_error('mp_allgather_rdp4', ierror)

         return
     end subroutine mp_allgather_rdp4

! mp_allgather_rdp5: gather real data from all processes and then
! redistribute it to all processes
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
         endif

! barrier until all processes reach here
         call mp_barrier(group)

! setup element count
         isize = size(send)

! invoke related mpi subroutines
         call MPI_ALLGATHER(send, isize, mpi_dreal, data, isize, mpi_dreal, group, ierror)

! handler for return code
         call mp_error('mp_allgather_rdp5', ierror)

         return
     end subroutine mp_allgather_rdp5

! mp_allgather_cdp1: gather complex data from all processes and then
! redistribute it to all processes
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
         endif

! barrier until all processes reach here
         call mp_barrier(group)

! setup element count
         isize = size(send)

! invoke related mpi subroutines
         call MPI_ALLGATHER(send, isize, mpi_dcmplx, data, isize, mpi_dcmplx, group, ierror)

! handler for return code
         call mp_error('mp_allgather_cdp1', ierror)

         return
     end subroutine mp_allgather_cdp1

! mp_allgather_cdp2: gather complex data from all processes and then
! redistribute it to all processes
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
         endif

! barrier until all processes reach here
         call mp_barrier(group)

! setup element count
         isize = size(send)

! invoke related mpi subroutines
         call MPI_ALLGATHER(send, isize, mpi_dcmplx, data, isize, mpi_dcmplx, group, ierror)

! handler for return code
         call mp_error('mp_allgather_cdp2', ierror)

         return
     end subroutine mp_allgather_cdp2

! mp_allgather_cdp3: gather complex data from all processes and then
! redistribute it to all processes
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
         endif

! barrier until all processes reach here
         call mp_barrier(group)

! setup element count
         isize = size(send)

! invoke related mpi subroutines
         call MPI_ALLGATHER(send, isize, mpi_dcmplx, data, isize, mpi_dcmplx, group, ierror)

! handler for return code
         call mp_error('mp_allgather_cdp3', ierror)

         return
     end subroutine mp_allgather_cdp3

! mp_allgather_cdp4: gather complex data from all processes and then
! redistribute it to all processes
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
         endif

! barrier until all processes reach here
         call mp_barrier(group)

! setup element count
         isize = size(send)

! invoke related mpi subroutines
         call MPI_ALLGATHER(send, isize, mpi_dcmplx, data, isize, mpi_dcmplx, group, ierror)

! handler for return code
         call mp_error('mp_allgather_cdp4', ierror)

         return
     end subroutine mp_allgather_cdp4

! mp_allgather_cdp5: gather complex data from all processes and then
! redistribute it to all processes
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
         endif

! barrier until all processes reach here
         call mp_barrier(group)

! setup element count
         isize = size(send)

! invoke related mpi subroutines
         call MPI_ALLGATHER(send, isize, mpi_dcmplx, data, isize, mpi_dcmplx, group, ierror)

! handler for return code
         call mp_error('mp_allgather_cdp5', ierror)

         return
     end subroutine mp_allgather_cdp5

!-------------------------------------------------------------------------
!::: MPI collective operations: allgatherving                          :::
!-------------------------------------------------------------------------

! mp_allgatherv_int1: gather integer data from all processes and then
! redistribute it to all processes
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
         endif

! barrier until all processes reach here
         call mp_barrier(group)

! setup element count
         isize = size(send)

! invoke related mpi subroutines
         call MPI_ALLGATHERV(send, isize, mpi_mint, data, recv, disp, mpi_mint, group, ierror)

! handler for return code
         call mp_error('mp_allgatherv_int1', ierror)

         return
     end subroutine mp_allgatherv_int1

! mp_allgatherv_int2: gather integer data from all processes and then
! redistribute it to all processes
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
         endif

! barrier until all processes reach here
         call mp_barrier(group)

! setup element count
         isize = size(send)

! invoke related mpi subroutines
         call MPI_ALLGATHERV(send, isize, mpi_mint, data, recv, disp, mpi_mint, group, ierror)

! handler for return code
         call mp_error('mp_allgatherv_int2', ierror)

         return
     end subroutine mp_allgatherv_int2

! mp_allgatherv_int3: gather integer data from all processes and then
! redistribute it to all processes
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
         endif

! barrier until all processes reach here
         call mp_barrier(group)

! setup element count
         isize = size(send)

! invoke related mpi subroutines
         call MPI_ALLGATHERV(send, isize, mpi_mint, data, recv, disp, mpi_mint, group, ierror)

! handler for return code
         call mp_error('mp_allgatherv_int3', ierror)

         return
     end subroutine mp_allgatherv_int3

! mp_allgatherv_int4: gather integer data from all processes and then
! redistribute it to all processes
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
         endif

! barrier until all processes reach here
         call mp_barrier(group)

! setup element count
         isize = size(send)

! invoke related mpi subroutines
         call MPI_ALLGATHERV(send, isize, mpi_mint, data, recv, disp, mpi_mint, group, ierror)

! handler for return code
         call mp_error('mp_allgatherv_int4', ierror)

         return
     end subroutine mp_allgatherv_int4

! mp_allgatherv_int5: gather integer data from all processes and then
! redistribute it to all processes
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
         endif

! barrier until all processes reach here
         call mp_barrier(group)

! setup element count
         isize = size(send)

! invoke related mpi subroutines
         call MPI_ALLGATHERV(send, isize, mpi_mint, data, recv, disp, mpi_mint, group, ierror)

! handler for return code
         call mp_error('mp_allgatherv_int5', ierror)

         return
     end subroutine mp_allgatherv_int5

! mp_allgatherv_rdp1: gather real data from all processes and then
! redistribute it to all processes
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
         endif

! barrier until all processes reach here
         call mp_barrier(group)

! setup element count
         isize = size(send)

! invoke related mpi subroutines
         call MPI_ALLGATHERV(send, isize, mpi_dreal, data, recv, disp, mpi_dreal, group, ierror)

! handler for return code
         call mp_error('mp_allgatherv_rdp1', ierror)

         return
     end subroutine mp_allgatherv_rdp1

! mp_allgatherv_rdp2: gather real data from all processes and then
! redistribute it to all processes
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
         endif

! barrier until all processes reach here
         call mp_barrier(group)

! setup element count
         isize = size(send)

! invoke related mpi subroutines
         call MPI_ALLGATHERV(send, isize, mpi_dreal, data, recv, disp, mpi_dreal, group, ierror)

! handler for return code
         call mp_error('mp_allgatherv_rdp2', ierror)

         return
     end subroutine mp_allgatherv_rdp2

! mp_allgatherv_rdp3: gather real data from all processes and then
! redistribute it to all processes
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
         endif

! barrier until all processes reach here
         call mp_barrier(group)

! setup element count
         isize = size(send)

! invoke related mpi subroutines
         call MPI_ALLGATHERV(send, isize, mpi_dreal, data, recv, disp, mpi_dreal, group, ierror)

! handler for return code
         call mp_error('mp_allgatherv_rdp3', ierror)

         return
     end subroutine mp_allgatherv_rdp3

! mp_allgatherv_rdp4: gather real data from all processes and then
! redistribute it to all processes
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
         endif

! barrier until all processes reach here
         call mp_barrier(group)

! setup element count
         isize = size(send)

! invoke related mpi subroutines
         call MPI_ALLGATHERV(send, isize, mpi_dreal, data, recv, disp, mpi_dreal, group, ierror)

! handler for return code
         call mp_error('mp_allgatherv_rdp4', ierror)

         return
     end subroutine mp_allgatherv_rdp4

! mp_allgatherv_rdp5: gather real data from all processes and then
! redistribute it to all processes
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
         endif

! barrier until all processes reach here
         call mp_barrier(group)

! setup element count
         isize = size(send)

! invoke related mpi subroutines
         call MPI_ALLGATHERV(send, isize, mpi_dreal, data, recv, disp, mpi_dreal, group, ierror)

! handler for return code
         call mp_error('mp_allgatherv_rdp5', ierror)

         return
     end subroutine mp_allgatherv_rdp5

! mp_allgatherv_cdp1: gather complex data from all processes and then
! redistribute it to all processes
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
         endif

! barrier until all processes reach here
         call mp_barrier(group)

! setup element count
         isize = size(send)

! invoke related mpi subroutines
         call MPI_ALLGATHERV(send, isize, mpi_dcmplx, data, recv, disp, mpi_dcmplx, group, ierror)

! handler for return code
         call mp_error('mp_allgatherv_cdp1', ierror)

         return
     end subroutine mp_allgatherv_cdp1

! mp_allgatherv_cdp2: gather complex data from all processes and then
! redistribute it to all processes
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
         endif

! barrier until all processes reach here
         call mp_barrier(group)

! setup element count
         isize = size(send)

! invoke related mpi subroutines
         call MPI_ALLGATHERV(send, isize, mpi_dcmplx, data, recv, disp, mpi_dcmplx, group, ierror)

! handler for return code
         call mp_error('mp_allgatherv_cdp2', ierror)

         return
     end subroutine mp_allgatherv_cdp2

! mp_allgatherv_cdp3: gather complex data from all processes and then
! redistribute it to all processes
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
         endif

! barrier until all processes reach here
         call mp_barrier(group)

! setup element count
         isize = size(send)

! invoke related mpi subroutines
         call MPI_ALLGATHERV(send, isize, mpi_dcmplx, data, recv, disp, mpi_dcmplx, group, ierror)

! handler for return code
         call mp_error('mp_allgatherv_cdp3', ierror)

         return
     end subroutine mp_allgatherv_cdp3

! mp_allgatherv_cdp4: gather complex data from all processes and then
! redistribute it to all processes
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
         endif

! barrier until all processes reach here
         call mp_barrier(group)

! setup element count
         isize = size(send)

! invoke related mpi subroutines
         call MPI_ALLGATHERV(send, isize, mpi_dcmplx, data, recv, disp, mpi_dcmplx, group, ierror)

! handler for return code
         call mp_error('mp_allgatherv_cdp4', ierror)

         return
     end subroutine mp_allgatherv_cdp4

! mp_allgatherv_cdp5: gather complex data from all processes and then
! redistribute it to all processes
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
         endif

! barrier until all processes reach here
         call mp_barrier(group)

! setup element count
         isize = size(send)

! invoke related mpi subroutines
         call MPI_ALLGATHERV(send, isize, mpi_dcmplx, data, recv, disp, mpi_dcmplx, group, ierror)

! handler for return code
         call mp_error('mp_allgatherv_cdp5', ierror)

         return
     end subroutine mp_allgatherv_cdp5

!-------------------------------------------------------------------------
!::: MPI collective operations: reducing                               :::
!-------------------------------------------------------------------------

! mp_reduce_int: reduce 1 integer from all processes
     subroutine mp_reduce_int(source, data, root, gid)
         implicit none

! external arguments
         integer, intent(in) :: source
         integer, intent(inout) :: data
         integer, intent(in) :: root
         integer, optional, intent(in) :: gid

! set current communicator
         if ( present(gid) .eqv. .true. ) then
             group = gid
         else
             group = MPI_COMM_WORLD
         endif

! barrier until all processes reach here
         call mp_barrier(group)

! invoke related mpi subroutines
         call MPI_REDUCE(source, data, 1, mpi_mint, mpi_sum, root, group, ierror)

! handler for return code
         call mp_error('mp_reduce_int', ierror)

         return
     end subroutine mp_reduce_int

! mp_reduce_int1: reduce integer vector from all processes
     subroutine mp_reduce_int1(source, data, root, gid)
         implicit none

! external arguments
         integer, intent(in) :: source(:)
         integer, intent(inout) :: data(:)
         integer, intent(in) :: root
         integer, optional, intent(in) :: gid

! set current communicator
         if ( present(gid) .eqv. .true. ) then
             group = gid
         else
             group = MPI_COMM_WORLD
         endif

! barrier until all processes reach here
         call mp_barrier(group)

! setup element count
         isize = size(source)

! invoke related mpi subroutines
         call MPI_REDUCE(source, data, isize, mpi_mint, mpi_sum, root, group, ierror)

! handler for return code
         call mp_error('mp_reduce_int1', ierror)

         return
     end subroutine mp_reduce_int1

! mp_reduce_int2: reduce integer matrix from all processes
     subroutine mp_reduce_int2(source, data, root, gid)
         implicit none

! external arguments
         integer, intent(in) :: source(:,:)
         integer, intent(inout) :: data(:,:)
         integer, intent(in) :: root
         integer, optional, intent(in) :: gid

! set current communicator
         if ( present(gid) .eqv. .true. ) then
             group = gid
         else
             group = MPI_COMM_WORLD
         endif

! barrier until all processes reach here
         call mp_barrier(group)

! setup element count
         isize = size(source)

! invoke related mpi subroutines
         call MPI_REDUCE(source, data, isize, mpi_mint, mpi_sum, root, group, ierror)

! handler for return code
         call mp_error('mp_reduce_int2', ierror)

         return
     end subroutine mp_reduce_int2

! mp_reduce_int3: reduce integer matrix from all processes
     subroutine mp_reduce_int3(source, data, root, gid)
         implicit none

! external arguments
         integer, intent(in) :: source(:,:,:)
         integer, intent(inout) :: data(:,:,:)
         integer, intent(in) :: root
         integer, optional, intent(in) :: gid

! set current communicator
         if ( present(gid) .eqv. .true. ) then
             group = gid
         else
             group = MPI_COMM_WORLD
         endif

! barrier until all processes reach here
         call mp_barrier(group)

! setup element count
         isize = size(source)

! invoke related mpi subroutines
         call MPI_REDUCE(source, data, isize, mpi_mint, mpi_sum, root, group, ierror)

! handler for return code
         call mp_error('mp_reduce_int3', ierror)

         return
     end subroutine mp_reduce_int3

! mp_reduce_int4: reduce integer matrix from all processes
     subroutine mp_reduce_int4(source, data, root, gid)
         implicit none

! external arguments
         integer, intent(in) :: source(:,:,:,:)
         integer, intent(inout) :: data(:,:,:,:)
         integer, intent(in) :: root
         integer, optional, intent(in) :: gid

! set current communicator
         if ( present(gid) .eqv. .true. ) then
             group = gid
         else
             group = MPI_COMM_WORLD
         endif

! barrier until all processes reach here
         call mp_barrier(group)

! setup element count
         isize = size(source)

! invoke related mpi subroutines
         call MPI_REDUCE(source, data, isize, mpi_mint, mpi_sum, root, group, ierror)

! handler for return code
         call mp_error('mp_reduce_int4', ierror)

         return
     end subroutine mp_reduce_int4

! mp_reduce_int5: reduce integer matrix from all processes
     subroutine mp_reduce_int5(source, data, root, gid)
         implicit none

! external arguments
         integer, intent(in) :: source(:,:,:,:,:)
         integer, intent(inout) :: data(:,:,:,:,:)
         integer, intent(in) :: root
         integer, optional, intent(in) :: gid

! set current communicator
         if ( present(gid) .eqv. .true. ) then
             group = gid
         else
             group = MPI_COMM_WORLD
         endif

! barrier until all processes reach here
         call mp_barrier(group)

! setup element count
         isize = size(source)

! invoke related mpi subroutines
         call MPI_REDUCE(source, data, isize, mpi_mint, mpi_sum, root, group, ierror)

! handler for return code
         call mp_error('mp_reduce_int5', ierror)

         return
     end subroutine mp_reduce_int5

! mp_reduce_rdp: reduce 1 real from all processes
     subroutine mp_reduce_rdp(source, data, root, gid)
         implicit none

! external arguments
         real(dp), intent(in) :: source
         real(dp), intent(inout) :: data
         integer, intent(in) :: root
         integer, optional, intent(in) :: gid

! set current communicator
         if ( present(gid) .eqv. .true. ) then
             group = gid
         else
             group = MPI_COMM_WORLD
         endif

! barrier until all processes reach here
         call mp_barrier(group)

! invoke related mpi subroutines
         call MPI_REDUCE(source, data, 1, mpi_dreal, mpi_sum, root, group, ierror)

! handler for return code
         call mp_error('mp_reduce_rdp', ierror)

         return
     end subroutine mp_reduce_rdp

! mp_reduce_rdp1: reduce real vector from all processes
     subroutine mp_reduce_rdp1(source, data, root, gid)
         implicit none

! external arguments
         real(dp), intent(in) :: source(:)
         real(dp), intent(inout) :: data(:)
         integer, intent(in) :: root
         integer, optional, intent(in) :: gid

! set current communicator
         if ( present(gid) .eqv. .true. ) then
             group = gid
         else
             group = MPI_COMM_WORLD
         endif

! barrier until all processes reach here
         call mp_barrier(group)

! setup element count
         isize = size(source)

! invoke related mpi subroutines
         call MPI_REDUCE(source, data, isize, mpi_dreal, mpi_sum, root, group, ierror)

! handler for return code
         call mp_error('mp_reduce_rdp1', ierror)

         return
     end subroutine mp_reduce_rdp1

! mp_reduce_rdp2: reduce real matrix from all processes
     subroutine mp_reduce_rdp2(source, data, root, gid)
         implicit none

! external arguments
         real(dp), intent(in) :: source(:,:)
         real(dp), intent(inout) :: data(:,:)
         integer, intent(in) :: root
         integer, optional, intent(in) :: gid

! set current communicator
         if ( present(gid) .eqv. .true. ) then
             group = gid
         else
             group = MPI_COMM_WORLD
         endif

! barrier until all processes reach here
         call mp_barrier(group)

! setup element count
         isize = size(source)

! invoke related mpi subroutines
         call MPI_REDUCE(source, data, isize, mpi_dreal, mpi_sum, root, group, ierror)

! handler for return code
         call mp_error('mp_reduce_rdp2', ierror)

         return
     end subroutine mp_reduce_rdp2

! mp_reduce_rdp3: reduce real matrix from all processes
     subroutine mp_reduce_rdp3(source, data, root, gid)
         implicit none

! external arguments
         real(dp), intent(in) :: source(:,:,:)
         real(dp), intent(inout) :: data(:,:,:)
         integer, intent(in) :: root
         integer, optional, intent(in) :: gid

! set current communicator
         if ( present(gid) .eqv. .true. ) then
             group = gid
         else
             group = MPI_COMM_WORLD
         endif

! barrier until all processes reach here
         call mp_barrier(group)

! setup element count
         isize = size(source)

! invoke related mpi subroutines
         call MPI_REDUCE(source, data, isize, mpi_dreal, mpi_sum, root, group, ierror)

! handler for return code
         call mp_error('mp_reduce_rdp3', ierror)

         return
     end subroutine mp_reduce_rdp3

! mp_reduce_rdp4: reduce real matrix from all processes
     subroutine mp_reduce_rdp4(source, data, root, gid)
         implicit none

! external arguments
         real(dp), intent(in) :: source(:,:,:,:)
         real(dp), intent(inout) :: data(:,:,:,:)
         integer, intent(in) :: root
         integer, optional, intent(in) :: gid

! set current communicator
         if ( present(gid) .eqv. .true. ) then
             group = gid
         else
             group = MPI_COMM_WORLD
         endif

! barrier until all processes reach here
         call mp_barrier(group)

! setup element count
         isize = size(source)

! invoke related mpi subroutines
         call MPI_REDUCE(source, data, isize, mpi_dreal, mpi_sum, root, group, ierror)

! handler for return code
         call mp_error('mp_reduce_rdp4', ierror)

         return
     end subroutine mp_reduce_rdp4

! mp_reduce_rdp5: reduce real matrix from all processes
     subroutine mp_reduce_rdp5(source, data, root, gid)
         implicit none

! external arguments
         real(dp), intent(in) :: source(:,:,:,:,:)
         real(dp), intent(inout) :: data(:,:,:,:,:)
         integer, intent(in) :: root
         integer, optional, intent(in) :: gid

! set current communicator
         if ( present(gid) .eqv. .true. ) then
             group = gid
         else
             group = MPI_COMM_WORLD
         endif

! barrier until all processes reach here
         call mp_barrier(group)

! setup element count
         isize = size(source)

! invoke related mpi subroutines
         call MPI_REDUCE(source, data, isize, mpi_dreal, mpi_sum, root, group, ierror)

! handler for return code
         call mp_error('mp_reduce_rdp5', ierror)

         return
     end subroutine mp_reduce_rdp5

! mp_reduce_cdp: reduce 1 complex from all processes
     subroutine mp_reduce_cdp(source, data, root, gid)
         implicit none

! external arguments
         complex(dp), intent(in) :: source
         complex(dp), intent(inout) :: data
         integer, intent(in) :: root
         integer, optional, intent(in) :: gid

! set current communicator
         if ( present(gid) .eqv. .true. ) then
             group = gid
         else
             group = MPI_COMM_WORLD
         endif

! barrier until all processes reach here
         call mp_barrier(group)

! invoke related mpi subroutines
         call MPI_REDUCE(source, data, 1, mpi_dcmplx, mpi_sum, root, group, ierror)

! handler for return code
         call mp_error('mp_reduce_cdp', ierror)

         return
     end subroutine mp_reduce_cdp

! mp_reduce_cdp1: reduce complex vector from all processes
     subroutine mp_reduce_cdp1(source, data, root, gid)
         implicit none

! external arguments
         complex(dp), intent(in) :: source(:)
         complex(dp), intent(inout) :: data(:)
         integer, intent(in) :: root
         integer, optional, intent(in) :: gid

! set current communicator
         if ( present(gid) .eqv. .true. ) then
             group = gid
         else
             group = MPI_COMM_WORLD
         endif

! barrier until all processes reach here
         call mp_barrier(group)

! setup element count
         isize = size(source)

! invoke related mpi subroutines
         call MPI_REDUCE(source, data, isize, mpi_dcmplx, mpi_sum, root, group, ierror)

! handler for return code
         call mp_error('mp_reduce_cdp1', ierror)

         return
     end subroutine mp_reduce_cdp1

! mp_reduce_cdp2: reduce complex matrix from all processes
     subroutine mp_reduce_cdp2(source, data, root, gid)
         implicit none

! external arguments
         complex(dp), intent(in) :: source(:,:)
         complex(dp), intent(inout) :: data(:,:)
         integer, intent(in) :: root
         integer, optional, intent(in) :: gid

! set current communicator
         if ( present(gid) .eqv. .true. ) then
             group = gid
         else
             group = MPI_COMM_WORLD
         endif

! barrier until all processes reach here
         call mp_barrier(group)

! setup element count
         isize = size(source)

! invoke related mpi subroutines
         call MPI_REDUCE(source, data, isize, mpi_dcmplx, mpi_sum, root, group, ierror)

! handler for return code
         call mp_error('mp_reduce_cdp2', ierror)

         return
     end subroutine mp_reduce_cdp2

! mp_reduce_cdp3: reduce complex matrix from all processes
     subroutine mp_reduce_cdp3(source, data, root, gid)
         implicit none

! external arguments
         complex(dp), intent(in) :: source(:,:,:)
         complex(dp), intent(inout) :: data(:,:,:)
         integer, intent(in) :: root
         integer, optional, intent(in) :: gid

! set current communicator
         if ( present(gid) .eqv. .true. ) then
             group = gid
         else
             group = MPI_COMM_WORLD
         endif

! barrier until all processes reach here
         call mp_barrier(group)

! setup element count
         isize = size(source)

! invoke related mpi subroutines
         call MPI_REDUCE(source, data, isize, mpi_dcmplx, mpi_sum, root, group, ierror)

! handler for return code
         call mp_error('mp_reduce_cdp3', ierror)

         return
     end subroutine mp_reduce_cdp3

! mp_reduce_cdp4: reduce complex matrix from all processes
     subroutine mp_reduce_cdp4(source, data, root, gid)
         implicit none

! external arguments
         complex(dp), intent(in) :: source(:,:,:,:)
         complex(dp), intent(inout) :: data(:,:,:,:)
         integer, intent(in) :: root
         integer, optional, intent(in) :: gid

! set current communicator
         if ( present(gid) .eqv. .true. ) then
             group = gid
         else
             group = MPI_COMM_WORLD
         endif

! barrier until all processes reach here
         call mp_barrier(group)

! setup element count
         isize = size(source)

! invoke related mpi subroutines
         call MPI_REDUCE(source, data, isize, mpi_dcmplx, mpi_sum, root, group, ierror)

! handler for return code
         call mp_error('mp_reduce_cdp4', ierror)

         return
     end subroutine mp_reduce_cdp4

! mp_reduce_cdp5: reduce complex matrix from all processes
     subroutine mp_reduce_cdp5(source, data, root, gid)
         implicit none

! external arguments
         complex(dp), intent(in) :: source(:,:,:,:,:)
         complex(dp), intent(inout) :: data(:,:,:,:,:)
         integer, intent(in) :: root
         integer, optional, intent(in) :: gid

! set current communicator
         if ( present(gid) .eqv. .true. ) then
             group = gid
         else
             group = MPI_COMM_WORLD
         endif

! barrier until all processes reach here
         call mp_barrier(group)

! setup element count
         isize = size(source)

! invoke related mpi subroutines
         call MPI_REDUCE(source, data, isize, mpi_dcmplx, mpi_sum, root, group, ierror)

! handler for return code
         call mp_error('mp_reduce_cdp5', ierror)

         return
     end subroutine mp_reduce_cdp5

!-------------------------------------------------------------------------
!::: MPI collective operations: allreducing                            :::
!-------------------------------------------------------------------------

! mp_allreduce_int: reduce 1 integer from all processes
     subroutine mp_allreduce_int(source, data, gid)
         implicit none

! external arguments
         integer, intent(in) :: source
         integer, intent(inout) :: data
         integer, optional, intent(in) :: gid

! set current communicator
         if ( present(gid) .eqv. .true. ) then
             group = gid
         else
             group = MPI_COMM_WORLD
         endif

! barrier until all processes reach here
         call mp_barrier(group)

! invoke related mpi subroutines
         call MPI_ALLREDUCE(source, data, 1, mpi_mint, mpi_sum, group, ierror)

! handler for return code
         call mp_error('mp_allreduce_int', ierror)

         return
     end subroutine mp_allreduce_int

! mp_allreduce_int1: reduce integer vector from all processes
     subroutine mp_allreduce_int1(source, data, gid)
         implicit none

! external arguments
         integer, intent(in) :: source(:)
         integer, intent(inout) :: data(:)
         integer, optional, intent(in) :: gid

! set current communicator
         if ( present(gid) .eqv. .true. ) then
             group = gid
         else
             group = MPI_COMM_WORLD
         endif

! barrier until all processes reach here
         call mp_barrier(group)

! setup element count
         isize = size(source)

! invoke related mpi subroutines
         call MPI_ALLREDUCE(source, data, isize, mpi_mint, mpi_sum, group, ierror)

! handler for return code
         call mp_error('mp_allreduce_int1', ierror)

         return
     end subroutine mp_allreduce_int1

! mp_allreduce_int2: reduce integer matrix from all processes
     subroutine mp_allreduce_int2(source, data, gid)
         implicit none

! external arguments
         integer, intent(in) :: source(:,:)
         integer, intent(inout) :: data(:,:)
         integer, optional, intent(in) :: gid

! set current communicator
         if ( present(gid) .eqv. .true. ) then
             group = gid
         else
             group = MPI_COMM_WORLD
         endif

! barrier until all processes reach here
         call mp_barrier(group)

! setup element count
         isize = size(source)

! invoke related mpi subroutines
         call MPI_ALLREDUCE(source, data, isize, mpi_mint, mpi_sum, group, ierror)

! handler for return code
         call mp_error('mp_allreduce_int2', ierror)

         return
     end subroutine mp_allreduce_int2

! mp_allreduce_int3: reduce integer matrix from all processes
     subroutine mp_allreduce_int3(source, data, gid)
         implicit none

! external arguments
         integer, intent(in) :: source(:,:,:)
         integer, intent(inout) :: data(:,:,:)
         integer, optional, intent(in) :: gid

! set current communicator
         if ( present(gid) .eqv. .true. ) then
             group = gid
         else
             group = MPI_COMM_WORLD
         endif

! barrier until all processes reach here
         call mp_barrier(group)

! setup element count
         isize = size(source)

! invoke related mpi subroutines
         call MPI_ALLREDUCE(source, data, isize, mpi_mint, mpi_sum, group, ierror)

! handler for return code
         call mp_error('mp_allreduce_int3', ierror)

         return
     end subroutine mp_allreduce_int3

! mp_allreduce_int4: reduce integer matrix from all processes
     subroutine mp_allreduce_int4(source, data, gid)
         implicit none

! external arguments
         integer, intent(in) :: source(:,:,:,:)
         integer, intent(inout) :: data(:,:,:,:)
         integer, optional, intent(in) :: gid

! set current communicator
         if ( present(gid) .eqv. .true. ) then
             group = gid
         else
             group = MPI_COMM_WORLD
         endif

! barrier until all processes reach here
         call mp_barrier(group)

! setup element count
         isize = size(source)

! invoke related mpi subroutines
         call MPI_ALLREDUCE(source, data, isize, mpi_mint, mpi_sum, group, ierror)

! handler for return code
         call mp_error('mp_allreduce_int4', ierror)

         return
     end subroutine mp_allreduce_int4

! mp_allreduce_int5: reduce integer matrix from all processes
     subroutine mp_allreduce_int5(source, data, gid)
         implicit none

! external arguments
         integer, intent(in) :: source(:,:,:,:,:)
         integer, intent(inout) :: data(:,:,:,:,:)
         integer, optional, intent(in) :: gid

! set current communicator
         if ( present(gid) .eqv. .true. ) then
             group = gid
         else
             group = MPI_COMM_WORLD
         endif

! barrier until all processes reach here
         call mp_barrier(group)

! setup element count
         isize = size(source)

! invoke related mpi subroutines
         call MPI_ALLREDUCE(source, data, isize, mpi_mint, mpi_sum, group, ierror)

! handler for return code
         call mp_error('mp_allreduce_int5', ierror)

         return
     end subroutine mp_allreduce_int5

! mp_allreduce_rdp: reduce 1 real from all processes
     subroutine mp_allreduce_rdp(source, data, gid)
         implicit none

! external arguments
         real(dp), intent(in) :: source
         real(dp), intent(inout) :: data
         integer, optional, intent(in) :: gid

! set current communicator
         if ( present(gid) .eqv. .true. ) then
             group = gid
         else
             group = MPI_COMM_WORLD
         endif

! barrier until all processes reach here
         call mp_barrier(group)

! invoke related mpi subroutines
         call MPI_ALLREDUCE(source, data, 1, mpi_dreal, mpi_sum, group, ierror)

! handler for return code
         call mp_error('mp_allreduce_rdp', ierror)

         return
     end subroutine mp_allreduce_rdp

! mp_allreduce_rdp1: reduce real vector from all processes
     subroutine mp_allreduce_rdp1(source, data, gid)
         implicit none

! external arguments
         real(dp), intent(in) :: source(:)
         real(dp), intent(inout) :: data(:)
         integer, optional, intent(in) :: gid

! set current communicator
         if ( present(gid) .eqv. .true. ) then
             group = gid
         else
             group = MPI_COMM_WORLD
         endif

! barrier until all processes reach here
         call mp_barrier(group)

! setup element count
         isize = size(source)

! invoke related mpi subroutines
         call MPI_ALLREDUCE(source, data, isize, mpi_dreal, mpi_sum, group, ierror)

! handler for return code
         call mp_error('mp_allreduce_rdp1', ierror)

         return
     end subroutine mp_allreduce_rdp1

! mp_allreduce_rdp2: reduce real matrix from all processes
     subroutine mp_allreduce_rdp2(source, data, gid)
         implicit none

! external arguments
         real(dp), intent(in) :: source(:,:)
         real(dp), intent(inout) :: data(:,:)
         integer, optional, intent(in) :: gid

! set current communicator
         if ( present(gid) .eqv. .true. ) then
             group = gid
         else
             group = MPI_COMM_WORLD
         endif

! barrier until all processes reach here
         call mp_barrier(group)

! setup element count
         isize = size(source)

! invoke related mpi subroutines
         call MPI_ALLREDUCE(source, data, isize, mpi_dreal, mpi_sum, group, ierror)

! handler for return code
         call mp_error('mp_allreduce_rdp2', ierror)

         return
     end subroutine mp_allreduce_rdp2

! mp_allreduce_rdp3: reduce real matrix from all processes
     subroutine mp_allreduce_rdp3(source, data, gid)
         implicit none

! external arguments
         real(dp), intent(in) :: source(:,:,:)
         real(dp), intent(inout) :: data(:,:,:)
         integer, optional, intent(in) :: gid

! set current communicator
         if ( present(gid) .eqv. .true. ) then
             group = gid
         else
             group = MPI_COMM_WORLD
         endif

! barrier until all processes reach here
         call mp_barrier(group)

! setup element count
         isize = size(source)

! invoke related mpi subroutines
         call MPI_ALLREDUCE(source, data, isize, mpi_dreal, mpi_sum, group, ierror)

! handler for return code
         call mp_error('mp_allreduce_rdp3', ierror)

         return
     end subroutine mp_allreduce_rdp3

! mp_allreduce_rdp4: reduce real matrix from all processes
     subroutine mp_allreduce_rdp4(source, data, gid)
         implicit none

! external arguments
         real(dp), intent(in) :: source(:,:,:,:)
         real(dp), intent(inout) :: data(:,:,:,:)
         integer, optional, intent(in) :: gid

! set current communicator
         if ( present(gid) .eqv. .true. ) then
             group = gid
         else
             group = MPI_COMM_WORLD
         endif

! barrier until all processes reach here
         call mp_barrier(group)

! setup element count
         isize = size(source)

! invoke related mpi subroutines
         call MPI_ALLREDUCE(source, data, isize, mpi_dreal, mpi_sum, group, ierror)

! handler for return code
         call mp_error('mp_allreduce_rdp4', ierror)

         return
     end subroutine mp_allreduce_rdp4

! mp_allreduce_rdp5: reduce real matrix from all processes
     subroutine mp_allreduce_rdp5(source, data, gid)
         implicit none

! external arguments
         real(dp), intent(in) :: source(:,:,:,:,:)
         real(dp), intent(inout) :: data(:,:,:,:,:)
         integer, optional, intent(in) :: gid

! set current communicator
         if ( present(gid) .eqv. .true. ) then
             group = gid
         else
             group = MPI_COMM_WORLD
         endif

! barrier until all processes reach here
         call mp_barrier(group)

! setup element count
         isize = size(source)

! invoke related mpi subroutines
         call MPI_ALLREDUCE(source, data, isize, mpi_dreal, mpi_sum, group, ierror)

! handler for return code
         call mp_error('mp_allreduce_rdp5', ierror)

         return
     end subroutine mp_allreduce_rdp5

! mp_allreduce_cdp: reduce 1 complex from all processes
     subroutine mp_allreduce_cdp(source, data, gid)
         implicit none

! external arguments
         complex(dp), intent(in) :: source
         complex(dp), intent(inout) :: data
         integer, optional, intent(in) :: gid

! set current communicator
         if ( present(gid) .eqv. .true. ) then
             group = gid
         else
             group = MPI_COMM_WORLD
         endif

! barrier until all processes reach here
         call mp_barrier(group)

! invoke related mpi subroutines
         call MPI_ALLREDUCE(source, data, 1, mpi_dcmplx, mpi_sum, group, ierror)

! handler for return code
         call mp_error('mp_allreduce_cdp', ierror)

         return
     end subroutine mp_allreduce_cdp

! mp_allreduce_cdp1: reduce complex vector from all processes
     subroutine mp_allreduce_cdp1(source, data, gid)
         implicit none

! external arguments
         complex(dp), intent(in) :: source(:)
         complex(dp), intent(inout) :: data(:)
         integer, optional, intent(in) :: gid

! set current communicator
         if ( present(gid) .eqv. .true. ) then
             group = gid
         else
             group = MPI_COMM_WORLD
         endif

! barrier until all processes reach here
         call mp_barrier(group)

! setup element count
         isize = size(source)

! invoke related mpi subroutines
         call MPI_ALLREDUCE(source, data, isize, mpi_dcmplx, mpi_sum, group, ierror)

! handler for return code
         call mp_error('mp_allreduce_cdp1', ierror)

         return
     end subroutine mp_allreduce_cdp1

! mp_allreduce_cdp2: reduce complex matrix from all processes
     subroutine mp_allreduce_cdp2(source, data, gid)
         implicit none

! external arguments
         complex(dp), intent(in) :: source(:,:)
         complex(dp), intent(inout) :: data(:,:)
         integer, optional, intent(in) :: gid

! set current communicator
         if ( present(gid) .eqv. .true. ) then
             group = gid
         else
             group = MPI_COMM_WORLD
         endif

! barrier until all processes reach here
         call mp_barrier(group)

! setup element count
         isize = size(source)

! invoke related mpi subroutines
         call MPI_ALLREDUCE(source, data, isize, mpi_dcmplx, mpi_sum, group, ierror)

! handler for return code
         call mp_error('mp_allreduce_cdp2', ierror)

         return
     end subroutine mp_allreduce_cdp2

! mp_allreduce_cdp3: reduce complex matrix from all processes
     subroutine mp_allreduce_cdp3(source, data, gid)
         implicit none

! external arguments
         complex(dp), intent(in) :: source(:,:,:)
         complex(dp), intent(inout) :: data(:,:,:)
         integer, optional, intent(in) :: gid

! set current communicator
         if ( present(gid) .eqv. .true. ) then
             group = gid
         else
             group = MPI_COMM_WORLD
         endif

! barrier until all processes reach here
         call mp_barrier(group)

! setup element count
         isize = size(source)

! invoke related mpi subroutines
         call MPI_ALLREDUCE(source, data, isize, mpi_dcmplx, mpi_sum, group, ierror)

! handler for return code
         call mp_error('mp_allreduce_cdp3', ierror)

         return
     end subroutine mp_allreduce_cdp3

! mp_allreduce_cdp4: reduce complex matrix from all processes
     subroutine mp_allreduce_cdp4(source, data, gid)
         implicit none

! external arguments
         complex(dp), intent(in) :: source(:,:,:,:)
         complex(dp), intent(inout) :: data(:,:,:,:)
         integer, optional, intent(in) :: gid

! set current communicator
         if ( present(gid) .eqv. .true. ) then
             group = gid
         else
             group = MPI_COMM_WORLD
         endif

! barrier until all processes reach here
         call mp_barrier(group)

! setup element count
         isize = size(source)

! invoke related mpi subroutines
         call MPI_ALLREDUCE(source, data, isize, mpi_dcmplx, mpi_sum, group, ierror)

! handler for return code
         call mp_error('mp_allreduce_cdp4', ierror)

         return
     end subroutine mp_allreduce_cdp4

! mp_allreduce_cdp5: reduce complex matrix from all processes
     subroutine mp_allreduce_cdp5(source, data, gid)
         implicit none

! external arguments
         complex(dp), intent(in) :: source(:,:,:,:,:)
         complex(dp), intent(inout) :: data(:,:,:,:,:)
         integer, optional, intent(in) :: gid

! set current communicator
         if ( present(gid) .eqv. .true. ) then
             group = gid
         else
             group = MPI_COMM_WORLD
         endif

! barrier until all processes reach here
         call mp_barrier(group)

! setup element count
         isize = size(source)

! invoke related mpi subroutines
         call MPI_ALLREDUCE(source, data, isize, mpi_dcmplx, mpi_sum, group, ierror)

! handler for return code
         call mp_error('mp_allreduce_cdp5', ierror)

         return
     end subroutine mp_allreduce_cdp5

!-------------------------------------------------------------------------
!::: MPI handler for return code                                       :::
!-------------------------------------------------------------------------

! mp_error: deal with the return code of MPI subroutine
     subroutine mp_error(sub, err)
         implicit none

! external arguments
! subroutine name
         character(len=*), intent(in) :: sub

! error no
         integer, intent(in) :: err

         select case (err)

             case (MPI_SUCCESS)
                 return

             case (MPI_ERR_COMM)
# define STR_ERR_COMM      'invalid communicator in mpi call.'
                 write(mystd,'(2a)') sub, STR_ERR_COMM

             case (MPI_ERR_COUNT)
# define STR_ERR_COUNT     'invalid count in mpi call.'
                 write(mystd,'(2a)') sub, STR_ERR_COUNT

             case (MPI_ERR_TYPE)
# define STR_ERR_TYPE      'invalid datatype in mpi call.'
                 write(mystd,'(2a)') sub, STR_ERR_TYPE

             case (MPI_ERR_BUFFER)
# define STR_ERR_BUFFER    'invalid buffer in mpi call.'
                 write(mystd,'(2a)') sub, STR_ERR_BUFFER

             case (MPI_ERR_ROOT)
# define STR_ERR_ROOT      'invalid root in mpi call.'
                 write(mystd,'(2a)') sub, STR_ERR_ROOT

             case (MPI_ERR_ARG)
# define STR_ERR_ARG       'invalid argument in mpi call.'
                 write(mystd,'(2a)') sub, STR_ERR_ARG

             case (MPI_ERR_TAG)
# define STR_ERR_TAG       'invalid tag in mpi call.'
                 write(mystd,'(2a)') sub, STR_ERR_TAG

             case (MPI_ERR_RANK)
# define STR_ERR_RANK      'invalid rank in mpi call.'
                 write(mystd,'(2a)') sub, STR_ERR_RANK

             case (MPI_ERR_GROUP)
# define STR_ERR_GROUP     'null group passed to mpi call.'
                 write(mystd,'(2a)') sub, STR_ERR_GROUP

             case (MPI_ERR_OP)
# define STR_ERR_OP        'invalid operation in mpi call.'
                 write(mystd,'(2a)') sub, STR_ERR_OP

             case (MPI_ERR_TOPOLOGY)
# define STR_ERR_TOPOLOGY  'invalid topology in mpi call.'
                 write(mystd,'(2a)') sub, STR_ERR_TOPOLOGY

             case (MPI_ERR_DIMS)
# define STR_ERR_DIMS      'illegal dimension argument in mpi call.'
                 write(mystd,'(2a)') sub, STR_ERR_DIMS

             case (MPI_ERR_UNKNOWN)
# define STR_ERR_UNKNOWN   'unknown error in mpi call.'
                 write(mystd,'(2a)') sub, STR_ERR_UNKNOWN

             case (MPI_ERR_TRUNCATE)
# define STR_ERR_TRUNCATE  'message truncated on receive in mpi call.'
                 write(mystd,'(2a)') sub, STR_ERR_TRUNCATE

             case (MPI_ERR_OTHER)
# define STR_ERR_OTHER     'other error in mpi call.'
                 write(mystd,'(2a)') sub, STR_ERR_OTHER

             case (MPI_ERR_INTERN)
# define STR_ERR_INTERN    'internal error code in mpi call.'
                 write(mystd,'(2a)') sub, STR_ERR_INTERN

             case (MPI_ERR_IN_STATUS)
# define STR_ERR_IN_STATUS 'look in status for error value.'
                 write(mystd,'(2a)') sub, STR_ERR_IN_STATUS

             case (MPI_ERR_PENDING)
# define STR_ERR_PENDING   'pending request in mpi call.'
                 write(mystd,'(2a)') sub, STR_ERR_PENDING

             case (MPI_ERR_REQUEST)
# define STR_ERR_REQUEST   'illegal mpi_request handle in mpi call.'
                 write(mystd,'(2a)') sub, STR_ERR_REQUEST

             case (MPI_ERR_LASTCODE)
# define STR_ERR_LASTCODE  'last error code in mpi call.'
                 write(mystd,'(2a)') sub, STR_ERR_LASTCODE

             case default
                 return

         end select

         return
     end subroutine mp_error

  end module mmpi

!>>> current used compiler is not mpif90
# else   /* MPI */

  module mmpi
     implicit none

!-------------------------------------------------------------------------
!::: declare global constants                                          :::
!-------------------------------------------------------------------------

! dp: number precision, double precision for reals
     integer, private, parameter :: dp = kind(1.0d0)

! mystd: device descriptor, console output
     integer, private, parameter :: mystd = 6

!-------------------------------------------------------------------------
!::: declare accessibility for module routines                         :::
!-------------------------------------------------------------------------

!>>> mpi information operation
     public :: mp_info

  contains

!-------------------------------------------------------------------------
!::: MPI information operations                                        :::
!-------------------------------------------------------------------------

! mp_info: return the current information about mpi environment
     subroutine mp_info()
         implicit none

         write(mystd,'(a)') 'Bad news, your compiler is mpi-incompatible.'

         return
     end subroutine mp_info

  end module mmpi

# endif  /* MPI */
