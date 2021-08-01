### Explicit Fortran types

We define some Fortran modules/types, which are used to passed parameters to the quantum impurity solvers and the other components.

**> JASMINE** component (atomic eigenvalue problem solver)

```fortran
  module japi
     implicit none

!!========================================================================
!!>>> declare private parameters                                       <<<
!!========================================================================

! dp: number precision, double precision for reals
     integer, private, parameter :: dp = kind(1.0d0)

!!========================================================================
!!>>> declare global data structure                                    <<<
!!========================================================================

! define type T_jasmine, which is used to describe the control parameters
! for the jasmine code
     public :: T_jasmine
     type :: T_jasmine
         integer :: ibasis
         integer :: ictqmc
         integer :: icu
         integer :: icf
         integer :: isoc

         integer :: nband
         integer :: nspin
         integer :: norbs
         integer :: ncfgs

         integer :: nmini
         integer :: nmaxi

         real(dp) :: Uc
         real(dp) :: Uv
         real(dp) :: Jz
         real(dp) :: Js
         real(dp) :: Jp
         real(dp) :: Ud
         real(dp) :: Jh
         real(dp) :: mune
         real(dp) :: lambda
     end type T_jasmine

  end module japi
```

---

**> CT-HYB** quantum impurity solvers

```fortran
  module capi
     implicit none

!!========================================================================
!!>>> declare private parameters                                       <<<
!!========================================================================

! dp: number precision, double precision for reals
     integer, private, parameter :: dp = kind(1.0d0)

!!========================================================================
!!>>> declare global constants                                         <<<
!!========================================================================

! solver identity
     integer, public, parameter :: solver_id_azalea         = 101
     integer, public, parameter :: solver_id_gardenia       = 102
     integer, public, parameter :: solver_id_narcissus      = 103
     integer, public, parameter :: solver_id_begonia        = 201
     integer, public, parameter :: solver_id_lavender       = 202
     integer, public, parameter :: solver_id_camellia       = 301
     integer, public, parameter :: solver_id_pansy          = 401
     integer, public, parameter :: solver_id_manjushaka     = 402

! solver status, 1 means ready, 0 means not ready
     integer, public, parameter :: solver_is_ready_azalea     = 1
     integer, public, parameter :: solver_is_ready_gardenia   = 1
     integer, public, parameter :: solver_is_ready_narcissus  = 1
     integer, public, parameter :: solver_is_ready_begonia    = 1
     integer, public, parameter :: solver_is_ready_lavender   = 1
     integer, public, parameter :: solver_is_ready_camellia   = 1
     integer, public, parameter :: solver_is_ready_pansy      = 1
     integer, public, parameter :: solver_is_ready_manjushaka = 1

!!========================================================================
!!>>> declare global data structure                                    <<<
!!========================================================================

! define type T_mpi, which is used to describe the mpi environment
     public :: T_mpi
     type :: T_mpi
         integer :: nprocs
         integer :: myid
         integer :: master
         integer :: cid
         integer :: cx
         integer :: cy
     end type T_mpi

! define type T_generic_solver, which is used to describe the generic
! abstract ctqmc impurity solver
! note: it can not be used directly
     private :: T_generic_solver
     type :: T_generic_solver
         integer :: isscf
         integer :: issun
         integer :: isspn
         integer :: isbin
         integer :: nband
         integer :: nspin
         integer :: norbs
         integer :: ncfgs
         integer :: niter
         integer :: mkink
         integer :: mfreq
         integer :: nfreq
         integer :: ntime
         integer :: nflip
         integer :: ntherm
         integer :: nsweep
         integer :: nwrite
         integer :: nclean
         integer :: nmonte
         integer :: ncarlo

         real(dp) :: U
         real(dp) :: Uc
         real(dp) :: Uv
         real(dp) :: Jz
         real(dp) :: Js
         real(dp) :: Jp
         real(dp) :: mune
         real(dp) :: beta
         real(dp) :: part
         real(dp) :: alpha
     end type T_generic_solver

! define type T_segment_solver, which is used to describe the ctqmc
! impurity solver which based on segment representation
! note: it can not be used directly
     private :: T_segment_solver
     type, extends (T_generic_solver) :: T_segment_solver
         character(len=10) :: solver_type = 'SEGMENT'
     end type T_segment_solver

! define type T_general_solver, which is used to describe the ctqmc
! impurity solver which based on general matrix formulation
! note: it can not be used directly
     private :: T_general_solver
     type, extends (T_generic_solver) :: T_general_solver
         character(len=10) :: solver_type = 'GENERAL'
     end type T_general_solver

! define type T_segment_azalea, which is used to describe the ctqmc
! impurity solver code azalea
     public :: T_segment_azalea
     type, extends (T_segment_solver) :: T_segment_azalea
         character(len=10) :: solver_name = 'AZALEA'
     end type T_segment_azalea

! define type T_segment_gardenia, which is used to describe the ctqmc
! impurity solver code gardenia
     public :: T_segment_gardenia
     type, extends (T_segment_solver) :: T_segment_gardenia
         character(len=10) :: solver_name = 'GARDENIA'

         integer :: isort
         integer :: issus
         integer :: isvrt
         integer :: lemax
         integer :: legrd
         integer :: chmax
         integer :: chgrd
         integer :: nffrq
         integer :: nbfrq
     end type T_segment_gardenia

! define type T_segment_narcissus, which is used to describe the ctqmc
! impurity solver code narcissus
     public :: T_segment_narcissus
     type, extends (T_segment_solver) :: T_segment_narcissus
         character(len=10) :: solver_name = 'NARCISSUS'

         integer :: isort
         integer :: issus
         integer :: isvrt
         integer :: isscr
         integer :: lemax
         integer :: legrd
         integer :: chmax
         integer :: chgrd
         integer :: nffrq
         integer :: nbfrq

         real(dp) :: lc
         real(dp) :: wc
     end type T_segment_narcissus

! define type T_general_begonia, which is used to describe the ctqmc
! impurity solver code begonia
     public :: T_general_begonia
     type, extends (T_general_solver) :: T_general_begonia
         character(len=10) :: solver_name = 'BEGONIA'

         integer :: nzero
         integer :: npart
     end type T_general_begonia

! define type T_general_lavender, which is used to describe the ctqmc
! impurity solver code lavender
     public :: T_general_lavender
     type, extends (T_general_solver) :: T_general_lavender
         character(len=10) :: solver_name = 'LAVENDER'

         integer :: isort
         integer :: issus
         integer :: isvrt
         integer :: nzero
         integer :: lemax
         integer :: legrd
         integer :: chmax
         integer :: chgrd
         integer :: nffrq
         integer :: nbfrq
         integer :: npart
     end type T_general_lavender

! define type T_general_camellia, which is used to describe the ctqmc
! impurity solver code camellia
     public :: T_general_camellia
     type, extends (T_general_solver) :: T_general_camellia
         character(len=10) :: solver_name = 'CAMELLIA'

         integer :: isort
         integer :: issus
         integer :: isvrt
         integer :: nzero
         integer :: lemax
         integer :: legrd
         integer :: chmax
         integer :: chgrd
         integer :: nffrq
         integer :: nbfrq
         integer :: nvect
         integer :: nleja
     end type T_general_camellia

! define type T_general_pansy, which is used to describe the ctqmc
! impurity solver code pansy
     public :: T_general_pansy
     type, extends (T_general_solver) :: T_general_pansy
         character(len=10) :: solver_name = 'PANSY'

         integer :: npart
     end type T_general_pansy

! define type T_general_manjushaka, which is used to describe the ctqmc
! impurity solver code manjushaka
     public :: T_general_manjushaka
     type, extends (T_general_solver) :: T_general_manjushaka
         character(len=10) :: solver_name = 'MANJUSHAKA'

         integer :: isort
         integer :: issus
         integer :: isvrt
         integer :: ifast
         integer :: itrun
         integer :: lemax
         integer :: legrd
         integer :: chmax
         integer :: chgrd
         integer :: nffrq
         integer :: nbfrq
         integer :: npart
     end type T_general_manjushaka

  end module capi
```

---

**> HF-QMC** quantum impurity solvers

```fortran
  module dapi
     implicit none

!!========================================================================
!!>>> declare private parameters                                       <<<
!!========================================================================

! dp: number precision, double precision for reals
     integer, private, parameter :: dp = kind(1.0d0)

!!========================================================================
!!>>> declare global constants                                         <<<
!!========================================================================

! solver identity
     integer, public, parameter :: solver_id_daisy = 901

! solver status, 1 means ready, 0 means not ready
     integer, public, parameter :: solver_is_ready_daisy = 1

!!========================================================================
!!>>> declare global data structure                                    <<<
!!========================================================================

! define type T_mpi, which is used to describe the mpi environment
     public :: T_mpi
     type :: T_mpi
         integer :: nprocs
         integer :: myid
         integer :: master
         integer :: cid
         integer :: cx
         integer :: cy
     end type T_mpi

! define type T_daisy, which is used to describe the control parameters
! for the daisy code
     public :: T_daisy
     type :: T_daisy
         integer :: isscf
         integer :: issun
         integer :: isspn
         integer :: isbin
         integer :: nband
         integer :: nspin
         integer :: norbs
         integer :: niter
         integer :: mstep
         integer :: mfreq
         integer :: nsing
         integer :: ntime
         integer :: ntherm
         integer :: nsweep
         integer :: nclean
         integer :: ncarlo

         real(dp) :: Uc
         real(dp) :: Jz
         real(dp) :: mune
         real(dp) :: beta
         real(dp) :: part
         real(dp) :: alpha
     end type T_daisy

  end module dapi
```