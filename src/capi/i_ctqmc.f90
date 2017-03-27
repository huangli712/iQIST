!!!-----------------------------------------------------------------------
!!! project : CAPI (Common Application Programming Interface)
!!! program : capi
!!! source  : i_ctqmc.f90
!!! type    : module
!!! author  : li huang (email:lihuang.dmft@gmail.com)
!!! history : 01/07/2014 by li huang (created)
!!!           03/27/2017 by li huang (last modified)
!!! purpose : the purpose of this module is to define a generic and robust
!!!           application programming interface (API) for continuous-time
!!!           quantum Monte Carlo impurity solver.
!!! status  : unstable
!!! comment :
!!!-----------------------------------------------------------------------

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
     integer, public, parameter :: solver_id_gardenia       = 101
     integer, public, parameter :: solver_id_narcissus      = 102
     integer, public, parameter :: solver_id_manjushaka     = 201

! solver status, 1 means ready, 0 means not ready
     integer, public, parameter :: solver_is_ready_gardenia   = 1
     integer, public, parameter :: solver_is_ready_narcissus  = 1
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
