!!!-----------------------------------------------------------------------
!!! project : CAPI (Common Application Programming Interface)
!!! program : dapi
!!! source  : hfqmc_api.f90
!!! type    : module
!!! author  : li huang (email:lihuang.dmft@gmail.com)
!!! history : 12/06/2014 by li huang (created)
!!!           08/17/2015 by li huang (last modified)
!!! purpose : the purpose of this module is to define a generic and robust
!!!           application programming interface (API) for the Hirsch-Fye
!!!           quantum Monte Carlo impurity solver.
!!! status  : unstable
!!! comment :
!!!-----------------------------------------------------------------------

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
