!!!-----------------------------------------------------------------------
!!! project : lilac
!!! program : dapi
!!! source  : hfqmc_api.f90
!!! type    : module
!!! author  : li huang (email:huangli712@gmail.com)
!!! history : 12/06/2014 by li huang
!!!           12/08/2014 by li huang
!!! purpose : the purpose of this module is to define a generic and robust
!!!           application programming interface (API) for the Hirsch-Fye
!!!           quantum Monte Carlo impurity solver
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
     integer, public, parameter :: solver_id_daisy     = 901

! solver status, 1 means ready, 0 means not ready
     integer, public, parameter :: solver_is_ready_daisy = 1

!!========================================================================
!!>>> declare global data structure                                    <<<
!!========================================================================

! note: now f2py does not support derived types, so we have to comment
! out them when f2py is used.

# if !defined (F2PY)

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

# endif  /* F2PY */

!!========================================================================
!!>>> declare accessibility for module routines                        <<<
!!========================================================================

     public :: solver_id
     public :: solver_status

     public :: init_hfqmc
     public :: exec_hfqmc
     public :: stop_hfqmc

     public :: set_wssf
     public :: set_symm
     public :: set_eimp

     public :: get_grnf
     public :: get_sigf
     public :: get_nmat
     public :: get_nnmat

  contains ! encapsulated functionality

!!>>> solver_id: return the solver identity
  subroutine solver_id(I_solver_id)
     implicit none

! external arguments
! solver identity
     integer, intent(out) :: I_solver_id

! declare f2py directives
!F2PY intent(out) I_solver_id
     call cat_solver_id(I_solver_id)

     return
  end subroutine solver_id

!!>>> solver_status: return the solver status
  subroutine solver_status(I_solver_status)
     implicit none

! external arguments
! solver status
     integer, intent(out) :: I_solver_status

! declare f2py directives
!F2PY intent(out) I_solver_status
     call cat_solver_status(I_solver_status)

     return
  end subroutine solver_status

  subroutine init_hfqmc()
     implicit none

     return
  end subroutine init_hfqmc

  subroutine exec_hfqmc()
     implicit none

     return
  end subroutine exec_hfqmc

  subroutine stop_hfqmc()
     implicit none

     return
  end subroutine stop_hfqmc

  end module dapi
