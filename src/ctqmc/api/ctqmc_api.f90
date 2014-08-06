!!!-----------------------------------------------------------------------
!!! project : lilac
!!! program : api
!!! source  : ctqmc_api.f90
!!! type    : module
!!! author  : li huang (email:huangli712@gmail.com)
!!! history : 01/07/2014 by li huang
!!!           01/11/2014 by li huang
!!!           01/13/2014 by li huang
!!! purpose : the purpose of this module is to define a generic and robust
!!!           application programming interface (CAPI) for continuous-time
!!!           quantum Monte Carlo impurity solver
!!! status  : unstable
!!! comment :
!!!-----------------------------------------------------------------------

!!
!!
!! Introduction
!! ============
!!
!!
!!
!!
!!
!!
!!
!!
!!
!!
!!
!!

  module api
     implicit none

!!========================================================================
!!>>> declare global parameters                                        <<<
!!========================================================================

! dp: number precision, double precision for reals
     integer, private, parameter :: dp = kind(1.0d0)

!!========================================================================
!!>>> declare global data structure                                    <<<
!!========================================================================

! define type T_mpi, which is used to describe the mpi environment
     type :: T_mpi
         integer :: nprocs
         integer :: myid
         integer :: master
         integer :: cid
         integer :: cx
         integer :: cy
     end type T_mpi

! define type T_solver, which is used to describe the generic abstract
! ctqmc impurity solver
! note: it can not be used directly
     type :: T_solver
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
     end type T_solver

! define type T_segment_solver, which is used to describe the ctqmc
! impurity solver which based on segment representation
! note: it can not be used directly
     type, extends (T_solver) :: T_segment_solver
         character(len=10) :: solver_type = 'SEGMENT'
     end type T_segment_solver

! define type T_general_solver, which is used to describe the ctqmc
! impurity solver which based on general matrix formulation 
! note: it can not be used directly
     type, extends (T_solver) :: T_general_solver
         character(len=10) :: solver_type = 'GENERAL'
     end type T_general_solver

! define type T_segment_azalea, which is used to describe the ctqmc
! impurity solver code azalea
     type, extends (T_segment_solver) :: T_segment_azalea
         character(len=10) :: solver_name = 'AZALEA'
         integer :: solver_id = 101
         integer :: solver_ready = 1
     end type T_segment_azalea

! define type T_segment_gardenia, which is used to describe the ctqmc
! impurity solver code gardenia
     type, extends (T_segment_solver) :: T_segment_gardenia
         character(len=10) :: solver_name = 'GARDENIA'
         integer :: solver_id = 102
         integer :: solver_ready = 1

         integer :: isort
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
     type, extends (T_segment_solver) :: T_segment_narcissus
         character(len=10) :: solver_name = 'NARCISSUS'
         integer :: solver_id = 103
         integer :: solver_ready = 1

         integer :: isort
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
     type, extends (T_general_solver) :: T_general_begonia
         character(len=10) :: solver_name = 'BEGONIA'
         integer :: solver_id = 201
         integer :: solver_ready = 1

         integer :: nzero
         integer :: npart
     end type T_general_begonia

! define type T_general_lavender, which is used to describe the ctqmc
! impurity solver code lavender
     type, extends (T_general_solver) :: T_general_lavender
         character(len=10) :: solver_name = 'LAVENDER'
         integer :: solver_id = 202
         integer :: solver_ready = 1

         integer :: isort
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

! define type T_general_pansy, which is used to describe the ctqmc
! impurity solver code pansy 
! TODO
     type, extends (T_general_solver) :: T_general_pansy
         character(len=10) :: solver_name = 'PANSY'
         integer :: solver_id = 301
         integer :: solver_ready = 1
     end type T_general_pansy

! define type T_general_manjushaka, which is used to describe the ctqmc
! impurity solver code manjushaka
! TODO
     type, extends (T_general_solver) :: T_general_manjushaka
         character(len=10) :: solver_name = 'MANJUSHAKA'
         integer :: solver_id = 302
         integer :: solver_ready = 1
     end type T_general_manjushaka

!-------------------------------------------------------------------------
!::: declare accessibility for module routines                         :::
!-------------------------------------------------------------------------

     private :: T_solver
     private :: T_segment_solver
     private :: T_general_solver

     public  :: T_segment_azalea
     public  :: T_segment_gardenia
     public  :: T_segment_narcissus

     public  :: T_general_begonia
     public  :: T_general_lavender

     public  :: T_mpi

     public  :: init_ctqmc
     public  :: exec_ctqmc
     public  :: stop_ctqmc

     public  :: set_hybf
     public  :: set_symm
     public  :: set_eimp

     public  :: get_grnf
     public  :: get_sigf

  contains ! encapsulated functionality

!>>> initialize the ctqmc quantum impurity solver
     subroutine init_ctqmc(I_mpi, I_solver)
         implicit none

! type structure of mpi
         class(*) :: I_mpi

! type structure of generic solver
         class(*) :: I_solver

         call cat_init_ctqmc(I_mpi, I_solver)

         return
     end subroutine init_ctqmc

!>>> execute the ctqmc quantum impurity solver
     subroutine exec_ctqmc(iter)
         implicit none

! current iteration number
         integer :: iter

         call cat_exec_ctqmc(iter)

         return
     end subroutine exec_ctqmc

!>>> stop the ctqmc quantum impurity solver
     subroutine stop_ctqmc()
         implicit none

         call cat_stop_ctqmc()

         return
     end subroutine stop_ctqmc

!>>> setup the hybridization function
     subroutine set_hybf(size_t, hybf_t)
         implicit none

! size of hybf
         integer :: size_t

! hybridization function
         complex(dp) :: hybf_t(size_t)

         call cat_set_hybf(size_t, hybf_t)

         return
     end subroutine set_hybf

!>>> setup the symmetry vector
     subroutine set_symm(size_t, symm_t)
         implicit none

! size of symm
         integer :: size_t

! symmetry vector
         integer :: symm_t(size_t)

         call cat_set_symm(size_t, symm_t)

         return
     end subroutine set_symm

!>>> setup the impurity level
     subroutine set_eimp(size_t, eimp_t)
         implicit none

! size of eimp
         integer :: size_t

! impurity level
         real(dp) :: eimp_t(size_t)

         call cat_set_eimp(size_t, eimp_t)

         return
     end subroutine set_eimp

!>>> extract the impurity green's function
     subroutine get_grnf(size_t, grnf_t)
         implicit none

! size of grnf
         integer :: size_t

! impurity green's function
         complex(dp) :: grnf_t(size_t)

         call cat_get_grnf(size_t, grnf_t)

         return
     end subroutine get_grnf

!>>> extract the self-energy function
     subroutine get_sigf(size_t, sigf_t)
         implicit none

! size of sigf
         integer :: size_t

! self-energy function
         complex(dp) :: sigf_t(size_t)

         call cat_get_sigf(size_t, sigf_t)

         return
     end subroutine get_sigf

  end module api
