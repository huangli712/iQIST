!!!-----------------------------------------------------------------------
!!! project : lilac
!!! program : api
!!! source  : ctqmc_api.f90
!!! type    : module
!!! author  : li huang (email:huangli712@gmail.com)
!!! history : 01/07/2014 by li huang
!!!           01/11/2014 by li huang
!!!           11/11/2014 by li huang
!!! purpose : the purpose of this module is to define a generic and robust
!!!           application programming interface (API) for continuous-time
!!!           quantum Monte Carlo impurity solver
!!! status  : unstable
!!! comment :
!!!-----------------------------------------------------------------------

!!
!!
!! Introduction
!! ============
!!
!! This module can provide a light weight interface (i.e., application
!! programming interface, API) for Fortran/Python language to the ctqmc
!! quantum impurity solver. The user can use it to access the azalea,
!! gardenia, narcissus, begonia, lavender, pansy, and manjushaka codes.
!!
!! How to build the Fortran API
!! ============================
!!
!! 1. edit src/build/make.sys
!! --------------------------
!!
!! Activate the API macro (keep F2PY macro disable).
!!
!! 2. compile api
!! --------------
!!
!! Please compile api (this directory) again. You can use the 'make api'
!! command in the src/build directory.
!!
!! 3. compile the ctqmc component
!! ------------------------------
!!
!! Please compile the desired ctqmc component again. You have to clean it
!! at first, and then compile it. Noted that you have to compile it in the
!! library mode, i.e., you must use 'make lib' (in the src/ctqmc/azalea
!! directory) or 'make azalea-lib' (in the src/build directory), etc.
!!
!! 4. get what you need
!! --------------------
!!
!! If everything is OK, you will find the libctqmc.a file in the ctqmc
!! component folder (for example, src/ctqmc/azalea directory). Please copy
!! it (together with the api.mod) to your own directory. That's all.
!!
!! How to build the Python API
!! ===========================
!!
!! 1. edit src/build/make.sys
!! --------------------------
!!
!! Activate the API macro and F2PY macro at the same time.
!!
!! 2. compile api
!! --------------
!!
!! Please compile api (this directory) again. You can use the 'make api'
!! command in the src/build directory. This step is mandatory.
!!
!! 3. compile the ctqmc component
!! ------------------------------
!!
!! Please compile the desired ctqmc component again. You have to clean it
!! at first, and then compile it. Noted that you have to compile it in the
!! library mode, i.e., you must use 'make lib' (in the src/ctqmc/azalea
!! directory) or 'make azalea-lib' (in the src/build directory), etc.
!!
!! 4. edit src/ctqmc/api/Makefile
!! ------------------------------
!!
!! check the target 'ctqmc', the original action is as follows:
!!     cp ../azalea/libctqmc.a .
!! If you want to use the other ctqmc components, instead of azalea, you
!! have to change the directory. BE CAREFUL!
!!
!! 5. generate pyiqist.so
!! ----------------------
!!
!! In the src/ctqmc/api directory, just input 'make pyiqist' command and
!! wait. At last you will get the pyiqist.so which is what you need.
!!
!! Usage (Fortran version)
!! =======================
!!
!! In the following, we will use azalea code as an example to show how to
!! use api to control it. When you want to compile your code, you have to
!! ensure that api.mod and libctqmc.a are in correct PATH. Or else the
!! compiler will complain that it can not find them.
!!
!! 1. import api support
!! ---------------------
!!
!! use api
!!
!! 2. create T_mpi
!! ---------------
!!
!! type (T_mpi) :: I_mpi           ! define I_mpi
!! ...
!! call mp_init()                  ! init mpi environment
!! call mp_comm_rank(I_mpi%myid)   ! init I_mpi structure
!! call mp_comm_size(I_mpi%nprocs) ! init I_mpi structure
!!
!! Note: The above codes need MPI support. Namely, you have to import the
!! mpi support explicitly.
!!
!! use mmpi ! import mpi support
!!
!! 3. create T_segment_azalea
!! --------------------------
!!
!! type (T_segment_azalea) :: I_solver ! define I_solver
!! ...
!! I_solver%isscf  = 1  ! setup I_solver
!! I_solver%issun  = 1
!! I_solver%isspn  = 2
!! I_solver%isbin  = 1
!! I_solver%nband  = 1
!! I_solver%nspin  = 2
!! I_solver%norbs  = 2
!! I_solver%ncfgs  = 4
!! I_solver%niter  = 20
!! I_solver%mkink  = 1024
!! I_solver%mfreq  = 8193
!! I_solver%nfreq  = 128
!! I_solver%ntime  = 1024
!! I_solver%nflip  = 10000
!! I_solver%ntherm = 20000
!! I_solver%nsweep = 20000000
!! I_solver%nwrite = 2000000
!! I_solver%nclean = 20000
!! I_solver%nmonte = 100
!! I_solver%ncarlo = 100
!!
!! I_solver%U     = 4.0
!! I_solver%Uc    = 4.0
!! I_solver%Uv    = 4.0
!! I_solver%Jz    = 0.0
!! I_solver%Js    = 0.0
!! I_solver%Jp    = 0.0
!! I_solver%mune  = 2.0
!! I_solver%beta  = 10.0
!! I_solver%part  = 0.50
!! I_solver%alpha = 0.50
!!
!! Note: If you want to use the other solvers, instead of the azalea code,
!! please choose suitable solver type.
!!
!! Note: Every parameter for quantum impurity solver must be initialized
!! here, or else the solver will not work properly.
!!
!! 4. init the ctqmc impurity solver
!! ---------------------------------
!!
!! call init_ctqmc(I_mpi, I_solver)
!!
!! 5. setup hybf, symm, eimp, and ktau
!! -----------------------------------
!!
!! For examples:
!!
!! integer :: size_t
!! complex(dp) :: hybf(size_t)
!! ...
!! call set_hybf(size_t, hybf) ! setup hybridization function: hybf
!!
!! Note: This step is optional, because the ctqmc will provide default
!! values for hybf, symm, eimp, and ktau or read them from external
!! disk files.
!!
!! 6. start the ctqmc impurity solver
!! ----------------------------------
!!
!! call exec_ctqmc(i)
!!
!! Here i is the current iteration number.
!!
!! 7. retrieve the calculated results
!! ----------------------------------
!!
!! Through this api, the user can only access the sigf (i.e., self-energy
!! function), grnf (i.e., impurity Green's function) directly, nmat (i.e.,
!! impurity occupation number), and nnmat (i.e., double occupation number)
!! via this api. As for the other physical observables, the user should
!! check the other output files generated by iqist.
!!
!! integer :: size_t
!! complex(dp) :: grnf(size_t)
!! call get_grnf(size_t, grnf)
!!
!! 8. close the ctqmc impurity solver
!! ----------------------------------
!!
!! call stop_ctqmc()
!!
!! 9. finalize the mpi environment
!! --------------------------------
!!
!! call mp_barrier()
!! call mp_finalize()
!!
!! Note: This step is also optional.
!!
!! Usage (Python version)
!! ======================
!!
!! In the following, we will use azalea code as an example to show how to
!! use api to control it. When you want to run your Python code, you have
!! to ensure that pyiqist.so is in correct PATH. Or else the Python will
!! complain that it can not find iqist.
!!
!! 1. import mpi support
!! ---------------------
!!
!! import pyalps.mpi
!!
!! We are not sure whether mpi4py will work. But pyalps.mpi always works.
!! This code will also start the mpi running environment implicitly.
!!
!! 2. import pyiqist
!! -----------------
!!
!! import pyiqist
!!
!! 3. configure the ctqmc impurity solver
!! --------------------------------------
!!
!! You have to setup the parameters for the ctqmc impurity solver, and
!! write them down to the 'solver.ctqmc.in' file. Now you must do that
!! manually. In the future we will provide a Python module to facilitate
!! this work (see src/tools/hibiscus/script/p_ctqmc.py).
!!
!! 4. init the ctqmc impurity solver
!! ---------------------------------
!!
!! pyiqist.api.init_ctqmc(my_id, num_procs)
!!
!! Here my_id means the rank for current process, and num_procs means
!! number of processes.
!!
!! 5. setup hybf, symm, eimp, and ktau
!! -----------------------------------
!!
!! This step has not been tested yet. I am sorry.
!!
!! 6. start the ctqmc impurity solver
!! ----------------------------------
!!
!! pyiqist.api.exec_ctqmc(i)
!!
!! Here i is the current iteration number.
!!
!! 7. retrieve the calculated results
!! ----------------------------------
!!
!! This step has not been tested yet. I am sorry.
!!
!! 8. close the ctqmc impurity solver
!! ----------------------------------
!!
!! pyiqist.api.stop_ctqmc()
!!
!! FAQ
!! ===
!!
!! Question:
!! ---------
!!
!! Can we change the ctqmc impurity solver at runtime?
!!
!! Answer:
!! -------
!!
!! No. You can not change the ctqmc impurity solver dynamically. Once
!! the pyiqist.so is generated, the ctqmc impurity solver is determined.
!! If you want to change the ctqmc impurity solver, you must regenerate
!! the pyiqist.so file at first.
!!
!!

  module api
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
     integer, public, parameter :: solver_id_pansy          = 301
     integer, public, parameter :: solver_id_manjushaka     = 302

! solver status, 1 means ready, 0 means not ready
     integer, public, parameter :: solver_is_ready_azalea     = 1
     integer, public, parameter :: solver_is_ready_gardenia   = 1
     integer, public, parameter :: solver_is_ready_narcissus  = 1
     integer, public, parameter :: solver_is_ready_begonia    = 1
     integer, public, parameter :: solver_is_ready_lavender   = 1
     integer, public, parameter :: solver_is_ready_pansy      = 1
     integer, public, parameter :: solver_is_ready_manjushaka = 1

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

! define type T_solver, which is used to describe the generic abstract
! ctqmc impurity solver
! note: it can not be used directly
     private :: T_solver
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
     private :: T_segment_solver
     type, extends (T_solver) :: T_segment_solver
         character(len=10) :: solver_type = 'SEGMENT'
     end type T_segment_solver

! define type T_general_solver, which is used to describe the ctqmc
! impurity solver which based on general matrix formulation
! note: it can not be used directly
     private :: T_general_solver
     type, extends (T_solver) :: T_general_solver
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
         integer :: isvrt
         integer :: ifast
         integer :: itrun
         integer :: ifast
         integer :: lemax
         integer :: legrd
         integer :: chmax
         integer :: chgrd
         integer :: nffrq
         integer :: nbfrq
         integer :: npart
     end type T_general_manjushaka

# endif  /* F2PY */

!!========================================================================
!!>>> declare accessibility for module routines                        <<<
!!========================================================================

     public :: solver_id
     public :: solver_status

     public :: init_ctqmc
     public :: exec_ctqmc
     public :: stop_ctqmc

     public :: set_hybf
     public :: set_symm
     public :: set_eimp
     public :: set_ktau

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

# if !defined (F2PY)

!!>>> init_ctqmc: initialize the ctqmc quantum impurity solver
!!>>> fortran version
  subroutine init_ctqmc(I_mpi, I_solver)
     implicit none

! external arguments
! type structure of mpi
     class(*), intent(in) :: I_mpi

! type structure of generic solver
     class(*), intent(in) :: I_solver

     call cat_init_ctqmc(I_mpi, I_solver)

     return
  end subroutine init_ctqmc

# else   /* F2PY */

!!>>> init_ctqmc: initialize the ctqmc quantum impurity solver
!!>>> python version
  subroutine init_ctqmc(my_id, num_procs)
     implicit none

! external arguments
! id for current process
     integer, intent(in) :: my_id

! number of processors
     integer, intent(in) :: num_procs

! declare f2py directives
!F2PY intent(in) my_id
!F2PY intent(in) num_procs

     call cat_init_ctqmc(my_id, num_procs)

     return
  end subroutine init_ctqmc

# endif  /* F2PY */

!!>>> exec_ctqmc: execute the ctqmc quantum impurity solver
  subroutine exec_ctqmc(iter)
     implicit none

! external arguments
! current iteration number
     integer, intent(in) :: iter

! declare f2py directives
!F2PY intent(in) iter

     call cat_exec_ctqmc(iter)

     return
  end subroutine exec_ctqmc

!!>>> stop_ctqmc: stop the ctqmc quantum impurity solver
  subroutine stop_ctqmc()
     implicit none

     call cat_stop_ctqmc()

     return
  end subroutine stop_ctqmc

!!>>> set_hybf: setup the impurity hybridization function
  subroutine set_hybf(size_t, hybf_t)
     implicit none

! external arguments
! size of hybf
     integer, intent(in)     :: size_t

! hybridization function
     complex(dp), intent(in) :: hybf_t(size_t)

! declare f2py directives
!F2PY intent(in) size_t
!F2PY intent(in) hybf_t
!F2PY depend(size_t) hybf_t

     call cat_set_hybf(size_t, hybf_t)

     return
  end subroutine set_hybf

!!>>> set_symm: setup the symmetry vector
  subroutine set_symm(size_t, symm_t)
     implicit none

! external arguments
! size of symm
     integer, intent(in) :: size_t

! symmetry vector
     integer, intent(in) :: symm_t(size_t)

! declare f2py directives
!F2PY intent(in) size_t
!F2PY intent(in) symm_t
!F2PY depend(size_t) symm_t

     call cat_set_symm(size_t, symm_t)

     return
  end subroutine set_symm

!!>>> set_eimp: setup the impurity energy level
  subroutine set_eimp(size_t, eimp_t)
     implicit none

! external arguments
! size of eimp
     integer, intent(in)  :: size_t

! impurity energy level
     real(dp), intent(in) :: eimp_t(size_t)

! declare f2py directives
!F2PY intent(in) size_t
!F2PY intent(in) eimp_t
!F2PY depend(size_t) eimp_t

     call cat_set_eimp(size_t, eimp_t)

     return
  end subroutine set_eimp

!!>>> set_ktau: setup the screening function and its first derivates for
!!>>> the dynamical screening effect
!!>>> note: only the narcissus code will implement the cat_set_ktau()
  subroutine set_ktau(size_t, ktau_t, ptau_t)
     implicit none

! external arguments
! size of ktau
     integer, intent(in)  :: size_t

! screening function K(\tau)
     real(dp), intent(in) :: ktau_t(size_t)

! first derivates of screening function K'(\tau)
     real(dp), intent(in) :: ptau_t(size_t)

! declare f2py directives
!F2PY intent(in) size_t
!F2PY intent(in) ktau_t
!F2PY intent(in) ptau_t
!F2PY depend(size_t) ktau_t
!F2PY depend(size_t) ptau_t

     call cat_set_ktau(size_t, ktau_t, ptau_t)

     return
  end subroutine set_ktau

!!>>> get_grnf: extract the impurity green's function
  subroutine get_grnf(size_t, grnf_t)
     implicit none

! external arguments
! size of grnf
     integer, intent(in)      :: size_t

! impurity green's function
     complex(dp), intent(out) :: grnf_t(size_t)

! declare f2py directives
!F2PY intent(in) size_t
!F2PY intent(out) grnf_t
!F2PY depend(size_t) grnf_t

     call cat_get_grnf(size_t, grnf_t)

     return
  end subroutine get_grnf

!!>>> get_sigf: extract the self-energy function
  subroutine get_sigf(size_t, sigf_t)
     implicit none

! external arguments
! size of sigf
     integer, intent(in)      :: size_t

! self-energy function
     complex(dp), intent(out) :: sigf_t(size_t)

! declare f2py directives
!F2PY intent(in) size_t
!F2PY intent(out) sigf_t
!F2PY depend(size_t) sigf_t

     call cat_get_sigf(size_t, sigf_t)

     return
  end subroutine get_sigf

!!>>> get_nmat: extract the occupation number
  subroutine get_nmat(size_t, nmat_t)
     implicit none

! external arguments
! size of nmat
     integer, intent(in)   :: size_t

! occupation number
     real(dp), intent(out) :: nmat_t(size_t)

! declare f2py directives
!F2PY intent(in) size_t
!F2PY intent(out) nmat_t
!F2PY depend(size_t) nmat_t

     call cat_get_nmat(size_t, nmat_t)

     return
  end subroutine get_nmat

!!>>> get_nnmat: extract the double occupation number
  subroutine get_nnmat(size_t, nnmat_t)
     implicit none

! external arguments
! size of nnmat
     integer, intent(in)   :: size_t

! double occupation number
     real(dp), intent(out) :: nnmat_t(size_t)

! declare f2py directives
!F2PY intent(in) size_t
!F2PY intent(out) nnmat_t
!F2PY depend(size_t) nnmat_t

     call cat_get_nnmat(size_t, nnmat_t)

     return
  end subroutine get_nnmat

  end module api
