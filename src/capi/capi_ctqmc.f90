!!!-----------------------------------------------------------------------
!!! project : CAPI (Common Application Programming Interface)
!!! program : capi
!!! source  : capi_ctqmc.f90
!!! type    : module
!!! author  : li huang (email:lihuang.dmft@gmail.com)
!!! history : 01/07/2014 by li huang (created)
!!!           01/27/2017 by li huang (last modified)
!!! purpose : the purpose of this module is to define a generic and robust
!!!           application programming interface (API) for continuous-time
!!!           quantum Monte Carlo impurity solver.
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
!! quantum impurity solver. The user can use it to access the gardenia,
!! narcissus, lavender, and manjushaka codes.
!!
!! How to build the Fortran API
!! ============================
!!
!! 1. edit iqist/build/make.sys
!! ----------------------------
!!
!! Setup the compiling environment correctly.
!!
!! 2. compile base and capi
!! ------------------------
!!
!! Please compile base and capi at first. You can type the 'make base'
!! and 'make capi' commands in the iqist/build directory, or use the
!! 'make' command in the src/base and src/capi directories.
!!
!! 3. compile the ctqmc component
!! ------------------------------
!!
!! Noted that you have to compile it in the library mode, i.e., you must
!! use 'make lib' in the src/ctqmc/gardenia directory or 'make gardenia-lib'
!! in the iqist/build directory.
!!
!! 4. get what you need
!! --------------------
!!
!! If everything is OK, you will find the libctqmc.a file in the ctqmc
!! component folder (for example, src/ctqmc/gardenia directory). Please
!! copy it (together with src/capi/capi.mod and src/base/libMM.a) to your
!! own directory. That's all.
!!
!! How to build the Python API
!! ===========================
!!
!! 1. edit iqist/build/make.sys
!! ----------------------------
!!
!! Setup the compiling environment correctly.
!!
!! 2. compile base and capi
!! ------------------------
!!
!! Please compile base and capi at first. You can type the 'make base'
!! and 'make capi' commands in the iqist/build directory, or use the
!! 'make' command in the src/base and src/capi directories.
!!
!! 3. compile the ctqmc component
!! ------------------------------
!!
!! Noted that you have to compile it in the python library mode, i.e.,
!! you have to type 'make pylib' in the src/ctqmc/gardenia directory or
!! 'make gardenia-pylib' in the iqist/build directory.
!!
!! 4. get what you need
!! --------------------
!!
!! If everything is OK, you will find the pyiqist.so file in the ctqmc
!! component folder (for example, src/ctqmc/gardenia directory). Please
!! copy it to your own directory. That's all.
!!
!! Usage (Fortran version)
!! =======================
!!
!! In the following, we will use gardenia code as an example to show how
!! to use api to control it. When you want to compile your code, you have
!! to ensure that capi.mod, libMM.a and libctqmc.a are in correct PATH. Or
!! else the compiler will complain that it can not find them.
!!
!! 1. import api support
!! ---------------------
!!
!! use capi
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
!! mpi support explicitly, such as,
!!
!! use mmpi ! import mpi support
!!
!! 3. create T_segment_gardenia
!! ----------------------------
!!
!! type (T_segment_gardenia) :: I_solver ! define I_solver
!! ...
!! I_solver%isscf  = 1  ! setup I_solver
!! I_solver%issun  = 2
!! I_solver%isspn  = 1
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
!! ...
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
!! ...
!!
!! Note: If you want to use the other solvers, instead of the gardenia
!! code, please choose suitable solver type.
!!
!! Note: Every parameter for quantum impurity solver must be initialized
!! here, or else the solver will not work properly.
!!
!! 4. init the ctqmc impurity solver
!! ---------------------------------
!!
!! call cat_init_ctqmc(I_mpi, I_solver)
!!
!! 5. setup hybf, symm, eimp, and ktau
!! -----------------------------------
!!
!! For examples:
!!
!! integer :: size_t
!! complex(dp) :: hybf(size_t)
!! ...
!! call cat_set_hybf(size_t, hybf) ! setup hybridization function: hybf
!!
!! Note: This step is optional, because the ctqmc will provide default
!! values for hybf, symm, eimp, and ktau or read them from external
!! disk files.
!!
!! 6. start the ctqmc impurity solver
!! ----------------------------------
!!
!! call cat_exec_ctqmc(i)
!!
!! Here i is the current iteration number.
!!
!! 7. retrieve the calculated results
!! ----------------------------------
!!
!! Through this api, the user can only access the sigf (i.e., self-energy
!! function), grnf (i.e., impurity Green's function), nmat (i.e., impurity
!! occupation number), and nnmat (i.e., double occupation number). As for
!! the other physical observables, the user should parse the other output
!! files generated by iqist.
!!
!! integer :: size_t
!! complex(dp) :: grnf(size_t)
!! call cat_get_grnf(size_t, grnf)
!!
!! 8. close the ctqmc impurity solver
!! ----------------------------------
!!
!! call cat_stop_ctqmc()
!!
!! 9. finalize the mpi environment
!! -------------------------------
!!
!! call mp_barrier()
!! call mp_finalize()
!!
!! Note: This step is also optional.
!!
!! Usage (Python version)
!! ======================
!!
!! In the following, we will use gardenia code as an example to show how
!! to use api to control it. When you want to run your Python code, you
!! have to ensure that pyiqist.so is in correct PATH. Or else the Python
!! will complain that it can not find iqist.
!!
!! 1. import mpi support
!! ---------------------
!!
!! Here we can use the pyalps.mpi package or mpi4py package to provide
!! the mpi support, such as,
!!
!! import pyalps.mpi
!!
!! or
!!
!! from mpi4py import MPI
!!
!! We recommend to use the mpi4py package since it is included in the
!! scipy package already. The above code will also start the mpi running
!! environment implicitly.
!!
!! 2. import pyiqist
!! -----------------
!!
!! import pyiqist
!!
!! You have to ensure that the pyiqist package is in the sys.path. For
!! example, you can use the following code to modify sys.path
!!
!! sys.path.append('../../src/ctqmc/gardenia/')
!!
!! 3. configure the ctqmc impurity solver
!! --------------------------------------
!!
!! You have to setup the parameters for the ctqmc impurity solver, and
!! write them down to the 'solver.ctqmc.in' file. Now you can do that
!! manually. On the other hand, we provide a powerful Python module to
!! facilitate this work (see src/tools/hibiscus/script/u_ctqmc.py).
!!
!! 4. init the ctqmc impurity solver
!! ---------------------------------
!!
!! pyiqist.cat_init_ctqmc(my_id, num_procs)
!!
!! Here my_id means the rank for current process, and num_procs means
!! number of processes. If you are using the mpi4py package to provide
!! mpi support, then you can use the following code to init the ctqmc
!! impurity solver.
!!
!! comm = MPI.COMM_WORLD
!! pyiqist.cat_init_ctqmc(comm.rank, comm.size)
!!
!! 5. setup hybf, symm, eimp, and ktau
!! -----------------------------------
!!
!! For examples:
!!
!! mfreq = 8193
!! norbs = 2
!! size_t = mfreq * norbs * norbs
!! hybf = numpy.zeros(size_t, dtype = numpy.complex, order = 'F')
!! ...
!! pyiqist.cat_set_hybf(size_t, hybf)
!!
!! Note: We strongly recommend to select the fortran rule to manage the
!! numpy array memory.
!!
!! Note: This step is optional, because the ctqmc will provide default
!! values for hybf, symm, eimp, and ktau or read them from external
!! disk files.
!!
!! 6. start the ctqmc impurity solver
!! ----------------------------------
!!
!! pyiqist.cat_exec_ctqmc(i)
!!
!! Here i is the current iteration number.
!!
!! 7. retrieve the calculated results
!! ----------------------------------
!!
!! For examples:
!!
!! size_t = norbs
!! nmat = pyiqist.cat_get_nmat(size_t)
!! print nmat
!!
!! size_t = mfreq * norbs * norbs
!! grnf = pyiqist.cat_get_grnf(size_t)
!! grnf = numpy.reshape(grnf, (mfreq, norbs, norbs), order = 'F')
!!
!! Note: You have to pay attention to the array order when you try to use
!! numpy.reshape() to convert a 1-D array to a 3-D array.
!!
!! 8. close the ctqmc impurity solver
!! ----------------------------------
!!
!! pyiqist.cat_stop_ctqmc()
!!
!! Examples
!! ========
!!
!! Fortran version
!! ---------------
!!
!! see iqist/tutor/t961/template.f90.
!!
!! Python version
!! --------------
!!
!! see iqist/tutor/t962/template.py.
!!
!! FAQ
!! ===
!!
!! Question
!! --------
!!
!! Can we change the ctqmc impurity solver at runtime?
!!
!! Answer
!! ------
!!
!! No. You can not change the ctqmc impurity solver dynamically. Once
!! the pyiqist.so is generated, the ctqmc impurity solver is determined.
!! If you want to change the ctqmc impurity solver, you must regenerate
!! the pyiqist.so file from scratch.
!!
!!

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
     integer, public, parameter :: solver_id_lavender       = 201
     integer, public, parameter :: solver_id_manjushaka     = 301

! solver status, 1 means ready, 0 means not ready
     integer, public, parameter :: solver_is_ready_gardenia   = 1
     integer, public, parameter :: solver_is_ready_narcissus  = 1
     integer, public, parameter :: solver_is_ready_lavender   = 1
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
