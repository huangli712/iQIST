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

!!
!!
!! Introduction
!! ============
!!
!! This module can provide a light weight interface (i.e., application
!! programming interface, API) for Fortran/Python language to the hfqmc
!! quantum impurity solver. The user can use it to access the daisy code.
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
!! 3. compile the hfqmc component
!! ------------------------------
!!
!! Noted that you have to compile it in the library mode, i.e., you must
!! use 'make lib' in the src/hfqmc/daisy directory or 'make daisy-lib'
!! in the iqist/build directory.
!!
!! 4. get what you need
!! --------------------
!!
!! If everything is OK, you will find the libhfqmc.a file in the hfqmc
!! component folder (for example, src/hfqmc/daisy directory). Please copy
!! it (together with src/capi/dapi.mod and src/base/libMM.a) to your own
!! directory. That's all.
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
!! 3. compile the hfqmc component
!! ------------------------------
!!
!! Noted that you have to compile it in the python library mode, i.e.,
!! you have to type 'make pylib' in the src/hfqmc/daisy directory or
!! 'make daisy-pylib' in the iqist/build directory.
!!
!! 4. get what you need
!! --------------------
!!
!! If everything is OK, you will find the pydaisy.so file in the hfqmc
!! component folder (for example, src/hfqmc/daisy directory). Please
!! copy it to your own directory. That's all.
!!
!! Usage (Fortran version)
!! =======================
!!
!! In the following, we will use daisy code as an example to show how to
!! use api to control it. When you want to compile your code, you have to
!! ensure that dapi.mod, libMM.a and libhfqmc.a are in correct PATH. Or
!! else the compiler will complain that it can not find them.
!!
!! 1. import api support
!! ---------------------
!!
!! use dapi
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
!! 3. create T_daisy
!! -----------------
!!
!! type (T_daisy) :: I_solver ! define I_solver
!! ...
!! I_solver%isscf  = 1  ! setup I_solver
!! I_solver%issun  = 2
!! I_solver%isspn  = 1
!! I_solver%isbin  = 1
!! I_solver%nband  = 1
!! I_solver%nspin  = 2
!! I_solver%norbs  = 2
!! I_solver%niter  = 20
!! I_solver%mstep  = 16
!! I_solver%mfreq  = 8193
!! I_solver%nsing  = 1
!! I_solver%ntime  = 128
!! I_solver%ntherm = 100
!! I_solver%nsweep = 240000
!! I_solver%nclean = 100
!! I_solver%ncarlo = 10
!!
!! I_solver%Uc    = 4.0
!! I_solver%Jz    = 0.0
!! I_solver%mune  = 2.0
!! I_solver%beta  = 10.0
!! I_solver%part  = 0.50
!! I_solver%alpha = 0.50
!!
!! Note: Every parameter for quantum impurity solver must be initialized
!! here, or else the solver will not work properly.
!!
!! 4. init the hfqmc impurity solver
!! ---------------------------------
!!
!! call cat_init_hfqmc(I_mpi, I_solver)
!!
!! 5. setup wssf, symm, eimp, and ktau
!! -----------------------------------
!!
!! For examples:
!!
!! integer :: size_t
!! complex(dp) :: wssf(size_t)
!! ...
!! call cat_set_wssf(size_t, wssf) ! setup bath weiss's function: wssf
!!
!! Note: This step is optional, because the hfqmc will provide default
!! values for wssf, symm, eimp, and ktau or read them from external
!! disk files.
!!
!! 6. start the hfqmc impurity solver
!! ----------------------------------
!!
!! call cat_exec_hfqmc(i)
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
!! 8. close the hfqmc impurity solver
!! ----------------------------------
!!
!! call cat_stop_hfqmc()
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
!! In the following, we will use daisy code as an example to show how to
!! use api to control it. When you want to run your Python code, you have
!! to ensure that pydaisy.so is in correct PATH. Or else the Python will
!! complain that it can not find iqist.
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
!! 2. import pydaisy
!! -----------------
!!
!! import pydaisy
!!
!! You have to ensure that the pydaisy package is in the sys.path. For
!! example, you can use the following code to modify sys.path
!!
!! sys.path.append('../../src/hfqmc/daisy/')
!!
!! 3. configure the hfqmc impurity solver
!! --------------------------------------
!!
!! You have to setup the parameters for the hfqmc impurity solver, and
!! write them down to the 'solver.hfqmc.in' file. Now you can do that
!! manually. On the other hand, we provide a powerful Python module to
!! facilitate this work (see src/tools/hibiscus/script/u_hfqmc.py).
!!
!! 4. init the hfqmc impurity solver
!! ---------------------------------
!!
!! pydaisy.cat_init_hfqmc(my_id, num_procs)
!!
!! Here my_id means the rank for current process, and num_procs means
!! number of processes. If you are using the mpi4py package to provide
!! mpi support, then you can use the following code to init the hfqmc
!! impurity solver.
!!
!! comm = MPI.COMM_WORLD
!! pydaisy.cat_init_hfqmc(comm.rank, comm.size)
!!
!! 5. setup wssf, symm, eimp, and ktau
!! -----------------------------------
!!
!! For examples:
!!
!! mfreq = 8193
!! norbs = 2
!! size_t = mfreq * norbs
!! wssf = numpy.zeros(size_t, dtype = numpy.complex, order = 'F')
!! ...
!! pydaisy.cat_set_wssf(size_t, wssf)
!!
!! Note: We strongly recommend to select the fortran rule to manage the
!! numpy array memory.
!!
!! Note: This step is optional, because the hfqmc will provide default
!! values for wssf, symm, eimp, and ktau or read them from external
!! disk files.
!!
!! 6. start the hfqmc impurity solver
!! ----------------------------------
!!
!! pydaisy.cat_exec_hfqmc(i)
!!
!! Here i is the current iteration number.
!!
!! 7. retrieve the calculated results
!! ----------------------------------
!!
!! For examples:
!!
!! size_t = norbs
!! nmat = pydaisy.cat_get_nmat(size_t)
!! print nmat
!!
!! size_t = mfreq * norbs
!! grnf = pydaisy.cat_get_grnf(size_t)
!! grnf = numpy.reshape(grnf, (mfreq, norbs), order = 'F')
!!
!! Note: You have to pay attention to the array order when you try to use
!! numpy.reshape() to convert a 1-D array to a 2-D array.
!!
!! 8. close the hfqmc impurity solver
!! ----------------------------------
!!
!! pydaisy.cat_stop_hfqmc()
!!
!! Examples
!! ========
!!
!! Fortran version
!! ---------------
!!
!! N/A.
!!
!! Python version
!! --------------
!!
!! see iqist/tutor/t963/template.py.
!!
!! FAQ
!! ===
!!
!! Question
!! --------
!!
!! Can we change the hfqmc impurity solver at runtime?
!!
!! Answer
!! ------
!!
!! No. You can not change the hfqmc impurity solver dynamically. Once
!! the pydaisy.so is generated, the hfqmc impurity solver is determined.
!! If you want to change the hfqmc impurity solver, you must regenerate
!! the pydaisy.so file from scratch.
!!
!!

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
