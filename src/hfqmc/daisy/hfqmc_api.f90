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
!! 1. edit src/build/make.sys
!! --------------------------
!!
!! Activate the API macro (keep MPY macro disable).
!!
!! 2. compile the hfqmc component
!! ------------------------------
!!
!! Please compile the desired hfqmc component again. You have to clean it
!! at first, and then compile it. Noted that you have to compile it in the
!! library mode, i.e., you must use 'make lib' (in the src/hfqmc/daisy
!! directory) or 'make daisy-lib' (in the src/build directory), etc.
!!
!! 3. get what you need
!! --------------------
!!
!! If everything is OK, you will find the libhfqmc.a file in the hfqmc
!! component folder (for example, src/hfqmc/daisy directory). Please copy
!! it (together with the dapi.mod) to your own directory. That's all.
!!
!! How to build the Python API
!! ===========================
!!
!! 1. edit src/build/make.sys
!! --------------------------
!!
!! Activate the API macro and MPY macro at the same time.
!!
!! 2. compile the hfqmc component
!! ------------------------------
!!
!! Please compile the desired hfqmc component again. You have to clean it
!! at first, and then compile it. Noted that you have to compile it in the
!! library mode, i.e., you must use 'make lib' (in the src/hfqmc/daisy
!! directory) or 'make daisy-lib' (in the src/build directory), etc.
!!
!! 3. generate pydaisy.so
!! ----------------------
!!
!! In the src/hfqmc/daisy directory, just input 'make pydaisy' command and
!! wait. At last you will get the pydaisy.so which is what you need.
!!
!! Usage (Fortran version)
!! =======================
!!
!! In the following, we will use daisy code as an example to show how to
!! use api to control it. When you want to compile your code, you have to
!! ensure that dapi.mod and libhfqmc.a are in correct PATH. Or else the
!! compiler will complain that it can not find them.
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
!! mpi support explicitly.
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
!! call init_hfqmc(I_mpi, I_solver)
!!
!! 5. setup wssf, symm, eimp, and ktau
!! -----------------------------------
!!
!! For examples:
!!
!! integer :: size_t
!! complex(dp) :: wssf(size_t)
!! ...
!! call set_wssf(size_t, wssf) ! setup bath weiss's function: hybf
!!
!! Note: This step is optional, because the hfqmc will provide default
!! values for wssf, symm, eimp, and ktau or read them from external
!! disk files.
!!
!! 6. start the hfqmc impurity solver
!! ----------------------------------
!!
!! call exec_hfqmc(i)
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
!! 8. close the hfqmc impurity solver
!! ----------------------------------
!!
!! call stop_hfqmc()
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
!! from pydaisy import dapi as hfqmc
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
!! hfqmc.init_hfqmc(my_id, num_procs)
!!
!! Here my_id means the rank for current process, and num_procs means
!! number of processes. If you are using the mpi4py package to provide
!! mpi support, then you can use the following code to init the hfqmc
!! impurity solver.
!!
!! comm = MPI.COMM_WORLD
!! hfqmc.init_hfqmc(comm.rank, comm.size)
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
!! hfqmc.set_wssf(size_t, wssf)
!!
!! Note: We strongly recommend to select the fortran rule to manage the
!! array memory.
!!
!! Note: This step is optional, because the hfqmc will provide default
!! values for wssf, symm, eimp, and ktau or read them from external
!! disk files.
!!
!! 6. start the hfqmc impurity solver
!! ----------------------------------
!!
!! hfqmc.exec_hfqmc(i)
!!
!! Here i is the current iteration number.
!!
!! 7. retrieve the calculated results
!! ----------------------------------
!!
!! For examples:
!!
!! size_t = norbs
!! nmat = hfqmc.get_nmat(size_t)
!! print nmat
!!
!! size_t = mfreq * norbs
!! grnf = hfqmc.get_grnf(size_t)
!! grnf = numpy.reshape(grnf, (mfreq, norbs), order = 'F')
!!
!! Note: You have to pay attention to the array order when you try to use
!! numpy.reshape() to convert a 1-D array to a 2-D array.
!!
!! 8. close the hfqmc impurity solver
!! ----------------------------------
!!
!! hfqmc.stop_hfqmc()
!!
!! FAQ
!! ===
!!
!! Question:
!! ---------
!!
!! Can we change the hfqmc impurity solver at runtime?
!!
!! Answer:
!! -------
!!
!! No. You can not change the hfqmc impurity solver dynamically. Once
!! the pydaisy.so is generated, the hfqmc impurity solver is determined.
!! If you want to change the hfqmc impurity solver, you must regenerate
!! the pydaisy.so file at first.
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
     integer, public, parameter :: solver_id_daisy     = 901

! solver status, 1 means ready, 0 means not ready
     integer, public, parameter :: solver_is_ready_daisy = 1

!!========================================================================
!!>>> declare global data structure                                    <<<
!!========================================================================

! note: now f2py does not support derived types, so we have to comment
! out them when f2py is used.

# if !defined (MPY)

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

# endif  /* MPY */

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

# if !defined (MPY)

!!>>> init_hfqmc: initialize the hfqmc quantum impurity solver
!!>>> fortran version
  subroutine init_hfqmc(I_mpi, I_solver)
     implicit none

! external arguments
! type structure of mpi
     class(*), intent(in) :: I_mpi

! type structure of generic solver
     class(*), intent(in) :: I_solver

     call cat_init_hfqmc(I_mpi, I_solver)

     return
  end subroutine init_hfqmc

# else   /* MPY */

!!>>> init_hfqmc: initialize the hfqmc quantum impurity solver
!!>>> python version
  subroutine init_hfqmc(my_id, num_procs)
     implicit none

! external arguments
! id for current process
     integer, intent(in) :: my_id

! number of processors
     integer, intent(in) :: num_procs

! declare f2py directives
!F2PY intent(in) my_id
!F2PY intent(in) num_procs

     call cat_init_hfqmc(my_id, num_procs)

     return
  end subroutine init_hfqmc

# endif  /* MPY */

!!>>> exec_hfqmc: execute the hfqmc quantum impurity solver
  subroutine exec_hfqmc(iter)
     implicit none

! external arguments
! current iteration number
     integer, intent(in) :: iter

! declare f2py directives
!F2PY intent(in) iter

     call cat_exec_hfqmc(iter)

     return
  end subroutine exec_hfqmc

!!>>> stop_hfqmc: stop the hfqmc quantum impurity solver
  subroutine stop_hfqmc()
     implicit none

     call cat_stop_hfqmc()

     return
  end subroutine stop_hfqmc

!!>>> set_wssf: setup the bath weiss's function
  subroutine set_wssf(size_t, wssf_t)
     implicit none

! external arguments
! size of wssf
     integer, intent(in)    :: size_t

! bath weiss's function
     complex(8), intent(in) :: wssf_t(size_t)

! declare f2py directives
!F2PY intent(in) size_t
!F2PY intent(in) wssf_t
!F2PY depend(size_t) wssf_t

     call cat_set_wssf(size_t, wssf_t)

     return
  end subroutine set_wssf

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
     integer, intent(in) :: size_t

! impurity energy level
     real(8), intent(in) :: eimp_t(size_t)

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
     integer, intent(in) :: size_t

! screening function K(\tau)
     real(8), intent(in) :: ktau_t(size_t)

! first derivates of screening function K'(\tau)
     real(8), intent(in) :: ptau_t(size_t)

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
     integer, intent(in)     :: size_t

! impurity green's function
     complex(8), intent(out) :: grnf_t(size_t)

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
     integer, intent(in)     :: size_t

! self-energy function
     complex(8), intent(out) :: sigf_t(size_t)

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
     integer, intent(in)  :: size_t

! occupation number
     real(8), intent(out) :: nmat_t(size_t)

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
     integer, intent(in)  :: size_t

! double occupation number
     real(8), intent(out) :: nnmat_t(size_t)

! declare f2py directives
!F2PY intent(in) size_t
!F2PY intent(out) nnmat_t
!F2PY depend(size_t) nnmat_t

     call cat_get_nnmat(size_t, nnmat_t)

     return
  end subroutine get_nnmat

  end module dapi
