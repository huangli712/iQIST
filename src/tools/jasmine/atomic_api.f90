!!!-----------------------------------------------------------------------
!!! project : lilac
!!! program : japi
!!!           japi@T_jasmine
!!!           japi@init_atomic
!!!           japi@exec_atomic
!!!           japi@stop_atomic
!!! source  : atomic_api.f90
!!! type    : module
!!! author  : li huang (email:huangli712@gmail.com)
!!! history : 10/29/2014 by li huang
!!! purpose : the purpose of this module is to define a generic and robust
!!!           application programming interface (API) for atomic eigenvalue
!!!           problem solver
!!! status  : unstable
!!! comment :
!!!-----------------------------------------------------------------------

!!
!!
!! Introduction
!! ============
!!
!! This module can provide a light weight interface (i.e., application
!! programming interface, API) for Fortran/Python language to the atomic
!! eigenvalue problem solver. The user can use it to access the jasmine
!! code.
!!
!! How to build the Fortran API
!! ============================
!!
!! 1. edit src/build/make.sys
!! --------------------------
!!
!! Activate the API macro (keep F2PY macro disable).
!!
!! 2. compile the jasmine component
!! --------------------------------
!!
!! Please compile the desired jasmine component again. You have to clean it
!! at first, and then compile it. Noted that you have to compile it in the
!! library mode, i.e., you must use 'make lib' (in the src/tools/jasmine
!! directory) or 'make jasmine-lib' (in the src/build directory), etc.
!!
!! 3. get what you need
!! --------------------
!!
!! If everything is OK, you will find the libatomic.a file in the jasmine
!! component folder (for example, src/tools/jasmine directory). Please copy
!! it (together with the japi.mod) to your own directory. That's all.
!!
!! How to build the Python API
!! ===========================
!!
!! 1. edit src/build/make.sys
!! --------------------------
!!
!! Activate the API macro and F2PY macro at the same time.
!!
!! 2. generate pyjasmine.so
!! ------------------------
!!
!! In the src/tools/jasmine directory, just input 'make pyjasmine' command
!! and wait. At last you will get the pyjasmine.so which is what you need.
!!
!! Usage (Fortran version)
!! =======================
!!
!! In the following, we will use an example to show you how to use api to
!! control the jasmine code. When you want to compile your code, you have
!! to ensure that japi.mod and libatomic.a are in correct PATH. Or else
!! the compiler will complain that it can not find them.
!!
!! 1. import api support
!! ---------------------
!!
!! use japi
!!
!! 2. create T_jasmine
!! -------------------
!!
!! type (T_jasmine) :: I_solver ! define I_solver
!! ...
!! I_solver%ibasis = 1  ! setup I_solver
!! I_solver%ictqmc = 1
!! ...
!! I_solver%Uc     = 4.0
!! I_solver%Uv     = 4.0
!! I_solver%Jz     = 0.0
!! I_solver%Js     = 0.0
!! I_solver%Jp     = 0.0
!! ...
!!
!! Note: Every parameter for the atomic eigenvalue problem solver must be
!! initialized here.
!!
!! 3. init the atomic eigenvalue problem solver
!! --------------------------------------------
!!
!! call init_atomic(I_solver)
!!
!! 4. start the atomic eigenvalue problem solver
!! ---------------------------------------------
!!
!! call exec_atomic()
!!
!! 5. close the atomic eigenvalue problem solver
!! ---------------------------------------------
!!
!! call stop_atomic()
!!
!! 6. access the computational results
!! -----------------------------------
!!
!! You have to write your fortran codes manually to access the results.
!!
!! Usage (Python version)
!! ======================
!!
!! In the following, we will use an example to show you how to use api to
!! control the jasmine code. When you want to run your Python code, you
!! have to ensure that pyjasmine.so is in correct PATH. Or else the Python
!! will complain that it can not find iqist.
!!
!! 1. import pyjasmine
!! -------------------
!!
!! import pyjasmine
!!
!! 2. configure the atomic eigenvalue problem solver
!! -------------------------------------------------
!!
!! You have to setup the key parameters for the atomic eigenvalue problem
!! solver, and write them down to the 'atom.config.in' file. Now you must
!! do that manually. In the future we will provide a Python module to
!! facilitate this work.
!!
!! 3. init the atomic eigenvalue problem solver
!! --------------------------------------------
!!
!! pyjasmine.japi.init_atomic() # there is no parameter for init_atomic()
!!
!! 4. start the atomic eigenvalue problem solver
!! ---------------------------------------------
!!
!! pyjasmine.japi.exec_atomic()
!!
!! 5. close the atomic eigenvalue problem solver
!! ---------------------------------------------
!!
!! pyjasmine.japi.stop_atomic()
!!
!! 6. access the computational results
!! -----------------------------------
!!
!! You have to write your Python codes manually to access the results.
!!
!!

  module japi
     implicit none

!!========================================================================
!!>>> declare global parameters                                        <<<
!!========================================================================

! dp: number precision, double precision for reals
     integer, private, parameter :: dp = kind(1.0d0)

!!========================================================================
!!>>> declare global data structure                                    <<<
!!========================================================================

! note: now f2py does not support derived types, so we have to comment
! out them when f2py is used.

# if !defined (F2PY)

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

# endif  /* F2PY */

!!========================================================================
!!>>> declare accessibility for module routines                        <<<
!!========================================================================

     public  :: init_atomic
     public  :: exec_atomic
     public  :: stop_atomic

  contains ! encapsulated functionality

# if !defined (F2PY)

!!>>> init_atomic: initialize the atomic eigenvalue problem solver
!!>>> fortran version
  subroutine init_atomic(I_solver)
     implicit none

! external arguments
! type structure of generic atomic eigenvalue problem solver
     class(*), intent(in) :: I_solver

     call cat_init_atomic(I_solver)

     return
  end subroutine init_atomic

# else   /* F2PY */

!!>>> init_atomic: initialize the atomic eigenvalue problem solver
!!>>> python version
  subroutine init_atomic()
     implicit none

     call cat_init_atomic()

     return
  end subroutine init_atomic

# endif  /* F2PY */

!!>>> exec_atomic: execute the atomic eigenvalue problem solver
  subroutine exec_atomic()
     implicit none

     call cat_exec_atomic()

     return
  end subroutine exec_atomic

!!>>> stop_atomic: stop the atomic eigenvalue problem solver
  subroutine stop_atomic()
     implicit none

     call cat_stop_atomic()

     return
  end subroutine stop_atomic

  end module japi
