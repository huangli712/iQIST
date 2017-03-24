!!!-----------------------------------------------------------------------
!!! project : CAPI (Common Application Programming Interface)
!!! program : japi
!!! source  : capi_atomic.f90
!!! type    : module
!!! author  : li huang (email:lihuang.dmft@gmail.com)
!!! history : 10/29/2014 by li huang (created)
!!!           01/27/2017 by li huang (last modified)
!!! purpose : the purpose of this module is to define a generic and robust
!!!           application programming interface (API) for atomic eigenvalue
!!!           problem solver.
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
!! 3. compile the jasmine component
!! --------------------------------
!!
!! Noted that you have to compile it in the library mode, i.e., you must
!! use 'make lib' in the src/tools/jasmine directory or 'make jasmine-lib'
!! in the iqist/build directory.
!!
!! 4. get what you need
!! --------------------
!!
!! If everything is OK, you will find the libatomic.a file in the jasmine
!! component folder (for example, src/tools/jasmine directory). Please copy
!! it (together with src/capi/japi.mod and src/base/libMM.a) to your own
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
!! 3. compile the jasmine component
!! --------------------------------
!!
!! Noted that you have to compile it in the python library mode, i.e.,
!! you have to type 'make pylib' in the src/tools/jasmine directory or
!! 'make jasmine-pylib' in the iqist/build directory.
!!
!! 4. get what you need
!! --------------------
!!
!! If everything is OK, you will find the pyjasmine.so file in the jasmine
!! component folder (for example, src/tools/jasmine directory). Please
!! copy it to your own directory. That's all.
!!
!! Usage (Fortran version)
!! =======================
!!
!! In the following, we will use an example to show you how to use api to
!! control the jasmine code. When you want to compile your code, you have
!! to ensure that japi.mod, libMM.a and libatomic.a are in correct PATH. Or
!! else the compiler will complain that it can not find them.
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
!! call cat_init_atomic(I_solver)
!!
!! 4. start the atomic eigenvalue problem solver
!! ---------------------------------------------
!!
!! call cat_exec_atomic()
!!
!! 5. close the atomic eigenvalue problem solver
!! ---------------------------------------------
!!
!! call cat_stop_atomic()
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
!! You have to ensure that the pyjasmine package is in the sys.path. For
!! example, you can use the following code to modify sys.path
!!
!! sys.path.append('../../src/tools/jasmine/')
!!
!! 2. configure the atomic eigenvalue problem solver
!! -------------------------------------------------
!!
!! You have to setup the key parameters for the atomic eigenvalue problem
!! solver, and write them down to the 'atom.config.in' file. Now you can
!! do that manually. On the other hand, we provide a powerful Python module
!! to facilitate this work (see src/tools/hibiscus/script/u_atomic.py).
!!
!! 3. init the atomic eigenvalue problem solver
!! --------------------------------------------
!!
!! pyjasmine.cat_init_atomic() # there is no parameter for cat_init_atomic()
!!
!! 4. start the atomic eigenvalue problem solver
!! ---------------------------------------------
!!
!! pyjasmine.cat_exec_atomic()
!!
!! 5. close the atomic eigenvalue problem solver
!! ---------------------------------------------
!!
!! pyjasmine.cat_stop_atomic()
!!
!! 6. access the computational results
!! -----------------------------------
!!
!! You have to write your own Python codes to access the results.
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
!! see iqist/tutor/t63/template.py.
!!
!!

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
