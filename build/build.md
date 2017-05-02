# iQIST (Interacting Quantum Impurity Solver Toolkit)

## Introduction

The make.sys file is the key component of the building system. You have to modify it to fulfill your requirements. If it is not configured correctly, the building system won't work correctly as well. So in the following we would like to provide a detailed explanations for it.

## Prerequisites

### Operation system
* Linux
* macOS

### Fortran compiler
* Intel Fortran Compiler

### Linear algebra library
* Apple Accelerate framework
* Reference implementations for BLAS and LAPACK at Netlib
* Intel Math Kernel Library

### MPI environment
* MPICH
* Openmpi

### OpenMP environment
* (Optional)

### Python environment
* (Optional) numpy, scipy, matplotlib

Though it is not mandatory, we still strongly recommend to update the above software components on your systems to the latest versions. The OpenMP and Python environments are optional.

## Explanations

### F90

The Fortran compiler. Both the parallel and sequential fortran compilers are supported. Note that only the Intel fortran compiler was tested. And we recommend to use the very latest version of Intel fortran compiler (i.e, Intel Parallel Studio 2017).

Possible options:

* mpif90
* mpifort
* mpiifort
* ifort

The internal compiler used by mpif90, mpifort, or mpiifort must be ifort.

### LINKER

Linker. Here it should be the same with the fortran compiler. Do not change it.

Possible options:

* $(F90)

### ARCHIVER

Archiver. It is used to pack the binary objects into a library. Do not modify it for ever.

Possible options:

* ar -ruv

### MPI

Specify whether MPI is enable. If you want to compile a sequential code, please comment it out with '#' symbol and then setup F90 to 'ifort'. We strongly suggest to compile the MPI parallelized codes.

Possible options:

* -DMPI

### OMP

Specify whether OpenMP is enable. If you want to disable it, please comment it out. In default it is disabled. So far the OpenMP was used by some ctqmc components to speedup the measurements of some selected two-particle quantities.

Possible options:

* -qopenmp

If you are using old version Intel fortran compiler, this option may be '-openmp'.

### FPP

Specify whether the fortran preprocessor (FPP) is used. It has to be enabled or else the iQIST can not be compiled correctly.

Possible options:

* -fpp

### CPP

Collection of preprocessor directives. Do not modify it unless you are an expert of iQIST.

Possible options:

* $(FPP)
* $(MPI)
* $(OMP)

Please make sure that the '$(FPP)' option is present.

### CHECK

Used to specify what types of check should be done.

Possible options:
* -nogen-interfaces
* -warn all
* -check all
* -traceback
* -g

The '-nogen-interfaces' option ask the compiler to do not generate an interface block for each routine defined in the source file. The '-warn all' option means the check is done in compiling. The '-check all' option means the check will be done in running. The '-traceback' option enables us to track the exact position (line number and file name) where the error occurs. The '-g' option enables the compiler to generate debug information and embed them into the final program. **Note that all of the '-check all', '-traceback', and '-g' options will decrease the efficiency greatly**.

### MTUNE

Collection of optimization options.

Possible options:

* -O3
* -xHost

The '-O3' option means the highest optimization. The '-xHost' option enables the compiler to try to generate the most suitable code for the current computer architecture.

### FFLAGS

Collection of Fortran compiler options. Do not modify them for ever.

Possible options:

* -c
* $(CPP)
* $(CHECK)
* $(MTUNE)

### LFLAGS

Collection of linker options. Do not modify them unless you know what you are doing.

Possible options:

* $(OMP)
* -Wl,-no_pie

The '-Wl,-no_pie' option is useful when you are using the macOS system and want to traceback the code (the -traceback option is applied). If you are using the Linux system, you can skip it.

### LIBS

Specify the external libraries. Now the iQIST software package depends on LAPACK and BLAS heavily. To achieve good performance, the highly optimized LAPACK and BLAS implementations are essential. Here we want to recommend the OpenBLAS and Intel MKL.

Possible options:

* -framework Accelerate
* -L/home/lihuang/lapack -llapack -lblas
* -L/opt/intel/mkl/lib -lmkl_core -lmkl_sequential -lmkl_rt

Here we provide three typical choices. (1) In the macOS system, we can use the Apple Accelerate framework. (2) We use the home-built BLAS and LAPACK libraries. Please pay attention to the path. You have to modify it to meet your software environment. (3) We link the iQIST code with the Intel MKL. Please pay attention to the path and the library's name. You have to modify them to meet your software environment. Please see the documentation about Intel MKL for more details.
