# iQIST (Interacting Quantum Impurity Solver Toolkit)

## Introduction

The make.sys file is the key component of the building system. You have to modify it to fulfill your requirements. If it is not configured correctly, the building system won't work correctly as well. So in the following we would like to provide a detailed explanations for it.

## Prerequisites

### Operation system
* Linux
* Mac OS X

### Fortran compiler
* Intel Fortran Compiler
* GNU gfortran

### Linear algebra library
* OpenBLAS
* Intel Math Kernel Library
* Reference implementations for BLAS and LAPACK at Netlib
* Apple Accelerate framework

### MPI environment
* MPICH
* Openmpi

### OpenMP environment
* Optional

### Python environment
* numpy, scipy, f2py, mpi4py, matplotlib

> NOTE:
>
> Though it is not mandatory, we still strongly recommend to update the above software components on your systems to the latest versions.

## Explanations

### F90

The Fortran compiler. Both the parallel and sequential fortran compilers are supported. Note that only the Intel fortran compiler and the GNU gfortran compiler were tested. And we recommend to use the very latest version of Intel fortran compiler (i.e, Intel Fortran Composer XE 2013 or Intel Parallel Studio 2015).

Possible options:

* mpif90
* mpifort
* mpiifort
* ifort
* gfortran

> NOTE:

> The internal compiler used by mpif90, mpifort, or mpiifort must be ifort or gfortran.

### LINKER

Linker. Here it should be the same with the fortran compiler. Do not change it.

Possible options:

* $(F90)

### ARCHIVER

Archiver. It is used to pack the binary objects into a library. Do not modify it for ever.

Possible options:

* ar -ruv

### MPI

Specify whether MPI is enable. If you want to compile a sequential code, please comment it out with '#' symbol and then setup F90 to 'ifort' or 'gfortran'. We strongly suggest to compile the MPI parallelized codes.

Possible options:

* -DMPI

### OMP

Specify whether OpenMP is enable. If you want to disable it, please comment it out. In default it is disabled. So far the OpenMP was used by some ctqmc components to speedup the measurements of some selected two-particle quantities.

Possible options:

* -openmp
* -fopenmp

> NOTE:

> If you are using the Intel fortran compiler, this option should be '-openmp'. On the other hand, if the GNU gfortran compiler is used, this option should be '-fopenmp'.

### FPP

Specify whether the fortran preprocessor (FPP) is used. It has to be enabled or else the iQIST can not be compiled correctly.

Possible options:

* -fpp
* -cpp

> NOTE:

> If you are using the Intel fortran compiler, this option should be '-fpp'. On the other hand, if the GNU gfortran compiler is used, this option should be '-cpp'.

### CPP

Collection of preprocessor directives. Do not modify it unless you are an expert of iQIST.

Possible options:

* $(FPP)
* $(MPI)
* $(OMP)

> NOTE:

> The '$(FPP)' option is necessary.

### GPROF

Specify whether the code profiling should be done. If it is enabled, then after the code is finished, a gmon.out file will be generated. Then you can use the gprof tool to analyze the runtime information and figure out the hotspot of the code. It is not a good idea to enable it to build a production code, because it will decrease the efficiency greatly.

Possible options:

* -p
* -pg

> NOTE:

> The '-p' option is useful for the GNU gfortran compiler, while the '-pg' option is useful for the Intel fortran compiler.

### CHECK

Used to specify what types of check should be done.

Possible options:

* -warn all
* -check all
* -traceback
* -g
* -Wall -Wunused -Wextra
* -fcheck=all
* -fbacktrace

If you are using the Intel fortran compiler, the '-warn all' option means the check is done in compiling. The '-check all' option means the check will be done in running. The '-traceback' option enables us to track the exact position (line number and file name) where the error occurs. The '-g' option enables the compiler to generate debug information and embed them into the final program. **Note that all of the '-check all', '-traceback', and '-g' options will decrease the efficiency greatly**.

If you are using the GNU gfortran compiler, the '-Wall' option will enable most warning messages. The '-Wunused' option will enable all -Wunused- warnings, such as the '-Wunused-parameter' option, etc. The '-Wextra' option will print some extra warnings (sometimes they are unwanted). The '-fcheck=all' option will specify that all of the runtime checks are to be performed. The '-fbacktrace' option will produce a backtrace when a runtime error is encountered. Finally, the '-g' option will enable the compiler to generate debug information and embed them into the final program.

### CDUMP

Specify whether the fortran compiler will output useful optimization information during the compiling process.

Possible options:

* -vec-report2
* -openmp-report2
* -nogen-interfaces
* -fopt-info

If you are using the Intel fortran compiler, the '-vec-report2' option tells the vectorizer to report on vectorized and non-vectorized loops. The '-openmp-report2' option controls the openmp parallelizer's level of diagnostic messages. Here number 2 is the level of openmp diagnostic messages to display. The '-nogen-interfaces' option ask the compiler to do not generate an interface block for each routine defined in the source file.

If you are using the GNU gfortran compiler, the '-fopt-info' option enables all optimization info dumps on stderr.

### LEVEL

Collection of optimization options.

Possible options:

* -O3
* -xHost
* -unroll-aggressive
* -align all
* -fPIC
* -Ofast
* -faggressive-loop-optimizations
* -fno-tree-pre

If you are using the Intel fortran compiler, the '-O3' option means the highest optimization. The '-xHost' option enables the compiler to try to generate the most suitable code for the current computer architecture. The '-unroll-aggressive' option means using aggressive method to unroll the loop structures. The '-align all' option means to align the arrays, structures, etc. The '-fPIC' option means to generate position independent code for the purpose of dynamic link. It should be enabled to compile the python API. But if you want to debug the code, it has to be commented out. Please modify them only if you are an expert of the Intel fortran compiler and you know what you are doing.

If you are using the GNU gfortran compiler, the '-Ofast' option means to optimize the code for speed disregarding exact standards compliance. The '-faggressive-loop-optimizations' option aggressively optimizes loops using language constraints. The '-fno-tree-pre' option will disable the SSA-PRE optimization on trees. In general, if the '-O2', '-O3', or '-Ofast' options are actived, the '-ftree-pre' option will be included automatically. However, based on some benchmark results, we found that the SSA-PRE optimization on trees would lead to strange runtime behaviors (we used the GNU gfortran 5.1.0). So, we decide to close the SSA-PRE optimization on trees explicitly, that is the reason why we need the '-fno-tree-pre' option. The '-fPIC' option means to generate position independent code for the purpose of dynamic link. It should be enabled to compile the python API. However, if you want to debug the code, it has to be commented out. Please modify them only if you are an expert of the GNU gfortran compiler and you know what you are doing.

### MARCH

Used to specify the instruction sets that the current system supports.

Possible options:

* -march=core2
* -march=corei7
* -march=corei7-avx
* -march=core-avx-i
* -march=core-avx2
* -mavx
* -mavx2
* -msse2
* -msse3
* -msse4

If you are using the Intel fortran compiler, the 'core2' option is the safest choice and it works always. But it may be not the best. Please modify it only when you understand what you are doing.

The 'core2' option will generate code for the Intel Core 2 processor family.

The 'corei7' option generates code for processors that support Intel SSE4 efficient accelerated string and text processing instructions. It may also generate code for Intel SSE4 vectorizing compiler and media accelerator, Intel SSE3, SSE2, SSE, and SSSE3 instructions.

The 'corei7-avx' option generates code for processors that support Intel advanced vector extensions (Intel AVX), Intel SSE4.2, SSE4.1, SSE3, SSE2, SSE, and SSSE3 instructions.

The 'core-avx-i' option generates code for processors that support the RDRND instruction, Intel advanced vector extensions (Intel AVX), Intel SSE4.2, SSE4.1, SSE3, SSE2, SSE, and SSSE3 instructions.

The 'core-avx2' option generates code for processors that support Intel advanced vector extensions 2 (Intel AVX2), Intel AVX, SSE4.2, SSE4.1, SSE3, SSE2, SSE, and SSSE3 instructions.

If you are using the GNU gfortran compiler, you can not use the above options to specify the instruction sets. The possible options are listed as follows:

The '-mavx' option supports the Intel MMX, SSE, SSE2, SSE3, SSSE3, SSE4.1, SSE4.2 and AVX built-in functions and code generation.

The '-mavx2' option supports the Intel MMX, SSE, SSE2, SSE3, SSSE3, SSE4.1, SSE4.2, AVX and AVX2 built-in functions and code generation.

The '-msse2' option supports the Intel MMX, SSE and SSE2 built-in functions and code generation.

The '-msse3' option supports the Intel MMX, SSE, SSE2 and SSE3 built -in functions and code generation.

The '-msse4' option supports the Intel MMX, SSE, SSE2, SSE3, SSSE3, SSE4.1 and SSE4.2 built-in functions and code generation.

> NOTE:

> You could simply cat a file (/proc/cpuinfo) on Linux and then glean a lot of information about the CPU. On Mac OS X, you can look at "About This Mac", or do a "sysctl -a hw" and try to look at relevant information in the output.

### FFLAGS

Collection of Fortran compiler options. Do not modify them for ever.

Possible options:

* -c
* $(CPP)
* $(CHECK)
* $(CDUMP)
* $(LEVEL)
* $(MARCH)
* $(GPROF)

### LFLAGS

Collection of linker options. Do not modify them unless you know what you are doing.

Possible options:

* $(OMP)
* $(GPROF)
* -Wl,-no_pie

The '-Wl,-no_pie' option is useful when you are using the Mac Os X system and want to traceback the code (-fbacktrace or -traceback is applied). If you are using the Linux system, you can skip it.

### LIBS

Specify the external libraries. Now the iQIST software package depends on LAPACK and BLAS heavily. To achieve good performance, the highly optimized LAPACK and BLAS implementations are essential. Here we want to recommend the OpenBLAS and Intel MKL.

Possible options:

* -framework Accelerate
* -L/home/lihuang/lapack -llapack -lblas
* -L/opt/intel/mkl/lib -lmkl_core -lmkl_sequential -lmkl_rt

Here we provide three typical choices. (1) In the Mac OS X system, we can use the Apple Accelerate framework. (2) We use the home-built BLAS and LAPACK libraries. Please pay attention to the path. You have to modify it to fit your software environment. (3) We link the iQIST code with the Intel MKL. Please pay attention to the path and the library's name. You have to modify them to fit your software environment. Please see the documentation about Intel MKL for more details.

### F2PY

Specify the command tool for the f2py package which is used to generate the python interface for the iQIST code.

Possible options:

* f2py

> NOTE:

> Sometimes you have to specify correct path for it. In the Mac OS X system, the path for f2py may be:

>
```
/System/Library/Frameworks/Python.framework/Versions/Current/Extras/bin/f2py
```

> You need to specify it here.

### F2PYC

Specify the fortran compiler used by the f2py tool. In general, it is the same with the fortran compiler we used to compile the iQIST code. Do not modify them unless you know what you are doing.

Possible options:

* $(F90)

### F2PYL

Specify the linear algebra library used by the f2py tool to generate the python modules.

Possible options:

* $(LIBS)
* --link-lapack_opt

In general, the '$(LIBS)' option works quite well. But when you are working on the Mac OS X system and using the Apple Accelerate framework as the linear algebra library, you have to use the '--link-lapack_opt' option. As for the details, please check the f2py's documentation.

## Templates

Here we provide a summary list for several make.sys templates. Please choose a suitable one and modify it to satisfy your requirements. Don't forget to use it to overwrite the default make.sys.

### make.sys

See template/macos/make.sys.gfortran. **PLEASE REPLACE IT WITH YOUR OWN VERSION**.

### template/macos/make.sys.standard

Machine   : MacBook Pro 2012

Processor : 2.3 GHz Intel Core i7

Memory    : 8 GB 1600 MHz DDR3

Software  :

```
Mac OS X 10.8.5
MPICH 3.0.3
Intel Fortran Compiler 13.0.0
Intel Math Kernel Library 11.0
Python 2.7.4
```

### template/macos/make.sys.veclib

Machine   : MacBook Pro 2012

Processor : 2.3 GHz Intel Core i7

Memory    : 8 GB 1600 MHz DDR3

Software  :

```
Mac OS X 10.8.5
MPICH 3.0.3
Intel Fortran Compiler 13.0.0
Apple Accelerate/vecLib framework
Python 2.7.4
```

### template/macos/make.sys.openmpi

Machine   : MacBook Pro 2012

Processor : 2.3 GHz Intel Core i7

Memory    : 8 GB 1600 MHz DDR3

Software  :

```
Mac OS X 10.8.5
OpenMPI 1.7.1
Intel Fortran Compiler 13.0.0
Intel Math Kernel Library 11.0
Python 2.7.4
```

### template/macos/make.sys.nompi

Machine   : MacBook Pro 2012

Processor : 2.3 GHz Intel Core i7

Memory    : 8 GB 1600 MHz DDR3

Software  :

```
Mac OS X 10.8.5
Intel Fortran Compiler 13.0.0
Intel Math Kernel Library 11.0
Python 2.7.4
```

### template/macos/make.sys.gfortran

Machine   : MacBook Pro Retina 2013

Processor : 2.4 GHz Intel Core i7

Memory    : 8 GB 1600 MHz DDR3

Software  :

```
Mac OS X 10.10.3
MPICH 3.1.4
GNU gfortran 5.1.0
Apple Accelerate/vecLib framework
Python 2.7.6
```

### template/linux/make.sys.standard

Machine   : Linux cluster (master node)

Processor : 2.0 GHz Intel Xeon E5-2620

Memory    : 32 GB RAM

Software  :

```
Linux kernel 2.6.32-431
OpenMPI 1.6.4
Intel Fortran Compiler 13.1.1
Netlib LAPACK/BLAS 3.4.2
Python 2.7.3
```

### template/linux/make.sys.gfortran

Machine   : Linux cluster (master node)

Processor : 2.0 GHz Intel Xeon E5-2620

Memory    : 32 GB RAM

Software  :

```
Linux kernel 2.6.32-431
OpenMPI 1.8.3
GNU gfortran 4.8.2
Netlib LAPACK/BLAS 3.4.2
Python 2.7.3
```
