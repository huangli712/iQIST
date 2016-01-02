# iQIST (Interacting Quantum Impurity Solver Toolkit)

Here we provide a summary list for several make.sys templates. Please choose a suitable one and modify it to satisfy your requirements. Don't forget to use it to override the default make.sys.

## make.sys

See template/macos/make.sys.gfortran. **PLEASE REPLACE IT WITH YOUR OWN VERSION**.

## template/macos/make.sys.standard

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

## template/macos/make.sys.veclib

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

## template/macos/make.sys.openmpi

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

## template/macos/make.sys.nompi

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

## template/macos/make.sys.gfortran

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

## template/linux/make.sys.standard

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

## template/linux/make.sys.gfortran

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

## template/tianhe/make.sys.standard

Machine   : TianHe-1 supercomputer (login node)

Processor : 2.93 GHz Intel Xeon X5670

Memory    : 23 GB RAM

Software  : 

```
Linux kernel 2.6.32-358
MPICH 3.0.4
Intel Fortran Compiler 13.0.0
Intel Math Kernel Library 11.0
Python 2.6.6
```

## template/tianhe/make.sys.openblas

Machine   : TianHe-1 supercomputer (login node)

Processor : 2.93 GHz Intel Xeon X5670

Memory    : 23 GB RAM

Software  : 

```
Linux kernel 2.6.32-358
MPICH 3.0.4
Intel Fortran Compiler 13.0.0
OpenBLAS 0.2.3
Python 2.6.6
```
