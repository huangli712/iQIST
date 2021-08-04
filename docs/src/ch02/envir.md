## Compiling environment

!!! note 

    The URLs/links provided in this page may be broken sometimes. If you meet the broken links, please contact us.

In order to compile/execute the iQIST software package successfully, please confirm whether the following software/libraries are correctly installed and configured in your system.

### Fortran Compiler

A modern Fortran compiler supporting ISO Fortran 2003 standard is necessary. We strongly advise you to choose the newest version of Intel Fortran Compiler (ifort) or GNU Fortran Compiler (gfortran). During the development of the iQIST software package, these two compilers are extensively used. 

Check the following websites for more details:

* [Intel Fortran compiler](https://software.intel.com/en-us/fortran-compilers)
* [GNU Fortran compiler](http://gcc.gnu.org/fortran/)

!!! note 

    Besides ifort and gfortran, in principle the iQIST software package can be compiled by any other Fortran 90 compilers, but we can not guarantee it.

### Parallel environment

We use the message passing interface (MPI) to implement the parallelism in the iQIST software package. Therefore, in order to accelerate the quantum impurity solvers, you have to install a MPI implementation in your system, and setup it to work with the Fortran compiler you choose correctly. There are many MPI implementations in the market provided by various vendors. But we recommend to use the newest versions of MPICH or Openmpi.

Check the following websites for more details:

* [Message passing interface standard](http://mpi-forum.org)
* [MPICH](http://www.mpich.org)
* [Openmpi](http://www.open-mpi.org)

Some features in the iQIST software package have been optimized using the OpenMP multi-thread technology. The chosen Fortran compilers should support the OpenMP standard version 3.0 or later version.

Check the following websites for more details:

* [OpenMP](http://openmp.org/wp/)
* [OpenMP support in GNU GCC](https://gcc.gnu.org/projects/gomp/)
* [OpenMP support in Intel compilers](https://software.intel.com/en-us/intel-parallel-studio-xe/details)

### Linear algebra library

The iQIST software package depends on the BLAS and LAPACK libraries heavily. Then you have to ensure that a compatible math library is installed in your system. We recommend to use the Intel Math Kernel Library (Intel MKL). If you are using the Mac OS X system, the Accelerate Framework may be an alternate choice. Keep in mind that the running efficiency of iQIST is always benefited from the highly optimized linear algebra library.

Check the following websites for more details:

* [BLAS](http://www.netlib.org/blas/)
* [LAPACK](http://www.netlib.org/lapack/)
* [Intel MKL](https://software.intel.com/en-us/intel-mkl/)
* [OpenBLAS](http://www.openblas.net)
* [Accelerate Framework](https://developer.apple.com/library/mac/documentation/Performance/Conceptual/vecLib/index.html)

### Python environment

There are a few python scripts in the iQIST software package. In order to run them, you need a Python interpreter. We recommend to install Python 2.6+, and the numpy, scipy, f2py, matplotlib packages must be installed as well.

Check the following websites for more details:

* [Python](https://www.python.org)
* [numpy](http://www.numpy.org)
* [scipy](http://www.scipy.org)
* [f2py](http://www.f2py.com)
* [matplotlib](http://matplotlib.org)