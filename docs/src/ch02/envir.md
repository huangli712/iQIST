# Compiling Environment

!!! note

    The URLs/links provided in this page may be broken sometimes. If you meet the broken links, please contact us.

In order to compile/execute the iQIST software package successfully, please confirm whether the following software/libraries are correctly installed and configured in your system.

**Fortran Compiler**

A modern Fortran compiler supporting ISO Fortran 2003 standard is necessary. We strongly advise you to choose the newest version of Intel Fortran Compiler (ifort) or GNU Fortran Compiler (gfortran). During the development of the iQIST software package, these two compilers are extensively used.

Check the following websites for more details:

* [Intel Fortran Compiler](https://software.intel.com/en-us/fortran-compilers)
* [GNU Fortran Compiler](http://gcc.gnu.org/fortran/)

!!! note

    Besides ifort and gfortran, in principle the iQIST software package can be compiled by any other Fortran 90 compilers, but we can not guarantee it.

---

**Parallel Environment**

We use the message passing interface (MPI) to implement the parallelism in the iQIST software package. Therefore, in order to accelerate the quantum impurity solvers, you have to install a MPI implementation in your system, and setup it to work with the Fortran compiler you choose correctly. There are many MPI implementations in the market provided by various vendors. But we recommend to use the newest versions of MPICH or Openmpi.

Check the following websites for more details:

* [Message Passing Interface Standard](http://mpi-forum.org)
* [MPICH](http://www.mpich.org)
* [Openmpi](http://www.open-mpi.org)

Some features in the iQIST software package have been optimized using the OpenMP multi-thread technology. The chosen Fortran compilers should support the OpenMP standard version 3.0 or later version.

Check the following websites for more details:

* [OpenMP](http://openmp.org/wp/)
* [OpenMP Support In GNU GCC](https://gcc.gnu.org/projects/gomp/)
* [OpenMP Support In Intel compilers](https://software.intel.com/en-us/intel-parallel-studio-xe/details)

---

**Linear Algebra Library**

The iQIST software package depends on the BLAS and LAPACK libraries heavily. Then you have to ensure that a compatible math library is installed in your system. We recommend to use the Intel Math Kernel Library (Intel MKL). If you are using the MacOS system, the Accelerate Framework may be an alternate choice. Keep in mind that the running efficiency of iQIST is always benefited from the highly optimized linear algebra library.

Check the following websites for more details:

* [BLAS](http://www.netlib.org/blas/)
* [LAPACK](http://www.netlib.org/lapack/)
* [Intel MKL](https://software.intel.com/en-us/intel-mkl/)
* [OpenBLAS](http://www.openblas.net)
* [Accelerate Framework](https://developer.apple.com/library/mac/documentation/Performance/Conceptual/vecLib/index.html)

---

**Python Environment**

There are a few python scripts in the iQIST software package. In order to run them, you need a Python interpreter. We recommend to install Python 3.8+, and the numpy, scipy, matplotlib packages must be installed as well.

Check the following websites for more details:

* [Python](https://www.python.org)
* [numpy](http://www.numpy.org)
* [scipy](http://www.scipy.org)
* [matplotlib](http://matplotlib.org)

!!! tip

    In the near future, all the scripts in the iQIST software package should be rewritten by using the Julia language. At that time, the Python environment is not necessary any more.

---

**Flink Library**

The Flink library is a collection of Fortran modules and subroutines for scientific computing. The iQIST software package depends on it.

Check the following websites for more details:

* [Source](https://github.com/huangli712/Flink)
* [Documentation](https://huangli712.github.io/projects/flink/index.html)

Please download the source codes of this library from github, and then uncompress it into the same directory with the iQIST software package. Next, compile it.

!!! tip

    The iQIST software package adopts an automatic algorithm to determine the directory that contains the Flink library. So, to make this algorithm works correctly, please make sure that the source codes of the Flink library are placed at *your_path/Flink* (**case sensitive**).

    If this automatic algorithm fails, please specify the environment variable **FLINK** manually. It should be associated with the *your_path/Flink/src* folder. For example:
    ```shell
    $ export FLINK=/home/lihuang/Flink/src
    ```
