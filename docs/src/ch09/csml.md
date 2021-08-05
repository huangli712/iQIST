### Common service module library

In order to facilitate the development of scientific software, we gather some numerical algorithms, datatypes, and constants, and implement them as common Fortran modules. They are widely used in the iQIST software package and the other private projects. The collection of these modules is the so-called common service module library, abbreviated **CSML** in this manual.

The source codes of the CSML are stored in the *iqist/src/base* directory. You can compile them via the following commands:

```
$ cd iqist/build
$ make base
```

or

```
$ cd iqist/src/base
$ make
```

The CSML contains many files (modules). Next we will introduce them one by one, so that the developers can be familiar with their functionality and limitations. Except where stated explicitly, these modules could be used individually.

!!! note

    This is not a API reference documentation. You should read the source codes by yourself to obtain detailed information.

**m_constants.f90**

```fortran
module constants
```
Some selected physical and numerical constants are provided in this modules, such as ``\pi`` and ``1i``.

**m_leja.f90**

```fortran
module leja
```
In this module, the Newton-Leja polynomial interpolation algorithm is implemented, which is used to evaluate 

```math
|\psi_f\rangle = e^{-H\tau} | \psi_i \rangle.
```

Now this module is only useful for the **CAMELLIA** component.

**m_linkedlist.f90**

```fortran
module linkedlist
```
As is well known, the linked list is a very useful datatype. In this module, a generic linked list is implemented. In the iQIST software package, this module is used to implement a parser for the configuration files.

**m_mpi.f90**

```fortran
module mmpi
```

In this module, some important MPI calls are decorated to simpler forms. 

**m_parser.f90**

```fortran
module parser
```

**m_skynet.f90**

```fortran
module skynet
```

**m_sparse.f90**

```fortran
module sparse
```

**m_spring.f90**

```fortran
module spring
```

**m_stack.f90**

```fortran
module stack
```