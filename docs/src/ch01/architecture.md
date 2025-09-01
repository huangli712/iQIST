# Software Architecture

To solve a quantum impurity model is not a straightforward job. Besides the necessary quantum impurity solvers, we need several auxiliary programs or tools. The iQIST is an *all-in-one* software package, which can be used to solve a broad range of quantum impurity problems. Thus, it is not surprising that iQIST is a collection of various codes and scripts. The core components contain about 50000 lines of code.

The software architecture of iQIST is slightly involved. We just use a layer model to illustrate it. Next, we will explain these layers one by one.

---

**Operating System Layer**

> The bottom layer is the operating system (OS). In principle, the iQIST is OS-independent. It can run properly on top of Unix/Linux, Mac OS X, FreeBSD, and Windows.

**System Layer**

> The second layer is the system layer, which contains highly optimized linear algebra math libraries (such as BLAS and LAPACK) and parallelism supports (such as MPI and OpenMP).

**Service Layer**

> The third layer is the service layer. In this layer, we implemented some commonly used modules and subroutines. They are provided in a separate library, namely [Flink](https://github.com/huangli712/Flink). The Flink library provides an useful interface between the system layer and the component layer and facilitate the development of core components. The features of the Flink library include basic data structures (stack and linked list), random number generators, spare matrix manipulations, linear algebra operations, string processing, linear interpolation, numerical integration, and fast Fourier transformation (FFT), etc.

**Component Layer**

> The core part of iQIST software package is in the fourth layer -- the component layer -- which contains various quantum impurity solvers as shown before. At present, iQIST contains three different components, including
>
> * **NARCISSUS** component
> * **MANJUSHAKA** component
> * **JASMINE** component
>
> Here, **NARCISSUS** and **MANJUSHAKA** are CT-HYB quantum impurity solvers, and **JASMINE** is an atomic eigenvalue solver. For more details about these components, please consult the following chapters.
