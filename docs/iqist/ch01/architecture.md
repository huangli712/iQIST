## Software architecture

To solve a quantum impurity model is not a straightforward job. Besides the necessary quantum impurity solvers, we need several auxiliary programs or tools. The iQIST is an *all-in-one* software package, which can be used to solve a broad range of quantum impurity problems. Thus, it is not surprising that iQIST is a collection of various codes and scripts. The core components contain about 50000 lines of code. 

![layer image](../figure/layer.png)

**Figure** | The hierarchical structure of the iQIST software package. Note that in the component layer, not all of the components are listed due to space limitations. See the main text for detailed explanations.

The software architecture of iQIST is slightly involved. In the above figure, we use a layer model to illustrate it. Next, we will explain these layers one by one.

### Operating system layer

The bottom layer is the operating system (OS). In principle, the iQIST is OS-independent. It can run properly on top of Unix/Linux, Mac OS X, FreeBSD, and Windows. 

### System layer

The second layer is the system layer, which contains highly optimized linear algebra math libraries (such as BLAS and LAPACK) and parallelism supports (such as MPI and OpenMP). 

### Service layer

The third layer is the service layer. In this layer, we implemented some commonly used modules and subroutines. They are named as common service subroutine library (CSSL) and common service module library (CSML), respectively. They provide a useful interface between the system layer and the component layer and facilitate the development of core components. The features of CSSL and CSML include basic data structures (stack and linked list), random number generators, spare matrix manipulations, linear algebra operations, string processing, linear interpolation, numerical integration, and fast Fourier transformation (FFT), etc. 

### Component layer

The core part of iQIST software package is in the fourth layer -- the component layer -- which contains various impurity solvers as shown before. At present, iQIST contains ten different components, including 

* **AZALEA**
* **GARDENIA**
* **NARCISSUS**
* **BEGONIA**
* **LAVENDER**
* **PANSY**
* **MANJUSHAKA**
* **DAISY**
* **JASMINE**
* **HIBISCUS**

Here, **AZALEA**, **GARDENIA**, **NARCISSUS**, **BEGONIA**, **LAVENDER**, **PANSY**, and **MANJUSHAKA** are all CT-HYB components, and **DAISY** is a HF-QMC impurity solver component. **JASMINE** is an atomic eigenvalue solver. **HIBISCUS** is a collection of several pre- and post-processing tools, including maximum entropy method, stochastic analytical continuation, Pade approximation, and Kramers-Kronig transformation, etc. For more details about these components, please consult the following chapters. 

### Interface layer

The top layer is the interface layer or user layer. On the one hand, users can execute the iQIST's components directly as usual. On the other hand, they can also invoke iQIST's components from other languages. The role of iQIST's components becomes a library or subroutine. To achieve this goal, in the interface layer, we offer the Fortran/Python language bindings for most of the iQIST's components, so that the users can develop their own codes on top of iQIST and consider it as a computational engine in black box.