# Atomic eigenvalue problem solver

The atomic eigenvalue problem solver is an important preprocessing tool for the general matrix version CT-HYB quantum impurity solvers. The purpose of it is to diagonalize the atomic Hamiltonian to obtain the eigenvalues and eigenvectors, which are essential input for some CT-HYB impurity solvers. In this chapter, we will introduce the **JASMINE** component, which is the only atomic eigenvalue problem solver in the iQIST software package. It should be used together with the **BEGONIA**, **LAVENDER**, **CAMELLIA**, **PANSY**, and **MANJUSHAKA** components.

The main topics are as follows:

* [Standard input files](input.md) // Input stuffs.
* [Standard output files](output.md) // Output stuffs.
* [Parameters](parameters.md) // Really a reference manual.

Maybe you are interested in the following related topics:

* [Getting started](../ch03/README.md) // To teach you how to setup and use the iQIST software package.
* [Quantum Monte Carlo impurity solvers](../ch04/README.md) // Another reference manual.
* [Auxiliary tools](../ch07/README.md) // The companion of quantum impurity solvers.
* [Application programming interfaces](../ch08/README.md) // Write your own Python/Fortran codes.
* [iQIST in action](../ch09/README.md) // Some lightweight tutorials.
* [Inside iQIST](../ch10/README.md) // Unleashing the secrets inside the iQIST.

**See also**:

Besides the above links, you can find many useful examples/cases under the *iqist/working/tools/jasmine* directory.