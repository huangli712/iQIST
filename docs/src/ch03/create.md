## Prepare input files

There are two ways to generate the necessary inpute files for the quantum impurity solvers in the iQIST software package. 

**Method 1**:

One is to prepare the input files manually. 

Usually, the necessary input file is *solver.ctqmc.in*. It is an ascii text file actually. You can use any text editor to create and edit it in principle. As for the file format of *solver.ctqmc.in*, please read:

* [solver.ctqmc.in](../ch04/in_ctqmc.md) // Configuration file for CT-QMC impurity solvers.
* [solver.hfqmc.in](../ch04/in_hfqmc.md) // Configuration file for HF-QMC impurity solvers.

Besides the *solver.ctqmc.in* file, some quantum impurity solvers also require the *atom.cix* file as input. You can use the **JASMINE** component to generate it. The necessary input file for the **JASMINE** component is *atom.config.in*. For more details about it, please read:

* [atom.config.in](../ch06/in_atom.md) // Configuration file for the **JASMINE** component.

**Method 2**:

Another approach is to use the *u\_ctqmc.py*, *u\_hfqmc.py*, and *u\_atomic.py* to generate *solver.ctqmc.in*, *solver.hfqmc.in*, and *atom.config.in* files, respectively. In this approach, you have to write some python scripts. Don't worry about it. It is a trivial task. Please see

* [Scripts](../ch07/script.md) // **HIBISCUS**/scripts codes.

for more details.

---

As for the other input files not mentioned here, please try to generate them by yourself. The detailed format descriptions can be seen in the following sections:

* [Standard input files](../ch04/input.md) // Input stuffs for CT-HYB/HF-QMC impurity solvers.
* [Standard input files](../ch06/input.md) // Input stuffs for atomic eigenvalue problem solver.