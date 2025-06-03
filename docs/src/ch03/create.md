# Prepare Input Files

There are three ways to generate the necessary inpute files for the quantum impurity solvers in the iQIST software package.

----

**Method 1**:

One is to prepare the input files manually.

Usually, the necessary input file is *solver.ctqmc.in*. It is an ascii text file actually. You can use any text editor to create and edit it in principle. As for the file format of *solver.ctqmc.in*, please read:

* [solver.ctqmc.in](../ch04/in_ctqmc.md) // Configuration file for CT-QMC quantum impurity solvers.

Besides the *solver.ctqmc.in* file, some quantum impurity solvers also require the *atom.cix* file as input. You can use the **JASMINE** component to generate it. The necessary input file for the **JASMINE** component is *solver.atomic.in*. For more details about it, please read:

* [solver.atomic.in](../ch05/in_atom.md) // Configuration file for the atomic eigenvalue problem solver.

---

**Method 2**:

Another approach is to use the *u\_ctqmc.py* and *u\_atomic.py* to generate the *solver.ctqmc.in* and *solver.atomic.in* files, respectively. In this approach, you have to write some python scripts. Don't worry about it. It is a trivial task. Please see

* [Scripts] // Python scripts.

for more details.

---

**Method 3**:

The [ZenGui](https://github.com/huangli712/ZenGui) code can be used to generate the *solver.ctqmc.in* and *solver.atomic.in* files.

---

As for the other input files not mentioned here, please try to generate them by yourself. The detailed format descriptions can be seen in the following sections:

* [Standard Input Files](../ch04/input.md) // Input stuffs for CT-HYB quantum impurity solvers.
* [Standard Input Files](../ch05/input.md) // Input stuffs for atomic eigenvalue problem solver.
