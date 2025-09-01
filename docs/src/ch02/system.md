# Compiling System

In this section, we will illustrate the compiling system of the iQIST software package. In fact, it is completely based on the well-known GNU GCC tool-chain.

The compiling system is in the *iqist/build* directory which includes the following files/folders:

* *Makefile*
* *make.inc*

Next, we will explain them in detail.

---

**Makefile**

It is the central file of the compiling system. Usually, when you type "make target" command in the terminal, the *make* utility will parse this file and then apply the rules defined in it to build the target. Now suppose that you are in the *iqist/build* directory, you can execute

```
$ make help
```

to see the available targets, and input

```
$ make help-more
```

for more details.

!!! warning

    **DO NOT** touch this file by yourself even you are very familiar with the iQIST software package.

---

**make.inc**

We design this file to configure the compiling system. In this file, we have to specify the Fortran compiler, the parallel environment, the linear algebra library, and the target hardware architecture, etc.

!!! warning

    The *make.inc* file is system-dependent, i.e., you have to modify it to fit your systems. Or else, the compiling will fail definitely.

We strongly recommend the users go through the *make.inc* file carefully and check whether the settings are correct before they start to compile the iQIST software package. See [About make.inc](inc.md) for more information about this file.
