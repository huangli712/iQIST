### Build libraries for Fortran

The quantum impurity solvers in the iQIST can not only run as standalone programs, but also can be called from the external Fortran or Python programs. In the latter case, the quantum impurity solvers can be considered as computational engines.

Thus, in this section we will show you how to compile the quantum impurity solvers as Fortran libraries, instead of executable programs.

!!! note

    Not only the quantum impurity solvers, but also the atomic eigenvalues solver can be compiled into Fortran libraries and Python modules.

Supposed that we would like to compile the **AZALEA** component.

**Method 1**:

```
$ cd iqist/build
$ editor make.sys
$ make base
$ make capi
$ make azalea-lib
```

!!! note

    1. Here **editor** means any ascii text editor you prefer.
    2. **azalea** can be any other component's name, such as **begonia**, **pansy**, **gardenia**, etc.
    3. If the base library and the application programming interfaces were already compiled successfully, then you can skip the following commands:
    ```sh
    $ make base
    $ make capi
    ```

After a few minutes, you will find *libctqmc.a* in the *iqist/src/ctqmc/azalea* directory. That is what you need, the Fortran library for the **AZALEA** component. Then you can add it to the system path, and start with your research. We are looking forward to your great finding!

**Method 2**:

In this approach, we go to the directory of the **AZALEA** component at first. Then 'make' it as usual.

```
$ cd iqist/src/ctqmc/azalea
$ make lib
```

Here, we assume that the compiling system was correctly configured, the base library and the application programming interfaces were already compiled successfully.

**Method 3**:

There is a trick to compile Fortran libraries for all of the quantum impurity solvers and atomic eigenvalues solvers at the same time:

```
$ cd iqist/build
$ editor make.sys
$ make lib
```

Here the base library and the application programming interfaces would be compiled implicitly.