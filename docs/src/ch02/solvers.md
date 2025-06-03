### Build Quantum Impurity Solvers: Integrally

>
> ehh...
>
> I want to compile all of the components of quantum impurity solvers. Is it possible?
>
> Of course. It is just a piece of cake. Please follow the following guidelines.
>

**Method**:

```shell
$ cd iqist/build
$ editor make.inc
$ make solver
```

!!! note

    Here **editor** means any ascii text editor you prefer.

Finally, you will see a lot of executable programs, such as *ctqmc*, etc., in the *iqist/src* directory. They are the CT-HYB quantum impurity solvers. You can execute them in the terminal directly.

---

### Build Quantum Impurity Solvers: Individually

Sometimes you only need to use some of the quantum impurity solver components of iQIST and you don't want to compile the whole iQIST software package. Are there any tricks to do? Yes! You will not be disappointed with iQIST.

Supposed that what you need is the **NARCISSUS** component:

**Method 1**:

```
$ cd iqist/build
$ editor make.inc
$ make ctseg
```

!!! note

    1. Here **editor** means any ascii text editor you prefer.
    2. If the user want to compile the **MANJUSHAKA** component, please execute
    ```shell
    $ make cthyb
    ```


After a few minutes, you will find an *ctqmc* file in the *iqist/src/ctseg* directory. That is what you need, the executable program for the **NARCISSUS** component. Then you can add it to the system path, and do your great research.

---

**Method 2**:

In this approach, we go to the directory of the **NARCISSUS** component at first. Then 'make' it as usual.

```
$ cd iqist/src/ctseg
$ make
```

Here, we assume that the compiling system was correctly configured.

!!! note

    Now the executable programs will not be copied into the *iqist/build* directory automatically.
