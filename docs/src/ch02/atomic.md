# Build Atomic Eigenvalue Problem Solver

The **JASMINE** component is the built-in atomic eigenvalues solver. If you want to use the **MANJUSHAKA** component, the **JASMINE** component is the necessary pre-processing tool. You can employ the following commands to compile it.

---

**Method 1**:

```shell
$ cd iqist/build
$ editor make.inc
$ make atomic
```

!!! note

    Here **editor** means any ascii text editor you prefer.

This is the most common and direct method. After a few minutes, you will find a *atomic* file in the *iqist/src/atomic* directory. That is what you need, the executable program for the **JASMINE** component. After that you can add it to the system path, and do some tests.

---

**Method 2**:

In this approach, we go to the directory of the **JASMINE** component at first. Then 'make' it as usual.

```shell
$ cd iqist/src/atomic
$ make
```

Here, we assume that the compiling system was correctly configured.

!!! note

    Now the executable program will not be copied into the *iqist/build* directory automatically.
