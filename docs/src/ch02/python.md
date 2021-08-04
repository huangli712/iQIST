### Build modules for Python

In the previous section, we introduce how to build Fortran library for a given component. In the section, we will introduce how to generate the required Python module.

!!! note 

    In order to generate the Python modules, the f2py package is of course necessary. Please check the following website for more details:
    ```
    www.f2py.com
    ```

Let's use the **AZALEA** component as example again.

**Method 1**:

```
$ cd iqist/build
$ editor make.sys
$ make base
$ make capi
$ make azalea-pylib
```

!!! note 

    1. Here **editor** means any ascii text editor you prefer.
    2. **azalea** can be any other component's name, such as **begonia**, **pansy**, **gardenia**, etc.
    3. If the base library and the application programming interfaces were already compiled successfully, then you can skip the following commands:
    ```
    $ make base
    $ make capi
    ```

After a few minutes, you will find *pyiqist.so* in the *iqist/src/ctqmc/azalea* directory. That is what you need, the Python module for the **AZALEA** component. Then you can add it to the system path, and do your great research.

**Method 2**:

In this approach, we go to the directory of the **AZALEA** component at first. Then 'make' it as usual.

```
$ cd iqist/src/ctqmc/azalea
$ make pylib
```

Here, we assume that the compiling system was correctly configured, the base library and the application programming interfaces were already compiled successfully.