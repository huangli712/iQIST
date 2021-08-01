### Build quantum impurity solvers: integrally

ehh...

I want to compile all of the components of quantum impurity solvers. Is it possible?

Of course. It is just a piece of cake. Please follow the following guidelines.

**Method**:

```
$ cd iqist/build
$ editor make.sys
$ make solver
$ ./x_setup.sh
```

> NOTE: 

> 1. Here **editor** means any ascii text editor you prefer.

> 2. The base library and the application programming interfaces should be compiled implicitly.

Finally, you will see a lot of executable programs, such as *azalea.x*, *begonia.x*, *daisy.x*, etc., in the *iqist/build* directory. They are the CT-HYB and HF-QMC impurity solvers. You can execute them in the terminal directly.

### Build quantum impurity solvers: individually

Sometimes you only need to use some of the quantum impurity solver components of iQIST and you don't want to compile the whole iQIST software package. Are there any tricks to do? Yes! You will not be disappointed with iQIST.

Supposed that what you need is the **AZALEA** component:

**Method 1**:

```
$ cd iqist/build
$ editor make.sys
$ make base
$ make capi
$ make azalea
$ ./x_setup.sh
```

> NOTE: 

> 1. Here **editor** means any ascii text editor you prefer.

> 2. **azalea** can be any other component's name, such as **begonia**, **pansy**, **gardenia**, etc.

> 3. If the base library and the application programming interfaces were already compiled successfully, then you can skip the following commands:
```
$ make base
$ make capi
```

After a few minutes, you will find an *azalea.x* file in the *iqist/build* directory. That is what you need, the executable program for the **AZALEA** component. Then you can add it to the system path, and do your great research.

**Method 2**:

In this approach, we go to the directory of the **AZALEA** component at first. Then 'make' it as usual.

```
$ cd iqist/src/ctqmc/azalea
$ make
```

Here, we assume that the compiling system was correctly configured, the base library and the application programming interfaces were already compiled successfully.


> NOTE: 

> Now the executable programs will not be copied into the *iqist/build* directory. You have to go to the *iqist/build* directory, and execute
```
$ ./x_setup.sh
```