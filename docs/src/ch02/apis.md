### Build application programming interfaces

The application programming interfaces are designed for the components of quantum impurity solvers and atomic eigenvalue problem solver. So if you would like to build these components, building application programming interfaces are prerequisite. There are two ways to compile it in the compiling system.

Once the compiling procedure is finished, you will find some `*.o` files in the *iqist/src/capi* directory.

**Method 1**:

```
$ cd iqist/build
$ editor make.sys
$ make capi
```

> NOTE: 

> Here **editor** means any ascii text editor you prefer, such as *vim*.

This is the most common and direct method.

**Method 2**:

You can also go to the directory where the application programming interfaces exist, and type 'make' command in the terminal.

```
$ cd iqist/src/capi
$ make
```