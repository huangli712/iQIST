### Build base library

In essence, the base library includes CSSL and CSML which form the underlying architecture of the iQIST software package. Actually, there are two ways to compile the base library in the compiling system.

Once the compiling is finished, you will find the *libMM.a* file in the *iqist/src/base* directory.

**Method 1**:

```
$ cd iqist/build
$ editor make.sys
$ make base
```

!!! tip

    Here **editor** means any ascii text editor you prefer, such as *vim*.

This is the most common and direct method.

**Method 2**:

You can also go to the directory where the base library exists, and type 'make' command in the terminal.

```
$ cd iqist/src/base
$ make
```

The shortcoming of this method is that it is not straightforward to edit the *make.sys* file once unexpected compiling errors take place.