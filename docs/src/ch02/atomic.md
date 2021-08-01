### Build atomic eigenvalue problem solver

The **JASMINE** component is the built-in atomic eigenvalues solver. If you want to use the **BEGONIA**, **LAVENDER**, **PANSY**, **MANJUSHAKA**, **CAMELLIA** components, the **JASMINE** component is the necessary pre-processing tool. You can employ the following commands to compile it.

**Method 1**:
```
$ cd iqist/build
$ editor make.sys
$ make base
$ make capi
$ make jasmine
$ ./x_setup.sh
```

> NOTE: 

> 1. Here **editor** means any ascii text editor you prefer.

> 2. If the base library and the application programming interfaces were already compiled successfully, then you can skip the following commands:
```
$ make base
$ make capi
```

This is the most common and direct method. After a few minutes, you will find a *jasmine.x* file in the *iqist/build* directory. That is what you need, the executable program for the **JASMINE** component. After that you can add it to the system path, and do some tests.

**Method 2**:

In this approach, we go to the directory of the **JASMINE** component at first. Then 'make' it as usual.

```
$ cd iqist/src/tools/jasmine
$ make
```

Here, we assume that the compiling system was correctly configured, the base library and the application programming interfaces were already compiled successfully.

> NOTE: 

> Now the executable program will not be copied into the *iqist/build* directory. You have to go to the *iqist/build* directory, and execute
```
$ ./x_setup.sh
```