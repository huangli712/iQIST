## Compiling system

In this section, we will illustrate the compiling system of the iQIST software package. In fact, it is completely based on the well-known GNU GCC tool-chain. 

The compiling system is in the *iqist/build* directory which includes the following files/folders:

* *Makefile*
* *make.sys*
* *build.md*
* *template/linux*
* *template/macos*
* *template/tianhe*
* *x_setup.sh*
* *x_clean.sh*

Next, we will explain them in detail.

### *Makefile*

It is the central file of the compiling system. Usually, when you type "make target" command in the terminal, the *make* utility will parse this file and then apply the rules defined in it to build the target. Now suppose that you are in the *iqist/build* directory, you can execute

```
$ make help
```

to see the available targets, and input

```
$ make help-more
```

for more details.

> NOTE: 

> **DO NOT** touch this file by yourself even you are very familiar with the iQIST software package.

### *make.sys*

We design this file to configure the compiling system. In this file, we have to specify the Fortran compiler, the parallel environment, the linear algebra library, and the target hardware architecture, etc. 

> NOTE: 

> The *make.sys* file is system-dependent, i.e., you have to modify it to fit your systems. Or else, the compiling will fail definitely.

In the *iqist/build/template* directory, we provide a lot of *make.sys* files for various systems we have tested. You can consider them as references. We strongly recommend the users go through the *make.sys* file carefully and check whether the settings are correct before they start to compile the iQIST.

### *build.md*

In this file, we provide some hints/explanations about the *make.sys* template files.

### *template*

In this folder, we provide a few *make.sys* files which are designed and tested for various systems. You can modify them to fit your requirement, and then copy it to override the default *make.sys* file (i.e., the *iqist/build/make.sys* file).

* *template/linux*, containing template *make.sys* files for Linux system.
* *template/macos*, containing template *make.sys* files for Mac OS X system.
* *template/tianhe*, containing template *make.sys* files for the TianHe-1 supercomputer system in Tianjin, China.

### *x_setup.sh*

When the compiling finishes, you can execute this script in the terminal to copy all of the (available) built targets into the *build* directory (i.e., the *iqist/build* directory).

```
$ ./x_setup.sh
```

### *x_clean.sh*

Sometimes, you come across the idea of deleting the executable programs in the *iqist/build* directory. The *x_clean.sh* script can help you to clean them and give you a fresh start.

```
$ ./x_clean.sh
```