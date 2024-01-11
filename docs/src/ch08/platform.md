# Development platform

**Programming languages**

The main part of the iQIST software package was developed with the modern Fortran 90 language. We extensively used advanced language features in the Fortran 2003/2008 standard such as an object oriented programming style (polymorphic, inheritance, and module, etc.) to improve the readability and re-usability of the source codes. The compilers are fixed to the Intel Fortran compiler and GNU GCC gfortran. We can not guarantee that the iQIST can be compiled successfully with other Fortran compilers. Some auxiliary scripts, pre- and post-processing tools are written using the Python language and Bash shell scripts. These scripts and tools act like a glue. They are very flexible and very easily extended or adapted to deal with various problems. In order to avoid incompatibilities, our Python codes only run on the Python 2.x runtime environment.

**Version control system**

We use *git* as the version control tool, and the source codes are hosted in some remote servers. Many thanks to Bitbucket and Github for providing free code repository services. The developers pull the source codes from the server into their local machines, and then try to improve them. Once the development is done, the source codes can be pushed back to the server and merged with the master branch. Then the other developers can access them and use them immediately to start further developments. The members of our developer team and the other users can access the public code repository anywhere and anytime.

The official code repositories are as follows:

* [Github](https://github.com/huangli712/iQIST)
