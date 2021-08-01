## Build iQIST at multiple steps

Well, maybe you don't like to compile the whole iQIST at one step. You want more controls about the compiling procedure, or you just don't want to waste time on compiling the components, or you even don't indend to compile/have. OK, the compiling system of iQIST provide an alternative way to fit your requirements.

Anyway, the compiling procedure can be split into three successive steps:

* [Build base library](base.md)

This step is mandatory, because all of the components of iQIST depend on it.

* [Build application programming interfaces](apis.md)

This step is optional. But all of the components about quantum impurity solvers and atomic eigenvalues solvers depend on it. So we strongly recommend to build it.

* Build selected components

Now you can feel free to compile what you want to have in the iQIST software package. More specially, in this step you can do the following jobs as you wish:

1. [Build quantum impurity solvers](solvers.md)
2. [Build applications](apps.md)
3. [Build atomic eigenvalues solver](atomic.md)
4. [Build auxiliary tools](tools.md)
5. [Build libraries for Fortran](fortran.md)
6. [Build modules for Python](python.md)

This step is also optional. You can do the above jobs in sequence or in a random manner.