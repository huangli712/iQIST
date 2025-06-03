# Summary

Well, maybe you don't like to compile the whole iQIST at one step. You want more controls about the compiling procedure, or you just don't want to waste time on compiling the components, or you even don't indend to compile/have. OK, the compiling system of iQIST provide an alternative way to meet your requirements.

Anyway, the compiling procedure can be split into three successive steps:

* [Build Quantum Impurity Solvers](solvers.md)
* [Build Atomic Eigenvalue Problem Solver](atomic.md)
* [Build Documentation](docs.md)

These steps are optional. You can do the above jobs in sequence or in a random manner.

!!! tip

    The iQIST software package depends on the [Flink](https://github.com/huangli712/Flink) library, so please make sure that the Flink library is correctly installed and the environment variable **FLINK** is correctly setup. See [Compiling Environment](envir.md) for more details.
