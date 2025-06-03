# Summary

The quantum impurity solvers requires a few input files to setup the impurity models, configure the impurity solvers, and control the DMFT iterations, etc. In principle, the quantum impurity solvers can run without any input files. But if the input files are available, the parameters read from the input files will **override** the default ones.

The standard input files supported by the CT-HYB quantum impurity solvers are as follows:

* [solver.ctqmc.in](in_ctqmc.md) // Configuration file for CT-QMC quantum impurity solvers.
* [solver.umat.in](in_umat.md) // Coulomb interaction matrix.
* [solver.eimp.in](in_eimp.md) // Impurity level and crystal-field splitting.
* [solver.anydos.in](in_anydos.md) // Density of states from various lattice or tight-binding models.
* [solver.ktau.in](in_ktau.md) // Screening interaction.
* [atom.cix](in_atom.md) // Atomic eigenstates data.
