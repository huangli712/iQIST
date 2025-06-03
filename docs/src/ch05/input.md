# Standard Input Files

The atomic eigenvalue problem solver, i.e., the **JASMINE** component, is designed to diagonalize the atomic Hamiltonian and generate the necessary *atom.cix* file for the CT-HYB quantum impurity solvers. Since the atomic Hamiltonian is a bit complex, which may include the chemical potential term, onsite Coulomb interaction matrix term, crystal field splitting term, and spin-orbit coupling term, and so on. Some of them are in matrix-form, and the others are vectors. So in order to define the atomic Hamiltonian conveniently, it is essential to introduce multiple input files.

The standard input files supported by the atomic eigenvalue problem solver are as follows:

* [solver.atomic.in](in_atom.md) // Primary configuration file for the **JASMINE** component.
* [atom.cmat.in](in_cmat.md) // Crystal field splitting.
* [atom.emat.in](in_emat.md) // On-site impurity energy level.
* [atom.tmat.in](in_tmat.md) // Transformation matrix.
