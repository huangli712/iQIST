## Standard Output Files

The atomic eigenvalue problem solvers (the **JASMINE** component) will generate a lot of data files at run time. In this section, we will depict their formats and usages in detail.

The data files could be classified as the following kinds:

* **Direct Output**
    * [Terminal output] // The runtime information of the code.
* **Eigensystem Data**
    * [atom.eigval.dat] // Eigenvalues.
    * [atom.eigvec.dat] // Eigenvectors.
* **Miscellaneous Data**
    * [atom.fock.dat] // Fock state.
    * [atom.tmat.dat] // Transformation matrix.
    * [atom.emat.dat] // On-site impurity level.
    * [atom.umat.dat] // On-site Coulomb interaction matrix.
    * [atom.sector.dat] // Configuration of subspace.
* **For The Quantum Impurity Solvers**
    * [solver.umat.in] // On-site Coulomb interaction matrix (only the density-density part).
    * [atom.cix] // All-in-one file for the atomic eigenstates.
