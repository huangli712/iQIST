## Standard output files

The atomic eigenvalue problem solvers (the **JASMINE** component) will generate a lot of data files at run time. In this section, we will depict their formats and usages in detail.

The data files could be classified as the following kinds:

* **Direct output**
    * [Terminal output](out_term.md) // The runtime information of the code.
* **Eigensystem data**
    * [atom.eigval.dat](out_val.md) // Eigenvalues.
    * [atom.eigvec.dat](out_vec.md) // Eigenvectors.
* **Miscellaneous data**
    * [atom.fock.dat](out_fock.md) // Fock state.
    * [atom.tmat.dat](out_tmat.md) // Transformation matrix.
    * [atom.emat.dat](out_emat.md) // On-site impurity level.
    * [atom.umat.dat](out_umat2.md) // On-site Coulomb interaction matrix.
    * [atom.sector.dat](out_sector.md) // Configuration of subspace.
* **For the quantum impurity solvers**
    * [solver.umat.in](out_umat1.md) // On-site Coulomb interaction matrix (only the density-density part).
    * [atom.cix](out_cix.md) // All-in-one file for the atomic eigenstates.