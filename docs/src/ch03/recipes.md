# iQIST Recipes

Using the iQIST software package is quite easy. Next we will show you the standard workflow for using the iQIST software package.

---

### Choose Suitable Component

At first, there are several CT-HYB quantum impurity solvers in the package. Their features and efficiency are somewhat different. Thus, it is the user's responsibility to choose suitable CT-HYB components to deal with the quantum impurity problem at hand.

**See also**:

* [Components](../ch01/components.md) // What does the iQIST software package include?
* [Features](../ch01/feature.md) // Is the required feature supported by the selected component?
* [How To Choose Suitable Quantum Impurity Solvers?](../ch04/choose.md) // A guideline.

---

### Design The Programs And Scripts

Second, the iQIST software package is in essence a computational engine, so users have to write scripts or programs to execute the selected CT-HYB quantum impurity solver. For example, if the users want to conduct CT-HYB/DMFT calculations, in principle they must implement the DMFT self-consistent equation by themselves.

There is a bonus. When the users want to study the Hubbard model on bethe lattice using the single-site dynamical mean-field theory, or solve the Anderson impurity models in one-shot mode, it is possible to execute the quantum impurity solvers directly without any additional scripts or programs.

---

### Prepare The Input Files

Third, an important task is to prepare proper input data for the selected CT-HYB quantum impurity solver. The optional inputs for the CT-HYB quantum impurity solver are the hybridization function [``\Delta(i\omega_n)``], impurity energy level (``E_{\alpha\beta}``), interaction parameters (``U``, ``J``, ``\lambda``, and ``\mu``), etc. If users do not provide them to the quantum impurity solver, it will use the default settings automatically. Specifically, if the Coulomb interaction matrix is general or the spin-orbit coupling is considered, users should use the **JASMINE** component to solve the local atomic Hamiltonian problem at first to generate the necessary eigenvalues and eigenvectors.

**See also**:

* [Prepare Input Files](create.md) // Are you ready?
* [Standard Input Files](../ch04/input.md) // Input stuffs for CT-HYB quantum impurity solvers.
* [Standard Input Files](../ch04/choose.md) // Input stuffs for atomic eigenvalue problem solver.

!!! tip

    We have developed a graphic user interface for the iQIST software package, namely ZenGui, which can be used to generate the essential input file (i.e., *solver.ctqmc.in*) for CT-HYB quantum impurity solver. See [ZenGui](https://github.com/huangli712/ZenGui) for more details.
---

### Let's Go.

Fourth, execute the CT-HYB quantum impurity solver directly or via some external scripts/programs.

**See also**:

* [Execute The Codes](execute.md) // MPI vs. OpenMP, paralleled or sequential.
* [Monitor The Codes](monitor.md) // What's the status of the code?

---

### Post-Processing

Finally, when the calculations are finished, users can use the tools contained in the *iqist/src/tools* directory to post-process the output data, such as the imaginary-time Green's function ``G(\tau)``, Matsubara self-energy function ``\Sigma(i\omega_n)``, and other physical observables.

**See also**:

* [Auxiliary Tools](../ch06/index.md) // Full descriptions about the auxiliary tools.

Good luck to you.

!!! tip

    A full-fledged analytic continuation toolkit, namely [ACFlow](https://github.com/huangli712/ACFlow), has been developed to deal with the imaginary-time or Matsubara Green's functions to extract the corresponding spectral functions.
