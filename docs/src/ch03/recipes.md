## iQIST recipes

Using the iQIST software package is quite easy. Next we will show you the standard workflow for using the iQIST software package.

### Choose suitable component

At first, there are several CT-HYB quantum impurity solvers in the package. Their features and efficiency are somewhat different. Thus, it is the user's responsibility to choose suitable CT-HYB components to deal with the impurity problem at hand. 

**See also**:

* [Components](../ch01/components.md) // What does the iQIST software package include?
* [Features](../ch01/feature.md) // Is the required feature supported by the selected component?
* [How to choose suitable quantum impurity solvers?](../ch04/choose.md) // A guideline.

### Design the programs and scripts

Second, the iQIST software package is in essence a computational engine, so users have to write scripts or programs to execute the selected CT-HYB impurity solver directly or to call it using the application programming interface. For example, if the users want to conduct CT-HYB/DMFT calculations, in principle they must implement the DMFT self-consistent equation by themselves.

There is a bonus. When the users want to study the Hubbard model on Bethe/cubic lattice using the single-site dynamical mean-field theory, or solve the Anderson impurity models in one-shot mode, it is possible to execute the quantum impurity solvers directly without any additional scripts or programs.

**See also**:

* [Application programming interfaces](../ch08/README.md) // How to use iQIST via external Python/Fortran programs.

### Prepare the input files

Third, an important task is to prepare proper input data for the selected CT-HYB impurity solver. The optional inputs for the CT-HYB impurity solver are the hybridization function [``\Delta(i\omega_n)``], impurity energy level (``E_{\alpha\beta}``), interaction parameters (``U``, ``J``, ``\lambda``, and ``\mu``), etc. If users do not provide them to the impurity solver, it will use the default settings automatically. Specifically, if the Coulomb interaction matrix is general or the spin-orbital coupling is considered, users should use the **JASMINE** component to solve the local atomic Hamiltonian problem at first to generate the necessary eigenvalues and eigenvectors. 

**See also**:

* [Prepare input files](create.md) // Are you ready?
* [Standard input files](../ch04/input.md) // Input stuffs for CT-HYB/HF-QMC impurity solvers.
* [Standard input files](../ch06/input.md) // Input stuffs for atomic eigenvalue problem solver.

### Let's go.

Fourth, execute the CT-HYB impurity solver directly or via some external scripts/programs.

**See also**:

* [Execute the codes](execute.md) // MPI vs. OpenMP, paralleled or sequential.
* [Monitor the codes](monitor.md) // What's the status of the code?

### Post-processing

Finally, when the calculations are finished, users can use the tools contained in the **HIBISCUS** component to post-process the output data, such as the imaginary-time Green's function ``G(\tau)``, Matsubara self-energy function ``\Sigma(i\omega_n)``, and other physical observables.

**See also**:

* [Auxiliary tools](../ch07/README.md) // Full descriptions about the auxiliary toolbox.

Good luck to you.