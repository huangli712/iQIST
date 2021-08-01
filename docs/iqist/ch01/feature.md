## Features

The main features of the core components (i.e., quantum impurity solvers) of the iQIST software package are as follows:

* **Model**
    * Density-density interaction
    * General interaction (Slater or Kanamori scheme)[^01]
    * Spin-orbital coupling and crystal field splitting[^02]
    * Hubbard-Holstein model[^03]
    * Frequency-dependent Coulomb interaction[^04]

---

* **Measurement tricks**
    * Orthogonal polynomial representation (Legendre and Chebyshev polynomials)[^05]
    * Kernel polynomial representation[^06]
    * Improved estimator for self-energy function[^07]

---

* **Observables**
    * Single-particle Green's function in imaginary time space
    * Single-particle Green's function in matsubara frequency space
    * Two-particle correlation function in matsubara frequency space (*experimental*)[^08]
    * Local irreducible vertex function in matsubara frequency space (*experimental*)[^09]
    * Pair susceptibility in matsubara frequency space (*experimental*)[^10]
    * Self-energy function in matsubara frequency space
    * Histogram of perturbation expansion order
    * Kurtosis and skewness of perturbation expansion order
    * Kinetic and potential energies
    * Orbital occupation numbers[^11]
    * Double occupation numbers[^12]
    * Magnetic moment
    * Atomic state probability
    * Spin-spin correlation function in imaginary time space[^13]
    * Spin-spin correlation function in matsubara frequency space[^13]
    * Orbital-orbital correlation function in imaginary time space[^14]
    * Orbital-orbital correlation function in matsubara frequency space[^14]
    * Fidelity susceptibility[^15]
    * Kinetic energy fluctuation $$\langle k^2\rangle - \langle k\rangle^2 - \langle k\rangle$$[^16]

---

* **Fast algorithms**
    * Segment algorithm for density-density interaction[^17]
    * Divide-and-conquer algorithm[^18]
    * Sparse matrix multiplication[^19]
    * Good quantum numbers ($$N, S_z, J_z$$, PS)[^20]
    * Lazy trace evaluation[^21]
    * Dynamical truncation approximation[^22]
    * Newton-Leja polynomial interpolation algorithm (*experimental*)[^23]

---

* **Parallelism**
    * MPI
    * OpenMP[^24]

---

* **API**
    * Python binding
    * Input file generator by Python
    * Fortran binding

---

* **Preprocessing**
    * Atomic eigenvalue problem solver[^25]

---

* **Postprocessing**
    * Maximum entropy method[^26]
    * Stochastic analytical continuation[^27]
    * Kramers-Kronig transformation[^28]
    * Pade approximation[^29]
    * Polynomial fitting for self-energy function[^30]
    * Many tools and scripts, etc.

---

> **Additional Limitations:**

> 
[^01] Only for the **BEGONIA**, **LAVENDER**, **CAMELLIA**, **PANSY**, and **MANJUSHAKA**.

> 
[^02] Only for the **BEGONIA**, **LAVENDER**, **PANSY**, and **MANJUSHAKA**.

> 
[^03] Only for the **NARCISSUS**.

> 
[^04] Only for the **NARCISSUS**.

> 
[^05] Only for the **GARDENIA**, **NARCISSUS**, **LAVENDER**, **CAMELLIA**, and **MANJUSHAKA**.

> 
[^06] Only for the **GARDENIA**, **NARCISSUS**, **LAVENDER**, **CAMELLIA**, and **MANJUSHAKA**.

> 
[^07] Only for the **GARDENIA** and **NARCISSUS**.

> 
[^08] Only for the **GARDENIA**, **NARCISSUS**, **LAVENDER**, **CAMELLIA**, and **MANJUSHAKA**.

> 
[^09] Only for the **GARDENIA**, **NARCISSUS**, **LAVENDER**, **CAMELLIA**, and **MANJUSHAKA**.

> 
[^10] Only for the **GARDENIA**, **NARCISSUS**, **LAVENDER**, **CAMELLIA**, and **MANJUSHAKA**.

> 
[^11] Only for the **AZALEA**, **GARDENIA**, and **NARCISSUS**.

> 
[^12] Only for the **AZALEA**, **GARDENIA**, and **NARCISSUS**.

> 
[^13] Only for the **GARDENIA** and **NARCISSUS**.

> 
[^14] Only for the **GARDENIA** and **NARCISSUS**.

> 
[^15] Only for the **GARDENIA**, **NARCISSUS**, **LAVENDER**, **CAMELLIA**, and **MANJUSHAKA**.

> 
[^16] Only for the **GARDENIA**, **NARCISSUS**, **LAVENDER**, **CAMELLIA**, and **MANJUSHAKA**.

> 
[^17] Only for the **AZALEA**, **GARDENIA**, and **NARCISSUS**.

> 
[^18] Only for the **BEGONIA**, **LAVENDER**, **PANSY**, and **MANJUSHAKA**.

> 
[^19] Only for the **BEGONIA**, **LAVENDER**, and **CAMELLIA**.

> 
[^20] Only for the **PANSY** and **MANJUSHAKA**.

> 
[^21] Only for the **MANJUSHAKA**.

> 
[^22] Only for the **MANJUSHAKA**.

> 
[^23] Only for the **CAMELLIA**.

> 
[^24] Only for the measurement of two-particle quantities.

> 
[^25] Only for the **JASMINE**.

> 
[^26] Only for the **HIBISCUS**.

> 
[^27] Only for the **HIBISCUS**.

> 
[^28] Only for the **HIBISCUS**.

> 
[^29] Only for the **HIBISCUS**.

> 
[^30] Only for the **HIBISCUS**.