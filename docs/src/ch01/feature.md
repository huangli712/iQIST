# Features

The main features of the core components (i.e., quantum impurity solvers) of the iQIST software package are summarized as follows:

---

* **Model**
    1. Density-density interaction
    2. General interaction (Slater or Kanamori scheme)
    3. Spin-orbital coupling and crystal field splitting
    4. Frequency-dependent Coulomb interaction

---

* **Measurement Tricks**
    1. Orthogonal polynomial representation (Legendre polynomials)
    2. Intermediate representation
    3. Improved estimator for self-energy function
    4. Improved estimator for two-particle green's function

---

* **Observables**
    1. Single-particle Green's function in imaginary time space
    2. Single-particle Green's function in matsubara frequency space
    3. Two-particle correlation function in matsubara frequency space
    4. Self-energy function in matsubara frequency space
    5. Histogram of perturbation expansion order
    6. Kurtosis and skewness of perturbation expansion order
    7. Kinetic and potential energies
    8. Orbital occupation numbers
    9. Double occupation numbers
    10. Magnetic moment
    11. Atomic state probability
    12. Spin-spin correlation function in imaginary time space
    13. Spin-spin correlation function in matsubara frequency space
    14. Orbital-orbital correlation function in imaginary time space
    15. Orbital-orbital correlation function in matsubara frequency space
    16. Fidelity susceptibility
    17. Kinetic energy fluctuation ``\langle k^2\rangle - \langle k\rangle^2 - \langle k\rangle``

---

* **Fast Algorithms**
    1. Segment algorithm for density-density interaction
    2. Divide-and-conquer algorithm
    3. Sparse matrix multiplication
    4. Good quantum numbers (``N, S_z, J_z``, PS, automatic partition)
    5. Lazy trace evaluation
    6. Dynamical truncation approximation

---

* **Parallelism**
    1. MPI
    2. OpenMP

---

* **Preprocessing**
    1. Atomic eigenvalue problem solver
