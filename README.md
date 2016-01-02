# iQIST (Interacting Quantum Impurity Solver Toolkit)

The iQIST software package includes several quantum impurity solvers which implement the hybridization expansion version continuous-time quantum Monte Carlo algorithm and Hirsch-Fye quantum Monte Carlo algorithm, and corresponding preprocessing and postprocessing tools.

## WARNING

The iQIST is still in heavy development. The codes are extremely unstable. Some features are still experimental. Everything could be changed in the future release. We can not guarantee that it is bug free. So be careful when you are using it and verify your data again and again before you submit your calculated results to any peer-reviewed journal.

Sometimes the latest commit will not be compiled correctly. So, please download the released version of iQIST which has an unique version tag.

## Version

v0.6.6 @ 2015.01.06T (devel)

## License

GNU General Public License Version 3

## Features

* Model
    * Density-density interaction
    * General interaction (Slater or Kanamori scheme)
    * SOC interaction and crystal field splitting
    * Hubbard-Holstein model
    * Frequency-dependent interaction

* Measurement tricks
    * Orthogonal polynomial representation (Legendre and Chebyshev polynomials)
    * Kernel polynomial representation
    * Improved estimator for self-energy function

* Observables
    * Single-particle Green's function in imaginary time space
    * Single-particle Green's function in matsubara frequency space
    * Two-particle correlation function in matsubara frequency space
    * Local irreducible vertex function in matsubara frequency space
    * Pair susceptibility in matsubara frequency space
    * Self-energy function in matsubara frequency space
    * Histogram of perturbation expansion order
    * Kurtosis and skewness of perturbation expansion order
    * Kinetic and potential energies
    * Orbital occupation numbers
    * Double occupation numbers
    * Magnetic moment
    * Atomic state probability
    * Spin-spin correlation function in imaginary time space
    * Orbital-orbital correlation function in imaginary time space
    * Fidelity susceptibility
    * kinetic energy fluctuation <k^2> - <k>^2 - <k>

* Fast algorithms
    * Segment algorithm for density-density interaction
    * Divide-and-conquer algorithm
    * Sparse matrix multiplication
    * Good quantum numbers (N, Sz, Jz, PS)
    * Lazy trace evaluation
    * Dynamical truncation approximation

* Parallelism
    * MPI
    * OpenMP (for the measurement of two-particle quantities)

* API
    * Python binding
    * Input file generator by Python
    * Fortran binding

* Preprocessing
    * Atomic eigenvalue problem solver

* Postprocessing
    * Maximum entropy method
    * Stochastic analytical continuation
    * Kramers-Kronig transformation
    * Pade approximation
    * Polynomial fitting for self-energy function
    * Many tools and scripts, etc.

## Installation

* Full Installation
```sh
$ cd iqist/build
$ editor make.sys
$ make all
$ ./x_setup.sh
```

* Partial Installation
```sh
$ cd iqist/build
$ editor make.sys
$ make base
$ make capi
$ make component (component could be azalea, gardenia, narcissus, etc.)
$ ./x_setup.sh
```

* Build Fortran Library
```sh
$ cd iqist/build
$ editor make.sys
$ make base
$ make capi
$ make component-lib (component could be azalea, gardenia, narcissus, etc.)
$ ./x_setup.sh
```

* Build Python Module
```sh
$ cd iqist/build
$ editor make.sys
$ make base
$ make capi
$ make component-pylib (component could be azalea, gardenia, narcissus, etc.)
$ ./x_setup.sh
```

Enjoy it!

If you want to know more about the compiling system implemented in the iQIST, please read the manual carefully.

## Documentation

see iQIST/doc/manual/ug.pdf (We are sorry. Currently this manual is far away from completeness, so we remove it temporally from the release).

## Development

The iQIST software package is developed and maintained by the iQIST Developer Team.

Find a bug? Want to contribute? Want new features? Great! Please contact with us as soon as possible.

## Reference

If you are using iQIST to do some studies and would like to publish your great works, it would be really appreciated if you can cite the following paper:

```sh
iQIST: An open source continuous-time quantum Monte Carlo impurity solver toolkit
Li Huang, Yilin Wang, Zi Yang Meng, Liang Du, Philipp Werner and Xi Dai
Computer Physics Communications 195, 140 (2015) or arXiv:1409.7573 (2014)
```

## Contact

```sh
Li Huang
Institute of Materials, China Academy of Engineering Physics, Sichuan, PRC
email: lihuang.dmft at gmail.com
```

or

```sh
Yilin Wang
Institute of Physics, Chinese Academy of Sciences, Beijing, PRC
email: qhwyl2006 at 126.com
```
