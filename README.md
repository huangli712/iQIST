# iQIST (Interacting Quantum Impurity Solver Toolkit)

The iQIST software package contains several state-of-the-art quantum impurity solvers, which implement the hybridization expansion version continuous-time quantum Monte Carlo algorithm (CT-HYB), and corresponding preprocessing and postprocessing tools.

## Version

v0.6.8 @ 2017.01.31D (devel)

## License

GNU General Public License Version 3

## Features

* The quantum impurity models could have the following terms
    * Density-density interaction
    * General interaction (Slater or Kanamori scheme)
    * Spin-orbital coupling
    * Crystal field splitting
    * Frequency-dependent interaction

* The following physical observables could be measured
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
    * Binder cumulant
    * Atomic state probability
    * Spin-spin correlation function in imaginary time space
    * Spin-spin correlation function in matsubara frequency space
    * Orbital-orbital correlation function in imaginary time space
    * Orbital-orbital correlation function in matsubara frequency space
    * Fidelity susceptibility
    * Kinetic energy fluctuation

* The following measurement tricks are supported
    * Orthogonal polynomial representation (Legendre and Chebyshev polynomials)
    * Kernel polynomial representation
    * Improved estimator for self-energy function

* The following optimized algorithms are adopted
    * Segment algorithm for density-density interaction
    * Divide-and-conquer algorithm
    * Good quantum numbers (N, Sz, Jz, PS)
    * Lazy trace evaluation
    * Dynamical truncation approximation

* The quantum impurity solvers are all parallelized
    * MPI
    * OpenMP (for the measurement of two-particle quantities)

* The application programming interfaces are provided
    * Fortran binding
    * Python binding
    * Input file generator by Python

* The preprocessing tools are provided
    * Atomic eigenvalue problem solver

* The postprocessing tools are provided
    * Many tools and scripts, etc.

> NOTE:
>
> The iQIST software package is still in heavy development. The codes are extremely unstable. Some features are still experimental. Everything could be changed or removed in the future release. We can not guarantee that it is bug free. So be careful when you are using it and verify your data again and again before you submit your calculated results to any peer-reviewed journal.

## Installation

* Full installation

```sh
$ cd iqist/build
$ editor make.sys
$ make all
```

* Partial installation

```sh
$ cd iqist/build
$ editor make.sys
$ make base
$ make capi
$ make component
```

* Build Fortran library

```sh
$ cd iqist/build
$ editor make.sys
$ make base
$ make capi
$ make component-lib
```

* Build Python module

```sh
$ cd iqist/build
$ editor make.sys
$ make base
$ make capi
$ make component-pylib
```

> NOTE:
>
> 1. 'iqist' is the directory where the iQIST software package is uncompressed.
>
> 2. 'editor' could be any ascii text editor which you prefer.
>
> 3. 'component' could be narcissus, manjushaka, etc.
>
> 4. Type 'make help-more' in the terminal for more details.
>
> 5. Sometimes the latest commit will not be compiled correctly. So, please download the released version of the iQIST software package which should have an unique version tag.

Enjoy it!

If you want to know more about the compiling system implemented in the iQIST software package, please read the reference manual carefully.

## Documentation

We provide a comprehensive reference manual for the iQIST software package via the Gitbook. Please go to:

```sh
https://www.gitbook.com/book/huangli712/iqist/
```

for more details.

## Development

The iQIST software package is developed and maintained by the iQIST Developer Team.

Find a bug? Want to contribute? Want new features? Great! Please contact with us as soon as possible.

## Reference

If you are using the iQIST software package to do some studies and would like to publish your great works, it would be really appreciated if you can cite the following paper:

```sh
iQIST: An open source continuous-time quantum Monte Carlo impurity solver toolkit
Li Huang, Yilin Wang, Zi Yang Meng, Liang Du, Philipp Werner and Xi Dai
Computer Physics Communications 195, 140 (2015) or arXiv:1409.7573 (2014)
```

## Contact

```sh
Li Huang
Institute of Materials, China Academy of Engineering Physics, Sichuan Jiangyou, PRC
email: lihuang.dmft at gmail.com
```

or

```sh
Yilin Wang
Institute of Physics, Chinese Academy of Sciences, Beijing, PRC
email: qhwyl2006 at 126.com
```
