# iQIST (Interacting Quantum Impurity Solver Toolkit)

The iQIST software package contains several state-of-the-art continuous-time quantum Monte Carlo impurity solvers (which implement the hybridization expansion algorithm), auxiliary tools, numerical libraries and a few typical applications.

## Version

v0.7.0 @ 2017.01.31D (devel)

## License

GNU General Public License Version 3

## Features

* The quantum impurity models could have the following terms
    * Dynamic or static density-density interaction
    * General interaction (Slater or Kanamori scheme)
    * Spin-orbital coupling
    * Crystal field splitting

* The following physical observables could be measured either directly or indirectly
    * One-particle Green's function in imaginary time space
    * One-particle Green's function in matsubara frequency space
    * Two-particle Green's function in matsubara frequency space
    * Self-energy function in matsubara frequency space
    * Spin-spin correlation function in imaginary time space
    * Spin-spin correlation function in matsubara frequency space
    * Orbital-orbital correlation function in imaginary time space
    * Orbital-orbital correlation function in matsubara frequency space
    * Histogram of perturbation expansion order
    * Kurtosis of perturbation expansion order
    * Skewness of perturbation expansion order
    * Kinetic energy
    * Potential energy
    * Orbital occupation numbers
    * Double occupation numbers
    * Local magnetic moment
    * Binder cumulant
    * Atomic state probability
    * Fidelity susceptibility
    * Kinetic energy fluctuation

* The following measurement tricks are supported
    * Legendre orthogonal polynomial representation
    * Singular value decomposition representation
    * Improved estimator for self-energy function
    * Improved estimator for two-particle vertex function

* The following optimized algorithms are adopted
    * Segment algorithm for density-density interaction
    * Divide-and-conquer algorithm
    * Good quantum numbers (N, Sz, Jz, PS)
    * Lazy trace evaluation
    * Dynamical truncation approximation

* The quantum impurity solvers are parallelized
    * MPI
    * OpenMP

* Many tools are provided
    * Atomic eigenvalue problem solver
    * Some users-oriented and developers-oriented scripts, etc.

The iQIST software package is still in heavy development. The codes are extremely unstable. Some features are still experimental. Everything could be changed or removed in the future release. We can not guarantee that it is bug free. So be careful when you are using it and verify your data again and again before you submit your calculated results to any peer-reviewed journal.

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
$ make component
```

Here '**iqist**' is the directory where the iQIST software package is uncompressed, '**editor**' could be any ascii text editor which you prefer, '**component**' could be **narcissus**, **manjushaka**, etc. Please type '**make help-more**' in the terminal for more details. Sometimes the latest commit will not be compiled correctly. So, please download the released version of the iQIST software package which should have an unique version tag.

## Documentation

We provide a comprehensive online [reference manual](https://www.gitbook.com/book/huangli712/iqist/) for the iQIST software package via the Gitbook.

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
