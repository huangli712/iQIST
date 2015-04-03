# iQIST (Interacting Quantum Impurity Solver Toolkit)

The iQIST software package includes several quantum impurity solvers which implement the hybridization expansion version continuous-time quantum Monte Carlo algorithm and Hirsch-Fye quantum Monte Carlo algorithm, and corresponding preprocessing and postprocessing tools.

### WARNING

The iQIST is still in heavy development. The codes are extremely unstable. Some features are still experimental. Everything could be changed in the future release. We can not guarantee that it is bug free. So be careful when you are using it and verify your data again and again before you submit your calculated results to any peer-reviewed journal.

### Version

v0.6.0 @ 2015.01.06T (beta)

### License

GNU General Public License Version 3

### Features

* Model
    * Density-density interaction
    * General interaction (Slater or Kanamori scheme)
    * SOC interaction and crystal field splitting
    * Hubbard-Holstein model
    * Frequency-dependent interaction

* Measurement tricks
    * Orthogonal polynomial representation (Legendre and Chebyshev polynomials)
    * Kernel polynomial representation
    * Improved estimator for self-energy

* Observables
    * Single-particle Green's function in imaginary time space
    * Single-particle Green's function in matsubara frequency space
    * Two-particle correlation function in matsubara frequency space
    * Local irreducible vertex function in matsubara frequency space
    * Pair susceptibility in matsubara frequency space
    * Self-energy function in matsubara frequency space
    * Histogram of perturbation expansion order
    * Kinetic and potential energies
    * (Double) occupation numbers, magnetic moment
    * Atomic state probability
    * Spin-spin correlation function
    * Orbital-orbital correlation function
    * Autocorrelation function and autocorrelation time

* Fast algorithms
    * Segment algorithm for density-density interaction
    * Divide-and-conquer algorithm
    * Sparse matrix multiplication
    * Good quantum numbers (N, Sz, Jz, PS)
    * Skip listing trick
    * Lazy trace evaluation
    * Dynamical truncation approximation

* Parallelism
    * MPI

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

### Installation
* Full Installation
```sh
$ cd iqist/src/build
$ editor make.sys
$ make all
$ cd ../../bin
$ ./setup.sh
```

* Partial Installation
```sh
$ cd iqist/src/build
$ editor make.sys
$ make common
$ make api
$ make component (component could be azalea, gardenia, narcissus, etc.)
$ cd ../../bin
$ ./setup.sh
```

Enjoy it!

If you want to know more about the compiling system implemented in the iQIST, please read the manual carefully.

### Documentation

see iQIST/doc/manual/ug.pdf (We are sorry. Currently this manual is far away from completeness).

### Development

The iQIST software package is developed and maintained by the iQIST Developer Team.

Find a bug? Want to contribute? Want new features? Great! Please contact with us as soon as possible.

### Reference

If you are using iQIST to do some studies and would like to publish your great works, it would be really appreciated if you can cite the following paper:

```sh
iQIST: An open source continuous-time quantum Monte Carlo impurity solver toolkit
Li Huang, Yilin Wang, Zi Yang Meng, Liang Du, Philipp Werner and Xi Dai
arXiv:1409.7573 (2014)
```

### Contact

```sh
Li Huang
Department of Physics, Fribourg University, Switzerland
email: huangli712 at gmail.com
```

or

```sh
Yilin Wang
Institute of Physics, Chinese Academy of Sciences, Beijing, PRC
email: qhwyl2006 at 126.com
```
