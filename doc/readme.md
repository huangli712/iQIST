

Interacting Quantum Impurity Solver Toolkit
===========================================


README
======



Introduction
------------

The Interacting Quantum Impurity Solver Toolkit (dubbed iQIST) is an open
source software package aiming to provide a full, reliable, flexible, and
powerful tool chain for various quantum impurity models. It contains some
continuous-time quantum Monte Carlo impurity solvers (hybridization expa-
nsion version), Hirsch-Fye quantum Monte Carlo impurity solver, and corr-
esponding prep-processed and post-processed tools.

Licence
-------

The iQIST software package is released under the General Public Licence (GPL)
version 3.

Versions
--------

The current release of iQIST is 0.2.x.

Features
--------

* Density-density interaction
* General interaction
* SOC interaction
* Hubbard model and Hubbard-Holstein model
* Frequency-dependent interaction
* Orthogonal polynomial representation
* Kernel polynomial representation
* Improved estimator for self-energy
* Single-particle Green’s function G(τ) 
* Single-particle Green’s function G(iωn) 
* Two-particle correlation function χ(ω, ω′, ν)
* Local irreducible vertex function Γ(ω, ω′, ν) 
* Self-energy function Σ(iωn)
* Histogram of perturbation expansion order 
* Kinetic and potential energies
* (Double) occupation numbers, magnetic moment 
* Atomic state probability
* Spin-spin correlation function
* Orbital-orbital correlation function 
* Autocorrelation function and autocorrelation time 
* Divide-and-conquer algorithm
* Sparse matrix multiplication
* Good quantum numbers
* Skip listing trick
* Lazy trace evaluation
* Dynamical truncation approximation

Obtain
------

The user can write a letter to one of the authors to request the newest copy
of iQIST or download it directly from the following website.

Prerequisite
------------

* Intel Fortran compiler
* MPICH2 or OpenMPI
* BLAS
* LAPACK
* Python 2.X
* scipy, numpy, and f2py


Installation
------------

The downloaded iQIST software package is likely a compressed file with zip or tar.gz suffix. The users should uncompress it at first. And then go to the iqist/src/build directory, edit the make.sys file to configure the compiling environment. The users must setup the Fortran compiler, MPI compiler, BLAS and LAPACK libraries manually. The components in iQIST can be successfully compiled using recent Intel Fortran compiler. Most of the MPI implementations, such as MPICH, MVAPICH, OpenMPI and Intel MPI
are compatible with iQIST. As for the BLAS implementation, we strongly recommend
the OpenBLAS. For the LAPACK, there is no doubt that the Intel Math Kernel Library is a good candidate. Of course, you can also use the linear algebra library provided by the operating system, for example, vecLib Framework in the Mac OS X. Some post-processed scripts contained in the HIBISCUS component are developed with Python language. In order to execute these scripts or use the Python language binding for iQIST, the users should ensure Python 2.x was installed successfully. Furthermore, the numpy, scipy, and f2py packages are also necessary. Once the compiling environment is configured, please run the make command in the top-level directory of iQIST. After a few minutes (depending on the performance of compiling platform), the iQIST is ready for you. Note that all of the executable programs will be copied into the iqist/bin directory automatically. Please add this directory into the system environment variable PATH.

Usage
-----

Documents
---------

Please see ~/doc/guide. We are sorry for it is only a Chinese guide and a
bit outdated. The English version of it will be released in the future.
So please be patient.

Support
-------

We are sorry. We DO NOT provide any technical support. If you meet some problems when
you are using iQIST. You can write a letter to us. But we can not guarantee we will
reply you.

Developers
----------

Li Huang
Yi-lin Wang

Contributors
------------

Zi Yang Meng
Liang Du


HISTORY
=======



v0.2.0 // Aug
-------------

* add bin/clean.sh
* add doc/guide support in src/build/Makefile


v0.1.9 // Aug 18, 2014
----------------------

* make new building/compiling system (src/build).
* make setup shell script (bin/).
* update CSSL code (src/common/s_vector.f90).
* refine azalea code (src/ctqmc/azalea).
* refine ctqmc api (src/ctqmc/api).


v0.1.8 // Aug 4, 2014
---------------------

* add pansy code (experimental).
* add manjushaka code (experimental).
* add jasmine code (experimental).
* implement CSSL and CSML codes (experimental).


v0.1.7 // Jul 2, 2014
---------------------

* add entropy code.
* add stochastic code.
* add swing code.
* add toolbox code.


v0.1.6 // Jul 2, 2014
---------------------

* add iris code.


v0.1.5 // Jul 2, 2014
---------------------

* add daisy code.


v0.1.4 // Jul 2, 2014
---------------------

* add lavender code.


v0.1.3 // Jul 2, 2014
---------------------

* add begonia code.


v0.1.2 // Jul 2, 2014
---------------------

* change the file mode for gardenia code.


v0.1.1 // Jul 2, 2014
---------------------

* change the file mode for azalea code.


v0.1.0 // Jul 2, 2014
---------------------

* init the whole directory structure.


v0.0.0 // Jul 2, 2014
---------------------

* init the project.
