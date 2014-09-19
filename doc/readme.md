

Interacting Quantum Impurity Solver Toolkit
===========================================


README
======



Introduction
------------

The Interacting Quantum Impurity Solver Toolkit (dubbed iQIST) is an open
source software package aiming to provide a full, reliable, flexible, and
powerful tool chain for various quantum impurity models. It contains a few
continuous-time quantum Monte Carlo impurity solvers (hybridization
expansion version), a Hirsch-Fye quantum Monte Carlo impurity solver, and
numerous prep-processed and post-processed tools. The iQIST is an all-in-one
package. With it you can solve quantum impurity models and anaylyze the
calculated results easily and efficiently.

Policy
------

The iQIST software package is released under the General Public Licence
 3.0 (GPL) or later version.

Versions
--------

The current stable release of iQIST is 0.2.x.

Features
--------

The iQIST is a powerful software package. It consists of many components 
(We just call the executable program as component in iQIST). The main
components of iQIST is the continuous-time quantum Monte Carlo impurity
solvers. So far these impurity solvers support the following features:

* Density-density interaction
* General interaction
* SOC interaction
* Hubbard model and Hubbard-Holstein model
* Frequency-dependent interaction
* Orthogonal polynomial representation
* Kernel polynomial representation
* Improved estimator for self-energy
* Single-particle Green’s function G(\tau) 
* Single-particle Green’s function G(i\omega_n) 
* Two-particle correlation function \chi(\omega, \omega', \nu)
* Local irreducible vertex function \Gamma(\omega, \omega', \nu) 
* Self-energy function \Sigma(i\omega_n)
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

The readers who are interested in it should write a letter to the authors
to request an electronic copy of the newest version of iQIST, or they can
download it directly from the public code repository. Please check the 
following website regularly:
    http://bitbucket.org

Prerequisite
------------

In order to compile and install iQIST correctly, you should ensure the
following softwares are correctly installed and configured in your OS.

* Intel Fortran compiler
* MPICH2 or OpenMPI
* BLAS
* LAPACK
* Python 2.X
* scipy, numpy, and f2py

Installation
------------

The downloaded iQIST software package is likely a compressed file with zip
or tar.gz suffix. The users should uncompress it at first. And then go to
the iqist/src/build directory, edit the make.sys file to configure the
compiling environment. Once the compiling environment is configured,
please run the make command in the top-level directory of iQIST. After a
few minutes (depending on the performance of compiling platform), the
iQIST is ready for you. Note that all of the executable programs will be
copied into the iqist/bin directory automatically. Please add this
directory into the system environment variable PATH.

Usage
-----

Documents
---------

Please see ~/doc/guide. We are sorry for it is only a Chinese guide and a
bit outdated. The English version of it will be released in the future. So
please be patient.

Support
-------

We are sorry. We DO NOT provide any technical support now. If you meet
some problems when you are using iQIST. You can write a letter to us. But
we can not guarantee we will reply you.

Developers
----------

Li Huang
Department of Physics, University of Fribourg, 1700 Fribourg, Switzerland

Yi-lin Wang
Beijing National Laboratory for Condensed Matter Physics and
Institute of Physics, Chinese Academy of Sciences, Beijing 100190, China

Contributors
------------

Zi Yang Meng
Beijing National Laboratory for Condensed Matter Physics and
Institute of Physics, Chinese Academy of Sciences, Beijing 100190, China
Department of Physics, University of Toronto, Toronto, Ontario M5S 1A7, Canada

Liang Du
Department of Physics, The University of Texas at Austin, Austin, Texas 78712, USA



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
