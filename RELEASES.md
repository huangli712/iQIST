# iQIST (Interacting Quantum Impurity Solver Toolkit)

## v0.6.8 (in progress)

* Remove azalea component.
* Remove begonia component.
* Remove camellia component.
* Remove pansy component.
* Remove lavender component.
* Update the compiling system.

## v0.6.7 @ c14fdc4 // Jan 13, 2017

* Remove ROADMAP.md.
* Remove iqist/doc directory.
* The building system is updated, remove support for doc target.
* The official reference manual is released at Gitbook.
* Fix the file format of solver.hist.dat.
* Add script/u_animator.py which can produce movie using the solver.diag.dat file.
* Now the CT-QMC impurity solvers can output the solver.diag.dat file.
* Now the HF-QMC impurity solver can output the solver.diag.dat file.
* Optimize the algorithm for the calculation of spin-spin correlation function (time space)
* Optimize the algorithm for the calculation of orbital-orbital correlation function (time space)
* Implement the measurement of spin-spin correlation function (matsubara space)
* Implement the measurement of orbital-orbital correlation function (matsubara space)
* Improve the code styles, fix many typos in comments.
* Add camellia component. **BUT IT DOES NOT WORK PROPERLY. DON'T USE IT!**
* Add kurtosis and skewness features in azalea.
* Add kurtosis and skewness features in begonia.
* Add kurtosis and skewness features in camellia.
* Add kurtosis and skewness features in gardenia.
* Add kurtosis and skewness features in lavender.
* Add kurtosis and skewness features in manjushaka.
* Add kurtosis and skewness features in narcissus.
* Add kurtosis and skewness features in pansy.
* Remove doc/manual. In the future, the documents will be published in Gitbook.
* The building system is modified. And the build directory is moved to the iqist directory.
* Remove iqist/bin directory, and the shell scripts are copied into iqist/build and renamed.

## v0.6.5 @ 90793f1 // Sep 06, 2015

* Add error bars support in azalea.
* Add error bars support in begonia.
* Add error bars support in gardenia.
* Add error bars support in lavender.
* Add error bars support in manjushaka.
* Add error bars support in narcissus.
* Add error bars support in pansy.
* Add error bars support in daisy.
* Re-design the Python and Fortran APIs.
* Re-design the output file format.
* Improve hibiscus to be compatible with the new file formats.
* Fix a serious bug in hibiscus/toolbox/maketau.f90
* Refine the examples/tutorials.
* Improve the code style
* Add <k^2> - <k>^2 support for gardenia.
* Add <k^2> - <k>^2 support for narcissus.
* Add <k^2> - <k>^2 support for lavender.
* Add <k^2> - <k>^2 support for manjushaka.
* Add <k^2> - <k>^2 support for script/u\_reader.py
* Add fidelity susceptibility support for gardenia.
* Add fidelity susceptibility support for narcissus.
* Add fidelity susceptibility support for lavender.
* Add fidelity susceptibility support for manjushaka.
* Add fidelity susceptibility support for script/u\_reader.py
* Now the iQIST is compatible with the GNU gfortran compiler.
* Add several templates for make.sys.

> WARNING:
>
> The file formats, Python and Fortran APIS are **NOT** compatible with the previous. Especially, the solver.hyb.in files generated using the previous iQIST are not valid for the current version.

## v0.6.0 @ aa24237 // Apr 03, 2015

* Update working/ctqmc/validation directory, add more examples.
* Fix the sign problems in the measurement.
* Redefine the behavior of control variable: isvrt.
* Introduce new control variable: issus.
* Fix a bug in the calculation of self-energy function for narcissus code.
* Add working/ctqmc/validation directory.
* Improve build/make.sys, add 'fPIC' option.
* Add OpenMP support for the measurements of two-particle quantities.
* Adjust the file format for solver.kernel.dat.
* Adjust the output for hybridization function.
* Fix the u\_reader.py to support new file format of solver.kernel.dat.
* Fix typo in makescr.f90 and comment in make.sys.
* Modify jasmine/atomic\_util.f90, remove unused variables.
* Add script/d\_archive.sh shell script.
* Fix a format bug in script/u\_writer.py.
* Modify ctqmc\_record.f90 in gardenia, narcissus, lavender, and manjushaka.
* Improve the computational efficiency for two-particle quantities.
* Improve the compiling system, add 'lib' target.
* Move src/ctqmc/api to src/api.
* Move jasmine/atomic\_api.f90 to src/api.
* Move daisy/hfqmc\_api.f90 to src/api.
* Update the compiling system to support new api directory.
* Add src/app directory.
* Update the tutor/t961-t964, fix path bug.

## v0.5.5 @ 797ace8 // Feb 05, 2015

* Rename check.py to d\_check.py, fix the comment.
* Rename clean.py to d\_clean.py, fix the comment.
* Reanme cmp.py to d\_cmp.py, fix the comment.
* Rename sar.sh to d\_sar.sh, fix the comment.
* Reanme trailing.sh to d\_trailing.sh, fix the comment.
* Add and implement u\_writer.py.
* Fix some deadly bugs in u\_reader.py.
* Fix the comments in ctqmc\_main.f90, include solver.anydos.in.
* Add icu = 3 in the jasmine code. See comments in the atomic\_control.f90 file.
* Add ROADMAP.md file.
* Remove the chinese directory in doc.
* Update the input files in tutor directory.
* Add t961, t962, t963, t964 directories, and tutor.md file in the tutor directory.
* Rename F2PY flag to MPY flag to avoid potential conflict.
* Correct some typos in working/advanced/begonia and working/advanced/lavender directories.
* Fix bugs/comments in ctqmc\_api.f90, hfqmc\_api.f90, and atomic\_api.f90.
* Now the Python API and Fortran API work.
* Update the manual.
* Add solver.umat.in support for the azalea, gardenia, and narcissus codes.
* Reconstruct the ctqmc\_make\_uumat() for the narcissus code.
* Add ctqmc\_make\_shift() in ctqmc\_util.f90 for the narcissus code.
* Modify the iqist/src/ctqmc/api/ctqmc\_api.f90, add set\_uumat().
* Modify all of the ctqmc\_main.f90 files, implement the cat\_set\_uumat() subroutine.
* Improve the comments.

## v0.5.0 @ e182bbf // Jan 16, 2015

* Add cmp.py file in the iqist/src/tools/hibiscus/script directory.
* Fix small bug in the jasmine code.
* Add solver.umat.in output in the jasmine code.
* Add u\_ready.py file in the iqist/src/tools/hibiscus/script directory.
* The benchmark files in the iqist/working directory are ready.
