# iQIST (Interacting Quantum Impurity Solver Toolkit)

### v0.5.5 @ 797ace8 // Feb 05, 2015

* Rename check.py to d\_check.py, fix the comment.
* Rename clean.py to d\_clean.py, fix the comment.
* Reanme cmp.py to d\_cmp.py, fix the comment.
* Rename sar.sh to d\_sar.sh, fix the comment.
* Reanme trailing.sh to d\_trailing.sh, fix the comment.
* Add and implement u\_writer.py.
* Fix some deadly bugs in u\_reader.py
* Fix the comments in ctqmc\_main.f90, include solver.anydos.in

### v0.5.4 @ 3ac57c2 // Feb 04, 2015

* Add icu = 3 in the jasmine code. See comments in the atomic\_control.f90 file.
* Add ROADMAP.md file.

### v0.5.3 @ ab08c42 // Feb 02, 2015

* Remove the chinese directory in doc.

### v0.5.2 @ 67b5561 // Jan 31, 2015

* Update the input files in tutor directory.
* Add t961, t962, t963, t964 directories, and tutor.md file in the tutor directory.
* Rename F2PY flag to MPY flag to avoid potential conflict.
* Correct some typos in working/advanced/begonia and working/advanced/lavender directories.
* Fix bugs/comments in ctqmc\_api.f90, hfqmc\_api.f90, and atomic\_api.f90.
* Now the Python API and Fortran API work.
* Update the manual.

### v0.5.1 @ 222982c // Jan 16, 2015

* Add solver.umat.in support for the azalea, gardenia, and narcissus codes.
* Reconstruct the ctqmc\_make\_uumat() for the narcissus code.
* Add ctqmc\_make\_shift() in ctqmc\_util.f90 for the narcissus code.
* Modify the iqist/src/ctqmc/api/ctqmc\_api.f90, add set\_uumat().
* Modify all of the ctqmc\_main.f90 files, implement the cat\_set\_uumat() subroutine.
* Improve the comments.

### v0.5.0 @ e182bbf // Jan 16, 2015

* Add cmp.py file in the iqist/src/tools/hibiscus/script directory.
* Fix small bug in the jasmine code.
* Add solver.umat.in output in the jasmine code.

### v0.4.9 @ 7a31b3c // Jan 09, 2015

* Add u\_ready.py file in the iqist/src/tools/hibiscus/script directory.

### v0.4.8 @ 66a3e72 // Jan 08, 2015

* The benchmark files in the iqist/working directory are ready.
