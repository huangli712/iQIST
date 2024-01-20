# solver.hyb.dat

**Introduction**

This file is used to store the Matsubara hybridization function ``\Delta(i\omega_n)``. It will be output by the quantum impurity solvers when they are **shut down**.

**Format**

The *solver.hyb.dat* file contains *norbs* block. Each block is appended by two blank lines. The format of each block is as follows:

---

*column 1*: orbital index ``i``, integer

*column 2*: Matsubara frequency point, ``\omega_n``, double precision

*column 3*: Matsubara hybridization function, ``\Re \Delta(i\omega_n)``, double precision

*column 4*: Matsubara hybridization function, ``\Im \Delta(i\omega_n)``, double precision

*column 5*: error bar, ``\Re [\delta \Delta(i\omega_n)]``, double precision

*column 6*: error bar, ``\Im [\delta \Delta(i\omega_n)]``, double precision

---

!!! note

    In the *solver.hyb.dat* file, we adopt the following orbital sequence:
    $$1\uparrow$$, $$2\uparrow$$, $$3\uparrow$$, $$\cdots$$, $$1\downarrow$$, $$2\downarrow$$, $$3\downarrow$$, $$\cdots$$
    In other words, the spin up part is always before the spin down part.

**Code**

The corresponding Fortran code block for the writing of *solver.hyb.dat* file is as follows:

```fortran
! open data file: solver.hyb.dat
open(mytmp, file='solver.hyb.dat', form='formatted', status='unknown')

! write it
do i=1,norbs
    do j=1,mfreq
        write(mytmp,'(i6,5f16.8)') i, rmesh(j), hybf(j,i,i), czero
    enddo ! over j={1,mfreq} loop
    write(mytmp,*) ! write empty lines
    write(mytmp,*)
enddo ! over i={1,norbs} loop

! close data file
close(mytmp)
```

!!! note

    The columns for the error bar are always zero in this file.

In the iQIST software package, we provide a Python module to read the output files of quantum impurity solvers. You can use it to read the *solver.hyb.dat* file. See [src/tools/u_reader.py] for more details.
