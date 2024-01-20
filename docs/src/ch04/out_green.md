# solver.green.dat

**Introduction**

This file is used to store the imaginary-time Green's function ``G(\tau)``. It will be output by the quantum impurity solvers when they are **shut down**.

**Format**

The *solver.green.dat* file contains *norbs* block. Each block is appended by two blank lines. The format of each block is as follows:

---

*column 1*: orbital index ``i``, integer

*column 2*: imaginary-time index ``j``, integer

*column 3*: imaginary-time point, ``\tau``, double precision

*column 4*: imaginary-time Green's function, ``G(\tau)``, double precision

*column 5*: error bar, ``\delta G(\tau)``, double precision

---

!!! note

    In the *solver.green.dat* file, we adopt the following orbital sequence:
    $$1\uparrow$$, $$2\uparrow$$, $$3\uparrow$$, $$\cdots$$, $$1\downarrow$$, $$2\downarrow$$, $$3\downarrow$$, $$\cdots$$
    In other words, the spin up part is always before the spin down part.

**Code**

The corresponding Fortran code block for the writing of *solver.green.dat* file is as follows:

```fortran
! open data file: solver.green.dat
open(mytmp, file='solver.green.dat', form='formatted', status='unknown')

! write it
do i=1,norbs
    do j=1,ntime
        write(mytmp,'(2i6,3f12.6)') i, j, tmesh(j), gtau(j,i,i), gerr(j,i,i)
    enddo ! over j={1,ntime} loop
    write(mytmp,*) ! write empty lines
    write(mytmp,*)
enddo ! over i={1,norbs} loop

! close data file
close(mytmp)
```

In the iQIST software package, we provide a Python module to read the output files of quantum impurity solvers. You can use it to read the *solver.green.dat* file. See [src/tools/u_reader.py] for more details.
