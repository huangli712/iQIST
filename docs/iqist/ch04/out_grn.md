### solver.grn.dat

**Introduction**

This file is used to store the Matsubara Green's function $$G(i\omega_n)$$. It will be output by the quantum impurity solvers when they are **shut down**.

**Format**

The *solver.grn.dat* file contains *norbs* block. Each block is appended by two blank lines. The format of each block is as follows:

---

*column 1*: orbital index $$i$$, integer

*column 2*: Matsubara frequency point, $$\omega_n$$, double precision

*column 3*: Matsubara Green's function, $$\Re G(i\omega_n)$$, double precision

*column 4*: Matsubara Green's function, $$\Im G(i\omega_n)$$, double precision

*column 5*: error bar, $$\Re [\delta G(i\omega_n)]$$, double precision

*column 6*: error bar, $$\Im [\delta G(i\omega_n)]$$, double precision

---

> NOTE:

> In the *solver.grn.dat* file, we adopt the following orbital sequence:

> $$1\uparrow$$, $$2\uparrow$$, $$3\uparrow$$, $$\cdots$$, $$1\downarrow$$, $$2\downarrow$$, $$3\downarrow$$, $$\cdots$$

> In other words, the spin up part is always before the spin down part.

**Code**

The corresponding Fortran code block for the writing of *solver.grn.dat* file is as follows:

```fortran
! open data file: solver.grn.dat
     open(mytmp, file='solver.grn.dat', form='formatted', status='unknown')

! write it
     do i=1,norbs
         do j=1,mfreq
             write(mytmp,'(i6,5f16.8)') i, rmesh(j), &
              real(grnf(j,i,i)), aimag(grnf(j,i,i)), &
              real(gerr(j,i,i)), aimag(gerr(j,i,i))
         enddo ! over j={1,mfreq} loop
         write(mytmp,*) ! write empty lines
         write(mytmp,*)
     enddo ! over i={1,norbs} loop

! close data file
     close(mytmp)
```

In the **HIBISCUS** component, we provide a Python module to read the output files of quantum impurity solvers. You can use it to read the *solver.grn.dat* file. See [script/u_reader.py](../ch07/reader.md) for more details.