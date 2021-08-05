### solver.sgm.dat

**Introduction**

This file is used to store the Matsubara self-energy function ``\Sigma(i\omega_n)``. It will be output by the quantum impurity solvers when they are **shut down**.

**Format**

The *solver.sgm.dat* file contains *norbs* block. Each block is appended by two blank lines. The format of each block is as follows:

---

*column 1*: orbital index ``i``, integer

*column 2*: Matsubara frequency point, ``\omega_n``, double precision

*column 3*: Matsubara self-energy function, ``\Re \Sigma(i\omega_n)``, double precision

*column 4*: Matsubara self-energy function, ``\Im \Sigma(i\omega_n)``, double precision

*column 5*: error bar, ``\Re [\delta \Sigma(i\omega_n)]``, double precision

*column 6*: error bar, ``\Im [\delta \Sigma(i\omega_n)]``, double precision

---

!!! note

    In the *solver.sgm.dat* file, we adopt the following orbital sequence:
    $$1\uparrow$$, $$2\uparrow$$, $$3\uparrow$$, $$\cdots$$, $$1\downarrow$$, $$2\downarrow$$, $$3\downarrow$$, $$\cdots$$
    In other words, the spin up part is always before the spin down part.

**Code**

The corresponding Fortran code block for the writing of *solver.sgm.dat* file is as follows:

```fortran
! open data file: solver.sgm.dat
     open(mytmp, file='solver.sgm.dat', form='formatted', status='unknown')

! write it
     do i=1,norbs
         do j=1,mfreq
             write(mytmp,'(i6,5f16.8)') i, rmesh(j), &
              real(sigf(j,i,i)), aimag(sigf(j,i,i)), &
                                         zero, zero
         enddo ! over j={1,mfreq} loop
         write(mytmp,*) ! write empty lines
         write(mytmp,*)
     enddo ! over i={1,norbs} loop

! close data file
     close(mytmp)
```

!!! note

    The columns for the error bar are always zero in this file.

In the **HIBISCUS** component, we provide a Python module to read the output files of quantum impurity solvers. You can use it to read the *solver.sgm.dat* file. See [script/u_reader.py](../ch07/reader.md) for more details.