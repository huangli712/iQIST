### solver.hub.dat

**Introduction**

This file is used to save the atomic self-energy function ``\Sigma_{\text{atomic}}(i\omega_n)`` and Green's function ``G_{\text{atomic}}(i\omega_n)``, which are calculated using the Hubbard-I approximation. It will be output by the quantum impurity solvers when they are **shut down**.

**Format**

The *solver.hub.dat* file contains *norbs* block. Each block is appended by two blank lines. The format of each block is as follows:

---

*column 1*: orbital index ``i``, integer

*column 2*: Matsubara frequency point, ``\omega_n``, double precision

*column 3*: atomic Green's function, ``\Re G_{\text{atomic}}(i\omega_n)``, double precision

*column 4*: atomic Green's function, ``\Im G_{\text{atomic}}(i\omega_n)``, double precision

*column 5*: atomic self-energy function, ``\Re \Sigma_{\text{atomic}}(i\omega_n)``, double precision

*column 6*: atomic self-energy function, ``\Im \Sigma_{\text{atomic}}(i\omega_n)``, double precision

---

!!! note

    In the *solver.hub.dat* file, we adopt the following orbital sequence:
    $$1\uparrow$$, $$2\uparrow$$, $$3\uparrow$$, $$\cdots$$, $$1\downarrow$$, $$2\downarrow$$, $$3\downarrow$$, $$\cdots$$
    That is to say, the spin up part is always before the spin down part.

**Code**

The corresponding Fortran code block for the writing of *solver.hub.dat* file is as follows:

```fortran
! open data file: solver.hub.dat
     open(mytmp, file='solver.hub.dat', form='formatted', status='unknown')

! write it
     do i=1,norbs
         do j=1,mfreq
             write(mytmp,'(i6,5f16.8)') i, rmesh(j), &
                  real(ghub(j,i)), aimag(ghub(j,i)), &
                  real(shub(j,i)), aimag(shub(j,i))
         enddo ! over j={1,mfreq} loop
         write(mytmp,*) ! write empty lines
         write(mytmp,*)
     enddo ! over i={1,norbs} loop

! close data file
     close(mytmp)
```

In the **HIBISCUS** component, we provide a Python module to read the output files of quantum impurity solvers. You can use it to read the *solver.hub.dat* file. Refer to [script/u_reader.py](../ch07/reader.md) for more details.