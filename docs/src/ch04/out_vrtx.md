### solver.vrtx.dat

**Introduction**

The *solver.vrtx.dat* file is designed to store the two-particle Green's function $$\chi(i\omega_n, i\omega'_n, i\nu_n)$$ and vertex function $$\mathcal{F}(i\omega_n, i\omega'_n, i\nu_n)$$. It will be output by the quantum impurity solvers when they are **shut down**.

> NOTE:

> 1. Only the **GARDENIA** and **NARCISSUS** components can generate the *solver.vrtx.dat* file.

> 2. The data stored in the *solver.twop.dat* file are generated using standard algorithm. However, those data in the *solver.vrtx.dat* file are generated using the improved estimator algorithm. Unfortunately, as for the **CAMELLIA**, **LAVENDER**, and **MANJUSHAKA** components, they don't support the improved estimator algorithm, so they couldn't generate the *solver.vrtx.dat* file.

**Format**

The *solver.vrtx.dat* file contains 

$$ \left(\sum\limits^{\text{norbs}}_{n=1} \frac{n(n+1)}{2}\right)\times \text{nbfrq}$$

blocks. Each block is appended by two blank lines. The format of each block is as follows:

---

*column 1*: Matsubara frequency point (fermionic type), $$\omega_n$$, integer, the unit is $$\pi/\beta$$

*column 2*: Matsubara frequency point (fermionic type), $$\omega'_n$$, integer, the unit is $$\pi/\beta$$

*column 3*: two-particle Green's function, $$\chi(i\omega_n, i\omega'_n, i\nu_n)$$, double precision

*column 4*: bubble function, $$\chi_0(i\omega_n, i\omega'_n, i\nu_n)$$, double precision

*column 5*: irreducible part of the two-particle Green's function, $$\chi_\text{irr}(i\omega_n, i\omega'_n, i\nu_n)$$, double precision

*column 6*: full vertex function, $$\mathcal{F}(i\omega_n, i\omega'_n, i\nu_n)$$, double precision

---

> NOTE:

> In the *solver.vrtx.dat* file, we adopt the following orbital sequence:

> $$1\uparrow$$, $$2\uparrow$$, $$3\uparrow$$, $$\cdots$$, $$1\downarrow$$, $$2\downarrow$$, $$3\downarrow$$, $$\cdots$$

> In other words, the spin up part is always before the spin down part.

**Code**

The corresponding Fortran code block for the writing of *solver.vrtx.dat* file is as follows:

```fortran
! open data file: solver.vrtx.dat
     open(mytmp, file='solver.vrtx.dat', form='formatted', status='unknown')

! write it
     do m=1,norbs
         do n=1,m
             do k=1,nbfrq
                 write(mytmp,'(a,i6)') '# flvr1:', m
                 write(mytmp,'(a,i6)') '# flvr2:', n
                 write(mytmp,'(a,i6)') '# nbfrq:', k
                 do j=1,nffrq

                     ......

                     do i=1,nffrq

                         ......

                         write(mytmp,'(2i6,8f16.8)') jt, it, chit, chi0, chii, chii/(g1*g2*g3*g4)
                     enddo ! over i={1,nffrq} loop
                 enddo ! over j={1,nffrq} loop
                 write(mytmp,*) ! write empty lines
                 write(mytmp,*)
             enddo ! over k={1,nbfrq} loop
         enddo ! over n={1,m} loop
     enddo ! over m={1,norbs} loop

! close data file
     close(mytmp)
```

In the **HIBISCUS** component, we provide a Python module to read the output files of quantum impurity solvers. You can use it to read the *solver.vrtx.dat* file. See [script/u_reader.py](../ch07/reader.md) for more details.