### solver.pair.dat

**Introduction**

The *solver.pair.dat* file is designed to store the particle-particle pairing susceptibility. It will be output by the quantum impurity solvers when they are **shut down**.

> NOTE:

> Only the **GARDENIA**, **NARCISSUS**, **CAMELLIA**, **LAVENDER**, and **MANJUSHAKA** components can generate the *solver.pair.dat* file.

**Format**

The *solver.pair.dat* file contains 

$$ \left(\sum\limits^{\text{norbs}}_{n=1} \frac{n(n+1)}{2}\right)\times \text{nbfrq}$$

blocks. Each block is appended by two blank lines. The format of each block is as follows:

---

*column 1*: Matsubara frequency point (fermionic type), $$\omega_n$$, integer, the unit is $$\pi/\beta$$

*column 2*: Matsubara frequency point (fermionic type), $$\omega'_n$$, integer, the unit is $$\pi/\beta$$

*column 3*: real part of the pairing susceptibility, $$\Re \chi_{\text{pair}}(i\omega_n, i\omega'_n, i\nu_n)$$, double precision

*column 4*: imaginary part of the pairing susceptibility, $$\Im \chi_{\text{pair}}(i\omega_n, i\omega'_n, i\nu_n)$$, double precision

---

> NOTE:

> In the *solver.pair.dat* file, we adopt the following orbital sequence:

> $$1\uparrow$$, $$2\uparrow$$, $$3\uparrow$$, $$\cdots$$, $$1\downarrow$$, $$2\downarrow$$, $$3\downarrow$$, $$\cdots$$

> In other words, the spin up part is always before the spin down part.

**Code**

The corresponding Fortran code block for the writing of *solver.pair.dat* file is as follows:

```fortran
! open data file: solver.pair.dat
     open(mytmp, file='solver.pair.dat', form='formatted', status='unknown')

! write it
     do m=1,norbs
         do n=1,m
             do k=1,nbfrq
                 write(mytmp,'(a,i6)') '# flvr1:', m
                 write(mytmp,'(a,i6)') '# flvr2:', n
                 write(mytmp,'(a,i6)') '# nbfrq:', k
                 do j=1,nffrq
                     do i=1,nffrq
                         ......
                         write(mytmp,'(2i6,2f16.8)') jt, it, ps_re(i,j,k,n,m), ps_im(i,j,k,n,m)
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

In the **HIBISCUS** component, we provide a Python module to read the output files of quantum impurity solvers. You can use it to read the *solver.pair.dat* file. See [script/u_reader.py](../ch07/reader.md) for more details.