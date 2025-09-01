# solver.ch_w.dat

**Introduction**

The *solver.ch_w.dat* file is used to store the orbital-orbital correlation function in frequency space, ``\chi_{\text{charge}}(i\nu_n)``. It will be output by the quantum impurity solvers when they are **shut down**.

!!! note

    Only the **NARCISSUS** component can generate the *solver.ch_w.dat* file.

**Format**

The *solver.ch_w.dat* file contains *norbs*``\times``*norbs* blocks. Each block is appended by two blank lines. The format of each block is as follows:

---

*column 1*: Matsubara frequency point (bosonic type), ``\nu_n``, double precision

*column 2*: orbital-orbital correlation function, ``\Re \chi_{\text{charge}}(i\nu_n)``, double precision

*column 3*: error bar, ``\Re [\delta\chi_{\text{charge}}(i\nu_n)]``, double precision

---

!!! note

    The imaginary part of ``\chi_{\text{charge}}(i\nu_n)`` should keeps zero.

!!! note

    In the *solver.ch_w.dat* file, we adopt the following orbital sequence:
    $$1\uparrow$$, $$2\uparrow$$, $$3\uparrow$$, $$\cdots$$, $$1\downarrow$$, $$2\downarrow$$, $$3\downarrow$$, $$\cdots$$
    In other words, the spin up part is always before the spin down part.

**Code**

The corresponding Fortran code block for the writing of *solver.ch_w.dat* file is as follows:

```fortran
! open data file: solver.ch_w.dat
     open(mytmp, file='solver.ch_w.dat', form='formatted', status='unknown')

! write it
     do k=1,norbs
         do j=1,norbs
             write(mytmp,'(2(a,i6))') '# flvr:', j, '  flvr:', k
             do i=1,nbfrq
                 write(mytmp,'(3f12.6)') bmesh(i), oofom(i,j,k), ooerr(i,j,k)
             enddo ! over i={1,nbfrq} loop
             write(mytmp,*) ! write empty lines
             write(mytmp,*)
         enddo ! over j={1,norbs} loop
     enddo ! over k={1,norbs} loop

! close data file
     close(mytmp)
```

In the iQIST software package, we provide a Python module to read the output files of quantum impurity solvers. You can use it to read the *solver.ch_w.dat* file. See [u_reader.py] for more details.
