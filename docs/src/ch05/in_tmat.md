# atom.tmat.in

**Introduction**

The purpose of the *atom.tmat.in* file is to supply the transformation matrix ``\mathcal{T}_{\alpha\beta}``, which transforms an operator from its original basis to the natural basis. Only when *ibasis* = 2, the *atom.tmat.in* file is used. See [ibasis](p_ibasis.md) for more details.

**Format**

The format of the *atom.tmat.in* file is as follows:

---

*column 1*: orbital index ``\alpha``, integer

*column 2*: orbital index ``\beta``, integer

*column 3*: elements of the transformation matrix ``\mathcal{T}_{\alpha\beta}``, double precision

---

!!! note

    In the *atom.tmat.in* file, we adopt the following orbital sequence:
    $$1\uparrow$$, $$2\uparrow$$, $$3\uparrow$$, $$\cdots$$, $$1\downarrow$$, $$2\downarrow$$, $$3\downarrow$$, $$\cdots$$
    In other words, the spin up part is always before the spin down part.

**Code**

The corresponding Fortran code block for the reading of the *atom.tmat.in* file is as follows:

```fortran
! open file atom.tmat.in
open(mytmp, file='atom.tmat.in', form='formatted', status='unknown')

! read the data file
do i=1,norbs
    do j=1,norbs
        read(mytmp,*) i1, i2, raux
        ! tmat is actually real
        tmat(i,j) = dcmplx(raux, zero)
    enddo ! over j={1,norbs} loop
enddo ! over i={1,norbs} loop

! close data file
close(mytmp)
```

In principle, the transformation matrix is complex. Here we think that its imaginary part is zero and only take the real part into consideration.
