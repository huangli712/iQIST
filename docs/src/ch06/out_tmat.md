### atom.tmat.dat

**Introduction**

The *atom.tmat.dat* file contains the transformation matrix $$\mathcal{T}_{\alpha,\beta}$$ from the original basis to natural basis.

**Format**

The format of the *atom.tmat.dat* file is as follows:

---

*column 1*: orbital index $$\alpha$$, integer

*column 2*: orbital index $$\beta$$, integer

*column 3*: Elements of the transformation matrix $$\mathcal{T}_{\alpha,\beta}$$, real part, double precision

*column 4*: Elements of the transformation matrix $$\mathcal{T}_{\alpha,\beta}$$, imaginary part, double precision

---

!!! note

In the *atom.tmat.dat* file, we adopt the following orbital sequence:
$$1\uparrow$$, $$2\uparrow$$, $$3\uparrow$$, $$\cdots$$, $$1\downarrow$$, $$2\downarrow$$, $$3\downarrow$$, $$\cdots$$
In other words, the spin up part is always before the spin down part.

**Code**

The corresponding Fortran code block for the writing of the *atom.tmat.dat* file is as follows:

```fortran
! open file atom.tmat.dat to write
     open(mytmp, file='atom.tmat.dat', form='formatted', status='unknown')

! write the header
     write(mytmp,'(75a1)') dash ! dashed line
     write(mytmp,'(a)') '# i | j | tmat_real | tmat_imag'
     write(mytmp,'(75a1)') dash ! dashed line

! write the data
     do i=1,norbs
         do j=1,norbs
             write(mytmp,'(2i6,2f16.8)') i, j, tmat(i,j)
         enddo ! over j={1,norbs} loop
     enddo ! over i={1,norbs} loop

! close data file
     close(mytmp)
```

In principle, the transformation matrix is complex.

See also [atom.tmat.in](in_tmat.md) for more details.