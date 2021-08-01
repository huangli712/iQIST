### atom.umat.dat

**Introduction**

The *atom.umat.dat* file contains the general (four-fermions) Coulomb interaction matrix. Only the non-zero matrix elements (> 1.0E-10) are outputted.

**Format**

The format of the *atom.umat.dat* file is as follows:

---

*column 1*: orbital index $$\alpha$$, integer

*column 2*: orbital index $$\beta$$, integer

*column 3*: orbital index $$\gamma$$, integer

*column 4*: orbital index $$\delta$$, integer

*column 5*: Elements of general Coulomb matrix $$U_{\alpha\beta\gamma\delta}$$, real part, double precision

*column 6*: Elements of general Coulomb matrix $$U_{\alpha\beta\gamma\delta}$$, imaginary part, double precision

---

> NOTE:

> In the *atom.umat.dat* file, we adopt the following orbital sequence:

> $$1\uparrow$$, $$2\uparrow$$, $$3\uparrow$$, $$\cdots$$, $$1\downarrow$$, $$2\downarrow$$, $$3\downarrow$$, $$\cdots$$

> In other words, the spin up part is always before the spin down part.

**Code**

The corresponding Fortran code block for the writing of the *atom.umat.dat* file is as follows:

```fortran
! open file atom.umat.dat to write
     open(mytmp, file='atom.umat.dat', form='formatted', status='unknown')

! write the header
     write(mytmp,'(75a1)') dash ! dashed line
     write(mytmp,'(a)') '# i | j | k | l | umat_real | umat_imag'
     write(mytmp,'(75a1)') dash ! dashed line

! write the data, only the non-zero elements are outputed
! note: we do not change the spin sequence here
     do i=1,norbs
         do j=1,norbs
             do k=1,norbs
                 do l=1,norbs
                     if ( abs( umat(i,j,k,l) ) > epst ) then
                         write(mytmp,'(4i6,2f16.8)') i, j, k, l, umat(i,j,k,l)
                     endif ! back if ( abs( umat(i,j,k,l) ) > epst ) block
                 enddo ! over l={1,norbs} loop
             enddo ! over k={1,norbs} loop
         enddo ! over j={1,norbs} loop
     enddo ! over i={1,norbs} loop

! close data file
     close(mytmp)
```

In general, the Coulomb interaction matrix is complex.

See also [solver.umat.in](out_umat1.md) for more details.