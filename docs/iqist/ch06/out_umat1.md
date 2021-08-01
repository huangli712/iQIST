### solver.umat.in

**Introduction**

The *solver.umat.in* file contains the density-density part of the general Coulomb interaction matrix. Some CT-HYB quantum impurity solvers can read it. It is compatible with the **AZALEA**, **GARDENIA**, **NARCISSUS** components. These components may utilize it to redefine the internal (default) Coulomb interaction matrix.

**Format**

The format of the *solver.umat.in* file is as follows:

---

*column 1*: orbital index $$i$$, integer

*column 2*: orbital index $$j$$, integer

*column 3*: Coulomb interaction matrix element $$U(i,j)$$, double precision

---

> NOTE:

> In the *solver.umat.in* file, we adopt the following orbital sequence:

> $$1\uparrow$$, $$2\uparrow$$, $$3\uparrow$$, $$\cdots$$, $$1\downarrow$$, $$2\downarrow$$, $$3\downarrow$$, $$\cdots$$

> In other words, the spin up part is always before the spin down part.

**Code**

The corresponding Fortran code block for the writing of the *solver.umat.in* file is as follows:

```fortran
! open file atom.umat.dat to write
     open(mytmp, file='solver.umat.in', form='formatted', status='unknown')

! write the data, all of the elements are outputed
! note: we have to change the spin sequence here
     do i=1,norbs
         if ( i <= nband ) then
             k = 2*i-1
         else
             k = 2*(i-nband)
         endif ! back if ( i <= nband ) block

         do j=1,norbs
             if ( j <= nband ) then
                 l = 2*j-1
             else
                 l = 2*(j-nband)
             endif ! back if ( j <= nband ) block

             write(mytmp,'(2i6,f16.8)') i, j, umat_t(k,l)
         enddo ! over j={1,norbs} loop
     enddo ! over i={1,norbs} loop

! close data file
     close(mytmp)
```

See also [solver.umat.in](../ch04/in_umat.md) for more details.