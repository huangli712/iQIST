### atom.emat.dat

**Introduction**

The *atom.emat.dat* file contains the on-site impurity level energy on natural basis.

**Format**

The format of the *atom.emat.dat* file is as follows:

---

*column 1*: orbital index $$\alpha$$, integer

*column 2*: On-site impurity energy level $$E_{\alpha,\alpha}$$, real part, double precision

*column 3*: On-site impurity energy level $$E_{\alpha,\alpha}$$, imaginary part, double precision

---

> NOTE:

> In the *atom.emat.dat* file, we adopt the following orbital sequence:

> $$1\uparrow$$, $$2\uparrow$$, $$3\uparrow$$, $$\cdots$$, $$1\downarrow$$, $$2\downarrow$$, $$3\downarrow$$, $$\cdots$$

> In other words, the spin up part is always before the spin down part.

**Code**

The corresponding Fortran code block for the writing of the *atom.emat.dat* file is as follows:

```fortran
! open file atom.emat.dat to write
     open(mytmp, file='atom.emat.dat', form='formatted', status='unknown')

! write the header
     write(mytmp,'(75a1)') dash ! dashed line
     write(mytmp,'(a)') '# i | emat_real | emat_imag'
     write(mytmp,'(75a1)') dash ! dashed line

! write the data
     do i=1,norbs
         if ( isoc == 0 ) then
             if ( i <= nband ) then
                 s_order = 2*i-1
             else
                 s_order = 2*(i-nband)
             endif ! back if ( i <= nband ) block
         else
             s_order = i
         endif ! back if ( isoc == 0 ) block
         write(mytmp,'(i6,2f16.8)') i, emat(s_order,s_order)
     enddo ! over i={1,norbs} loop

! close data file
     close(mytmp)
```

In principle, it is a complex matrix. Here we only write its diagonal part.

See also [atom.emat.in](in_emat.md) for more details.