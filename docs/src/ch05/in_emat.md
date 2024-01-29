# atom.emat.in

**Introduction**

The purpose of the *atom.emat.in* file is to supply the on-site impurity energy level ``E_{\alpha\beta}``, which is in fact the *CFS* + *SOC*. It is a diagonal matrix in the natural basis. Only when *ibasis* = 2, the *atom.emat.in* file is used. See [ibasis](p_ibasis.md) for more details.

!!! note

    1. The so-called natural basis is the eigenstates of *CFS* + *SOC* matrix.
    2. Only the diagonal elements of ``E_{\alpha\beta}`` are included in the *atom.emat.in* file.

**Format**

The format of the *atom.emat.in* file is as follows:

---

*column 1*: orbital index ``\alpha``, integer

*column 2*: orbital index ``\beta``, integer

*column 3*: on-site impurity energy level ``E_{\alpha,\beta}``, double precision

---

!!! note

    In the *atom.emat.in* file, we adopt the following orbital sequence:
    ``1\uparrow``, ``2\uparrow``, ``3\uparrow``, ``\cdots``, ``1\downarrow``, ``2\downarrow``, ``3\downarrow``, ``\cdots``
    In other words, the spin up part is always before the spin down part.

**Code**

The corresponding Fortran code block for the reading of the *atom.emat.in* file is as follows:

```fortran
! open file atom.emat.in
open(mytmp, file='atom.emat.in', form='formatted', status='unknown')

! read the data file
do i=1,norbs
    read(mytmp,*) i1, i2, raux
    ! emat is actually real and diagonal in natural basis
    emat(i,i) = dcmplx(raux, zero)
enddo ! over i={1,norbs} loop

! close data file
close(mytmp)
```
