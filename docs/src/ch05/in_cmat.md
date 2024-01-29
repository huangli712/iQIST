# atom.cmat.in

**Introduction**

The purpose of the *atom.cmat.in* file is to define the crystal field splitting of the atomic Hamiltonian. The **JASMINE** component supports both diagonal and non-diagonal crystal field matrix. Only when *ibasis* = 1, the *atom.cmat.in* file is read. See [ibasis](p_ibasis.md) parameter for more details.

**Format**

The format of the *atom.cmat.in* file is as follows:

---

*column 1*: orbital index ``\alpha``, integer

*column 2*: orbital index ``\beta``, integer

*column 3*: crystal field splitting matrix element ``\Delta_{\alpha\beta}``, double precision

---

!!! note

    In the *atom.cmat.in* file, we adopt the following orbital sequence:
    ``1\uparrow``, ``2\uparrow``, ``3\uparrow``, ``\cdots``, ``1\downarrow``, ``2\downarrow``, ``3\downarrow``, ``\cdots``
    In other words, the spin up part is always before the spin down part.

**Code**

The corresponding Fortran code block for the reading of the *atom.cmat.in* file is as follows:

```fortran
! open file atom.cmat.in
open(mytmp, file='atom.cmat.in', form='formatted', status='unknown')

! read the data until EOF
do
    read(mytmp,*,iostat = ierr) i1, i2, raux
    if ( ierr == iostat_end ) EXIT
    !
    ! crystal field splitting is actually real
    call s_assert( i1 <= norbs .and. i2 <= norbs )
    cmat(i1,i2) = dcmplx(raux, zero)
enddo ! over do while loop

! close data file
close(mytmp)
```
