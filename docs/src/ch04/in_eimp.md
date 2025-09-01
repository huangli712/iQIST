# solver.eimp.in

**Introduction**

The *solver.eimp.in* file is used to define impurity level matrix ``E_{\alpha\beta}`` and symmetry matrix *symm(``\alpha``,``\beta``)*. Unfortunately, now only the diagonal elements of the two matrices are supported.

If the orbitals are degenerated, the impurity levels should be a constant and can be absorbed into the chemical potential ``\mu``. If the *isbnd* parameter is set to 2, the orbital-resolved observables are then symmetrized according the symmetry matrix. See [isbnd](p_isbnd.md) for more details.

---

**Format**

The format of the *solver.eimp.in* file is as follows:

>
> *column 1*: orbital index, ``\alpha``, integer
>
> *column 2*: diagonal element of the impurity level, ``E_{\alpha\alpha}``, double precision
>
> *column 3*: symmetry vector, *symm(``\alpha``)*, integer
>

!!! tip

    In the *solver.eimp.in* file, we adopt the following orbital sequence:
    $$1\uparrow$$, $$2\uparrow$$, $$3\uparrow$$, $$\cdots$$, $$1\downarrow$$, $$2\downarrow$$, $$3\downarrow$$, $$\cdots$$
    In other words, the spin up part is always before the spin down part.

---

**Code**

The corresponding Fortran code block for the reading of *solver.eimp.in* file is as follows:

```fortran
open(mytmp, file='solver.eimp.in', form='formatted', status='unknown')
do i=1,norbs
    read(mytmp,*) k, eimp(i), symm(i)
enddo ! over i={1,norbs} loop
close(mytmp)
```

Usually, you have to edit the *solver.eimp.in* file by yourself.
