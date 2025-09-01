# solver.umat.in

**Introduction**

The quantum impurity solvers will automatically generate the Coulomb interaction matrix using the ``U`` and ``J`` parameters provided by you. But sometimes you may want to customize the Coulomb interaction matrix by yourself.

Is it possible within the iQIST code?

Yes, of course.

You can define your Coulomb interaction matrix in the *solver.umat.in* file. And the iQIST codes will read data from it if it is available. Then the default Coulomb interaction matrix will be replaced with the new one. That's all.

!!! warning

    The continuous-time quantum Monte Carlo impurity solvers in the general matrix representation, i.e., the **MANJUSHAKA** component (in *iqist/src/cthyb*), does not support this file/feature. All the information about the interaction matrix is already encapsulated in the *atom.cix* file.

---

**Format**

The format of the *solver.umat.in* file is as follows:

>
> *column 1*: orbital index ``i``, integer
>
> *column 2*: orbital index ``j``, integer
>
> *column 3*: Coulomb interaction matrix element ``U(i,j)``, double precision
>

!!! tip

    In the *solver.umat.in* file, we employed the following orbital sequence:
    $$1\uparrow$$, $$2\uparrow$$, $$3\uparrow$$, $$\cdots$$, $$1\downarrow$$, $$2\downarrow$$, $$3\downarrow$$, $$\cdots$$
    Namely, we put the spin up part before the spin down part.

---

**Code**

The corresponding Fortran code block for the reading of *solver.umat.in* file is as follows:

```fortran
open(mytmp, file='solver.umat.in', form='formatted', status='unknown')
do i=1,norbs
    do j=1,norbs
        read(mytmp,*) k, l, rtmp
        uumat(k,l) = rtmp
    enddo ! over j={1,norbs} loop
enddo ! over i={1,norbs} loop
close(mytmp)
```

Usually, you have to edit the *solver.umat.in* file by yourself.
