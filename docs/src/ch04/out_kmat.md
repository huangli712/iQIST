# solver.kmat.dat

**Introduction**

The *solver.kmat.dat* is used to store the orbital-dependent kinetic energy data (Here we use the expectation value of the perturbation expansion order ``\langle k \rangle`` to represent the kinetic energy). It will be output by the quantum impurity solvers when they are **shut down**.

!!! note

    Only the **NARCISSUS** and **MANJUSHAKA** components can generate the *solver.kmat.dat* file.

**Format**

The *solver.kmat.dat* file contains two blocks. One is for the orbital-dependent kinetic energy ``\langle k_i\rangle``, and another one is for the double kinetic energy matrix ``\langle k_i k_j\rangle``. In each block, the error bar data are always shown in the rightmost column.

!!! note

    In the *solver.kmat.dat* file, we adopt the following orbital sequence:
    $$1\uparrow$$, $$2\uparrow$$, $$3\uparrow$$, $$\cdots$$, $$1\downarrow$$, $$2\downarrow$$, $$3\downarrow$$, $$\cdots$$
    In other words, the spin up part is always before the spin down part.

**Code**

The corresponding Fortran code block for the writing of *solver.kmat.dat* file is as follows:

```fortran
! open data file: solver.kmat.dat
     open(mytmp, file='solver.kmat.dat', form='formatted', status='unknown')

! write it
     write(mytmp,'(a)') '# <  k  > data:'
     do i=1,norbs
         write(mytmp,'(i6,2f12.6)') i, kmat(i), kerr(i)
     enddo ! over i={1,norbs} loop
     write(mytmp,'(a6,2f12.6)') 'k_sum', sum( kmat ), sum( kerr )

     write(mytmp,'(a)') '# < k^2 > data:'
     do i=1,norbs
         do j=1,norbs
             write(mytmp,'(2i6,2f12.6)') i, j, kkmat(i,j), kkerr(i,j)
         enddo ! over j={1,norbs} loop
     enddo ! over i={1,norbs} loop
     write(mytmp,'(a6,2f12.6)') 'kksum', sum( kkmat ), sum( kkerr )
     write(mytmp,'(a6,2f12.6)') 'final', f_val, f_err

! close data file
     close(mytmp)
```

In the iQIST software package, we provide a Python module to read the output files of quantum impurity solvers. You can use it to read the *solver.kmat.dat* file. See [u_reader.py] for more details.
