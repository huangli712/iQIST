### solver.lmat.dat

**Introduction**

The *solver.lmat.dat* file is used to store the necessary data used to build the generalized fidelity susceptibility ``\chi_{\text{FS}}``, which can be calculated via ``\chi_{\text{FS}} = \langle k_L k_R \rangle - \langle k_L \rangle \langle k_R \rangle``. It will be output by the quantum impurity solvers when they are **shut down**.

!!! note

    Only the **GARDENIA**, **NARCISSUS**, **CAMELLIA**, **LAVENDER**, and **MANJUSHAKA** components can generate the *solver.lmat.dat* file.

**Format**

The *solver.lmat.dat* file contains two blocks. One is for the ``\langle k_L \rangle`` and ``\langle k_R \rangle`` data, and another one is for the ``\langle k_L k_R \rangle``. In each block, the error bar data are always shown in the rightmost column.

!!! note

    In the *solver.lmat.dat* file, we adopt the following orbital sequence:
    $$1\uparrow$$, $$2\uparrow$$, $$3\uparrow$$, $$\cdots$$, $$1\downarrow$$, $$2\downarrow$$, $$3\downarrow$$, $$\cdots$$
    In other words, the spin up part is always before the spin down part.

**Code**

The corresponding Fortran code block for the writing of *solver.lmat.dat* file is as follows:

```fortran
! open data file: solver.lmat.dat
     open(mytmp, file='solver.lmat.dat', form='formatted', status='unknown')

! write it
     write(mytmp,'(a)') '# < k_l > < k_r > data:'
     do i=1,norbs
         write(mytmp,'(i6,4f12.6)') i, lmat(i), rmat(i), lerr(i), rerr(i)
     enddo ! over i={1,norbs} loop
     write(mytmp,'(a6,2f12.6)') 'l_sum', sum( lmat ), sum( lerr )
     write(mytmp,'(a6,2f12.6)') 'r_sum', sum( rmat ), sum( rerr )

     write(mytmp,'(a)') '# < k_l k_r > data:'
     do i=1,norbs
         do j=1,norbs
             write(mytmp,'(2i6,2f12.6)') i, j, lrmat(i,j), lrerr(i,j)
         enddo ! over j={1,norbs} loop
     enddo ! over i={1,norbs} loop
     write(mytmp,'(a6,2f12.6)') 'lrsum', sum( lrmat ), sum( lrerr )
     write(mytmp,'(a6,2f12.6)') 'fidel', f_val, f_err

! close data file
     close(mytmp)
```

In the **HIBISCUS** component, we provide a Python module to read the output files of quantum impurity solvers. You can use it to read the *solver.lmat.dat* file. See [script/u_reader.py](../ch07/reader.md) for more details.