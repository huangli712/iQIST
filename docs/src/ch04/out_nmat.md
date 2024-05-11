# solver.nmat.dat

**Introduction**

The *solver.nmat.dat* file is used to memory the impurity occupancy and double occupancy matrix. It will be output by the quantum impurity solvers when they **stop working**.

**Format**

The *solver.nmat.dat* file contains two blocks. One is for the impurity occupancy ``\langle n_i\rangle``, and another one is for the double occupancy matrix ``\langle n_i n_j\rangle``. In each block, the error bar data are always shown in the rightmost column.

!!! note

    In the *solver.nmat.dat* file, we adopt the following orbital sequence:
    $$1\uparrow$$, $$2\uparrow$$, $$3\uparrow$$, $$\cdots$$, $$1\downarrow$$, $$2\downarrow$$, $$3\downarrow$$, $$\cdots$$
    In other words, the spin up part is always before the spin down part.

**Code**

The corresponding Fortran code block for the writing of *solver.nmat.dat* file is as follows:

```fortran
! open data file: solver.nmat.dat
     open(mytmp, file='solver.nmat.dat', form='formatted', status='unknown')

! write it
     write(mytmp,'(a)') '#   < n_i >   data:'
     do i=1,norbs
         write(mytmp,'(i6,2f12.6)') i, nmat(i), nerr(i)
     enddo ! over i={1,norbs} loop
     write(mytmp,'(a6,2f12.6)') 'sup', sum( nmat(1:nband) ), sum( nerr(1:nband) )
     write(mytmp,'(a6,2f12.6)') 'sdn', sum( nmat(nband+1:norbs) ), sum( nerr(nband+1:norbs) )
     write(mytmp,'(a6,2f12.6)') 'sum', sum( nmat(1:norbs) ), sum( nerr(1:norbs) )

     write(mytmp,'(a)') '# < n_i n_j > data:'
     do i=1,norbs
         do j=1,norbs
             write(mytmp,'(2i6,2f12.6)') i, j, nnmat(i,j), nnerr(i,j)
         enddo ! over j={1,norbs} loop
     enddo ! over i={1,norbs} loop

! close data file
     close(mytmp)
```

Now we provide a Python module to read the output files of quantum impurity solvers. You can use it to read the *solver.nmat.dat* file. See [script/u_reader.py](../ch06/reader.md) for more details.
