# solver.hist.dat

**Introduction**

The *solver.hist.dat* file is used to store the histogram for the perturbation expansion orders of the continuous-time quantum Monte Carlo impurity solvers. It will be output by the quantum impurity solvers periodically (every *nwrite* Monte Carlo sampling steps, see [nwrite](p_nwrite.md) for more details).

**Format**

The *solver.hist.dat* file only contains one block. The format of the block is as follows:

---

*column 1*: index of perturbation expansion order, integer

*column 2*: count of the perturbation expansion order, integer

*column 3*: percent of the perturbation expansion order, double precision

*column 4*: error bar of the histogram

---

**Code**

The corresponding Fortran code block for the writing of *solver.hist.dat* file is as follows:

```fortran
! open data file: solver.hist.dat
     open(mytmp, file='solver.hist.dat', form='formatted', status='unknown')

! write it
     write(mytmp,'(a)') '# histogram: order | count | percent'
     do i=1,mkink
         write(mytmp,'(i6,i12,2f12.6)') i-1, int( hint(i) ), haux(i), htmp(i)
     enddo ! over i={1,mkink} loop

! close data file
     close(mytmp)
```

Now we provide a Python module to read the output files of quantum impurity solvers. You can use it to read the *solver.hist.dat* file. See [script/u_reader.py](../ch06/reader.md) for more details.
