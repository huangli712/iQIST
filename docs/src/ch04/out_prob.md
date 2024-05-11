# solver.prob.dat

**Introduction**

The *solver.prob.dat* file is used to store the atomic state probability and related information about the quantum impurity solvers. It will be output by the quantum impurity solvers when they are **shut down**. We can use the data in it to analyze the valence fluctuation.

**Format**

The format of the *solver.prob.dat* file is somewhat complex. It contains three blocks. They show the probabilities for atomic states, orbitals, and spins, respectively. In the first block, the fifth column is used to represent the error bar.

!!! note

    For the **MANJUSHAKA** component, there is an additional block for showing the probabilities for sectors/superstates.

**Code**

The corresponding Fortran code block for the writing of *solver.prob.dat* file is as follows:

```fortran
! open data file: solver.prob.dat
open(mytmp, file='solver.prob.dat', form='formatted', status='unknown')

! write it
write(mytmp,'(a)') '# state probability: index | prob | occupy | spin'
do i=1,ncfgs
    write(mytmp,'(i6,4f12.6)') i, prob(i), noccs(i), soccs(i) * half, perr(i)
enddo ! over i={1,ncfgs} loop

write(mytmp,'(a)') '# orbital probability: index | occupy | prob'
do i=0,norbs
    write(mytmp,'(i6,3f12.6)') i + 1, real(i), oprob(i), operr(i)
enddo ! over i={0,norbs} loop
write(mytmp,'(a6,12X,f12.6)') 'sum', sum(oprob)

write(mytmp,'(a)') '# spin probability: index | spin | prob'
do i=-nband,nband
    write(mytmp,'(i6,3f12.6)') i + nband + 1, i * half, sprob(i), sperr(i)
enddo ! over i={-nband,nband} loop
write(mytmp,'(a6,12X,f12.6)') 'sum', sum(sprob)

! close data file
close(mytmp)
```

Now we provide a Python module to read the output files of quantum impurity solvers. You can use it to read the *solver.prob.dat* file. See [script/u_reader.py](../ch06/reader.md) for more details.
