### solver.schi.dat

**Introduction**

The *solver.schi.dat* file is used to store the spin-spin correlation function in time space, ``\chi_{\text{spin}}(\tau) = \langle S_z(0)S_z(\tau)\rangle``. It will be output by the quantum impurity solvers when they are **shut down**.

!!! note

    Only the **GARDENIA** and **NARCISSUS** components can generate the *solver.schi.dat* file.

**Format**

The *solver.schi.dat* file contains a few (*nband* + 2) blocks. The first *nband* blocks are orbital-resolved spin-spin correlation functions. The next block is the total spin-spin correlation function. The final block is the sum of orbital-resolved spin-spin correlation functions. Each block is appended by two blank lines. In each block, the error bar data are always shown in the rightmost column.

**Code**

The corresponding Fortran code block for the writing of *solver.schi.dat* file is as follows:

```fortran
! open data file: solver.schi.dat
     open(mytmp, file='solver.schi.dat', form='formatted', status='unknown')

! write it
     do j=1,nband
         write(mytmp,'(a,i6)') '# flvr:', j
         do i=1,ntime
             write(mytmp,'(3f12.6)') tmesh(i), sschi(i,j), sserr(i,j)
         enddo ! over i={1,ntime} loop
         write(mytmp,*) ! write empty lines
         write(mytmp,*)
     enddo ! over j={1,nband} loop

     write(mytmp,'(a,i6)') '# flvr:', 8888
     do i=1,ntime
         write(mytmp,'(3f12.6)') tmesh(i), schi(i), serr(i) 
     enddo ! over i={1,ntime} loop
     write(mytmp,*) ! write empty lines
     write(mytmp,*)

     write(mytmp,'(a,i6)') '# flvr:', 9999
     do i=1,ntime
         write(mytmp,'(3f12.6)') tmesh(i), sum( sschi(i,:) ), sum( sserr(i,:) )
     enddo ! over i={1,ntime} loop
     write(mytmp,*) ! write empty lines
     write(mytmp,*)

! close data file
     close(mytmp)
```

In the **HIBISCUS** component, we provide a Python module to read the output files of quantum impurity solvers. You can use it to read the *solver.schi.dat* file. See [script/u_reader.py](../ch07/reader.md) for more details.