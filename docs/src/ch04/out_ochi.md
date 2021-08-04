### solver.ochi.dat

**Introduction**

The *solver.ochi.dat* file is used to store the orbital-orbital correlation function in time space, $$\chi_{\text{charge}}(\tau) = \langle n(0)n(\tau)\rangle$$. It will be output by the quantum impurity solvers when they are **shut down**.

!!! note

    Only the **GARDENIA** and **NARCISSUS** components can generate the *solver.ochi.dat* file.

**Format**

The *solver.ochi.dat* file contains a few (*norbs*$$\times$$*norbs* + 2) blocks. The first *norbs*$$\times$$*norbs* blocks are orbital-resolved orbital-orbital correlation functions. The next block is the total orbital-orbital correlation function. The final block is the sum of orbital-resolved orbital-orbital correlation functions. Each block is appended by two blank lines. In each block, the error bar data are always shown in the rightmost column.

!!! note

    In the *solver.ochi.dat* file, we adopt the following orbital sequence:
    $$1\uparrow$$, $$2\uparrow$$, $$3\uparrow$$, $$\cdots$$, $$1\downarrow$$, $$2\downarrow$$, $$3\downarrow$$, $$\cdots$$
    In other words, the spin up part is always before the spin down part.

**Code**

The corresponding Fortran code block for the writing of *solver.ochi.dat* file is as follows:

```fortran
! open data file: solver.ochi.dat
     open(mytmp, file='solver.ochi.dat', form='formatted', status='unknown')

! write it
     do k=1,norbs
         do j=1,norbs
             write(mytmp,'(2(a,i6))') '# flvr:', j, '  flvr:', k
             do i=1,ntime
                 write(mytmp,'(3f12.6)') tmesh(i), oochi(i,j,k), ooerr(i,j,k)
             enddo ! over i={1,ntime} loop
             write(mytmp,*) ! write empty lines
             write(mytmp,*)
         enddo ! over j={1,norbs} loop
     enddo ! over k={1,norbs} loop

     write(mytmp,'(a,i6)') '# flvr:', 8888
     do i=1,ntime
         write(mytmp,'(3f12.6)') tmesh(i), ochi(i), oerr(i) 
     enddo ! over i={1,ntime} loop
     write(mytmp,*) ! write empty lines
     write(mytmp,*)

     write(mytmp,'(a,i6)') '# flvr:', 9999
     do i=1,ntime
         write(mytmp,'(3f12.6)') tmesh(i), sum( oochi(i,:,:) ), sum( ooerr(i,:,:) ) 
     enddo ! over i={1,ntime} loop
     write(mytmp,*) ! write empty lines
     write(mytmp,*)

! close data file
     close(mytmp)
```

In the **HIBISCUS** component, we provide a Python module to read the output files of quantum impurity solvers. You can use it to read the *solver.ochi.dat* file. See [script/u_reader.py](../ch07/reader.md) for more details.