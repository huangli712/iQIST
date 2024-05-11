# solver.diag.dat

**Introduction**

The *solver.diag.dat* file is used to store the diagram configurations of the perturbation expansion series. Owing to the limitation of disk capacity and computational efficiency, it is impossible to save all of the visited diagram configurations. Consequently we only save the current diagram configurations every *nwrite* Monte Carlo step. Namely, the frequency for writing diagram configurations is controlled by the *nwrite* parameter. See also [nwrite](p_nwrite.md) for more details.

**Format**

The *solver.diag.dat* file consists many blocks. The blocks are separated by two blank lines. A typical block looks like as follows:

```
>> cur_iter:   1 tot_iter:  20
>> cur_diag:   4 tot_diag:1000
# flvr:   1 rank:   2
   1      6.49617550      0.35585709
   1      7.30348648      7.28301304
# flvr:   2 rank:   1
   2      0.43805150      6.89040598
```

Here *cur_iter* means current iteration number, *tot_iter* the total iteration number, *cur_diag* the index of the current diagram configuration, *tot_diag* the total number of diagram configurations. Actually, we have

```math
\text{tot\_diag} = \frac{\text{nsweep}}{\text{nwrite}}.
```

Next, the positions for creator/destroy operators are given for all orbitals (flavors).

---

*column 1*: the index of orbitals.

*column 2*: ``\tau_s``, the imaginary-time points for creator operators

*column 3*: ``\tau_e``, the imaginary-time points for destroy operators

---

We can utilize the data in the *solver.diag.dat* file to produce animation movie. A block can be used to generate one frame in the movie.

**Code**

The corresponding Fortran code block for the writing of *solver.diag.dat* file is as follows:

```fortran
! write the snapshot
! open data file: solver.diag.dat
     open(mytmp, file='solver.diag.dat', form='formatted', status='unknown', position='append')

! write diagram info
     write(mytmp,'(2(a,i4))') '>> cur_iter:', iter, ' tot_iter:', niter
     write(mytmp,'(2(a,i4))') '>> cur_diag:', cstep/nwrite, ' tot_diag:', nsweep/nwrite

! write the position of operators
     do i=1,norbs
         write(mytmp,'(2(a,i4))') '# flvr:', i, ' rank:', rank(i)
         do j=1,rank(i)
             write(mytmp,'(i4,2f16.8)') i, time_s( index_s(j, i), i ), time_e( index_e(j, i), i )
         enddo ! over j={1,rank(i)} loop
     enddo ! over i={1,norbs} loop

! write two blank lines
     write(mytmp,*)
     write(mytmp,*)

! close data file
     close(mytmp)
```

You can use the *u_movie.py* script to visualize the diagram configurations in the *solver.diag.dat* file. See also [script/u_movie.py](../ch06/movie.md) for more details.
