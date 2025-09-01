# solver.sp_w.dat

**Introduction**

The *solver.sp_w.dat* file is used to store the spin-spin correlation function in frequency space, ``\chi_{\text{spin}}(i\nu_n)``. It will be output by the quantum impurity solvers when they are **shut down**.

!!! note

    Only the **NARCISSUS** component can generate the *solver.sp_w.dat* file.

**Format**

The *solver.sp_w.dat* file contains *nband* blocks. Each block is appended by two blank lines. The format of each block is as follows:

---

*column 1*: Matsubara frequency point (bosonic type), ``\nu_n``, double precision

*column 2*: spin-spin correlation function, ``\Re \chi_{\text{spin}}(i\nu_n)``, double precision

*column 3*: error bar, ``\Re [\delta\chi_{\text{spin}}(i\nu_n)]``, double precision

---

!!! note

    The imaginary part of ``\chi_{\text{spin}}(i\nu_n)`` should be kept zero.

**Code**

The corresponding Fortran code block for the writing of *solver.sp_w.dat* file is as follows:

```fortran
! open data file: solver.sp_w.dat
     open(mytmp, file='solver.sp_w.dat', form='formatted', status='unknown')

! write it
     do j=1,nband
         write(mytmp,'(a,i6)') '# flvr:', j
         do i=1,nbfrq
             write(mytmp,'(3f12.6)') bmesh(i), ssfom(i,j), sserr(i,j)
         enddo ! over i={1,nbfrq} loop
         write(mytmp,*) ! write empty lines
         write(mytmp,*)
     enddo ! over j={1,nband} loop

! close data file
     close(mytmp)
```

In the iQIST software package, we provide a Python module to read the output files of quantum impurity solvers. You can use it to read the *solver.sp_w.dat* file. See [u_reader.py] for more details.
