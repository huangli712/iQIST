# solver.kernel.dat

**Introduction**

The *solver.kernel.dat* is designed to store the screening function ``K(\tau)`` and it derivates ``K'(\tau)`` in imaginary-time space. It will be output by the quantum impurity solvers when they are **shut down**.

!!! note

    1. The screening function ``K(\tau)`` and its derivates ``K'(\tau)`` are the input to the quantum impurity solvers (see [solver.ktau.in](in_ktau.md) for more details). In other words, they won't be changed by the quantum impurity solvers.
    2. Only the **NARCISSUS** component can generate the *solver.kernel.dat* file.

**Format**

The *solver.kernel.dat* file only contains one block. The format of the block is as follows:

---

*column 1*: index of imaginary-time point, integer

*column 2*: imaginary-time point, ``\tau``, double precision

*column 3*: screening function, ``K(\tau)``, double precision

*column 4*: the first derivate of screening function, ``K'(\tau)``, double precision

*column 5*: the second derivate of screening function, ``K''(\tau)``, double precision

*column 6*: the third derivate of screening function, ``K'''(\tau)``, double precision

---

!!! note

    In the **NARCISSUS** component, the ``K(\tau)`` and ``K'(\tau)`` are assumed to be orbital-independent.

**Code**

The corresponding Fortran code block for the writing of *solver.kernel.dat* file is as follows:

```fortran
! open data file: solver.kernel.dat
     open(mytmp, file='solver.kernel.dat', form='formatted', status='unknown')

! write it
     do i=1,ntime
         write(mytmp,'(i6,5f12.6)') i, tmesh(i), ktau(i), ptau(i), ksed(i), psed(i)
     enddo ! over i={1,ntime} loop

! close data file
     close(mytmp)
```

In the iQIST software package, we provide a Python module to read the output files of quantum impurity solvers. You can use it to read the *solver.kernel.dat* file. See [u_reader.py] for more details.
