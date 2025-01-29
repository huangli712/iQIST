### solver.anydos.in

**Introduction**

In the CT-QMC impurity solvers contained in the iQIST software package, we provide a mini dynamical mean-field theory engine. This engine implements a self-consistent condition for the Bethe lattice which has a semi-circular density of states with bandwith ``4t``. Sometimes you may want to try the other models with arbitrary density of states.

Is it possible?

Yes, of course. You can define your density of states in the *solver.anydos.in* file. And then you have to hack the *ctqmc\_dmft.f90* file.

Change the following codes

```fortran
......
     use context, only : grnf
......
! calculate new hybridization function using self-consistent condition
     call ctqmc_dmft_bethe(hybf, grnf)
```

to

```fortran
......
     use context, only : grnf, sig2
......
! calculate new hybridization function using self-consistent condition
     call ctqmc_dmft_anydos(hybf, grnf, sig2)
```

You have to be very careful. Finally, recompiling the CT-QMC impurity solvers is necessary.

!!! note

    1. Don't forget to set the *isscf* parameter to 2, or else the CT-QMC impurity solvers will skip the dynamical mean-field theory engine and perform one-shot calculation only.
    2. The HF-QMC impurity solver (the **DAISY** component) does not support this file/feature.

**Format**

The format of the *solver.anydos.in* file is as follows:

---

*column 1*: frequency point, ``\epsilon``, double precision

*column 2*: density of states, ``\rho(\epsilon)``, double precision

---

!!! note

    We assume that:
    1. The orbitals are degenerated.
    2. The number of frequency points is 801.
    If you are not satisfied with these assumptions, you have to hack the corresponding *ctqmc\_dmft.f90* by yourself.

**Code**

The corresponding Fortran code block for the reading of *solver.anydos.in* file is as follows:

```fortran
open(mytmp, file='solver.anydos.in', form='formatted', status='unknown')
do i=1,ngrid
    read(mytmp,*) epsi(i), pdos(i,1)
    do j=2,norbs
        pdos(:,j) = pdos(:,1)
    enddo ! over j={2,norbs} loop
enddo ! over i={1,ngrid} loop
close(mytmp)
```

For some typical models, the *solver.anydos.in* file can be generated using the toolbox/makedos tool which is included in the **HIBISCUS** component. See [toolbox/makedos](../ch07/dos.md) for more details.