# solver.anydos.in

!!! warning

    This feature is already disable. And this page will be removed from this user manual soon.

**Introduction**

In the CT-QMC impurity solvers contained in the iQIST software package, we provide a mini dynamical mean-field theory engine. This engine implements a self-consistent condition for the bethe lattice which has a semi-circular density of states with bandwith ``4t``. Sometimes you may want to try the other models with arbitrary density of states.

Is it possible?

Yes, of course. You can define your density of states in the *solver.anydos.in* file. And then you have to hack the *iqist/src/ctseg/ctqmc\_dmft.f90* or *iqist/src/cthyb/ctqmc\_dmft.f90* file.

Please change the following codes

```fortran
......
     use context, only : grnf
......
     ! apply the self-consistent condition. here we consider a Hubbard model
     ! on a bethe lattice. of course you can replace it with your implements
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

You have to be very careful. Finally, recompiling the CT-QMC quantum impurity solvers is necessary.

!!! note

    Don't forget to set the *isscf* parameter to 2, or else the CT-QMC impurity solvers will skip the dynamical mean-field theory engine and perform one-shot calculation only.

---

**Format**

The format of the *solver.anydos.in* file is as follows:

>
> *column 1*: frequency point, ``\epsilon``, double precision
>
> *column 2*: density of states, ``\rho(\epsilon)``, double precision
>

!!! note

    We assume that:
    1. The orbitals are degenerated.
    2. The number of frequency points is 801.
    If you are not satisfied with these assumptions, you have to hack the corresponding *ctqmc\_dmft.f90* by yourself.

---

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
