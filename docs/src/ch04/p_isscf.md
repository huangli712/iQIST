# Parameter: isscf

**Definition**

> It is the most important control flag, which is used to control the running mode of the quantum impurity solvers. As mentioned before, the quantum impurity solvers in the iQIST software package contain a mini dynamical mean-field theory self-consistent engine, concerning the Hubbard model on a Bethe lattice. You can use the *isscf* parameter to control whether to active this self-consistent engine or not.

**Type**

> Integer

**Default value**

> 1

**Component**

> ALL

**Behavior**

> There are two possible values for the *isscf* parameter so far:
>
> * *isscf* = 1, one-shot non-self-consistent scheme, usually used in local density approximation plus dynamical mean-field theory scenario. The internal dynamical mean-field theory self-consistent is disabled.
>
> * *isscf* = 2, self-consistent scheme, used in standard Hubbard model (on Bethe lattice) plus dynamical mean-field theory scenario. The internal dynamical mean-field theory self-consistent engine is used.

**Comment**

> The internal dynamical mean-field theory self-consistent engine is implemented in the *ctqmc\_dmft.f90*. The experienced users can customize their own versions of it.
