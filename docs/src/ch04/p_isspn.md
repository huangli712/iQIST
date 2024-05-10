# Parameter: isspn

**Definition**

> Key control flag, determine the symmetry of spin orientation freedom of degree.

**Type**

> Integer

**Default value**

> 1

**Component**

> ALL

**Behavior**

> There are two possible values for the *isspn* parameter so far:
>
>* *isspn* = 1, let spin up and spin down states evolve independently.
>* *isspn* = 2, enforce spin up = spin down.
>
> The quantum impurity solvers will symmetrize the relevant physical quantities according to the *isspn* parameter.

**Comment**

> Usually, the *isspn* and *isbnd* are used at the same time to setup the magnetic order of the system. For example, if you want to simulate a model with anti-ferromagnetic order, then *isspn* should be 1, and *isbnd* can be 1 or 2 (a *solver.eimp.in* file must be available for the latter case).
>
> See [isbnd](p_isbnd.md) parameter for more details.
