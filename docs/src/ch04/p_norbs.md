# Parameter: norbs

**Definition**

> Number of correlated orbitals considered in the calculations.

**Type**

> Integer

**Default value**

> 2

**Component**

> ALL

**Behavior**

> Determine the size of involved arrays.

**Comment**

> The following relations always hold.
>
> ```math
> \text{norbs} = \text{npsin} * \text{nband}
> ```
>
> ```math
> \text{ncfgs} = 2^{\text{norbs}}
> ```
>
> You have to ensure the value of *norbs* is compatible with *nband*, *nspin* and *ncfgs*. The quantum impurity solvers will not check and correct them automatically. So you have to setup them in the *solver.ctqmc.in* and *solver.hfqmc.in* files explicitly.
>
> See also [nspin](p_nspin.md), [nband](p_nband.md), and [ncfgs](p_ncfgs.md) parameters for more details.
