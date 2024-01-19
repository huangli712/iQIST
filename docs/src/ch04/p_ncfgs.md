# Parameter: ncfgs

**Definition**

> Number of atomic states considered in the calculations.

**Type**

> Integer

**Default value**

> 4

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
> You have to ensure the value of *norbs* is compatible with *nband*, *nspin* and *ncfgs*. The quantum impurity solvers will not check and correct them automatically. So you have to setup them in the *solver.ctqmc.in* file explicitly.
>
> See also [nspin](p_nspin.md), [nband](p_nband.md), and [norbs](p_norbs.md) parameters for more details.
