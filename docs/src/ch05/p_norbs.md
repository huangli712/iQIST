# Parameter: norbs

**Definition**

> Number of correlated orbitals.

**Type**

> Integer

**Default value**

> 2

**Component**

> Only for the **JASMINE** component.

**Behavior**

> Determine the size of involved arrays.

**Comment**:

> The following relations always hold:
>
> ```math
> \text{norbs} = \text{npsin} * \text{nband}
> ```
>
> ```math
> \text{ncfgs} = 2^{\text{norbs}}
> ```
>
> You have to ensure the value of *norbs* is compatible with *nband*, *nspin* and *ncfgs*. The **JASMINE** component will not check and correct them automatically. See also [nspin](p_nspin.md), [nband](p_nband.md), and [ncfgs](p_ncfgs.md) parameters for more details.
