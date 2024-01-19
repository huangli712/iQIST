# Parameter: ncfgs

**Definition**

> Number of atomic states or many-body configurations, the dimension of Hilbert space.

**Type**

> Integer

**Default value**

> 4

**Component**

> Only for the **JASMINE** component.

**Behavior**

> Determine the size of involved arrays.

**Comment**

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
> You have to ensure the value of *ncfgs* is compatible with *nband*, *nspin* and *norbs*. The **JASMINE** component will not check and correct them automatically. See also [nspin](p_nspin.md), [nband](p_nband.md), and [norbs](p_norbs.md) parameters for more details.
