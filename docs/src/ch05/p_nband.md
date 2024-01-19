# Parameter: nband

**Definition**

> Number of correlated bands.

**Type**

> Integer

**Default value**

> 1

**Component**

> Only for the **JASMINE** component.

**Behavior**

> Determine the size of involved arrays.

**Comment**

> In the iQIST software package, when we say *nband*, we always do not consider the spin degree of freedom. So for ``d``−electron system, it should be a five-band model (*nband* = 5), while for ``f``−electron, it should be a seven-band model (*nband* = 7).

> See [nspin](p_nspin.md), [norbs](p_norbs.md), and [ncfgs](p_ncfgs.md) parameters for more details.
