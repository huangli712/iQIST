# Parameter: nband

**Definition**

> Number of correlated bands in the model Hamiltonian.

**Type**

> Integer

**Default value**

> 1

**Component**

> ALL

**Behavior**

> Determine the size of involved arrays.

**Comment**

> In the iQIST software package, when we say *nband*, we always do not consider the spin degree of freedom. So for ``d``−electron system, it should be a five-band model (*nband* = 5), while for ``f``−electron, it should be a seven-band model (*nband* = 7).
>
> When you try to setup the *nband* parameter in the *solver.ctqmc.in* file, the corresponding *norbs* and *ncfgs* parameters must be set as well.
>
> See [nspin](p_nspin.md), [norbs](p_norbs.md), and [ncfgs](p_ncfgs.md) parameters for more details.
