# Parameter: nspin

**Definition**

> Number of spin projection considered in the calculations.

**Type**

> Integer

**Default value**

> 2

**Component**

> ALL

**Behavior**

> Determine the size of involved arrays.

**Comment**

> See also [nband](p_nband.md), [norbs](p_norbs.md), and [ncfgs](p_ncfgs.md) parameters for more details.

!!! warning

    Please DO NOT modify it. Now the quantum impurity solvers in the iQIST software package don't support *nspin* = 1.
