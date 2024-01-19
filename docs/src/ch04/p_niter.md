# Parameter: niter

**Definition**

> Maximum number of continuous time quantum Monte Carlo quantum impurity solver plus dynamical mean field theory self-consistent iterations.

**Type**

> Integer

**Default value**

> 20

**Component**

> ALL

**Behavior**

> To terminate the quantum impurity solvers when the *niter* DMFT iteration is reached.
>
> If *isscf* = 1, the *niter* becomes meaningless. The self-consistent iteration won't be carried out which is equivalent to *niter* = 1.

**Comment**

> See also [isscf](p_isscf.md) and [alpha](p_alpha.md) parameters for more details.

!!! note

    If the convergence is obtained and the minimal iteration number is reached (the default value is 16, you can modify it in the *ctqmc\_dmft.f90*), the quantum impurity solvers will exit as well.
