# Parameter: ncarlo

**Definition**

> How often to measure the physical observables.

**Type**

> Integer

**Default value**

> 10

**Component**

> All

**Behavior**

> Every *ncarlo* Monte Carlo sampling steps, the quantum impurity solvers try to measure the physical observables. The affected physical observables are as follows:
>
> * Atomic state probability, ``P_{\Gamma}``,
> * Imaginary-time Green's function, ``G(\tau)``,
> * Auxiliary imaginary-time correlation function, ``F(\tau)``.
>
> The measuring period for the other physical observables are controlled by the *nmonte* parameter. See also [nmonte](p_nmonte.md) parameter for more details.
>

!!! note

    The histogram for the perturbation expansion series is measured in each Monte Carlo sampling step.

**Comment**

> You can increase *ncarlo* to relieve the auto-correlation of physical observables > between two successive measurements. However, large *ncarlo* will waste the CPU times. What's the best choices? We guess a good *ncarlo* should satisfy the following relation:
>
> ```math
> \text{ncarlo} * \text{P}_{\text{insert}} >= 1.0
> ```
>
> Here, ``\text{P}_{\text{insert}}`` means the accepted ratio for the insert update actions.
