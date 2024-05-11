# Parameter: nsweep

**Definition**

> Maximum number of quantum Monte Carlo sampling steps.

**Type**

> Integer

**Default value**

> 20000000 for CT-HYB impurity solvers

**Component**

> ALL

**Behavior**

> This is the number of Monte Carlo sampling steps conducted by the current process. If you run the quantum impurity solvers parallelly with ``N_{\text{procs}}`` processes, then the total number of Monte Carlo sampling steps in the calculation is
>
> ```math
> N_{\text{procs}} \times \text{nsweep}
> ```
>
> Larger the *nsweep* is, more accurate and more time-consuming the calculation is. If you conduct the quantum impurity solver on many cores, you can decrease *nsweep*. But the minimal value of it should not be less than its default value.

**Comment**

> In order to improve the numerical quality and suppress the numerical noise, the most direct and practicable route is to increase *nsweep*.
>
> If the data binning mode is activated, *nsweep* should not be larger than 200000000.
>
> *nsweep* should be always larger than *ntherm*, *nwrite*, *nmonte*, and *ncarlo*.

> See also [ntherm](p_ntherm.md), [nwrite](p_nwrite.md), [nmonte](p_nmonte.md), [ncarlo](p_ncarlo.md) parameters for more details.
