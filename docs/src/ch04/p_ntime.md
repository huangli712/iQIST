# Parameter: ntime

**Definition**

> Number of imaginary time slices sampling by continuous time or Hirsch-Fye quantum Monte Carlo quantum impurity solver.
>
> In the iQIST software package, the imaginary time mesh is calculated as follows:
>
> ```math
> \tau_i = \frac{\beta (i-1)}{n-1}
> ```
>
> where ``i \in [1,n] `` and ``\tau_i \in [0,\beta]``. So, the value of ``n`` is just *ntime*.

**Type**

> Integer

**Default value**

> 1024 (for CT-HYB impurity solvers)

**Component**

> ALL

**Behavior**

> Determine the size of involved arrays, such as the imaginary-time Green's function ``G(\tau)``, the hybridization function ``\Delta(\tau)``, the auxiliary correlation function ``F(\tau)``, the spin-spin correlation function ``\langle S_z(0) S_z(\tau) \rangle``, etc.

**Comment**

> The *ntime = 1024* is an optimal value for most cases. But if the system temperature is too low, it is useful to increase *ntime* to obtain higher accurate. For example, if ``\beta = 100``, it is better to let *ntime = 2048 or 4096*.
>
> See [mfreq](p_mfreq.md) and [nfreq](p_nfreq.md) parameters for more details.
