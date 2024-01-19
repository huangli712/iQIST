# Parameter: nffrq

**Definition**

> Number of fermionic Matsubara frequency points for the two-particle green's function.
>
> The two-particle Green's function ``\chi(i\omega_n, i\omega'_n, i\nu_n)`` and vertex function ``\mathcal{F}(i\omega_n, i\omega'_n, i\nu_n)`` have three frequency indices where ``\omega_n`` and ``\omega'_n`` are fermionic frequencies, and ``\nu_n`` bosonic frequency:
>
> ```math
> \nu_n = \frac{2n\pi}{\beta}
> ```
>
> ```math
> \omega_n = \frac{(2n+1)\pi}{\beta}
> ```
>
> The *nbfrq* parameter is used to define and generate the bosonic mesh ``\nu_n``. The corresponding fermionic mesh ``\omega_n`` is defined by the *nffrq* parameter.

**Type**

> Integer

**Default value**

> 32

**Component**

> ALL

**Behavior**

> Determine the size of involved arrays. Only useful when we need to compute the two-particle quantities.

**Comment**

> If the *nffrq* is too large, the computational time is not bearable. So we suggest to set *nffrq* ``< 128``. The computation of two-particle quantities is extremely time-consuming, though we have tried our best optimizing it.
>
> See [nbfrq](p_nbfrq.md) and [nfreq](p_nfreq.md) for more details.

!!! tip

    The computation of two-particle quantities has been optimized using the OpenMP multi-thread technology. You can enable this feature in the compiling procedure. See also [Compiling environment](../ch02/envir.md) for more details.
