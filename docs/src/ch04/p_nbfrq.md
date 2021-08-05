### Parameter: nbfrq

**Definition**

Number of bosonic Matsubara frequency points for the two-particle green's function.

The two-particle Green's function ``\chi(i\omega_n, i\omega'_n, i\nu_n)`` and vertex function ``\mathcal{F}(i\omega_n, i\omega'_n, i\nu_n)`` have three frequency indices where ``\omega_n`` and ``\omega'_n`` are fermionic frequencies, and ``\nu_n`` bosonic frequency: 

```math
\nu_n = \frac{2n\pi}{\beta}
```

```math
\omega_n = \frac{(2n+1)\pi}{\beta}
```

The *nbfrq* parameter is used to define and generate ``\nu_n`` bosonic mesh. The corresponding ``\omega_n`` fermionic mesh is defined by the *nffrq* parameter.

**Type**

Integer

**Default value**

8

**Component**

Only for the **GARDENIA**, **NARCISSUS**, **CAMELLIA**, **LAVENDER**, and **MANJUSHAKA** components.

**Behavior**

Determine the size of involved arrays. Only useful when we need to compute the two-particle quantities.

**Comment**

If the *nbfrq* is too large, the computational time is not bearable. So we recommend to set *nbfrq* `` \leq 32`` The computation of two-particle quantities is extremely time-consuming, though we have tried our best optimizing it.

See [nffrq](p_nffrq.md) and [nfreq](p_nfreq.md) for more details.

!!! tip

    The computation of two-particle quantities has been optimized using the OpenMP multi-thread technology. You can enable this feature in the compiling procedure. See also [Compiling environment](../ch02/envir.md) for more details.