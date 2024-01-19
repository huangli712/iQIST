# Parameter: mfreq

**Definition**

> Maximum number of Matsubara frequency points ``n_{\text{max}}``. The Matsubara frequency mesh for fermions is defined as follows:
>
> ```math
> \omega_n = \frac{(2n + 1)\pi}{\beta}
> ```
>
> where ``n`` = 0, 1, 2, 3, ``\cdots``, ``n_{\text{max}}``.

**Type**

> Integer

**Default value**

> 8193 ``(\equiv 2^{13}+1)``

**Component**

> ALL

**Behavior**

> Determine the size of involved arrays, such as ``G(i\omega_n)``, ``G_0(i\omega_n)``, ``\Sigma(i\omega_n)``, and ``\Delta(i\omega_n)`` etc.

**Comment**

> *mfreq = 8193* is a safe setting. We don't recommend to change it.
