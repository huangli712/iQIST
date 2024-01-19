# Parameter: part

**Definition**

> Hopping (or coupling) parameter ``t`` for Hubbard model.

**Type**

> Float, double precision

**Default value**

> 0.5

**Component**

> ALL

**Behavior**

> It is used to initialize the default hybridization function ``\Delta(i\omega_n)``, and in the implementation of the dynamical mean-field theory self-consistent equation for the Bethe lattice.
>
> Para-magnetic phase:
>
> ```math
> \begin{equation}
> G_{\alpha\sigma}(i\omega_n) = t^2 \Delta_{\alpha\sigma}(i\omega_n)
> \end{equation}
> ```
>
> Anti-ferromagnetic phase:
>
> ```math
> \begin{equation}
> G_{\alpha\bar{\sigma}}(i\omega_n) = t^2 \Delta_{\alpha\bar{\sigma}}(i\omega_n)
> \end{equation}
> ```

**Comment**

> It is useful only when *isscf* = 2. See also [isscf](p_isscf.md) parameter for more details.
