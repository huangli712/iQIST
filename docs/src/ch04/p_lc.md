# Parameter: lc

**Definition**

> It is a model parameter for the Hubbard-Holstein model or dynamical screening effect. The exact definition of *lc* depends on the value of *isscr* parameter.

**Type**

> Float, double precision

**Default value**

> 1.0

**Component**

> Only for the **NARCISSUS** component.

**Behavior**

> There are five possible choices for the *lc* parameter so far:
>
> * *isscr* = 1, no meaning.
>
> * *isscr* = 2, the ``\lambda`` parameter in the plasmon-pole model for dynamical screening effect.
>
> * *isscr* = 3, the ``\alpha`` parameter in the ohmic model for dynamical screening effect.
>
> * *isscr* = 4, no meaning.

**Comment**

> Only when *isscr* > 1, it matters. Please see [isscr](p_isscr.md) parameter for more details. About CT-QMC algorithms for the Hubbard-Holstein model and dynamical screening effect (including plasmon-pole model and ohmic model), please refer to Philipp Werner's papers[^1][^2].

**Reference**

[^1]: Philipp Werner and Andrew J. Millis, Phys. Rev. Lett 99, 146404 (2007).

[^2]: Philipp Werner and Andrew J. Millis, Phys. Rev. Lett. 104, 146401 (2010).
