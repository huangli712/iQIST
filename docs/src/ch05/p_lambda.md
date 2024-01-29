# Parameter: lambda

**Definition**

> Strength of the spin-orbit coupling ``\lambda``.

**Type**

> Float, double precision

**Default value**

> 0.0

**Component**

> Only for the **JASMINE** component.

**Behavior**

> This parameter is used to build the spin-orbit coupling term of the local interaction Hamiltonian. It is valid only when *isoc* = 1.
>
> ```math
> \hat{H}_{\text{SOC}} = \lambda \hat{\textbf{l}} \cdot \hat{\textbf{s}}.
> ```

**Comment**

> See also [isoc](p_isoc.md) parameter for more details.
