### Parameter: lambda

**Definition**

Strength of the spin-orbital coupling $$\lambda$$.

**Type**

Float, double precision

**Default value**

0.0

**Component**

Only for the **JASMINE** component.

**Behavior**

This parameter is used to build the spin-orbital coupling term of the local interaction Hamiltonian. It is valid only when *isoc* = 1.

$$
H_{\text{SOC}} = \lambda \textbf{L} \cdot \textbf{S}
$$

**Comment**

See also [isoc](p_isoc.md) parameter for more details.