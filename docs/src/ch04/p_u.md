# Parameter: U

**Definition**

> Averaged Coulomb interaction ``U``.

**Type**

> Float, double precision

**Default value**

> 4.0

**Component**

> ALL

**Behavior**

> Actually, it is not used by the quantum impurity solvers to build the Coulomb interaction matrix. It will be output by the impurity solver as a reference. You can set it to any values. I forgot why I had designed this parameter.

**Comment**

> See also [Uc](p_uc.md), [Uv](p_uv.md), [Jz](p_jz.md), [Js](p_js.md) and [Jp](p_jp.md) parameters for more details.

!!! warning

    In the future, perhaps we will remove this parameter.
