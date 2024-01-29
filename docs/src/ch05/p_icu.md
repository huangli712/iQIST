# Parameter: icu

**Definition**

> Key control flag, determine the scheme to build the Coulomb interaction matrix.

**Type**

> Integer

**Default value**

> 1

**Component**

> Only for the **JASMINE** component.

**Behavior**

> There are three possible values for the *icu* parameter.
>
> * *icu* = 1, Kanamori scheme, using ``U_c``, ``U_v``, ``J_z``, ``J_s``, ``J_p`` parameters, isotropic Hund's rule coupling.
>
> * *icu* = 2, Slater-Cordon scheme, using ``U_d``, ``J_h`` parameters to build ``F_0``, ``F_2``, ``F_4``, ``F_6``.
>
> * *icu* = 3, Kanamori scheme, using ``U_c``, ``U_v``, ``J_z``, ``J_s``, ``J_p`` parameters, anisotropic Hund's rule coupling.
>
> When *icu* = 2, it is not compatible with *ictqmc* = 4.

**Comment**

> See [Uc](p_uc.md), [Uv](p_uv.md), [Jz](p_jz.md), [Js](p_js.md), and [Jp](p_jp.md) parameters for more details.
