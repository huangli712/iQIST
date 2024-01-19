# Parameter: Jz

**Definition**

> Coupling constant for the Hund's exchange interaction in z axis, ``J_z``.

**Type**:

> Float, double precision

**Default value**

> 0.0

**Component**

> ALL

**Behavior**

> It is used to determine the Coulomb interaction matrix.

**Comment**

> Actually, only the **NARCISSUS** component needs it. For the other quantum impurity solvers components (such as the **MANJUSHAKA** component), the information about the Coulomb interaction matrix is imported via the *atom.cix* file. So you can set it to any values for the **MANJUSHAKA** component.
>
> Usually, the ``U_c``, ``U_v`` and ``J_z`` should satisfy the following relation:
>
> ```math
> \begin{equation}
> U_c = U_v - 2J_z
> \end{equation}
> ```
>
> See also [Uc](p_uc.md), [Uv](p_uv.md) parameters for more details.
>
> Usually, ``J_z = J_s = J_p = J``. As for a single-band model, they are zero. See also [Js](p_js.md), [Jp](p_jp.md) parameters for more details.

!!! note

    For the **NARCISSUS** component, the Coulomb interaction matrix can be imported via the *solver.umat.in* file which has the highest priority. See [solver.umat.in](in_umat.md) for more details.
