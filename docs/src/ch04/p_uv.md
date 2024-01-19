# Parameter: Uv

**Definition**

> Inter-orbital Coulomb interaction ``U_v``.

**Type**

> Float, double precision

**Default value**

> 4.0

**Component**

> ALL

**Behavior**

> It is used to determine the Coulomb interaction matrix.

**Comment**

> Actually, only the **NARCISSUS** component need it. For the other quantum impurity solvers components (such as the **MANJUSHAKA** component), the information about the Coulomb interaction matrix is imported via the *atom.cix* file. So you can set it to any values for the **MANJUSHAKA** component.
>
> Usually, the ``U_c``, ``U_v`` and ``J_z`` should satisfy the following relation:
>
> ```math
> \begin{equation}
> U_c = U_v - 2J_z
> \end{equation}
> ```
>
> See also [Uc](p_uc.md), [Jz](p_jz.md) parameters for more details.

!!! note

    For the **NARCISSUS** component, the Coulomb interaction matrix can be imported via the *solver.umat.in* file which has the highest priority. See [solver.umat.in](in_umat.md) for more details.
