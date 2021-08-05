### Parameter: Uc

**Definition**

Intra-orbital Coulomb interaction ``U_c``.

**Type**

Float, double precision

**Default value**

4.0

**Component**

ALL

**Behavior**

It is used to determine the Coulomb interaction matrix.

**Comment**

Actually, only the **AZALEA**, **GARDENIA**, **NARCISSUS**, and **DAISY** components need it. For the other quantum impurity solvers components (i.e., the **CAMELLIA**, **BEGONIA**, **LAVENDER**, **MANJUSHAKA** and **PANSY** components), the information about the Coulomb interaction matrix is imported via the *atom.cix* file. So you can set it to any values for the latter five components.

Usually, the ``U_c``, ``U_v`` and ``J_z`` should meet the following relation:

```math
\begin{equation}
U_c = U_v - 2J_z
\end{equation}
```

See also [Uv](p_uv.md), [Jz](p_jz.md) parameters for more details.

!!! note

    For the **AZALEA**, **GARDENIA**, **NARCISSUS**, and **DAISY** components, the Coulomb interaction matrix can be imported via the *solver.umat.in* file which has the highest priority. See [solver.umat.in](in_umat.md) for more details.