# Parameter: Jp

**Definition**

> Strength of the pair-hopping term in the interaction term, ``J_p``.

**Type**

> Float, double precision

**Default value**

> 0.0

**Component**

> ALL

**Behavior**

> In principle, it is used to build the Coulomb interaction matrix.
>
> But the **NARCISSUS** component does not support general interaction. In other words, it does not support the ``J_p`` term in the local impurity Hamiltonian.
>
> As for the **MANJUSHAKA** component, the Coulomb interaction matrix is not built within the quantum impurity solvers. All of the information is encapsulated in the *atom.cix* file. So they don't need the ``J_p`` parameter as well.
>
> So, in summary, we don't use this parameter actually. You can set it to any values as you wish. But we strongly suggest to set it to the exact value.

**Comment**

> See also [Jz](p_jz.md) and [Js](p_js.md) for more details.
