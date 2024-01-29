# Parameter: ictqmc

**Definition**

> Key control flag, define which algorithm is used to diagonalize the atomic Hamiltonian matrix.

**Type**

> Integer

**Default value**

> 1

**Component**

> Only for the **JASMINE** component.

**Behavior**

> There are five possible values for the *ictqmc* parameter:
>
> * *ictqmc* = 1, direct diagonalization in full Hilbert space. The generated *atom.cix* file is only suitable for the **LAVENDER** code. actually, the **LAVENDER** code has been deprecated, so we keep this option only for reference.
>
> * *ictqmc* = 2, subspace diagonalization using good quantum number ``N``. The generated *atom.cix* file is only suitable for the **MANJUSHAKA** code.
>
> * *ictqmc* = 3, subspace diagonalization using good quantum numbers ``N``, ``S_z``. The generated *atom.cix* file is only suitable for the **MANJUSHAKA** code.
>
> * *ictqmc* = 4, subspace diagonalization using good quantum numbers ``N``, ``S_z``, PS. The generated *atom.cix* file is only suitable for the **MANJUSHAKA** code.
>
> * *ictqmc* = 5, subspace diagonalization using good quantum numbers ``N``, ``J_z``. The generated *atom.cix* file is only suitable for the **MANJUSHAKA** code.
>
> When *ictqmc* = 4, it is not compatible with *icu* = 2.

**Comment**

> See also [icu](p_icu.md) for more details.
