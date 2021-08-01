### Parameter: ictqmc

**Definition**

Key control flag, define which algorithm is used to diagonalize the atomic Hamiltonian matrix.

**Type**

Integer

**Default value**

1

**Component**

Only for the **JASMINE** component.

**Behavior**

There are six possible values for the *ictqmc* parameter:

* *ictqmc* = 0, direct diagonalization in full Hilbert space. The generated *atom.cix* file is only suitable for the **CAMELLIA** code.

* *ictqmc* = 1, direct diagonalization in full Hilbert space. The generated *atom.cix* file is only suitable for the **BEGONIA** and **LAVENDER** codes.

* *ictqmc* = 2, diagonalization in subspace using good quantum number $$N$$. The generated *atom.cix* file is only suitable for the **MANJUSHAKA** and **PANSY** codes.

* *ictqmc* = 3, diagonalization in subspace using good quantum numbers $$N$$, $$S_z$$. The generated *atom.cix* file is only suitable for the **MANJUSHAKA** and **PANSY** codes.

* *ictqmc* = 4, diagonalization in subspace using good quantum numbers $$N$$, $$S_z$$, PS. The generated *atom.cix* file is only suitable for the **MANJUSHAKA** and **PANSY** codes.

* *ictqmc* = 5, diagonalization in subspace using good quantum numbers $$N$$, $$J_z$$. The generated *atom.cix* file is only suitable for the **MANJUSHAKA** and **PANSY** codes.

When *ictqmc* = 4, it is not compatible with *icu* = 2.

**Comment**

See also [icu](p_icu.md) for more details.