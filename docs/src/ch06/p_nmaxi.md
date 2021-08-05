### Parameter: nmaxi

**Definition**

The maximal total occupancy ``N`` which will be kept in the construction of atomic eigenstates.

**Type**

Integer

**Default value**

2

**Component**

Only for the **JASMINE** component.

**Behavior**

Those atomic states in which the total occupancy ``N \in [\text{nmini}, \text{nmaxi}] `` will be kept. The other atomic states will be discarded. It is an aggressive truncation and may led to significant derivations.

**Comment**

See also [nmini](p_nmini.md) parameter for more details.