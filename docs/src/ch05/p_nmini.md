# Parameter: nmini

**Definition**

> The minimal total occupancy ``N`` which will be kept in the construction of atomic eigenstates.

**Type**

> Integer

**Default value**

> 0

**Component**

> Only for the **JASMINE** component.

**Behavior**

> Those atomic eigenstates in which the total occupancy ``N \in [\text{nmini}, \text{nmaxi}] `` will be kept. The other atomic eigenstates will be discarded. It is an aggressive truncation and may led to significant derivations.

**Comment**

> See also [nmaxi](p_nmaxi.md) parameter for more details.
