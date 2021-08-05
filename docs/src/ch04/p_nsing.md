### Parameter: nsing

**Definition**

The Hirsch-Fye algorithm is an auxiliary field quantum Monte Carlo algorithm. We implemented in the **DAISY** component in the iQIST software package. The *nsing* parameter is used to determine the number of auxiliary ising-like fields.

**Type**

Integer

**Default value**

1

**Component**

Only for the **DAISY** component.

**Behavior**

Determine the size of involved arrays.

The following relation always holds:

```math
\text{nsing} = \frac{\text{norbs} * (\text{norbs} - 1)} { 2 }
```

So for single-band model, we obtain *nband* = 1, *norbs* = 2, *nsing* = 1.

**Comment**

N/A