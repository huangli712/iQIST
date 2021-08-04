### Parameter: itrun

**Definition**

It is a control flag, used to control which scheme should be used to truncate the Hilbert space dynamically. 

**Type**

Integer

**Default value**

1

**Component**

Only for the **MANJUSHAKA** component.

**Behavior**

There are two possible values for *itrun* parameter so far:

* *itrun* = 1, don't truncate it.

* *itrun* = 2, those atomic states with low probability ($$P_{\Gamma} < 1.0E − 6$$) will be truncated for next DMFT iteration.

!!! note 

    The truncation according to the occupancy has been done in the **JASMINE** component. See [Atomic eigenvalue problem solver](../ch06/README.md) for more details.

If the truncation is applied, the quantum impurity solvers will discard the un-selected atomic states which can improve the computational efficiency. On the other hand, the accuracy will be deteriorated. So it is a trade-off.

**Comment**

For 7-band system, sometimes it is essential to consider severe truncation. Or else the computational time is too long to be bearable.

This feature is experimental. **BE CAREFULLY!** Please check your calculated results carefully if you used the truncation approximation. To perform the truncation, the **MANJUSHAKA** will read the *solver.prob.dat* file at first to get the probability data. So if this file is not available, the truncation will not be done. 

To perform truncation over occupation number, you should use the **JASMINE** component to generate suitable *atom.cix* file.