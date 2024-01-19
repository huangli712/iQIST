# Parameter: mkink

**Definition**

> Maximum allowable diagrammatic perturbation order in the continuous-time quantum Monte Carlo algorithm.

**Type**

> Integer

**Default value**

> 1024

**Component**

> ALL

**Behavior**

> Determine the size of involved arrays. If the perturbation expansion order exceeds *mkink*, fatal error/exceptions will happen.

**Comment**

> *mkink = 1024* is a safe setting on the condition that the temperature (controlled by the *beta* parameter) is not too low, the interaction (controlled by the *Uc*, *Uv*, *Jz*, etc,  parameters) is not too weak, the number of orbitals (controlled by the *nband* parameter) is not too large.
