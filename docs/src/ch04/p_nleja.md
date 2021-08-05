### Parameter: nleja

**Definition**

In the **CAMELLIA** component, we used the Newton-Leja interpolation algorithm to evaluate the local operator trace. Specially, the real leja points are used in the newton interpolation to evaluate 

```math
\exp( -\tau H ) |v\rangle
``` 

efficiently. Here ``H`` is the local Hamiltonian, ``|v\rangle`` the propagated vector, ``\tau`` imaginary-time point. The *nleja* parameter is a key parameter to control the Newton-Leja algorithm. It means the maximum number of real leja points. 

**Type**

Integer

**Default value**

64

**Component**

Only for the **CAMELLIA** component.

**Behavior**

Determine the size of involved arrays.

This Newton-Leja algorithm is unstable. Sometimes it is not easy to obtain convergence with the given leja points. At that time, you can increase *nleja* and redo the calculation.

**Comment**

The **CAMELLIA** component is not ready now.