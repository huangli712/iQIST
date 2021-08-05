### Parameter: nvect

**Definition**

Number of selected eigen-vectors (maximum value is *ncfgs*, minimum value is 1) for the calculation of local operator trace.

**Type**

Integer

**Default value**

4

**Component**

Only for the **CAMELLIA** component.

**Behavior**

In the **CAMELLIA** component, we used the Newton-Leja interpolation algorithm to calculate the local operator trace. In this algorithm, we have to consider the time evolution of each eigenvectors of the local Hamiltonian. In order to improve the computational efficiency, we applied the truncation approximation further. In other words, we only keep a few eigenvectors in the calculations. And *nvect* is the number of the selected eigenstates. 

If we only retain the ground states, it is the ``O(1)`` approximation. If both the ground states and the first low-lying states are retained, it is the ``O(2)`` approximation. According to our experience, for low temperature system, the ``O(1)`` approximation is reasonable. Obviously, the larger *nvect* is, the more accurate and slower the calculation will be.

**Comment**

The **CAMELLIA** component is not ready now.