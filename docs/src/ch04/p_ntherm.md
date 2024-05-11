# Parameter: ntherm

**Definition**

> Maximum number of Monte Carlo thermalization steps.

**Type**

> Integer

**Default value**

> 200000 for CT-HYB impurity solvers

**Component**

> ALL

**Behavior**

> Before starting to measure the physical observables, we have to ensure the system has reached thermal equilibrium state. So during the first *ntherm* Monte Carlo steps, the quantum impurity solvers just update the diagrammatic configuration, but don't record the physical quantities.

**Comment**

> According to our experience, the optimal value for *ntherm* can be one tenth of the value of *nsweep*. See [nsweep](p_nsweep.md) parameter for more details.
