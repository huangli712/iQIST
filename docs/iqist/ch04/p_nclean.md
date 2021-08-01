### Parameter: nclean

**Definition**

Clean update period for quantum impurity solvers.

**Type**

Integer

**Default value**

100000 (for CT-HYB impurity solvers) or 100 (for HF-QMC impurity solver)

**Component**

ALL

**Behavior**

After a few Monte Carlo sampling steps, the numerical accuracy may be deteriorated. In order to retain the numerical accuracy, the quantum impurity solvers will conduct a *clean update* (in other words, recalculate everything from scratch) every *nclean* Monte Carlo sampling step.

**Comment**

*nclean = 100000* (for CT-HYB) or *nclean = 100* (for HF-QMC) is an optimal setting. Do not modify it casually.