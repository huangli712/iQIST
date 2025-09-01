# Parameter: legrd

**Definition**

> Number of linear grid points in [-1,1] which is used to define the Legendre polynomials.

**Type**

> Integer

**Default value**

> 20001

**Component**

> Only for the **NARCISSUS** and **MANJUSHAKA** components.

**Behavior**

> The Legendre orthogonal polynomials are defined in a linear grid in [-1,1]. And the number of grid points are controlled by this parameter.

**Comment**

> Only when *isort* = 2 or *isort* = 5 this parameter is useful. See [isort](p_isort.md) for more details. The 20001 is an optimal value for *legrd*. It is not suggested to modify it. 
>
> See also [lemax](p_lemax.md) for more details.
>
> As for the applications of orthogonal polynomials in CT-QMC impurity solver, please refer to Lewin's[^1] and Hartmann's[^2] papers.

**Reference**

[^1]: Lewin Boehnke, Hartmut Hafermann, Michel Ferrero, Frank Lechermann, and Olivier Parcollet, *Phys. Rev. B* **84**, 075145 (2011).

[^2]: Hartmut Hafermann, Kelly R. Patton, and Philipp Werner, *Phys. Rev. B* **85**, 205106 (2012).
