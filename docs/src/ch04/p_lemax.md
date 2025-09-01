# Parameter: lemax

**Definition**

> The maximum allowable expansion order ``n`` for the Legendre orthogonal polynomials. It must be greater than 2.
>
> The recursive definition for the Legendre polynomials is as follows:
>
> ```math
> P_{0}(x) =1
> ```
>
> ```math
> P_{1}(x) = x
> ```
>
> ```math
> (n+1) P_{n+1}(x) = (2n+1)xP_{n}(x) - n P_{n-1}(x)
> ```
>
>An important property of the Legendre polynomials is that they are orthogonal with respect to the ``L^2`` inner product on the interval ``−1 \leq x\leq 1``:
>
> ```math
> \int^{1}_{-1} P_{m}(x) P_{n}(x) dx = \frac{2}{2n+1}\delta_{mn}
> ```

**Type**

> Integer

**Default value**

> 32

**Component**

> Only for the **NARCISSUS** and **MANJUSHAKA** components.

**Behavior**

> The parameter is used as a cutoff to limit the maximum expansion order for the Legendre orthogonal polynomials.

**Comment**

> Only when *isort* = 2 or *isort* = 5 this parameter is useful. See [isort](p_isort.md) parameter for more details. How to choose a suitable *lemax* parameters is a tricky job. If *lemax* is too small, the calculated results won't be accurate. If *lemax* is too large, the so-called Gibbs oscillation will occur dramatically. According to our experiences, 32 or 48 may be a reasonable choice. It is worthy to emphasis that due to the limitation of implementation, *lemax* must be less than 50.
>
> See also [legrd](p_legrd.md) for more details.
>
> As for the applications of orthogonal polynomials in CT-QMC impurity solver, please refer to Lewin's[^1] and Hartmann's[^2] papers.

**Reference**

[^1]: Lewin Boehnke, Hartmut Hafermann, Michel Ferrero, Frank Lechermann, and Olivier Parcollet, *Phys. Rev. B* **84**, 075145 (2011).

[^2]: Hartmut Hafermann, Kelly R. Patton, and Philipp Werner, *Phys. Rev. B* **85**, 205106 (2012).
