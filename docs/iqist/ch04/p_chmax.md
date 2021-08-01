### Parameter: chmax

**Definition**

The maximum allowable expansion order $$n$$ for the second kind Chebyshev orthogonal polynomials. It must be greater than 2.

The recursive definition for the second kind Chebyshev polynomials is as follows:

$$
U_{0} = 1
$$

$$
U_{1} = 2x
$$

$$
U_{n+1}(x) = 2xU_{n}(x) - U_{n-1}(x)
$$

The second kind Chebyshev polynomials are orthogonal with respect to the weight 

$$
\sqrt{1-x^2}
$$

on the interval [-1,1], i.e., we have:

$$
\begin{equation}
\int_{-1}^1 U_n(x)U_m(x)\sqrt{1-x^2}\,dx = 
\begin{cases} 0 &: n\ne m, \\ \pi/2 &: n=m. \end{cases} 
\end{equation}
$$
**Type**

Integer

**Default value**

32

**Component**

Only for the **GARDENIA**, **NARCISSUS**, **LAVENDER**, **CAMELLIA**, and **MANJUSHAKA** components.

**Behavior**

The parameter is used as a cutoff to limit the maximum expansion order for the second kind Chebyshev orthogonal polynomials.

**Comment**

Only when *isort* = 3 or *isort* = 6 this parameter is useful. See [isort](p_isort.md) parameter for more details. How to choose a suitable *chmax* parameters is a tricky job. If *chmax* is too small, the calculated results won't be accurate. If *chmax* is too large, the so-called Gibbs oscillation will occur significantly. According to our experiences, 32 or 48 may be a reasonable choice. Though there is no upper limit for *isort*, the larger the value of *chmax* is, the heavier the computational burden will be. 

See also [chgrd](p_chgrd.md) for more details.

As for the applications of orthogonal polynomials in CT-QMC impurity solver, please refer to Lewin's[^1] and Hartmann's[^2] papers.

**Reference**

[^1] Lewin Boehnke, Hartmut Hafermann, Michel Ferrero, Frank Lechermann, and Olivier Parcollet, *Phys. Rev. B* **84**, 075145 (2011).

[^2] Hartmut Hafermann, Kelly R. Patton, and Philipp Werner, *Phys. Rev. B* **85**, 205106 (2012).