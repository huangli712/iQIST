# Parameter: isort

**Definition**

> Key control flag. It is used to determine whether we should use the orthogonal polynomials trick to improve the accuracy and suppress the numerical noises. 

**Type**

> Integer

**Default value**

> 1

**Component**

> Only for the **NARCISSUS** and **MANJUSHAKA** components.

**Behavior**

> There are three possible values for *isort* parameter so far:
>
> * *isort* = 1, the standard method is used to measure ``G(\tau)`` and ``F(\tau)``.
>
> * *isort* = 2, the Legendre orthogonal polynomials trick is used to measure ``G(\tau)`` and ``F(\tau)``.
>
> * *isort* = 3, the Chebyshev orthogonal polynomials trick is used to measure ``G(\tau)`` and ``F(\tau)``.
>
> Here ``G(\tau)`` is the imaginary-time Green's function and ``F(\tau)`` auxiliary imaginary-time function which can be used to calculate ``\Sigma(i\omega_n)`` analytically. 
>
> ```math
> \Sigma(i\omega_n) = \frac{F(i\omega_n)}{G(i\omega_n)}
> ```
>
> If the orthogonal polynomials trick is used, the ``G(\tau)`` and ``F(\tau)`` are expanded using orthogonal polynomials. And then the expansion coefficients are measured during Monte Carlo sampling procedure. Finally, according the expansion coefficients the ``G(\tau)`` and ``F(\tau)`` can be calculated analytically.
>
> The expansion series are defined by the *svmax*, *svgrd*, *lemax*, and *legrd* parameters.
>
> The improved estimator for self-energy function was not implemented in the **MANJUSHAKA** component. In other words, the **MANJUSHAKA** component doesn't support measurement of ``F(\tau)``.

**Comment**

> See [svmax], [svgrd], [lemax](p_lemax.md), and [legrd](p_legrd.md) parameters for more details.
