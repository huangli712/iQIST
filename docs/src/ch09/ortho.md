### Orthogonal polynomial representation

**Legendre polynomial representation**

Boehnke *et al.*[^1] proposed that the Legendre polynomial can be used to improve the measurements of single-particle and two-particle Green's functions. Thanks to the Legendre polynomial representation, the numerical noise and memory space needed to store the Green's function are greatly reduced.

The imaginary time Green's function ``G(\tau)`` is expressed using the Legendre polynomial ``P_n(x)`` defined in [-1,1]:
```math
\begin{equation}
G(\tau) = \frac{1}{\beta} \sum_{n \leq 0} \sqrt{2n + 1} P_{n}[x(\tau)] G_n,
\end{equation}
```
here ``n`` is the order of Legendre polynomial, ``G_{n}`` is the expansion coefficient, ``x(\tau)``
maps ``\tau \in [0,\beta]`` to ``x \in [-1,1]``:
```math
\begin{equation}
x(\tau) = \frac{2\tau}{\beta} - 1.
\end{equation}
```
Using the orthogonal relations of Legendre polynomials, we obtain 
```math
\begin{equation}
G_n = \sqrt{2n + 1} \int^{\beta}_{0} d\tau P_{n}[x(\tau)] G(\tau).
\end{equation}
```
If we substitute ``G(\tau)`` into ``G_n``, we get 
```math
\begin{equation}
G_n = -\frac{\sqrt{2n + 1}}{\beta} \left\langle \sum^{k}_{i=1} \sum^{k}_{j=1}
\mathcal{M}_{ji} \tilde{P}_{n}(\tau^e_i - \tau^s_j) \right\rangle,
\end{equation}
```
where
```math
\begin{equation}
\tilde{P}_n(\tau) = 
\begin{cases}
P_n [x(\tau)], & \tau > 0, \\
-P_n [x(\tau + \beta)], & \tau < 0, \\
\end{cases}
\end{equation}
```
and ``\tau^{s}`` and ``\tau^{e}`` denote the positions of creation and annihilation operators on the imaginary time axis, respectively. We can also express the Matsubara Green's function ``G(i\omega_n)`` using Legendre polynomials:
```math
\begin{equation}
G(i\omega_m) = \sum_{n \leq 0} T_{mn} G_n.
\end{equation}
```
The transformation matrix ``T_{mn}`` is defined as
```math
\begin{equation}
T_{mn} = (-1)^m i^{n+1} \sqrt{2n + 1} j_n \left[\frac{(2m + 1)\pi}{2}\right],
\end{equation}
```
where ``j_n(z)`` is the spheric Bessel function. Actually, in the Monte Carlo simulation, only the expansion coefficients ``G_n`` are measured. When the calculation is finished, the final Green's function can be evaluated. It is worthwhile to note that the ``T_{mn}`` do not depend on the inverse temperature ``\beta``, so that we can calculate and store the matrix elements beforehand to save computer time.

**Chebyshev polynomial representation**

It is easily to extend this formalism to other orthogonal polynomials. For example, in the iQIST software package, we not only implemented the Legendre polynomial representation, but also the Chebyshev polynomial representation. In the Chebyshev polynomial representation, the imaginary time Green's function ``G(\tau)`` is expanded as follows:
```math
\begin{equation}
G(\tau) = \frac{2}{\beta} \sum_{n \leq 0} U_n [x({\tau})]G_{n},
\end{equation}
```
here the ``U_n(x)`` denote the *second kind* Chebyshev polynomials and ``x \in [-1,1]``. The equation for the expansion coefficients ``G_n`` is:
```math
\begin{equation}
G_n = -\frac{2}{\pi\beta} \left\langle  \sum^{k}_{i=1} \sum^{k}_{j=1} 
\mathcal{M}_{ji} 
\tilde{U}_{n}(\tau^e_i - \tau^s_j)
\sqrt{1 - \tilde{x}(\tau^e_i - \tau^s_j)^2}
\right\rangle,
\end{equation}
```
where
```math
\begin{equation}
\tilde{U}_n (x) = 
\begin{cases}
U_n[x(\tau)], & \tau > 0, \\
-U_n[x(\tau+\beta)], & \tau < 0, \\
\end{cases}
\end{equation}
```
and
```math
\begin{equation}
\tilde{x}(\tau) = 
\begin{cases}
x(\tau), & \tau > 0, \\
x(\tau + \beta), & \tau < 0. \\
\end{cases}
\end{equation}
```
Unfortunately, there is no explicit expression for ``G(i\omega_n)`` in the Chebyshev polynomial representation. So the Legendre polynomial representation is better.

**See also**

The orthogonal polynomial representation has been implemented in the **GARDENIA**, **NARCISSUS**, **LAVENDER**, **MANJUSHAKA**, **CAMELLIA**, and **HIBISCUS**/stoch components.

See [isort](../ch04/p_isort.md), [legrd](../ch04/p_legrd.md), [lemax](../ch04/p_lemax.md), [chgrd](../ch04/p_chgrd.md), [chmax](../ch04/p_chmax.md) parameters for more details.

**Reference**

[^1]: Lewin Boehnke, Hartmut Hafermann, Michel Ferrero, Frank Lechermann, and Olivier Parcollet, *Phys. Rev. B* **84**, 075145 (2011)