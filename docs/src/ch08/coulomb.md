# Coulomb interaction

The standard form of Coulomb interaction in second quantization form is:
```math
\begin{equation}
\hat{H}_{U}=\frac{1}{2}\sum_{\sigma,\sigma^{\prime}}\sum_{a,b,c,d}
\left\langle a\sigma,b\sigma^{\prime}\left|\frac{1}{r_{12}}\right|c\sigma,d\sigma^{\prime}\right\rangle
\hat{f}_{a\sigma}^{\dagger}\hat{f}_{b\sigma^{\prime}}^{\dagger}\hat{f}_{d\sigma^{\prime}}\hat{f}_{c\sigma}
\end{equation}
```

```math
\begin{equation}
\left\langle a\sigma,b\sigma^{\prime}\left|\frac{1}{r_{12}}\right|c\sigma,d\sigma^{\prime}\right\rangle=
\int d\vec{r}_{1}d\vec{r}_{2}\phi_{\alpha\sigma}^{*}(\vec{r}_{1})\phi_{b\sigma^{\prime}}^{*}(\vec{r}_{2})\frac{1}{r_{12}}\phi_{c\sigma}
(\vec{r}_{1})\phi_{d\sigma^{\prime}}(\vec{r}_{2})
\end{equation}
```

where, ``\frac{1}{r_{12}}`` is the Coulomb interaction, ``r_{12}=|\vec{r}_{1}-\vec{r}_{2}|``, ``a,b,c,d`` is orbital index and ``\sigma,\sigma^{\prime}=\uparrow,\downarrow`` is spin index. In the **JASMINE** component, we use a array *UMAT* to store the ``U`` tensor. **BE VERY CAREFUL WITH THE ORBITAL ORDER.** The indices order of *UMAT* is the same with that of the fermion operators, it is **NOT** the same with that of the ``U`` tensor.

**Slater Type**

Expand ``\frac{1}{r_{12}}`` in terms of complex spherical harmonics ``Y_{l}^{m}``,

```math
\begin{equation}
\frac{1}{r_{12}}=\sum_{k}\frac{r_{<}^{k}}{r_{>}^{k+1}}\sum_{m}(-1)^{m}C_{m}^{k}(\theta_{1}\phi_{1})C_{-m}^{k}(\theta_{2}\phi_{2})
\end{equation}
```
where,

```math
C_{m}^{k}(\theta\phi)=\sqrt{\frac{4\pi}{2k+1}}Y_{l}^{m}(\theta\phi).
```

Set ``\phi(\vec{r})=R_{nl}(r)Y_{l}^{m}(\theta\phi)``, then for fixed ``n,l``, we obtain,

```math
\begin{align}
\left\langle m_{1},m_{2}\left|\frac{1}{r_{12}}\right|m_{1}^{\prime},m_{2}^{\prime}\right\rangle & = & \sum_{km}(-1)^{m}c_{l}^{k}(m_{1},m_{1}^{\prime})c_{l}^{k}(m_{2},m_{2}^{\prime})\delta(m+m_{1}^{\prime},m_{1})\delta(-m+m_{2}^{\prime},m_{2})F_{nl}^{k} \\
 & = & \delta(m_{1}+m_{2},m_{1}^{\prime}+m_{2}^{\prime})(-1)^{m_{2}^{\prime}-m_{2}}\sum_{k}c_{l}^{k}(m_{1},m_{1}^{\prime})c_{l}^{k}(m_{2},m_{2}^{\prime})F_{nl}^{k} \\
 & = & \delta(m_{1}+m_{2},m_{1}^{\prime}+m_{2}^{\prime})\sum_{k}c_{l}^{k}(m_{1},m_{1}^{\prime})c_{l}^{k}(m_{2}^{\prime},m_{2})F_{nl}^{k}
\end{align}
```

where,

```math
\begin{equation}
c_{l}^{k}(m^{\prime},m^{\prime\prime})=\int d\phi d\theta \sin(\theta)Y_{l}^{m^{\prime}*}(\theta\phi)C_{m^{\prime}-m^{\prime\prime}}^{k}(\theta\phi)Y_{l}^{m^{\prime\prime}}(\theta\phi)
\end{equation}
```

is the Gaunt coefficient for fixed ``l``, and

```math
\begin{equation}
F_{nl}^{k}=\int_{0}^{\infty}r_{1}^{2}dr_{1}\int_{0}^{\infty}r_{2}^{2}dr_{2}R_{nl}^{2}(r_{1})R_{nl}^{2}(r_{2})\frac{r_{<}^{k}}{r_{>}^{k+1}}
\end{equation}
```
is the Slater integrals.

In this single particle basis, the Coulomb interaction Hamiltonian is:
```math
\begin{equation}
\hat{H}_{U}=\frac{1}{2}\sum_{m_{1},m_{2},m_{1}^{\prime},m_{2}^{\prime},\sigma,\sigma^{\prime}}\left\langle m_{1}\sigma,m_{2}\sigma^{\prime}\left|\frac{1}{r_{12}}\right|m_{1}^{\prime}\sigma,m_{2}^{\prime}\sigma^{\prime}\right\rangle \hat{f}_{m_{1}\sigma}^{\dagger}\hat{f}_{m_{2}\sigma^{\prime}}^{\dagger}\hat{f}_{m_{2}^{\prime}\sigma^{\prime}}\hat{f}_{m_{1}^{\prime}\sigma}
\end{equation}
```
where, ``\sigma,\sigma^{\prime}`` is spin index. Set

```math
\alpha=m_{1}\sigma,\ \beta=m_{2}\sigma^{\prime},\ \gamma=m_{1}^{\prime}\sigma,\ \delta=m_{2}^{\prime}\sigma^{\prime},
```

thus the *UMAT* in the **JASMINE** component reads:

```math
\begin{equation}
\textbf{UMAT}(\alpha,\beta,\delta,\gamma)=\frac{1}{2}\delta(m_{1}+m_{2},m_{1}^{\prime}+m_{2}^{\prime})\sum_{k}c_{l}^{k}(m_{1},m_{1}^{\prime})c_{l}^{k}(m_{2}^{\prime},m_{2})F_{nl}^{k}.
\end{equation}
```

**Kanamori Type**

The Kanomori type of Coulomb interaction Hamiltonian in the **JASMINE** component is defined as:

```math
\begin{align}
\hat{H}_{U} & = & U\sum_{a}\hat{f}_{a,\uparrow}^{\dagger}\hat{f}_{a,\uparrow}\hat{f}_{a,\downarrow}^{\dagger}\hat{f}_{a,\downarrow}\\
 & = & U^{\prime}\sum_{a\text{<}b,\sigma}\hat{f}_{a,\sigma}^{\dagger}\hat{f}_{a,\sigma}\hat{f}_{b,-\sigma}^{\dagger}\hat{f}_{b,-\sigma}
  + (U^{\prime}-J_{z})\sum_{a<b,\sigma}\hat{f}_{a,\sigma}^{\dagger}\hat{f}_{a,\sigma}\hat{f}_{b,\sigma}^{\dagger}\hat{f}_{b,\sigma}\\
 & - & J_{s}\sum_{a<b,\sigma}\hat{f}_{a,\sigma}^{\dagger}\hat{f}_{a,-\sigma}\hat{f}_{b,-\sigma}^{\dagger}\hat{f}_{b,\sigma}
  +  J_{p}\sum_{a\neq b}\hat{f}_{a,\uparrow}^{\dagger}\hat{f}_{a,\downarrow}^{\dagger}\hat{f}_{b,\downarrow}\hat{f}_{b,\uparrow}
\end{align}
```
where, ``a,b`` is orbital index, and ``\sigma=\uparrow,\downarrow`` is spin index.
