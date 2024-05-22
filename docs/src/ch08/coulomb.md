# Coulomb interaction

The standard form of Coulomb interaction in second quantization form is:

```math
\begin{equation}
\hat{H}_{\text{Coulomb}}=
\frac{1}{2}\sum_{\sigma\sigma^{\prime}}\sum_{abcd}
\left\langle a\sigma,b\sigma^{\prime}\left|
\frac{1}{r_{12}}
\right|c\sigma,d\sigma^{\prime}\right\rangle
\hat{f}_{a\sigma}^{\dagger}
\hat{f}_{b\sigma^{\prime}}^{\dagger}
\hat{f}_{d\sigma^{\prime}}
\hat{f}_{c\sigma},
\end{equation}
```

where ``\frac{1}{r_{12}}`` is the Coulomb interaction, ``r_{12}=|\vec{r}_{1}-\vec{r}_{2}|``, ``a,b,c,d`` is orbital index and ``\sigma,\sigma^{\prime}=\uparrow,\downarrow`` is spin index.

---

**Slater Type Interaction**

In this single particle basis, the Coulomb interaction Hamiltonian is:

```math
\begin{equation}
\hat{H}_{U}=\frac{1}{2}\sum_{m_{1},m_{2},m_{1}^{\prime},m_{2}^{\prime},\sigma,\sigma^{\prime}}\left\langle m_{1}\sigma,m_{2}\sigma^{\prime}\left|\frac{1}{r_{12}}\right|m_{1}^{\prime}\sigma,m_{2}^{\prime}\sigma^{\prime}\right\rangle \hat{f}_{m_{1}\sigma}^{\dagger}\hat{f}_{m_{2}\sigma^{\prime}}^{\dagger}\hat{f}_{m_{2}^{\prime}\sigma^{\prime}}\hat{f}_{m_{1}^{\prime}\sigma}.
\end{equation}
```

Set

```math
\alpha=m_{1}\sigma,\ \beta=m_{2}\sigma^{\prime},\ \gamma=m_{1}^{\prime}\sigma,\ \delta=m_{2}^{\prime}\sigma^{\prime},
```

thus the Coulomb *U*-tensor (*UMAT*) in the **JASMINE** component reads:

```math
\begin{equation}
\textbf{UMAT}(\alpha,\beta,\delta,\gamma)=\frac{1}{2}\delta(m_{1}+m_{2},m_{1}^{\prime}+m_{2}^{\prime})\sum_{k}c_{l}^{k}(m_{1},m_{1}^{\prime})c_{l}^{k}(m_{2}^{\prime},m_{2})F_{nl}^{k}.
\end{equation}
```

Here, ``F_{nl}^{k}`` is the Slater integrals:

```math
\begin{equation}
F_{nl}^{k}=\int_{0}^{\infty}r_{1}^{2}dr_{1}\int_{0}^{\infty}r_{2}^{2}dr_{2}R_{nl}^{2}(r_{1})R_{nl}^{2}(r_{2})\frac{r_{<}^{k}}{r_{>}^{k+1}}.
\end{equation}
```

``c_{l}^{k}(m^{\prime},m^{\prime\prime})`` is related to the Gaunt coefficient for fixed ``l``. It is defined as:

```math
\begin{equation}
c_{l}^{k}(m^{\prime},m^{\prime\prime}) = 
\sqrt{\frac{4\pi}{2k+1}} 
\int d\phi d\theta\sin(\theta)
Y_{l}^{m^{\prime}*}(\theta,\phi)
Y_{l}^{m^{\prime}-m^{\prime\prime}}(\theta,\phi)
Y_{l}^{m^{\prime\prime}}(\theta,\phi)
\end{equation}
```

---

**Kanamori Type Interaction**

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
