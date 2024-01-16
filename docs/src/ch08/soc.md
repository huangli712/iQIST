# Spin-orbit coupling

**Spin-orbit interaction Hamiltonian**

The spin-orbit coupling (SOC) is implemented at atomic level,

```math
\begin{equation}
\hat{H}_{\text{SOC}}=\lambda\sum_{i}\vec{\mathbf{l}}_{i}\cdot\vec{\mathbf{s}}_{i},
\end{equation}
```

where ``\vec{\mathbf{l}}`` is orbital angular momentum, and ``\vec{\mathbf{s}}`` is spin angular momentum. In second quantization form,

```math
\begin{equation}
\hat{H}_{\text{SOC}}=\lambda\sum_{\alpha\sigma,\beta\sigma^{\prime}}\left\langle\alpha\sigma\left|\vec{\mathbf{l}}\cdot\vec{\mathbf{s}}\right|\beta\sigma^{\prime}\right\rangle\hat{f}_{\alpha\sigma}^{\dagger}\hat{f}_{\beta\sigma^{\prime}},
\end{equation}
```

where ``\alpha`` is orbital index and ``\sigma`` is spin index. We note that

```math
\begin{equation}
\vec{\mathbf{l}}\cdot\vec{\mathbf{s}} = 
\frac{1}{2} \vec{\mathbf{l}} \cdot \vec{\mathbf{\sigma}}.
\end{equation}
```

where ``\vec{\mathbf{\sigma}}`` is the Pauli operator:

```math
\begin{equation}
\vec{\mathbf{\sigma}} = \hat{\sigma}_x \hat{x} + \hat{\sigma}_y \hat{y} + {\sigma}_z \hat{z}.
\end{equation}
```

```math
\begin{equation}
\hat{\sigma}_x = \left[
\begin{array}{cc}
0 & 1 \\
1 & 0 \\
\end{array}
\right],~
\hat{\sigma}_y = \left[
\begin{array}{cc}
0 & -i \\
i & 0 \\
\end{array}
\right],~
\hat{\sigma}_z = \left[
\begin{array}{cc}
1 & 0 \\
0 & -1 \\
\end{array}
\right].
\end{equation}
```

Now the question is how to write down the matrix elements for ``\vec{\mathbf{l}}\cdot\vec{\mathbf{\sigma}}``.

```math
\begin{equation}
\vec{\mathbf{l}}\cdot\vec{\mathbf{\sigma}} =
\hat{l}_{x}\hat{\sigma}_{x}+
\hat{l}_{y}\hat{\sigma}_{y}+
\hat{l}_{z}\hat{\sigma}_{z}
= \left[\begin{array}{cc}
\hat{l}_{z} & \hat{l}_{-}\\
\hat{l}_{+} & -\hat{l}_{z}
\end{array}\right],
\end{equation}
```

where ``\hat{l}_{\pm}=\hat{l}_{x}\pm\hat{l}_{y}``. It is easy to prove that Eq.(6) can be rewritten as

```math
\begin{equation}
\vec{\mathbf{l}}\cdot\vec{\mathbf{\sigma}} =
\hat{l}_{+}\hat{\sigma}_{-}+
\hat{l}_{-}\hat{\sigma}_{+}+
\hat{l}_{z}\hat{\sigma}_{z},
\end{equation}
```

where

```math
\begin{equation}
\hat{\sigma}_+ = \left[
\begin{array}{cc}
0 & 1 \\
0 & 0 \\
\end{array}
\right],~
\hat{\sigma}_- = \left[
\begin{array}{cc}
0 & 0 \\
1 & 0 \\
\end{array}
\right].
\end{equation}
```

**Spin-orbit interaction in the complex orbital basis**

We just write down ``\vec{\mathbf{l}}\cdot\vec{\mathbf{\sigma}}`` in the complex spherical harmonics basis. The following relations are used to derive the matrix elements of the spin-orbit interaction:

```math
\begin{equation}
\hat{l}_{\pm}Y_{l}^{m}=\sqrt{(l\mp m)(l\pm m+1)}Y_{l}^{m\pm1}.
\end{equation}
```

```math
\begin{equation}
\hat{\sigma}_+ |\uparrow\rangle = 0,~
\hat{\sigma}_+ |\downarrow\rangle = |\uparrow\rangle.
\end{equation}
```

```math
\begin{equation}
\hat{\sigma}_- |\uparrow\rangle = |\downarrow\rangle,~
\hat{\sigma}_- |\downarrow\rangle = 0.
\end{equation}
```

```math
\begin{equation}
\hat{\sigma}_z |\uparrow\rangle = 1,~
\hat{\sigma}_z |\downarrow\rangle = -1.
\end{equation}
```

For ``p`` system,
```math
\begin{equation}
\vec{\mathbf{l}}\cdot\vec{\mathbf{\sigma}}=\left[\begin{array}{ccc|ccc}
0 & 0 & 0 & 0 & 0 & 0\\
0 & 0 & 0 & 0 & 0 & 0\\
0 & 0 & 0 & 0 & 0 & 0\\
\hline
0 & 0 & 0 & 0 & 0 & 0\\
0 & 0 & 0 & 0 & 0 & 0\\
0 & 0 & 0 & 0 & 0 & 0\\
\end{array}\right]
\end{equation}
```

For ``d`` system,
```math
\begin{equation}
\vec{\mathbf{l}}\cdot\vec{\mathbf{\sigma}}=\left[\begin{array}{ccccc|ccccc}
0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0\\
0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0\\
0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0\\
0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0\\
0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0\\
\hline
0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0\\
0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0\\
0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0\\
0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0\\
0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0\\
\end{array}\right]
\end{equation}
```

For ``f`` system,
```math
\begin{equation}
\vec{\mathbf{l}}\cdot\vec{\mathbf{\sigma}}=\left[\begin{array}{ccccccc|ccccccc}
0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0\\
0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0\\
0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0\\
0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0\\
0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0\\
0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0\\
0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0\\
\hline
0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0\\
0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0\\
0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0\\
0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0\\
0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0\\
0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0\\
0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0\\
\end{array}\right]
\end{equation}
```
