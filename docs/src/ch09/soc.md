### Spin-orbit coupling

The spin-orbit coupling (SOC) is implemented at atomic level, 

```math
\begin{equation}
\hat{H}_{\text{SOC}}=\lambda\sum_{i}\vec{\mathbf{l}}_{i}\cdot\vec{\mathbf{s}}_{i},
\end{equation}
```

where, ``\vec{\mathbf{l}}`` is orbital angular momentum, and ``\vec{\mathbf{s}}`` is spin angular momentum. In second quantization form,

```math
\begin{equation}
\hat{H}_{\text{SOC}}=\lambda\sum_{\alpha\sigma,\beta\sigma^{\prime}}\left\langle\alpha\sigma\left|\vec{\mathbf{l}}\cdot\vec{\mathbf{s}}\right|\beta\sigma^{\prime}\right\rangle\hat{f}_{\alpha\sigma}^{\dagger}\hat{f}_{\beta\sigma^{\prime}},
\end{equation}
```

where, ``\alpha`` is orbital index and ``\sigma`` is spin index, and 

```math
\begin{equation}
\vec{\mathbf{l}}\cdot\vec{\mathbf{s}} =  \frac{1}{2}\vec{\mathbf{l}}\cdot\vec{\mathbf{\mathbf{\sigma},}}
\end{equation}
```

where, ``\vec{\mathbf{\sigma}}`` is Pauli operator.

Now the question is how to write down the matrix elements for ``\vec{\mathbf{l}}\cdot\vec{\mathbf{\sigma}}``

```math
\begin{align*}
\vec{\mathbf{l}}\cdot\vec{\mathbf{\sigma}} & = & \hat{l}_{x}\hat{\sigma}_{x}+\hat{l}_{y}\hat{\sigma}_{y}+\hat{l}_{z}\hat{\sigma}_{z}\\
 & = & \left[\begin{array}{cc}
0 & \hat{l}_{x}\\
\hat{l}_{x} & 0
\end{array}\right]+\left[\begin{array}{cc}
0 & -i\hat{l}_{y}\\
i\hat{l}_{y} & 0
\end{array}\right]+\left[\begin{array}{cc}
\hat{l}_{z} & 0\\
0 & -\hat{l}_{z}
\end{array}\right]\\
 & = & \left[\begin{array}{cc}
\hat{l}_{z} & \hat{l}_{x}-i\hat{l}_{y}\\
\hat{l}_{x}+i\hat{l}_{y} & -\hat{l}_{z}
\end{array}\right]\\
 & = & \left[\begin{array}{cc}
\hat{l}_{z} & \hat{l}_{-}\\
\hat{l}_{+} & -\hat{l}_{z}
\end{array}\right]
\end{align*}
```

where, ``\hat{l}_{\pm}=\hat{l}_{x}\pm\hat{l}_{y}``, and

```math
\begin{equation}
\hat{l}_{\pm}Y_{l}^{m}=\sqrt{(l\mp m)(l\pm m+1)}Y_{l}^{m\pm1}.
\end{equation}
```

We just write down ``\vec{\mathbf{l}}\cdot\vec{\mathbf{\sigma}}`` in the complex shperical harmonics basis, the orbital order is:
```math
\begin{equation}
Y_{l}^{-l}\uparrow,Y_{l}^{-l}\downarrow,Y_{l}^{-l+1}\uparrow,Y_{l}^{-l+1}\downarrow,\cdots,Y_{l}^{l}\uparrow,Y_{l}^{l}\downarrow.
\end{equation}
```

For ``p`` system,
```math
\begin{equation}
\vec{\mathbf{l}}\cdot\vec{\mathbf{\sigma}}=\left[\begin{array}{cccccc}
-1 & 0 & 0 & \sqrt{2} & 0 & 0\\
0 & 1 & 0 & 0 & 0 & 0\\
0 & 0 & 0 & 0 & 0 & \sqrt{2}\\
\sqrt{2} & 0 & 0 & 0 & 0 & 0\\
0 & 0 & 0 & 0 & 1 & 0\\
0 & 0 & \sqrt{2} & 0 & 0 & -1
\end{array}\right]
\end{equation}
```

For ``t_{2g}`` system, 
```math
\begin{equation}
\vec{\mathbf{l}}\cdot\vec{\mathbf{\sigma}}=-\left[\begin{array}{cccccc}
-1 & 0 & 0 & \sqrt{2} & 0 & 0\\
0 & 1 & 0 & 0 & 0 & 0\\
0 & 0 & 0 & 0 & 0 & \sqrt{2}\\
\sqrt{2} & 0 & 0 & 0 & 0 & 0\\
0 & 0 & 0 & 0 & 1 & 0\\
0 & 0 & \sqrt{2} & 0 & 0 & -1
\end{array}\right]
\end{equation}
```

For ``d`` system,
```math
\begin{equation}
\vec{\mathbf{l}}\cdot\vec{\mathbf{\sigma}}=\left[\begin{array}{cccccccccc}
-2 & 0 & 0 & 2 & 0 & 0 & 0 & 0 & 0 & 0\\
0 & 2 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0\\
0 & 0 & -1 & 0 & 0 & \sqrt{6} & 0 & 0 & 0 & 0\\
2 & 0 & 0 & 1 & 0 & 0 & 0 & 0 & 0 & 0\\
0 & 0 & 0 & 0 & 0 & 0 & 0 & \sqrt{6} & 0 & 0\\
0 & 0 & \sqrt{6} & 0 & 0 & 0 & 0 & 0 & 0 & 0\\
0 & 0 & 0 & 0 & 0 & 0 & 1 & 0 & 0 & 2\\
0 & 0 & 0 & 0 & \sqrt{6} & 0 & 0 & -1 & 0 & 0\\
0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 2 & 0\\
0 & 0 & 0 & 0 & 0 & 0 & 2 & 0 & 0 & -2
\end{array}\right]
\end{equation}
```

For ``f`` system,
```math
\begin{equation}
\vec{\mathbf{l}}\cdot\vec{\mathbf{\sigma}}=\left[\begin{array}{cccccccccccccc}
-3 & 0 & 0 & \sqrt{6} & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0\\
0 & 3 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0\\
0 & 0 & -2 & 0 & 0 & \sqrt{10} & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0\\
\sqrt{6} & 0 & 0 & 2 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0\\
0 & 0 & 0 & 0 & -1 & 0 & 0 & \sqrt{12} & 0 & 0 & 0 & 0 & 0 & 0\\
0 & 0 & \sqrt{10} & 0 & 0 & 1 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0\\
0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & \sqrt{12} & 0 & 0 & 0 & 0\\
0 & 0 & 0 & 0 & \sqrt{12} & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0\\
0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 1 & 0 & 0 & \sqrt{10} & 0 & 0\\
0 & 0 & 0 & 0 & 0 & 0 & \sqrt{12} & 0 & 0 & -1 & 0 & 0 & 0 & 0\\
0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 2 & 0 & 0 & \sqrt{6}\\
0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & \sqrt{10} & 0 & 0 & -2 & 0 & 0\\
0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 3 & 0\\
0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & \sqrt{6} & 0 & 0 & -3
\end{array}\right]
\end{equation}
```