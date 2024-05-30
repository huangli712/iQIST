# Single particle basis

In this section, we define some single particle basis used in the **JASMINE** component to write down the atomic Hamiltonian ``\hat{H}_{\text{atom}}``. We set ``\hbar=1`` in this note.

---

**Spherical harmonics**

The complex spherical harmonics ``Y_{l}^{m}(\theta,\phi)`` are the eigenstates of operators ``\hat{l}^{2}`` and ``\hat{l}_{z}``,

```math
\begin{equation}
\hat{l}^{2}Y_{l}^{m}=l(l+1)Y_{l}^{m},
\end{equation}
```

```math
\begin{equation}
\hat{l}_{z}Y_{l}^{m}=mY_{l}^{m},
\end{equation}
```

where ``l`` is the azimuthal quantum number (``l = 0,~1,~2,~\cdots,~n-1``), and ``m`` is the magnetic quantum number (``m=-l,~-l+1,~\cdots,~l``)[^1][^2]. They are defined as follows:

```math
\begin{equation}
Y^m_l(\theta,\phi) = \sqrt{\frac{2l+1}{4\pi}\frac{(l-m)!}{(l+m)!}}
P^m_l(\cos{\theta}) e^{im\phi},
\end{equation}
```

where ``\theta`` is taken as the polar (colatitudinal) coordinate with ``\theta \in [0,\pi]``, and ``\phi`` as the azimuthal (longitudinal) coordinate with ``\phi \in [0,2\pi]``, and ``P^m_l(z)`` is an associated Legendre polynomial.

The spherical harmonics are orthonormal

```math
\begin{equation}
\int^{\pi}_{\theta = 0} \int^{2\pi}_{\phi = 0}
Y^m_l(\theta,\phi) Y^{m'*}_{l'}(\theta,\phi)~d\Omega = \delta_{ll'} \delta_{mm'},
\end{equation}
```
where ``\delta_{ij}`` is the Kronecker delta and ``d\Omega = \sin(\theta) d\phi d\theta``.

---

**Real Spherical Harmonics**

The real spherical harmonics ``Y_{lm}`` are defined as[^3]

```math
\begin{gather}
Y_{lm}=\begin{cases}
  \frac{i}{\sqrt{2}}\left(Y_{l}^{-|m|}-(-1)^{m}Y_{l}^{|m|}\right) & \text{if}\ m<0,\\
  Y_{l}^{0} & \text{if}\ m=0,\\
  \frac{1}{\sqrt{2}}\left(Y_{l}^{-|m|}+(-1)^{m}Y_{l}^{|m|}\right) & \text{if}\ m>0.
\end{cases}
\end{gather}
```

The corresponding inverse equations defining the complex spherical harmonics ``Y^m_l`` in terms of the real spherical harmonics ``Y_{lm}`` read:

```math
\begin{gather}
Y_{l}^{m}=\begin{cases}
  \frac{1}{\sqrt{2}}\left(Y_{l|m|}-iY_{l,-|m|}\right) & \text{if}\ m<0,\\
  Y_{l0} & \text{if}\ m=0,\\
  \frac{(-1)^m}{\sqrt{2}}\left(Y_{l|m|}+iY_{l,-|m|}\right) & \text{if}\ m>0.
\end{cases}
\end{gather}
```

The real spherical harmonics ``Y_{lm}`` are sometimes known as tesseral spherical harmonics. These functions have the same orthonormality properties as the complex ones ``Y_{l}^{m}``.


---

**Spinor spherical harmonics**

The spinor spherical harmonics ``\Omega^l_{jm_j}(\theta,\phi)`` are eigenstates of the operators ``\hat{j}^2``, ``\hat{j}_z``, ``\hat{l}^2``, and ``\hat{s}^2``,

```math
\begin{equation}
\hat{j}^2 \Omega^l_{jm_j} = j (j + 1) \Omega^l_{jm_j},
\end{equation}
```

```math
\begin{equation}
\hat{j}_z \Omega^l_{jm_j} = m_j \Omega^l_{jm_j},
\end{equation}
```

```math
\begin{equation}
\hat{l}^2 \Omega^l_{jm_j} = l (l + 1) \Omega^l_{jm_j},
\end{equation}
```

```math
\begin{equation}
\hat{s}^2 \Omega^l_{jm_j} = s(s+1) \Omega^l_{jm_j}.
\end{equation}
```

For given ``j`` only two values of ``l`` are possible, ``l = j \pm \frac{1}{2}``, while ``m_j`` assumes ``2j + 1`` values (``m_j = -j,~-j + 1,~\cdots,~j - 1,~j``) [^1][^2].

For ``j = l + \frac{1}{2}``, ``m_j = m + \frac{1}{2}``,

```math
\begin{equation}
\Omega^l_{jm_j} = \sqrt{\frac{l+m+1}{2l+1}} Y^m_l \chi_{\uparrow}
                + \sqrt{\frac{l-m}{2l+1}} Y^{m+1}_l \chi_{\downarrow}.
\end{equation}
```

For ``j = l - \frac{1}{2} (l \neq 0)``, ``m_j = m + \frac{1}{2}``,

```math
\begin{equation}
\Omega^l_{jm_j} = -\sqrt{\frac{l-m}{2l+1}} Y^m_l \chi_{\uparrow}
                + \sqrt{\frac{l+m+1}{2l+1}} Y^{m+1}_l \chi_{\downarrow}.
\end{equation}
```

---

**Real orbital basis**

The basis functions are the real spherical harmonics ``Y_{lm}(\theta,\phi)``.

For ``p`` system, the basis order is[^4]

```math
\begin{equation}
\begin{split}
&
|p_y, \uparrow \rangle,~
|p_z, \uparrow \rangle,~
|p_x, \uparrow \rangle,~\\
&
|p_y, \downarrow \rangle,~
|p_z, \downarrow \rangle,~
|p_x, \downarrow \rangle.
\end{split}
\end{equation}
```

```math
\begin{equation}
|p_{y}\rangle = Y_{1,-1}=\frac{i}{\sqrt{2}}\left(Y_{1}^{-1}+Y_{1}^{1}\right),
\end{equation}
```

```math
\begin{equation}
|p_{z}\rangle =  Y_{10}=Y_{1}^{0},
\end{equation}
```

```math
\begin{equation}
|p_{x}\rangle  =  Y_{11}=\frac{1}{\sqrt{2}}\left(Y_{1}^{-1}-Y_{1}^{1}\right).
\end{equation}
```

For ``d`` system, the basis order is[^4]

```math
\begin{equation}
\begin{split}
&
|d_{xy}, \uparrow \rangle,~
|d_{yz}, \uparrow \rangle,~
|d_{z^2}, \uparrow \rangle,~
|d_{xz}, \uparrow \rangle,~
|d_{x^2-y^2}, \uparrow \rangle,~\\
&
|d_{xy}, \downarrow \rangle,~
|d_{yz}, \downarrow \rangle,~
|d_{z^2}, \downarrow \rangle,~
|d_{xz}, \downarrow \rangle,~
|d_{x^2-y^2}, \downarrow \rangle.
\end{split}
\end{equation}
```

```math
\begin{equation}
d_{xy} = Y_{2,-2}=\frac{i}{\sqrt{2}}\left(Y_{2}^{-2}-Y_{2}^{2}\right),
\end{equation}
```

```math
\begin{equation}
d_{yz} = Y_{2,-1}=\frac{i}{\sqrt{2}}\left(Y_{2}^{-1}+Y_{2}^{1}\right),
\end{equation}
```

```math
\begin{equation}
d_{z^{2}} = Y_{20}=Y_{2}^{0},
\end{equation}
```

```math
\begin{equation}
d_{xz} = Y_{21}=\frac{1}{\sqrt{2}}\left(Y_{2}^{-1}-Y_{2}^{1}\right),
\end{equation}
```

```math
\begin{equation}
d_{x^{2}-y^{2}} = Y_{22}=\frac{1}{\sqrt{2}}\left(Y_{2}^{-2}+Y_{2}^{2}\right).
\end{equation}
```

For ``f`` system, the basis order is[^4]

```math
\begin{equation}
\begin{split}
&
|f_{y(3x^2-y^2)}, \uparrow \rangle,~
|f_{xyz}, \uparrow \rangle,~
|f_{yz^2}, \uparrow \rangle,~
|f_{z^3}, \uparrow \rangle,~
|f_{xz^2}, \uparrow \rangle,~
|f_{z(x^2-y^2)}, \uparrow \rangle,~
|f_{x(x^2-3y^2)}, \uparrow \rangle,\\
&
|f_{y(3x^2-y^2)}, \downarrow \rangle,~
|f_{xyz}, \downarrow \rangle,~
|f_{yz^2}, \downarrow \rangle,~
|f_{z^3}, \downarrow \rangle,~
|f_{xz^2}, \downarrow \rangle,~
|f_{z(x^2-y^2)}, \downarrow \rangle,~
|f_{x(x^2-3y^2)}, \downarrow \rangle.
\end{split}
\end{equation}
```

```math
\begin{equation}
f_{y(3x^{2}-y^{2})} = Y_{3,-3}=\frac{i}{\sqrt{2}}\left(Y_{3}^{-3}+Y_{3}^{3}\right),
\end{equation}
```

```math
\begin{equation}
f_{xyz} = Y_{3,-2}=\frac{i}{\sqrt{2}}\left(Y_{3}^{-2}-Y_{3}^{2}\right),
\end{equation}
```

```math
\begin{equation}
f_{yz^{2}} = Y_{3,-1}=\frac{i}{\sqrt{2}}\left(Y_{3}^{-1}+Y_{3}^{1}\right),
\end{equation}
```

```math
\begin{equation}
f_{z^{3}} = Y_{30}=Y_{3}^{0},
\end{equation}
```

```math
\begin{equation}
f_{xz^{2}} = Y_{31}=\frac{1}{\sqrt{2}}\left(Y_{3}^{-1}-Y_{3}^{1}\right),
\end{equation}
```

```math
\begin{equation}
f_{z(x^{2}-y^{2})} = Y_{32}=\frac{1}{\sqrt{2}}\left(Y_{3}^{-2}+Y_{3}^{2}\right),
\end{equation}
```

```math
\begin{equation}
f_{x(x^{2}-3y^{2})} = Y_{33}=\frac{1}{\sqrt{2}}\left(Y_{3}^{-3}-Y_{3}^{3}\right).
\end{equation}
```

For ``t_{2g}`` system, we have a ``T-P`` equivalence,

```math
\begin{equation}
d_{xz} \rightarrow p_{y}=\frac{i}{\sqrt{2}}\left(Y_{1}^{-1}+Y_{1}^{1}\right),
\end{equation}
```

```math
\begin{equation}
d_{xy} \rightarrow p_{z}=Y_{1}^{0},
\end{equation}
```

```math
\begin{equation}
d_{yz} \rightarrow p_{x}=\frac{1}{\sqrt{2}}\left(Y_{1}^{-1}-Y_{1}^{1}\right).
\end{equation}
```

---

**Complex orbital basis**

It is also called the ``|l^2,l_z\rangle`` basis. The basis functions are the complex spherical harmonics ``Y^{m}_l(\theta,\phi)``. We just use ``l`` and ``m`` to label the basis functions ``|l,m\rangle``.

For ``p`` system (``l = 1,~m = \pm 1,~0``), the basis order is

```math
\begin{equation}
\begin{split}
&| 1,   -1, \uparrow \rangle, \\
&| 1, ~~~0, \uparrow \rangle, \\
&| 1, ~~~1, \uparrow \rangle, \\
&| 1,   -1, \downarrow \rangle, \\
&| 1, ~~~0, \downarrow \rangle, \\
&| 1, ~~~1, \downarrow \rangle.
\end{split}
\end{equation}
```

For ``d`` system (``l = 2,~m = \pm 2,~\pm 1,~0``), the basis order is

```math
\begin{equation}
\begin{split}
&| 2,   -2, \uparrow \rangle, \\
&| 2,   -1, \uparrow \rangle, \\
&| 2, ~~~0, \uparrow \rangle, \\
&| 2, ~~~1, \uparrow \rangle, \\
&| 2, ~~~2, \uparrow \rangle, \\
&| 2,   -2, \downarrow \rangle, \\
&| 2,   -1, \downarrow \rangle, \\
&| 2, ~~~0, \downarrow \rangle, \\
&| 2, ~~~1, \downarrow \rangle, \\
&| 2, ~~~2, \downarrow \rangle.
\end{split}
\end{equation}
```

For ``f`` system (``l = 3,~m = \pm 3,~\pm 2,~\pm 1,~0``), the basis order is

```math
\begin{equation}
\begin{split}
&| 3,   -3, \uparrow \rangle, \\
&| 3,   -2, \uparrow \rangle, \\
&| 3,   -1, \uparrow \rangle, \\
&| 3, ~~~0, \uparrow \rangle, \\
&| 3, ~~~1, \uparrow \rangle, \\
&| 3, ~~~2, \uparrow \rangle, \\
&| 3, ~~~3, \uparrow \rangle, \\
&| 3,   -3, \downarrow \rangle, \\
&| 3,   -2, \downarrow \rangle, \\
&| 3,   -1, \downarrow \rangle, \\
&| 3, ~~~0, \downarrow \rangle, \\
&| 3, ~~~1, \downarrow \rangle, \\
&| 3, ~~~2, \downarrow \rangle, \\
&| 3, ~~~3, \downarrow \rangle.
\end{split}
\end{equation}
```

---

**``\hat{j}^{2}-\hat{j}_{z}-\hat{l}^2-\hat{s}^2`` diagonal basis**

We just use ``j`` and ``m_j`` to label the eigenfunctions ``|j, m_j\rangle``, which are just the spinor spherical harmonics ``\Omega^l_{jm_j}(\theta,\phi)``.

For ``p`` system, ``l = 1``, ``j = \frac{1}{2}`` or ``\frac{3}{2}``, the basis order is

```math
\begin{equation}
\begin{split}
&\left|\frac{1}{2}, -\frac{1}{2}\right\rangle   = -\sqrt{\frac{2}{3}}Y^{-1}_{1}\chi_{\uparrow} + \sqrt{\frac{1}{3}}Y^{ 0}_{1}\chi_{\downarrow}, \\
&\left|\frac{1}{2}, ~~~\frac{1}{2}\right\rangle = -\sqrt{\frac{1}{3}}Y^{ 0}_{1}\chi_{\uparrow} + \sqrt{\frac{2}{3}}Y^{ 1}_{1}\chi_{\downarrow}, \\
&\left|\frac{3}{2}, -\frac{3}{2}\right\rangle   =  \sqrt{\frac{0}{3}}Y^{-2}_{1}\chi_{\uparrow} + \sqrt{\frac{3}{3}}Y^{-1}_{1}\chi_{\downarrow} = Y^{-1}_{1} \chi_{\downarrow}, \\
&\left|\frac{3}{2}, -\frac{1}{2}\right\rangle   =  \sqrt{\frac{1}{3}}Y^{-1}_{1}\chi_{\uparrow} + \sqrt{\frac{2}{3}}Y^{ 0}_{1}\chi_{\downarrow}, \\
&\left|\frac{3}{2}, ~~~\frac{1}{2}\right\rangle =  \sqrt{\frac{2}{3}}Y^{ 0}_{1}\chi_{\uparrow} + \sqrt{\frac{1}{3}}Y^{ 1}_{1}\chi_{\downarrow}, \\
&\left|\frac{3}{2}, ~~~\frac{3}{2}\right\rangle =  \sqrt{\frac{3}{3}}Y^{ 1}_{1}\chi_{\uparrow} + \sqrt{\frac{0}{3}}Y^{ 2}_{1}\chi_{\downarrow} = Y^{1}_{1}\chi_{\uparrow}.
\end{split}
\end{equation}
```

For ``d`` system, ``l = 2``, ``j = \frac{3}{2}`` or ``\frac{5}{2}``, the basis order is

```math
\begin{equation}
\begin{split}
&\left|\frac{3}{2}, -\frac{3}{2}\right\rangle   = -\sqrt{\frac{4}{5}}Y^{-2}_{2}\chi_{\uparrow} + \sqrt{\frac{1}{5}}Y^{-1}_{2}\chi_{\downarrow}, \\
&\left|\frac{3}{2}, -\frac{1}{2}\right\rangle   = -\sqrt{\frac{3}{5}}Y^{-1}_{2}\chi_{\uparrow} + \sqrt{\frac{2}{5}}Y^{ 0}_{2}\chi_{\downarrow}, \\
&\left|\frac{3}{2}, ~~~\frac{1}{2}\right\rangle = -\sqrt{\frac{2}{5}}Y^{ 0}_{2}\chi_{\uparrow} + \sqrt{\frac{3}{5}}Y^{ 1}_{2}\chi_{\downarrow}, \\
&\left|\frac{3}{2}, ~~~\frac{3}{2}\right\rangle = -\sqrt{\frac{1}{5}}Y^{ 1}_{2}\chi_{\uparrow} + \sqrt{\frac{4}{5}}Y^{ 2}_{2}\chi_{\downarrow}, \\
&\left|\frac{5}{2}, -\frac{5}{2}\right\rangle   =  \sqrt{\frac{0}{5}}Y^{-3}_{2}\chi_{\uparrow} + \sqrt{\frac{5}{5}}Y^{-2}_{2}\chi_{\downarrow} = Y^{-2}_{2}\chi_{\downarrow}, \\
&\left|\frac{5}{2}, -\frac{3}{2}\right\rangle   =  \sqrt{\frac{1}{5}}Y^{-2}_{2}\chi_{\uparrow} + \sqrt{\frac{4}{5}}Y^{-1}_{2}\chi_{\downarrow}, \\
&\left|\frac{5}{2}, -\frac{1}{2}\right\rangle   =  \sqrt{\frac{2}{5}}Y^{-1}_{2}\chi_{\uparrow} + \sqrt{\frac{3}{5}}Y^{ 0}_{2}\chi_{\downarrow}, \\
&\left|\frac{5}{2}, ~~~\frac{1}{2}\right\rangle =  \sqrt{\frac{3}{5}}Y^{ 0}_{2}\chi_{\uparrow} + \sqrt{\frac{2}{5}}Y^{ 1}_{2}\chi_{\downarrow}, \\
&\left|\frac{5}{2}, ~~~\frac{3}{2}\right\rangle =  \sqrt{\frac{4}{5}}Y^{ 1}_{2}\chi_{\uparrow} + \sqrt{\frac{1}{5}}Y^{ 2}_{2}\chi_{\downarrow}, \\
&\left|\frac{5}{2}, ~~~\frac{5}{2}\right\rangle =  \sqrt{\frac{5}{5}}Y^{ 2}_{2}\chi_{\uparrow} + \sqrt{\frac{0}{5}}Y^{ 3}_{2}\chi_{\downarrow} = Y^{2}_{2}\chi_{\uparrow}.
\end{split}
\end{equation}
```

For ``f`` system, ``l = 3``, ``j = \frac{5}{2}`` or ``\frac{7}{2}``, the basis order is

```math
\begin{equation}
\begin{split}
&\left|\frac{5}{2}, -\frac{5}{2}\right\rangle   = -\sqrt{\frac{6}{7}}Y^{-3}_{3}\chi_{\uparrow} + \sqrt{\frac{1}{7}}Y^{-2}_{3}\chi_{\downarrow}, \\
&\left|\frac{5}{2}, -\frac{3}{2}\right\rangle   = -\sqrt{\frac{5}{7}}Y^{-2}_{3}\chi_{\uparrow} + \sqrt{\frac{2}{7}}Y^{-1}_{3}\chi_{\downarrow}, \\
&\left|\frac{5}{2}, -\frac{1}{2}\right\rangle   = -\sqrt{\frac{4}{7}}Y^{-1}_{3}\chi_{\uparrow} + \sqrt{\frac{3}{7}}Y^{ 0}_{3}\chi_{\downarrow}, \\
&\left|\frac{5}{2}, ~~~\frac{1}{2}\right\rangle = -\sqrt{\frac{3}{7}}Y^{ 0}_{3}\chi_{\uparrow} + \sqrt{\frac{4}{7}}Y^{ 1}_{3}\chi_{\downarrow}, \\
&\left|\frac{5}{2}, ~~~\frac{3}{2}\right\rangle = -\sqrt{\frac{2}{7}}Y^{ 1}_{3}\chi_{\uparrow} + \sqrt{\frac{5}{7}}Y^{ 2}_{3}\chi_{\downarrow}, \\
&\left|\frac{5}{2}, ~~~\frac{5}{2}\right\rangle = -\sqrt{\frac{1}{7}}Y^{ 2}_{3}\chi_{\uparrow} + \sqrt{\frac{6}{7}}Y^{ 3}_{3}\chi_{\downarrow}, \\
&\left|\frac{7}{2}, -\frac{7}{2}\right\rangle   =  \sqrt{\frac{0}{7}}Y^{-4}_{3}\chi_{\uparrow} + \sqrt{\frac{7}{7}}Y^{-3}_{3}\chi_{\downarrow} = Y^{-3}_{3}\chi_{\downarrow}, \\
&\left|\frac{7}{2}, -\frac{5}{2}\right\rangle   =  \sqrt{\frac{1}{7}}Y^{-3}_{3}\chi_{\uparrow} + \sqrt{\frac{6}{7}}Y^{-2}_{3}\chi_{\downarrow}, \\
&\left|\frac{7}{2}, -\frac{3}{2}\right\rangle   =  \sqrt{\frac{2}{7}}Y^{-2}_{3}\chi_{\uparrow} + \sqrt{\frac{5}{7}}Y^{-1}_{3}\chi_{\downarrow}, \\
&\left|\frac{7}{2}, -\frac{1}{2}\right\rangle   =  \sqrt{\frac{3}{7}}Y^{-1}_{3}\chi_{\uparrow} + \sqrt{\frac{4}{7}}Y^{ 0}_{3}\chi_{\downarrow}, \\
&\left|\frac{7}{2}, ~~~\frac{1}{2}\right\rangle =  \sqrt{\frac{4}{7}}Y^{ 0}_{3}\chi_{\uparrow} + \sqrt{\frac{3}{7}}Y^{ 1}_{3}\chi_{\downarrow}, \\
&\left|\frac{7}{2}, ~~~\frac{3}{2}\right\rangle =  \sqrt{\frac{5}{7}}Y^{ 1}_{3}\chi_{\uparrow} + \sqrt{\frac{2}{7}}Y^{ 2}_{3}\chi_{\downarrow}, \\
&\left|\frac{7}{2}, ~~~\frac{5}{2}\right\rangle =  \sqrt{\frac{6}{7}}Y^{ 2}_{3}\chi_{\uparrow} + \sqrt{\frac{1}{7}}Y^{ 3}_{3}\chi_{\downarrow}, \\
&\left|\frac{7}{2}, ~~~\frac{7}{2}\right\rangle =  \sqrt{\frac{7}{7}}Y^{ 3}_{3}\chi_{\uparrow} + \sqrt{\frac{0}{7}}Y^{ 4}_{3}\chi_{\downarrow} = Y^{3}_{3}\chi_{\uparrow}. \\
\end{split}
\end{equation}
```

---

**Natural basis**

The natural basis is defined as the diagonal basis of on-site term ``E_{\alpha\beta}``.

---

**Transformation matrix from complex orbital basis to real orbital basis**

For *p* system, the transformation matrix reads

```math
\begin{equation}
T = \left[
\begin{array}{ccc|ccc}
\frac{i}{\sqrt{2}} & 0 & \frac{1}{\sqrt{2}} & 0 & 0 & 0 \\
0 & 1 & 0 & 0 & 0 & 0 \\
\frac{i}{\sqrt{2}} & 0 & -\frac{1}{\sqrt{2}} & 0 & 0 & 0 \\
\hline
0 & 0 & 0 & \frac{i}{\sqrt{2}} & 0 & \frac{1}{\sqrt{2}} \\
0 & 0 & 0 & 0 & 1 & 0 \\
0 & 0 & 0 & \frac{i}{\sqrt{2}} & 0 & -\frac{1}{\sqrt{2}} \\
\end{array}
\right]
\end{equation}
```

For *d* system, the transformation matrix reads

```math
T = \left[
\begin{array}{ccccc|ccccc}
\frac{i}{\sqrt{2}} & 0 & 0 & 0 & \frac{1}{\sqrt{2}} & 0 & 0 & 0 & 0 & 0 \\
0 & \frac{i}{\sqrt{2}} & 0 & \frac{1}{\sqrt{2}} & 0 & 0 & 0 & 0 & 0 & 0 \\
0 & 0 & 1 & 0 & 0 & 0 & 0 & 0 & 0 & 0 \\
0 & \frac{i}{\sqrt{2}} & 0 & -\frac{1}{\sqrt{2}} & 0 & 0 & 0 & 0 & 0 & 0 \\
-\frac{i}{\sqrt{2}} & 0 & 0 & 0 & \frac{1}{\sqrt{2}} & 0 & 0 & 0 & 0 & 0 \\
\hline
0 & 0 & 0 & 0 & 0 & \frac{i}{\sqrt{2}} & 0 & 0 & 0 & \frac{1}{\sqrt{2}} \\
0 & 0 & 0 & 0 & 0 & 0 & \frac{i}{\sqrt{2}} & 0 & \frac{1}{\sqrt{2}} & 0 \\
0 & 0 & 0 & 0 & 0 & 0 & 0 & 1 & 0 & 0 \\
0 & 0 & 0 & 0 & 0 & 0 & \frac{i}{\sqrt{2}} & 0 & -\frac{1}{\sqrt{2}} & 0 \\
0 & 0 & 0 & 0 & 0 & -\frac{i}{\sqrt{2}} & 0 & 0 & 0 & \frac{1}{\sqrt{2}} \\
\end{array}
\right]
```

For *f* system, the transformation matrix reads

```math
T= \left[
\begin{array}{ccccccc|ccccccc}
\frac{i}{\sqrt{2}} & 0 & 0 & 0 & 0 & 0 & \frac{1}{\sqrt{2}}     & 0 & 0 & 0 & 0 & 0 & 0 & 0 \\
0 & \frac{i}{\sqrt{2}} & 0 & 0 & 0 & \frac{1}{\sqrt{2}} & 0     & 0 & 0 & 0 & 0 & 0 & 0 & 0 \\
0 & 0 & \frac{i}{\sqrt{2}} & 0 & \frac{1}{\sqrt{2}} & 0 & 0     & 0 & 0 & 0 & 0 & 0 & 0 & 0 \\
0 & 0 & 0 & 1 & 0 & 0 & 0     & 0 & 0 & 0 & 0 & 0 & 0 & 0 \\
0 & 0 & \frac{i}{\sqrt{2}} & 0 & -\frac{1}{\sqrt{2}} & 0 & 0     & 0 & 0 & 0 & 0 & 0 & 0 & 0 \\
0 & -\frac{i}{\sqrt{2}} & 0 & 0 & 0 & \frac{1}{\sqrt{2}} & 0     & 0 & 0 & 0 & 0 & 0 & 0 & 0 \\
\frac{i}{\sqrt{2}} & 0 & 0 & 0 & 0 & 0 & -\frac{1}{\sqrt{2}}     & 0 & 0 & 0 & 0 & 0 & 0 & 0 \\
\hline
0 & 0 & 0 & 0 & 0 & 0 & 0     & \frac{i}{\sqrt{2}} & 0 & 0 & 0 & 0 & 0 & \frac{1}{\sqrt{2}} \\
0 & 0 & 0 & 0 & 0 & 0 & 0     & 0 & \frac{i}{\sqrt{2}} & 0 & 0 & 0 & \frac{1}{\sqrt{2}} & 0 \\
0 & 0 & 0 & 0 & 0 & 0 & 0     & 0 & 0 & \frac{i}{\sqrt{2}} & 0 & \frac{1}{\sqrt{2}} & 0 & 0 \\
0 & 0 & 0 & 0 & 0 & 0 & 0     & 0 & 0 & 0 & 1 & 0 & 0 & 0 \\
0 & 0 & 0 & 0 & 0 & 0 & 0     & 0 & 0 & \frac{i}{\sqrt{2}} & 0 & -\frac{1}{\sqrt{2}} & 0 & 0 \\
0 & 0 & 0 & 0 & 0 & 0 & 0     & 0 & -\frac{i}{\sqrt{2}} & 0 & 0 & 0 & \frac{1}{\sqrt{2}} & 0 \\
0 & 0 & 0 & 0 & 0 & 0 & 0     & \frac{i}{\sqrt{2}} & 0 & 0 & 0 & 0 & 0 & -\frac{1}{\sqrt{2}} \\
\end{array}
\right]
```

The following Julia script is used to construct the complex orbital basis and the real orbital basis, and the transformation matrix between them.

```julia
# To calculate the transformation matrix from the complex orbital basis
# to the real orbital basis.
function calc_matrix(l::Int64)
    println("Construct complex orbital basis for 𝑙 = $l")
    COB = [] # To save the complex orbital basis
    # m = -l, -l+1, ..., l-1, l
    mlist = collect(-l:1:l)
    for s in ("up", "down")
        for m in mlist
            push!(COB, [m, s])
        end
    end
    #
    for i in eachindex(COB)
        m = COB[i][1]
        s = COB[i][2] == "up" ? "↑" : "↓"
        println("$i -> | $l, $m, $s ⟩")
    end

    println("Construct real orbital basis for 𝑙 = $l")
    RO = [] # To save the real orbital basis basis
    ROB = [] # To save the detailed expressions for the real orbital basis
    for s in ("up", "down")
        for m in mlist
            if m < 0
                b = ["i/sqrt(2)", -abs(m), -(-1)^m, abs(m)]
            elseif m == 0
                b = [0]
            elseif m > 0
                b = ["1/sqrt(2)", -abs(m),  (-1)^m, abs(m)]
            end
            push!(RO, [m,s])
            push!(ROB, b)
        end
    end
    #
    for i in eachindex(RO)
        m = RO[i][1]
        s = RO[i][2] == "up" ? "↑" : "↓"
        println("$i -> Y_{$l,$m} χ$s")        
    end

    println("Evaluate transformation matrix for 𝑙 = $l")
    for m in eachindex(COB)
        for n in eachindex(RO)
            if COB[m][2] == RO[n][2] # their spins are the same
                # for Y_{l0} case
                if length(ROB[n]) == 1
                    if COB[m][1] == ROB[n][1]
                        println("T($m,$n) -> 1")
                    end
                # for Y_{lm} case where m /= 0
                else
                    if COB[m][1] == ROB[n][2]
                        println("T($m,$n) -> ", ROB[n][1])
                    end
                    if COB[m][1] == ROB[n][4]
                        s = ROB[n][3] < 0 ? "-" : ""
                        println("T($m,$n) -> $s", ROB[n][1])
                    end
                end
            end
        end
    end
end
```

---

**Transformation matrix from complex orbital basis to ``\hat{j}^{2}-\hat{j}_{z}-\hat{l}^2-\hat{s}^2`` diagonal basis**

For *p* system, the transformation matrix reads

```math
\begin{equation}
T = \left[
\begin{array}{ccc|ccc}
-\sqrt{\frac{2}{3}} & 0 & 0 & \sqrt{\frac{1}{3}} & 0 & 0 \\
0 & -\sqrt{\frac{1}{3}} & 0 & 0 & \sqrt{\frac{2}{3}} & 0 \\
0 & 0 & 0 & 0 & 0 & 1.0 \\
\hline
0 & 0 & 1.0 & 0 & 0 & 0 \\
\sqrt{\frac{1}{3}} & 0 & 0 & \sqrt{\frac{2}{3}} & 0 & 0 \\
0 & \sqrt{\frac{2}{3}} & 0 & 0 & \sqrt{\frac{1}{3}} & 0 \\
\end{array}
\right]
\end{equation}
```

For *d* system, the transformation matrix reads

```math
T = \left[
\begin{array}{ccccc|ccccc}
-\sqrt{\frac{4}{5}} & 0 & 0 & 0 & 0 & \sqrt{\frac{1}{5}} & 0 & 0 & 0 & 0 \\
0 & -\sqrt{\frac{3}{5}} & 0 & 0 & 0 & 0 & \sqrt{\frac{2}{5}} & 0 & 0 & 0 \\
0 & 0 & -\sqrt{\frac{2}{5}} & 0 & 0 & 0 & 0 & \sqrt{\frac{3}{5}} & 0 & 0 \\
0 & 0 & 0 & -\sqrt{\frac{1}{5}} & 0 & 0 & 0 & 0 & \sqrt{\frac{4}{5}} & 0 \\
0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 1.0 \\
\hline
0 & 0 & 0 & 0 & 1.0 & 0 & 0 & 0 & 0 & 0 \\
\sqrt{\frac{1}{5}} & 0 & 0 & 0 & 0 & \sqrt{\frac{4}{5}} & 0 & 0 & 0 & 0 \\
0 & \sqrt{\frac{2}{5}} & 0 & 0 & 0 & 0 & \sqrt{\frac{3}{5}} & 0 & 0 & 0 \\
0 & 0 & \sqrt{\frac{3}{5}} & 0 & 0 & 0 & 0 & \sqrt{\frac{2}{5}} & 0 & 0 \\
0 & 0 & 0 & \sqrt{\frac{4}{5}} & 0 & 0 & 0 & 0 & \sqrt{\frac{1}{5}} & 0 \\
\end{array}
\right]
```

For *f* system, the transformation matrix reads

```math
T = \left[
\begin{array}{ccccccc|ccccccc}
-\sqrt{\frac{6}{7}} & 0 & 0 & 0 & 0 & 0 & 0 & \sqrt{\frac{1}{7}} & 0 & 0 & 0 & 0 & 0 & 0 \\
0 & -\sqrt{\frac{5}{7}} & 0 & 0 & 0 & 0 & 0 & 0 & \sqrt{\frac{2}{7}} & 0 & 0 & 0 & 0 & 0 \\
0 & 0 & -\sqrt{\frac{4}{7}} & 0 & 0 & 0 & 0 & 0 & 0 & \sqrt{\frac{3}{7}} & 0 & 0 & 0 & 0 \\
0 & 0 & 0 & -\sqrt{\frac{3}{7}} & 0 & 0 & 0 & 0 & 0 & 0 & \sqrt{\frac{4}{7}} & 0 & 0 & 0 \\
0 & 0 & 0 & 0 & -\sqrt{\frac{2}{7}} & 0 & 0 & 0 & 0 & 0 & 0 & \sqrt{\frac{5}{7}} & 0 & 0 \\
0 & 0 & 0 & 0 & 0 & -\sqrt{\frac{1}{7}} & 0 & 0 & 0 & 0 & 0 & 0 & \sqrt{\frac{6}{7}} & 0 \\
0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 1.0 \\
\hline
0 & 0 & 0 & 0 & 0 & 0 & 1.0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 \\
\sqrt{\frac{1}{7}} & 0 & 0 & 0 & 0 & 0 & 0 & \sqrt{\frac{6}{7}} & 0 & 0 & 0 & 0 & 0 & 0 \\
0 & \sqrt{\frac{2}{7}} & 0 & 0 & 0 & 0 & 0 & 0 & \sqrt{\frac{5}{7}} & 0 & 0 & 0 & 0 & 0 \\
0 & 0 & \sqrt{\frac{3}{7}} & 0 & 0 & 0 & 0 & 0 & 0 & \sqrt{\frac{4}{7}} & 0 & 0 & 0 & 0 \\
0 & 0 & 0 & \sqrt{\frac{4}{7}} & 0 & 0 & 0 & 0 & 0 & 0 & \sqrt{\frac{3}{7}} & 0 & 0 & 0 \\
0 & 0 & 0 & 0 & \sqrt{\frac{5}{7}} & 0 & 0 & 0 & 0 & 0 & 0 & \sqrt{\frac{2}{7}} & 0 & 0 \\
0 & 0 & 0 & 0 & 0 & \sqrt{\frac{6}{7}} & 0 & 0 & 0 & 0 & 0 & 0 & \sqrt{\frac{1}{7}} & 0 \\
\end{array}
\right]
```

The following Julia script is used to construct the complex orbital basis and the ``j^2-j_z`` basis, and the transformation matrix between them.

```julia
# To calculate the transformation matrix from the complex orbital basis
# to the j²-jz basis.
function calc_matrix(l::Int64)
    println("Construct complex orbital basis for 𝑙 = $l")
    COB = [] # To save the complex orbital basis
    # m = -l, -l+1, ..., l-1, l
    mlist = collect(-l:1:l)
    for s in ("up", "down")
        for m in mlist
            push!(COB, [m, s])
        end
    end
    #
    for i in eachindex(COB)
        m = COB[i][1]
        s = COB[i][2] == "up" ? "↑" : "↓"
        println("$i -> | $l, $m, $s ⟩")
    end

    println("Construct j²-jz basis for 𝑙 = $l")
    JJ = []  # To save the j²-jz basis
    JJB = [] # To save the detailed expressions for the j²-jz basis
    jlist = [l-1//2, l+1//2]
    for j in jlist
        # mⱼ = -j, -j+1, ..., j-1, j
        mⱼlist = collect(-j:2//2:j)
        for mⱼ in mⱼlist
            push!(JJ, [j, mⱼ])
            m = mⱼ-1//2
            if j == l-1//2 # For j = l-1/2
                jj = ["-", (l-m)/(2*l+1), Int(m), (l+m+1)/(2*l+1), Int(m+1)]
            else           # For j = l+1/2
                jj = ["" , (l+m+1)/(2*l+1), Int(m), (l-m)/(2*l+1), Int(m+1)]
            end
            push!(JJB, jj)
        end
    end
    #
    for i in eachindex(JJ)
        j = JJ[i][1]
        mⱼ = JJ[i][2]
        print("$i -> | $j, $mⱼ ⟩ = ")

        jj = JJB[i]
        print(jj[1])
        print("sqrt(", jj[2], ")Y^{", jj[3],"}_{$l}χ↑ + ")
        print("sqrt(", jj[4], ")Y^{", jj[5],"}_{$l}χ↓\n")
    end

    println("Evaluate transformation matrix for 𝑙 = $l")
    for m in eachindex(COB)
        for n in eachindex(JJ)
            if COB[m][2] == "up"
                if COB[m][1] == JJB[n][3]
                    println("T($m,$n) -> ", JJB[n][1], "sqrt(", JJB[n][2],")")
                end
            else
                if COB[m][1] == JJB[n][5]
                    println("T($m,$n) -> ", "sqrt(", JJB[n][4],")")
                end
            end
        end
    end
end
```

---

**Reference**

[^1]: 曾谨言, 量子力学（卷1, 第四版）， 科学出版社， 2007。

[^2]: D. A. Varshalovich, A. N. Moskalev, and V. K. Khersonskii, Quantum Theory of Angular Momentum, World Scientific, 1988.

[^3]: See https://handwiki.org/wiki/Physics:Spherical_harmonics.

[^4]: The orbital orders are consistent with the definition of local basis used by VASP (see https://www.vasp.at/wiki/index.php/LOCPROJ), and the definition in HandWiki (see https://handwiki.org/wiki/Table\_of\_spherical\_harmonics).
