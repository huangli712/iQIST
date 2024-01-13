# Single particle basis

In this section, we define some single particle basis used in the **JASMINE** component to write down the atomic Hamiltonian ``\hat{H}_{\text{atom}}``. We set ``\hbar=1`` in this note.

---

**Spherical harmonics**

The complex spherical harmonics ``Y_{l}^{m}(\theta,\phi)`` is the eigenstate of operators ``\hat{l}^{2}`` and ``\hat{l}_{z}``,

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

where ``m=-l,~-l+1,~\cdots,~l``[^1][^2]. The real spherical harmonics ``Y_{lm}`` is defined as

```math
\begin{gather}
Y_{lm}=\begin{cases}
  \frac{i}{\sqrt{2}}\left(Y_{l}^{-|m|}-(-1)^{m}Y_{l}^{|m|}\right) & \text{if}\ m<0,\\
  Y_{l}^{0} & \text{if}\ m=0,\\
  \frac{1}{\sqrt{2}}\left(Y_{l}^{-|m|}+(-1)^{m}Y_{l}^{|m|}\right) & \text{if}\ m>0.
\end{cases}
\end{gather}
```

The spherical harmonics are orthonormal 

```math
\begin{equation}
\int^{\pi}_{\theta = 0} \int^{2\pi}_{\phi = 0}
Y^m_l Y^{m'*}_{l'} d\Omega = \delta_{ll'} \delta_{mm'},
\end{equation}
```
where ``\delta_{ij}`` is the Kronecker delta and ``d\Omega = \sin(\theta) d\phi d\theta``.

---

**Spinor spherical harmonics**

The spinor spherical harmonics ``\Omega^l_{jm_j}(\theta,\phi)`` is eigenstate of the operators ``\hat{j}^2``, ``\hat{j}_z``, ``\hat{l}^2``, and ``\hat{s}^2``,

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
\hat{s}^2 \Omega^l_{jm_j} = \frac{3}{4} \Omega^l_{jm_j}.
\end{equation}
```

For given ``j`` only two values of ``l`` are possible, ``l = j \pm \frac{1}{2}``, while ``m_j`` assumes ``2j + 1`` values: ``m_j = -j,~-j + 1,~\cdots,~j - 1,~j``[^1][^2].

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

The basis functions are the real spherical harmonics ``Y_{lm}``.

For ``p`` system, the orbital order is[^3]

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

For ``d`` system, the orbital order is[^3]

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

For ``f`` system, the orbital order is[^3]

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

It is also called the ``|l^2,l_z\rangle`` basis. The basis functions are the complex spherical harmonics ``Y^{m}_l(\theta,\phi)``.

For ``p`` system (``l = 1,~m = \pm 1,~0``), the orbital order is

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

For ``d`` system (``l = 2,~m = \pm 2,~\pm 1,~0``), the orbital order is

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

For ``f`` system (``l = 3,~m = \pm 3,~\pm 2,~\pm 1,~0``), the orbital order is

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

**``j^{2},j_{z}`` diagonal basis**

Define ``\phi_{ljm_{j}}`` as the eigenstate of ``j^{2},j_{z}``,

```math
\begin{equation}
j^{2}\phi_{ljm_{j}}=j(j+1)\phi_{ljm_{j}},
\end{equation}
```

```math
\begin{equation}
j_{z}\phi_{ljm_{j}}=m_{j}\phi_{ljm_{j}}.
\end{equation}
```

For ``j=l+\frac{1}{2},m_{j}=m+\frac{1}{2}``,

```math
\begin{equation}
\phi_{ljm_{j}}=\sqrt{\frac{l+m+1}{2l+1}}Y_{l}^{m}\uparrow+\sqrt{\frac{l-m}{2l+1}}Y_{l}^{m+1}\downarrow.
\end{equation}
```

For ``j=l-\frac{1}{2},m_{j}=m+\frac{1}{2}``,

```math
\begin{equation}
\phi_{ljm_{j}}=-\sqrt{\frac{l-m}{2l+1}}Y_{l}^{m}\uparrow+\sqrt{\frac{l+m+1}{2l+1}}Y_{l}^{m+1}\downarrow.
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

---

**Reference**

[^1]: 曾谨言, 量子力学（卷1）， 科学出版社， 第四版， 2007。

[^2]: D. A. Varshalovich, A. N. Moskalev, and V. K. Khersonskii, Quantum Theory of Angular Momentum, World Scientific, 1988.

[^3]: The orbital orders are consistent with the definition of local basis used by VASP (see https://www.vasp.at/wiki/index.php/LOCPROJ), and the definition in HandWiki (see https://handwiki.org/wiki/Table\_of\_spherical\_harmonics).
