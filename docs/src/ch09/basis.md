### Definition of single particle basis

In this section, we define some single particle basis used in the **JASMINE** component to write down the atomic Hamiltonian $$\hat{H}_{\text{atom}}$$. We set $$\hbar=1$$ in this note.

**Complex spherical harmonics basis**

The complex spherical harmonics $$Y_{l}^{m}(\theta,\phi)$$ is the eigenstate of $$l^{2},l_{z}$$,

$$
\begin{equation}
l^{2}Y_{l}^{m}=l(l+1)Y_{l}^{m},
\end{equation}
$$
$$
\begin{equation}
l_{z}Y_{l}^{m}=mY_{l}^{m},\ m=-l,-l+1,\cdots,l.
\end{equation}
$$

**Real spherical harmonics basis**

The real spherical harmonics $$Y_{lm}$$ is defined as

$$
\begin{gather}
Y_{lm}=\begin{cases}
\frac{i}{\sqrt{2}}\left(Y_{l}^{-|m|}-(-1)^{m}Y_{l}^{|m|}\right) & \text{if}\ m<0,\\
Y_{l}^{0} & \text{if}\ m=0,\\
\frac{1}{\sqrt{2}}\left(Y_{l}^{-|m|}+(-1)^{m}Y_{l}^{|m|}\right) & \text{if}\ m>0.
\end{cases}
\end{gather}
$$

For $$p$$ system,
$$
\begin{eqnarray}
p_{x} & = & Y_{11}=\frac{1}{\sqrt{2}}\left(Y_{1}^{-1}-Y_{1}^{1}\right),\nonumber \\
p_{y} & = & Y_{1,-1}=\frac{i}{\sqrt{2}}\left(Y_{1}^{-1}+Y_{1}^{1}\right),\\
p_{z} & = & Y_{10}=Y_{1}^{0}.\nonumber 
\end{eqnarray}
$$

For $$d$$ system,
$$
\begin{eqnarray}
d_{z^{2}} & = & Y_{20}=Y_{2}^{0},\nonumber \\
d_{zx} & = & Y_{21}=\frac{1}{\sqrt{2}}\left(Y_{2}^{-1}-Y_{2}^{1}\right),\nonumber \\
d_{zy} & = & Y_{2,-1}=\frac{i}{\sqrt{2}}\left(Y_{2}^{-1}+Y_{2}^{1}\right),\\
d_{x^{2}-y^{2}} & = & Y_{22}=\frac{1}{\sqrt{2}}\left(Y_{2}^{-2}+Y_{2}^{2}\right),\nonumber \\
d_{xy} & = & Y_{2,-2}=\frac{i}{\sqrt{2}}\left(Y_{2}^{-2}-Y_{2}^{2}\right).\nonumber 
\end{eqnarray}
$$

For $$t_{2g}$$ system $$(l\approx-1)$$, we have a $$T-P$$ equivalence, 
$$
\begin{eqnarray}
d_{zx} & \rightarrow & p_{y}=\frac{i}{\sqrt{2}}\left(Y_{1}^{-1}+Y_{1}^{1}\right),\nonumber \\
d_{zy} & \rightarrow & p_{x}=\frac{1}{\sqrt{2}}\left(Y_{1}^{-1}-Y_{1}^{1}\right),\\
d_{xy} & \rightarrow & p_{z}=Y_{1}^{0}.\nonumber 
\end{eqnarray}
$$

For $$f$$ system, 
$$
\begin{eqnarray}
f_{z^{3}} & = & Y_{30}=Y_{3}^{0},\nonumber \\
f_{xz^{2}} & = & Y_{31}=\frac{1}{\sqrt{2}}\left(Y_{3}^{-1}-Y_{3}^{1}\right),\nonumber \\
f_{yz^{2}} & = & Y_{3,-1}=\frac{i}{\sqrt{2}}\left(Y_{3}^{-1}+Y_{3}^{1}\right),\nonumber \\
f_{z(x^{2}-y^{2})} & = & Y_{32}=\frac{1}{\sqrt{2}}\left(Y_{3}^{-2}+Y_{3}^{2}\right),\\
f_{xyz} & = & Y_{3,-2}=\frac{i}{\sqrt{2}}\left(Y_{3}^{-2}-Y_{3}^{2}\right),\nonumber \\
f_{x(x^{2}-3y^{2})} & = & Y_{33}=\frac{1}{\sqrt{2}}\left(Y_{3}^{-3}-Y_{3}^{3}\right),\nonumber \\
f_{y(3x^{2}-y^{2})} & = & Y_{3,-3}=\frac{i}{\sqrt{2}}\left(Y_{3}^{-3}+Y_{3}^{3}\right).\nonumber 
\end{eqnarray}
$$

**Cubic spherical harmonics basis**

The cubic spherical harmonics is defined as the basis of the irreducible representation of cubic point group $$O_{h}$$.

For $$p$$ orbitals,
$$
\begin{equation}
T_{1u}:p_{x},p_{y},p_{z}.
\end{equation}
$$

For $$d$$ orbitals,
$$
\begin{gather}
\begin{cases}
E_{g}: & d_{z^{2}},d_{x^{2}-y^{2}},\\
T_{2g}: & d_{zx},d_{zy},d_{xy}.
\end{cases}
\end{gather}
$$

For $$f$$ orbitals,
$$
\begin{eqnarray}
f_{x^{3}} & = & -\frac{\sqrt{6}}{4}f_{xz^{2}}+\frac{\sqrt{10}}{4}f_{x(x^{2}-3y^{2})},\nonumber \\
f_{y^{3}} & = & -\frac{\sqrt{6}}{4}f_{yz^{2}}-\frac{\sqrt{10}}{4}f_{y(3x^{2}-y^{2})},\nonumber \\
f_{z^{3}} & = & f_{z^{3}},\nonumber \\
f_{x(y^{2}-z^{2})} & = & -\frac{\sqrt{10}}{4}f_{xz^{2}}-\frac{\sqrt{6}}{4}f_{x(x^{2}-3y^{2})},\nonumber \\
f_{y(z^{2}-x^{2})} & = & \frac{\sqrt{10}}{4}f_{yz^{2}}-\frac{\sqrt{6}}{4}f_{y(3x^{2}-y^{2})},\nonumber \\
f_{z(x^{2}-y^{2})} & = & f_{z(x^{2}-y^{2})},\nonumber \\
f_{xyz} & = & f_{xyz}.\nonumber 
\end{eqnarray}
$$

> NOTE:
$$
\begin{gather}
\begin{cases}
T_{1u}: & f_{x^{3}},f_{y^{3}},f_{z^{3}}\\
T_{2u}: & f_{x(y^{2}-z^{2})},f_{y(z^{2}-x^{2})},f_{z(x^{2}-y^{2})}\\
A_{2u}: & f_{xyz}
\end{cases}
\end{gather}
$$

**$$j^{2},j_{z}$$ diagonal basis**

Define $$\phi_{ljm_{j}}$$ as the eigenstate of $$j^{2},j_{z}$$,

$$
\begin{equation}
j^{2}\phi_{ljm_{j}}=j(j+1)\phi_{ljm_{j}},
\end{equation}
$$

$$
\begin{equation}
j_{z}\phi_{ljm_{j}}=m_{j}\phi_{ljm_{j}}.
\end{equation}
$$

For $$j=l+\frac{1}{2},m_{j}=m+\frac{1}{2}$$,
$$
\begin{equation}
\phi_{ljm_{j}}=\sqrt{\frac{l+m+1}{2l+1}}Y_{l}^{m}\uparrow+\sqrt{\frac{l-m}{2l+1}}Y_{l}^{m+1}\downarrow.
\end{equation}
$$

For $$j=l-\frac{1}{2},m_{j}=m+\frac{1}{2}$$,
$$
\begin{equation}
\phi_{ljm_{j}}=-\sqrt{\frac{l-m}{2l+1}}Y_{l}^{m}\uparrow+\sqrt{\frac{l+m+1}{2l+1}}Y_{l}^{m+1}\downarrow.
\end{equation}
$$

**Natural basis**

The natural basis is defined as the diagonal basis of on-site term $$E_{\alpha\beta}$$.