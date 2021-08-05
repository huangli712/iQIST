### Subspace and symmetry

For a local Hamiltonian ``H_{\text{loc}}`` with general interactions, the evaluation of local trace is heavily time-consuming,
```math
\begin{equation}
\omega_{d}(\mathcal{C}) = 
\text{Tr}_{\text{loc}} (T_{2k+1}F_{2k}T_{2k} \cdots F_{2}T_{2}F_{1}T_{1}),
\end{equation}  
```
where ``T=e^{-\tau H_{\text{loc}}}`` is time evolution operator, ``F`` is a fermion creation or annihilation operator, and ``k`` is the expansion order for the current diagrammatic configuration ``\mathcal{C}``. The straightforward method to evaluate this trace is to insert the complete eigenstates ``\{ \Gamma \}`` of ``H_{\text{loc}}`` into the RHS of the above equation, then 
```math
\begin{align}
\text{Tr}_{\text{loc}} &= &\sum_{\{\Gamma_{1}\Gamma_{2} \cdots \Gamma_{2k}\}} 
            \langle\Gamma_{1}|T_{2k+1}|\Gamma_{1}\rangle
            \langle\Gamma_{1}|F_{2k}|\Gamma_{2k}\rangle
            \langle\Gamma_{2k}|T_{2k}|\Gamma_{2k}\rangle \cdots \\  
          & &\times \langle\Gamma_{3}|F_{2}|\Gamma_{2}\rangle
            \langle\Gamma_{2}|T_{2}|\Gamma_{2}\rangle
            \langle\Gamma_{2}|F_{1}|\Gamma_{1}\rangle
            \langle\Gamma_{1}|T_{1}|\Gamma_{1}\rangle.
\end{align}
```
Thus, we must do ``4k+1`` matrix-matrix multiplications with the dimension of the Hilbert space of ``H_{\text{loc}}``. This method is robust but very slow for large multi-orbital impurity model as the dimension of the matrix is impractically large for 5- and 7-band systems, and the expansion order ``k`` is large as well.

Actually, the matrices of the fermion operators (``F``-matrix) are very sparse due to the symmetry of ``H_{\text{loc}}``. We can take advantage of this to speed up the matrix-matrix multiplications. We consider the symmetry of ``H_{\text{loc}}`` to find some good quantum numbers (GQNs) and divide the full Hilbert space of ``H_{\text{loc}}`` with very large dimension into much smaller subspaces labeled by these GQNs[^1]. We call a subspace ``|\alpha\rangle`` as a superstate[^2] which consists of all the ``n_{\alpha}`` eigenstates of this subspace, ``|\alpha\rangle=\{ \Gamma_{1}, \Gamma_{2}, \cdots, \Gamma_{n_{\alpha}}\}``. The ``F``-matrix element can only be nonzero between pairs of superstates with different values of GQNs. One fermion operator may bring one initial superstate ``|\alpha\rangle`` to some other final superstates ``|\beta\rangle``,
```math
\begin{equation}
F|\alpha\rangle= |\beta\rangle,
\end{equation}
```
or outside of the full Hilbert space. We have to carefully choose some GQNs to make sure that for a fixed initial superstate ``|\alpha\rangle`` and a fixed fermion operator, there is one and only one final superstate ``|\beta\rangle`` if it doesn't go outside of the full Hilbert space. Given an arbitrary diagrammatic configuration, starting with a superstate ``|\alpha_{1}\rangle``, there will be only one possible evolution path. That is,
```math
\begin{equation}
|\alpha_{1}\rangle \xrightarrow{F_{1}} 
|\alpha_{2}\rangle \xrightarrow{F_{2}} 
|\alpha_{3}\rangle \xrightarrow{F_{3}} 
|\alpha_{4}\rangle \cdots 
|\alpha_{2k-1}\rangle \xrightarrow{F_{2k-1}} 
|\alpha_{2k}\rangle   \xrightarrow{F_{2k}} 
|\alpha_{1}\rangle.
\end{equation}
```
The path may break at some point because it goes outside of the full Hilbert space or violates the Pauli principle. For a successful path starting with ``|\alpha_{1}\rangle``, its contribution to the local trace is
```math
\begin{align}
\text{Tr}_{\alpha_{1}} &=& 
\sum_{\{ \Gamma_{\alpha_{1}}, \Gamma_{\alpha_{2}}, \cdots, \Gamma_{\alpha_{2k}}\}}
\langle\Gamma_{\alpha_{1}}|T_{2k+1}|\Gamma_{\alpha_{1}}\rangle
\langle\Gamma_{\alpha_{1}}|F_{2k}|\Gamma_{\alpha_{2k}}\rangle
\langle\Gamma_{\alpha_{2k}}|T_{2k}|\Gamma_{\alpha_{2k}}\rangle \cdots \\ 
& & \times 
\langle\Gamma_{\alpha_{3}}|F_{2}|\Gamma_{\alpha_{2}}\rangle
\langle\Gamma_{\alpha_{2}}|T_{2}|\Gamma_{\alpha_{2}}\rangle
\langle\Gamma_{\alpha_{2}}|F_{1}|\Gamma_{\alpha_{1}}\rangle
\langle\Gamma_{\alpha_{1}}|T_{1}|\Gamma_{\alpha_{1}}\rangle,
\end{align}
```
where ``\{ \Gamma_{\alpha_{i}} \}`` are the eigenstates of subspace ``\alpha_{i}``. Thus, the final local trace should be
```math
\begin{equation}
\text{Tr}_{\text{loc}} = \sum_{i} \text{Tr}_{\alpha_{i}}.
\end{equation}
```
As a result, the original ``4k+1`` matrix-matrix multiplications with large dimension reduces to several ``4k+1`` matrix-matrix multiplications with much smaller dimensions, resulting in a huge speedup.

| GQNs | Kanamori-``U`` | Slater-``U`` | SOC |
| -- | -- | -- | -- |
|``N``, ``S_{z}``          | Yes          | Yes        | No   |
|``N``, ``S_{z}``, PS      | Yes          | No         | No   | 
|``N``, ``J_{z}``          | Yes          | Yes        | Yes  |
|``N``                     | Yes          | Yes        | Yes  |

**Table** | The GQNs supports for various types of local Hamiltonians ``H_{\text{loc}}``.

In our codes, we implemented several GQNs schemes for different types of local Hamiltonians ``H_{\text{loc}}``, which is summarized in the above table. For ``H_{\text{loc}}`` without spin-orbit coupling (SOC), we have two choices: (1) with Slater parameterized Coulomb interaction matrix, we use the total occupation number ``N``, the ``z`` component of total spin ``S_{z}`` as GQNs; (2) with Kanamori parameterized Coulomb interaction matrix, besides ``N`` and ``S_{z}``, we can use another powerful GQN, the so-called PS number[^3]. It is defined as,
```math
\begin{equation}
\text{PS} = \sum_{\alpha=1}^{N_{\text{orb}}} 
             (n_{\alpha\uparrow}-n_{\alpha\downarrow})^2 \times 2^{\alpha},
\end{equation}
```
where ``\alpha`` is the orbital index, ``\{\uparrow, \downarrow\}`` is spin index, ``n_{\alpha\uparrow}`` and ``n_{\alpha\downarrow}`` are the orbital occupancy numbers. The PS number labels the occupation number basis with the same singly occupied orbitals. With its help, the dimensions of the subspaces become very small, such that we can treat 5-band Kanamori systems efficiently without any approximations. For ``H_{\text{loc}}`` with SOC, we can use the total occupancy number ``N`` and the ``z`` component of total angular momentum ``J_{z}`` as GQNs. We summarize the total number of subspaces, maximum and mean dimension of subspaces for different GQNs schemes and multi-orbital impurity models in the below table. Obviously, using these GQNs can largely reduce the dimension of the ``F``-matrix, and make accurate DMFT calculations for complex electronic systems (such as ``d``- and ``f``-electron materials) possible. 

|               | 2-band       | 3-band       | 5-band       | 7-band      |
| -- | -- | -- | -- | -- |
|GQNs           | ``N``/max/mean | ``N``/max/mean | ``N``/max/mean | ``N``/max/mean   | 
|``N``, ``S_{z}``     |  9/4/1.78    | 16/9/4.00    | 36/100/28.44 | 64/1225/256.00 |
|``N``, ``S_{z}``, PS |  14/2/1.14   | 44/3/1.45    | 352/10/2.91  | 2368/35/6.92   |
|``N``, ``J_{z}``     |  -           | 26/5/2.46    | 96/37/10.67  | 246/327/66.60  |
|``N``            |  5/6/3.20    | 7/20/9.14    | 11/252/93.09 | 15/3432/1092.27|

**Table** | The total number of subspaces ``N``, maximum and mean dimension of subspaces for different GQNs schemes and multi-orbital models.

**Reference**

[^1]: Emanuel Gull, Andrew J. Millis, Alexander I. Lichtenstein, Alexey N. Rubtsov, Matthias Troyer, and Philipp Werner, *Rev. Mod. Phys.* **83**, 349 (2011)

[^2]: Kristjan Haule, *Phys. Rev. B* **75**, 155113 (2007)

[^3]: Nicolaus Parragh, Alessandro Toschi, Karsten Held, and Giorgio Sangiovanni, *Phys. Rev. B* **86**, 155158