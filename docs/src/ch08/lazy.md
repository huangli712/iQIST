# Lazy trace evaluation

The diagrammatic Monte Carlo sampling algorithm consists of the following steps: (1) Propose an update for the current diagrammatic configuration. (2) Calculate the acceptance probability ``p`` according to the Metropolis-Hasting algorithm,
```math
\begin{equation}
p = \text{min} (1, \frac{A^\prime}{A} \left| \frac{\omega_{c}}{\omega_{c}^{\prime}}\right|
     \left|\frac{\omega_{d}}{\omega_{d}^{\prime}} \right|),
\end{equation}
```
where, ``A`` is the proposal probability for the current update and ``A^\prime`` for the inverse update, ``\omega_{c}`` and ``\omega_{c}^{\prime}`` are the determinants for the new and old configurations, respectively, and ``\omega_{d}`` and ``\omega_{d}^{\prime}`` are the local traces for the new and old configurations, respectively. (3) Generate a random number ``r``. If ``p>r``, the proposed update is accepted, otherwise it is rejected. (4) Update the current diagrammatic configuration if the proposed update is accepted. It turns out that for CT-HYB, ``p`` is usually low (``1\% \sim 20\%``), especially in the low temperature region. On the other hand, the calculation of ``p`` involves a costly local trace evaluation. To avoid wasting computation time when the acceptance probability is very low, in the subspace algorithm, we implemented the so-called lazy trace evaluation proposed in Reference[^1].

The basic idea of the lazy trace evaluation is simple. For the proposed Monte Carlo move, we first generate a random number ``r``. Then, instead of calculating the local trace from scratch to determine ``p``, we calculate bounds for ``\left|\text{Tr}_{\text{loc}}\right|``,
```math
\begin{equation}
\left|\omega_{d}\right| = \left|\text{Tr}_{\text{loc}}\right| \leq \sum_{i} \left|\text{Tr}_{i}\right| \leq \sum_{i} B_{i},
\end{equation}
```
where ``B_i \geq \left|\text{Tr}_{i}\right|``. ``B_{i}`` is a product of some chosen matrix norms of ``T`` and ``F`` matrices:
```math
\begin{equation}
B_i = C  \left\| T_{2k+1}
\right\| \left\| F_{2k} \right\| \cdots \left\| F_{2} \right\| \left\| T_{2} \right\|
\left\| F_{1} \right\| \left\| T_{1} \right\| \geq
\left|\text{Tr}(T_{2k+1}F_{2k} \cdots F_{2}T_{2}F_{1}T_{1})\right|,
\end{equation}
```
where ``C`` is a parameter depending on the specific type of matrix norm, and ``\left\| \cdot \right\|`` denotes a matrix norm. If ``\text{Tr}_{i^\prime}`` denotes the exact trace of some subspaces, then we have
```math
\begin{equation}
\left| \left|\text{Tr}_{\text{loc}}\right| - \sum_{i^\prime}\left|\text{Tr}_{i^\prime}\right| \right|
\leq \sum_{i \neq i^\prime} B_{i}.
\end{equation}
```
Thus, we can determine the upper ``p_{\text{max}}`` and lower ``p_{\text{min}}`` bounds of ``p`` as
```math
\begin{equation}
\begin{aligned}
p_{\text{max}}=R \left(\sum_{i^\prime}\left|\text{Tr}_{i^\prime}\right| + \sum_{i \neq i^\prime} B_{i}\right),\\
p_{\text{min}}=R \left(\sum_{i^\prime}\left|\text{Tr}_{i^\prime}\right| - \sum_{i \neq i^\prime} B_{i}\right),
\end{aligned}
\end{equation}
```
where ``R=\frac{A^\prime}{A} \left| \frac{\omega_{c}}{\omega_{c}^{\prime}}\right| \left|\frac{1}{\omega_{d}^{\prime}} \right|``.

If ``r>p_{\text{max}}``, we reject this move immediately. If ``r<p_{\text{min}}``, we accept the move and calculate the determinant and local trace from scratch. If `` p_{\text{min}} < r < p_{\text{max}} ``, we refine the bounds by calculating the local trace of one more subspace ``\text{Tr}_{i}`` until we can reject or accept the move. The calculation of these bounds involves only simple linear algebra calculations of matrix norms which cost little computation time, and one refining operation involves only one subspace trace evaluation. On average, it saves a lot of computation time, as confirmed by our benchmarks.

**Reference**

[^1]: P. Sémon, Chuck-Hou Yee, Kristjan Haule, and A.-M. S. Tremblay, Phys. Rev. B 90, 075149
