# Improved estimator for the self-energy function

Recently, Hafermann *et al.* proposed efficient measurement procedures for the self-energy and vertex functions within the CT-HYB algorithm[^1,2]. In their method, some higher-order correlation functions (related to the quantities being sought through the equation of motion) are measured. For the case of interactions of density-density type, the segment algorithm is available[^3]. Thus, the additional correlators can be obtained essentially at no additional computational cost. When the calculations are completed, the required self-energy function and vertex function can be evaluated analytically.

The improved estimator for the self-energy function can be expressed in the following form:
```math
\begin{equation}
\Sigma_{ab}(i\omega_n) = \frac{1}{2}
\sum_{ij} G^{-1}_{ai}(i\omega_n) (U_{jb} + U_{bj}) F^{j}_{ib}(i\omega_n),
\end{equation}
```
where ``U_{ab}`` is the Coulomb interaction matrix element. The expression for the new two-particle correlator ``F^{j}_{ab}(\tau - \tau')`` reads
```math
\begin{equation}
F^{j}_{ab}(\tau-\tau')
= -\langle \mathcal{T} d_{a}(\tau) d^{\dagger}_{b}(\tau') n_{j}(\tau') \rangle,
\end{equation}
```
and ``F^{j}_{ab}(i\omega_n)`` is its Fourier transform. The actual measurement formula is
```math
\begin{equation}
F^{j}_{ab}(\tau - \tau') =
-\frac{1}{\beta}
\left\langle
\sum_{\alpha\beta = 1}^{k}
\mathcal{M}_{\beta\alpha}\delta^{-}(\tau-\tau', \tau^{e}_{\alpha} - \tau^{s}_{\beta})
n_{j}(\tau^s_\beta)\delta_{a,\alpha}\delta_{b,\beta}
\right\rangle.
\end{equation}
```
As one can see, this equation looks quite similar to the one for imaginary-time Green's function. Thus we use the same method to measure ``F^{j}_{ab}(\tau - \tau')`` and finally get the self-energy function via the first equation. Here, the matrix element ``n_{j}(\tau^s_\beta)`` (one or zero) denotes whether or not flavor ``j`` is occupied (whether or not a segment is present) at time ``\tau^s_\beta``.

This method can be combined with the orthogonal polynomial representation[^4] as introduced in the previous subsection to suppress fluctuations and filter out the Monte Carlo noise. Using this technique, we can obtain the self-energy and vertex functions with unprecedented accuracy, which leads to an enhanced stability in the analytical continuations of those quantities[^2]. In the iQIST software package, we only implemented the improved estimator for the self-energy function. Note that when the interaction matrix is frequency-dependent, the first equation should be modified slightly[^1].

**Reference**

[^1]: Hartmut Hafermann, *Phys. Rev. B* **89**, 235128 (2014)

[^2]: Hartmut Hafermann, Kelly R. Patton, and Philipp Werner, *Phys. Rev. B* **85**, 205106 (2012)

[^3]: Philipp Werner, Armin Comanac, Luca de’ Medici, Matthias Troyer, and Andrew J. Millis, *Phys. Rev. Lett.* **97**, 076405 (2006)

[^4]: Lewin Boehnke, Hartmut Hafermann, Michel Ferrero, Frank Lechermann, and Olivier Parcollet, *Phys. Rev. B* **84**, 075145 (2011)
