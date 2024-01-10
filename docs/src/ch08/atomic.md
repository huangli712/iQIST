# Atomic eigenvalue solver

The **JASMINE** component of iQIST software package is used to solve a local atomic Hamiltonian defined as,
```math
\begin{equation}
\hat{H}_{\text{atom}}=\hat{H}_{0}+\hat{H}_{\text{Coulomb}},
\end{equation}
```
where, ``\hat{H}_{0}`` is the on-site term, and ``\hat{H}_{\text{Coulomb}}`` is the Coulomb interaction term.

The ``\hat{H}_{0}`` term includes the crystal field (CF) splitting term ``\hat{H}_{\text{CF}}``, spin-orbit coupling (SOC) term ``\hat{H}_{\text{SOC}}``, and other terms such as a Zeeman term in case of the presence of external magnetic field,
```math
\begin{equation}
\hat{H}_{0}=\hat{H}_{\text{CF}}+\hat{H}_{\text{SOC}}+\text{other on-site terms},
\end{equation}
```

We can re-write the ``\hat{H}_{\text{atom}}`` in second quantization form,
```math
\begin{equation}
\hat{H}_{\text{atom}}=\sum_{\alpha\beta}E_{\alpha\beta}\hat{f}_{\alpha}^{\dagger}\hat{f}_{\beta}
                     +\sum_{\alpha\beta\gamma\delta}U_{\alpha\beta\gamma\delta}\hat{f}_{\alpha}^{\dagger}
                     \hat{f}_{\beta}^{\dagger}\hat{f}_{\delta}\hat{f}_{\gamma},
\end{equation}
```
where, ``\alpha,\beta,\gamma,\delta`` is the single particle orbital-spin index, the first term is the on-site term, the second one is the Coulomb interaction term.
