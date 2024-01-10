# Hybridization expansion

In the hybridization expansion algorithm (CT-HYB), due to the fact that ``H_1`` does not mix the impurity and bath states, the trace in ``\omega(\mathcal{C}_n)`` can be written as

```math
\begin{equation}
\text{Tr} = \text{Tr}_d \text{Tr}_c.
\end{equation}
```

As a result, we can split the weight of each configuration as

```math
\begin{equation}
\omega(\mathcal{C}_n) = \omega_{d}(\mathcal{C}_n) \omega_{c}(\mathcal{C}_n) \prod\limits_{i=1}^{n} d\tau_i.
\end{equation}
```

``\omega_{d}(\mathcal{C}_n)`` is the trace over impurity operators (``\text{Tr}_d``), ``\omega_{c}(\mathcal{C}_n)`` is the trace over bath operators (``\text{Tr}_c``). Further, since the Wick's theorem is applicable for the ``\omega_c(\mathcal{C}_n)`` part, we can represent it as a determinant of a matrix ``\mathcal{Z}_{\text{bath}}\mathcal{M}^{-1}`` with

```math
\begin{equation}
\mathcal{Z}_{\text{bath}}=\text{Tr}_c e^{-\beta H_{\text{bath}}}
\end{equation}
```

and

```math
\begin{equation}
(\mathcal{M}^{-1})_{ij} = \Delta(\tau_i - \tau_j).
\end{equation}
```

The ``\omega_{d}(\mathcal{C}_n)`` part can be expressed using segment representation when ``[n_{\alpha}, H_{\text{loc}}] = 0``[^1]. However, if this condition is not fulfilled, we have to calculate the trace explicitly, which is called the general matrix algorithm[^2][^3]. The explicit calculation of the trace for a large multi-orbital AIM with general interactions is computationally expensive. Many tricks and strategies have been implemented in the iQIST software package to address this challenge. Please refer to next section for more details.

In this package, we used importance sampling and the Metropolis algorithm to evaluate the partition function ``\mathcal{Z}``. The following four local update procedures, with which the ergodicity of Monte Carlo algorithm is guaranteed, are used to generate the Markov chain:

* Insert a pair of creation and annihilation operators in the time interval ``[0,\beta)``.
* Remove a pair of creation and annihilation operators from the current configuration.
* Select a creation operator randomly and shift its position on the imaginary time axis.
* Select a annihilation operator randomly and shift its position on the imaginary time axis.

In the Monte Carlo simulations, sometimes the system can be trapped by some (for example symmetry-broken) state. In order to avoid unphysical trapping, we also consider the following two global updates:
* Swap the operators of randomly selected spin up and spin down flavors.
* Swap the creation and annihilation operators globally.

**Reference**

[^1]: Philipp Werner, Armin Comanac, Luca de’ Medici, Matthias Troyer, and Andrew J. Millis, Phys. Rev. Lett. 97, 076405 (2006)

[^2]: Philipp Werner and Andrew J. Millis, *Phys. Rev. B* **74**, 155107 (2006)

[^3]: Kristjan Haule, *Phys. Rev. B* **75**, 155113 (2007)
