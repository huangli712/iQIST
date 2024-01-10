# Principles of continuous-time quantum Monte Carlo algorithm

Now we already have an impurity model Hamiltonian ``H_{\text{imp}}``, the question is how to solve it using the Monte Carlo algorithm?

We first split the impurity Hamiltonian ``H_{\text{imp}}`` into two separate parts,

```math
\begin{equation}
H_{\text{imp}} = H_1 + H_2,
\end{equation}
```

then treat ``H_2`` as a perturbation term, and expand the partition function ``\mathcal{Z}`` in powers of ``H_2``,
```math
\begin{equation}
\mathcal{Z} =  \text{Tr} e^{-\beta H} = \sum_{n=0}^{\infty} \int_{0}^{\beta} \cdots \int_{\tau_{n-1}}^\beta \omega(\mathcal{C}_n),
\end{equation}
```
with
```math
\begin{equation}
\omega(\mathcal{C}_n)=d\tau_1 \cdots d\tau_n \text{Tr}\left\{ e^{-\beta H_1}[-H_2(\tau_n)]\cdots [-H_2(\tau_1)]\right\},
\end{equation}
```
where ``H_2(\tau)`` is defined in the interaction picture with

```math
\begin{equation}
H_2(\tau) = e^{\tau H_1} H_2 e^{-\tau H_1}.
\end{equation}
```

Each term in the right side of the second equation can be regarded as a diagram or configuration (labeled by ``\mathcal{C}``), and ``\omega(\mathcal{C}_n)`` is the diagrammatic weight of a specific order-``n`` configuration. Next we use a stochastic Monte Carlo algorithm to sample the terms of this series. This is the core spirit of the continuous-time quantum Monte Carlo impurity solver. The idea is very simple, but the realization is not.

Depending on the different choices of ``H_{2}`` term, there are multiple variations for the continuous-time quantum Monte Carlo impurity solver. According to our knowledge, the variations at least include

* CT-INT
* CT-HYB
* CT-J
* CT-AUX

In the CT-INT and CT-AUX quantum impurity solvers[^1][^2], the interaction term is the perturbation term, namely, ``H_2 = H_{\text{int}}``, while ``H_2 = H_{\text{hyb}}`` is chosen for the CT-HYB quantum impurity solver[^3]. The CT-J quantum impurity solver is designed for the Kondo lattice model only[^4]. We won't discuss it at here. In the intermediate and strong interaction region, CT-HYB is much more efficient than CT-INT and CT-AUX. We could even say that it is the most powerful and efficient quantum impurity solver so far. This is also the main reason that we only implemented the CT-HYB quantum impurity solvers in the ``i``QIST software package.

**Reference**

[^1]: A. N. Rubtsov, V. V. Savkin, and A. I. Lichtenstein, *Phys. Rev. B* **72**, 035122 (2005)

[^2]: Emanuel Gull, Philipp Werner, Olivier Parcollet, Matthias Troyer, *EPL* **82**, 57003 (2008)

[^3]: Philipp Werner, Armin Comanac, Luca de’ Medici, Matthias Troyer, and Andrew J. Millis, *Phys. Rev. Lett.* **97**, 076405 (2006)

[^4]: Junya Otsuki, Hiroaki Kusunose, Philipp Werner, and Yoshio Kuramoto, *J. Phys. Soc. Jpn.* **76**, 114707 (2007)
