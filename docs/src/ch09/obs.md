### Physical observable

Many physical observables are measured in our CT-HYB quantum impurity solvers. Here we provide a list of them.

**Single-particle Green's function $$G(\tau)$$**

The most important observable is the single-particle Green's function $$G(\tau)$$, which is measured using the elements of the matrix $$\mathcal{M}$$, 
$$
\begin{align}
G(\tau) = \left\langle \frac{1}{\beta} \sum_{ij}\delta^{-}(\tau, \tau_i - \tau_j) \mathcal{M}_{ji}\right\rangle,
\end{align}
$$
with
$$
\begin{align}
\delta^{-}(\tau, \tau') = 
\begin{cases} 
\delta(\tau - \tau'), & \tau' > 0, \\
-\delta(\tau - \tau' + \beta), & \tau' < 0.
\end{cases}
\end{align}
$$

**Single-particle Green's function $$G(i\omega_n)$$**

Note that in the iQIST software package, the Matsubara Green's function $$G(i\omega_n)$$ is also measured directly, instead of being calculated from $$G(\tau)$$ using Fourier transformation.

$$
G(i\omega_n) = -\frac{1}{\beta} \sum_{ij} e^{i\omega_n\tau_i}\mathcal{M}_{ij}e^{-i\omega_n\tau_j}
$$

**Two-particle correlation function $$\chi_{\alpha\beta}(\omega,\omega',\nu)$$**

The two-particle correlation functions are often used to construct lattice susceptibilities within DMFT and diagrammatic extensions of DMFT. However, the measurements of two-particle correlation functions are a nontrivial task[^1] as it is very time-consuming to obtain good quality data, and most of the previous publications in this field are restricted to measurements of two-particle correlation functions in one-band models. Thanks to the development of efficient CT-HYB algorithms, the calculation of two-particle correlation functions for multi-orbital impurity models now become affordable[^2,3,4]. In the iQIST software package, we implemented the measurement for the two-particle correlation function $$\chi_{\alpha\beta}(\tau_a,\tau_b,\tau_c,\tau_d)$$, which is defined as follows:
$$
\begin{equation}
\chi_{\alpha\beta}(\tau_a,\tau_b,\tau_c,\tau_d)
= \langle c_{\alpha}(\tau_a)c^{\dagger}_{\alpha}(\tau_b)c_{\beta}(\tau_c)c^{\dagger}_{\beta}(\tau_d)\rangle.
\end{equation}
$$
Due to the memory restrictions, the actual measurement is performed in the frequency space, for which we use the following definition of the Fourier transform:
$$
\begin{align}
\chi_{\alpha\beta}(\omega,\omega',\nu) &= \frac{1}{\beta}
\int^{\beta}_{0}d\tau_a\int^{\beta}_{0}d\tau_b\int^{\beta}_{0}d\tau_c\int^{\beta}_{0}d\tau_d
\nonumber\\
&\times \chi_{\alpha\beta}(\tau_a,\tau_b,\tau_c,\tau_d) 
e^{i(\omega+\nu)\tau_a}e^{-i\omega\tau_b}e^{-i\omega'\tau_c}e^{-i(\omega'+\nu)\tau_d}.
\end{align}
$$
where $$\omega$$ and $$\omega'$$ [$$\equiv (2n+1)\pi\beta$$] are fermionic frequencies, and $$\nu$$ is bosonic ($$\equiv 2n\pi/\beta$$).

**Local irreducible vertex functions $$\Gamma_{\alpha\beta}(\omega,\omega',\nu)$$**

From the two-particle Green's function $$\chi_{\alpha\beta}(\omega,\omega',\nu)$$, the local irreducible vertex function $$\Gamma_{\alpha\beta}(\omega,\omega',\nu)$$ can be calculated easily, via the Bethe-Salpeter equation[^3,4,5]:
$$
\begin{equation}
\Gamma_{\alpha\beta}(\omega,\omega',\nu) = 
\frac{\chi_{\alpha\beta}(\omega,\omega',\nu) 
- \beta[G_\alpha(\omega+\nu)G_\beta(\omega')\delta_{\nu,0} 
- G_\alpha(\omega+\nu) G_\beta(\omega') \delta_{\alpha\beta}\delta_{\omega\omega'}]}
{G_\alpha(\omega+\nu)G_\alpha(\omega)G_\beta(\omega')G_\beta(\omega'+\nu)}.
\end{equation}
$$
The $$G(i\omega_n)$$ and $$\Gamma_{\alpha\beta}(\omega,\omega',\nu)$$ are essential inputs for the ladder dual fermion code **ROSEMARY**, see section [Ladder dual fermions](../ch05/ladder.md) for more details.

**Pair susceptibility**

**Impurity self-energy function $$\Sigma(i\omega_n)$$**

The self-energy $$\Sigma(i\omega_n)$$ is calculated using Dyson's equation directly
$$
\begin{equation}
\Sigma(i\omega_n) = G^{-1}_{0}(i\omega_n) - G^{-1}(i\omega_n),
\end{equation}
$$
 or measured using the so-called improved estimator[^3,4]. Noted that now the latter approach only works when the segment representation is used. 

**Histogram of the perturbation expansion order**

We record the histogram of the perturbation expansion order $$k$$, which can be used to evaluate the kinetic energy.

**Kurtosis and skewness of perturbation expansion order**

Skewness

$$
\gamma_1 = \frac{E[(k - \langle k \rangle)^3]}{(E[(k - \langle k \rangle)^2])^{3/2}}
$$

Kurtosis

$$
\frac{E[(k - \langle k \rangle)^4]}{(E[(k - \langle k \rangle)^2])^{2}}
$$

Actually, in the iQIST software package, only the $$\langle k \rangle$$, $$\langle k^2 \rangle$$, $$\langle k^3 \rangle$$, and $$\langle k^4 \rangle$$ are measured. And then they are used to evaluate the skewness and kurtosis.

**Kinetic energy**

The expression for the system kinetic energy reads
$$
E_{\text{kin}} = -\frac{1}{\beta} \langle k \rangle,
$$
where $$k$$ is the perturbation expansion order.

**Potential energy**

**Occupation number and double occupation number**

The orbital occupation number $$\langle n_\alpha\rangle$$ and double occupation number $$\langle n_\alpha n_\beta \rangle$$ are measured. From them we can calculate for example the charge fluctuation $$\sqrt{\langle N^2 \rangle - \langle N \rangle^2}$$, where $$N$$ is the total occupation number:
$$
\begin{equation}
N = \sum_{\alpha} n_{\alpha}.
\end{equation}
$$

**Magnetic moment**

Actually, we only measure $$\langle S_{z} \rangle$$.

**Spin-spin correlation function**

For a system with spin rotational symmetry, the expression for the spin-spin correlation function reads
$$
\begin{equation}
\chi_{ss}(\tau) = \langle S_{z}(\tau) S_{z}(0) \rangle,
\end{equation}
$$
where $$S_{z} = n_{\uparrow} - n_{\downarrow}$$. From it we can calculate the effective magnetic moment:
$$
\begin{equation}
\mu_{\text{eff}} = \int^{\beta}_{0}d\tau \chi_{ss}(\tau).
\end{equation}
$$

**Orbital-orbital correlation function**

The expression for the orbital-orbital correlation function reads 
$$
\begin{equation}
\chi^{nn}_{\alpha\beta}(\tau) = \langle n_{\alpha}(\tau) n_{\beta}(0) \rangle.
\end{equation}
$$

**Atomic state probability**

The expression for the atomic state probability is 
$$
\begin{equation}
p_{\Gamma} = \langle |\Gamma \rangle \langle \Gamma| \rangle,
\end{equation}
$$
where $$\Gamma$$ is the atomic state.

**Fidelity susceptibility**

$$
\chi_{\text{FS}} = \langle k_{L} k_{R} \rangle - \langle k_{L} \rangle \langle k_{R} \rangle
$$

**Kinetic energy fluctuation**

$$
\chi_{\text{k}} = \langle k^2 \rangle - \langle k \rangle^2 - \langle k \rangle
$$

**Reference**

[^1]: Jan Kuneš, *Phys. Rev. B* **83**, 085102 (2011)

[^2]: Lewin Boehnke, Hartmut Hafermann, Michel Ferrero, Frank Lechermann, and Olivier Parcollet, *Phys. Rev. B* **84**, 075145 (2011)

[^3]: Hartmut Hafermann, *Phys. Rev. B* **89**, 235128 (2014)

[^4]: Hartmut Hafermann, Kelly R. Patton, and Philipp Werner, *Phys. Rev. B* **85**, 205106 (2012)

[^5]: G. Rohringer, A. Valli, and A. Toschi, *Phys. Rev. B* **86**, 125114 (2012)