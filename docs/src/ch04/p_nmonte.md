# Parameter: nmonte

**Definition**

> How often to measure the physical observables.

**Type**

> Integer

**Default value**

> 10

**Component**

> ALL

**Behavior**

> Every *nmonte* Monte Carlo sampling steps, the quantum impurity solvers try to measure the physical observables. The affected physical observables are as follows:
>
> * Impurity occupancy, ``\langle n_{\alpha} \rangle``, ``\langle n_{\alpha}n_{\beta} \rangle``,
> * Impurity Matsubara Green's function, ``G(i\omega_n)``,
> * Spin-spin correlation function, ``\langle S_z(0) S_z(\tau) \rangle``, ``\chi_{\text{sp}}(i\nu_n)``,
> * Orbital-orbital correlation function, ``\langle n_{\alpha}(0) n_{\beta}(\tau)\rangle``, ``\chi_{\text{ch}}(i\nu_n)``,
> * Momentum of kinetic energy, ``\langle k^2 \rangle - \langle k \rangle^2``,
> * Fidelity susceptibility, ``\langle k_L k_R \rangle - \langle k_L \rangle \langle k_R \rangle``,
> * Two-particle Green's function, ``\chi(i\omega_n,i\omega_n,i\nu_n)``,
> * Two-particle vertex function, ``\mathcal{F}(i\omega_n,i\omega_n,i\nu_n)``,
> * Pairing susceptibility, ``\Gamma(i\omega_n,i\omega_n,i\nu_n)``,.
>
> The measuring period for the other physical observables are controlled by the *ncarlo* parameter. See also [ncarlo](p_ncarlo.md) parameter for more details.

!!! note

    The histogram for the perturbation expansion series is measured in each Monte Carlo sampling step.

**Comment**

> You can increase *nmonte* to relieve the auto-correlation of physical observables between two successive measurements. However, large *nmonte* will waste the CPU times. It is recommended that a good *nmonte* should satisfy the following relation:
>
> ```math
> \text{nmonte} * \text{P}_{\text{insert}} >= 1.0
> ```
>
> Here, ``\text{P}_{\text{insert}}`` means the accepted ratio for the insert update actions.
