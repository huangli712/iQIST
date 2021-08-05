## Standard output files

The quantum impurity solvers will generate a lot of data files at run time. In this section, we will depict their formats and usages in detail.

The data files could be classified as the following kinds:

* **Direct output**
    * [Terminal output](out_term.md) // Many useful stuffs in it
* **Green's function data**
    * [solver.green.dat](out_green.md) // ``G(\tau)``
    * [solver.weiss.dat](out_weiss.md) // ``G_0(\tau)``
    * [solver.grn.dat](out_grn.md) // ``G(i\omega_n)``
    * [solver.wss.dat](out_wss.md) // ``G_0(i\omega_n)``
* **Hybridization function data**
    * [solver.hybri.dat](out_hybri.md) // ``\Delta(\tau)``
    * [solver.hyb.dat](out_hyb.md) // ``\Delta(i\omega_n)``
* **Self-energy function data**
    * [solver.sgm.dat](out_sgm.md) // ``\Sigma(i\omega_n)``
    * [solver.hub.dat](out_hub.md) // ``\Sigma_{\text{atomic}}(i\omega_n)`` and ``G_{\text{atomic}}(i\omega_n)``
* **Two-particle Green's function data**
    * [solver.twop.dat](out_twop.md) // Two-particle Green's function and vertex function
    * [solver.vrtx.dat](out_vrtx.md) // Two-particle Green's function and vertex function
    * [solver.pair.dat](out_pair.md) // Pairing susceptibility
* **Susceptibility data**
    * [solver.schi.dat](out_schi.md) // Spin-spin correlation function in ``\tau`` space
    * [solver.sfom.dat](out_sfom.md) // Spin-spin correlation function in Matsubara frequency space
    * [solver.ochi.dat](out_ochi.md) // Orbital-orbital correlation in ``\tau`` space
    * [solver.ofom.dat](out_ofom.md) // Orbital-orbital correlation in Matsubara frequency space
* **Physical observables data**
    * [solver.hist.dat](out_hist.md) // Histogram for perturbation expansion order
    * [solver.prob.dat](out_prob.md) // Atomic state probability
    * [solver.nmat.dat](out_nmat.md) // Occupation number
    * [solver.kmat.dat](out_kmat.md) // ``\langle k \rangle`` and ``\langle k^2 \rangle``
    * [solver.lmat.dat](out_lmat.md) // ``\langle k_L \rangle`` and ``\langle k_R \rangle``
* **Miscellaneous data**
    * [solver.kernel.dat](out_kern.md) // Frequency-dependent interaction
    * [solver.status.dat](out_stat.md) // Status of the quantum impurity solvers
    * [solver.diag.dat](out_diag.md) // Snapshot of the current configuration for diagrammatic perturbation expansion