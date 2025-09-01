# Summary

The quantum impurity solvers will generate a lot of data files at run time. In this section, we will depict their formats and usages in detail.

The data files could be classified as the following kinds:

* **Direct Output**
    * [Terminal output](out_term.md) // Many useful stuffs in it
* **Green's Function Data**
    * [solver.green.dat](out_green.md) // ``G(\tau)``
    * [solver.weiss.dat](out_weiss.md) // ``G_0(\tau)``
    * [solver.grn.dat](out_grn.md) // ``G(i\omega_n)``
    * [solver.wss.dat](out_wss.md) // ``G_0(i\omega_n)``
* **Auxiliary Green's Function Data**
    * [solver.fcorr.dat]
    * [solver.frn.dat]
* **Hybridization Function Data**
    * [solver.hybri.dat](out_hybri.md) // ``\Delta(\tau)``
    * [solver.hyb.dat](out_hyb.md) // ``\Delta(i\omega_n)``
* **Self-energy Function Data**
    * [solver.sgm.dat](out_sgm.md) // ``\Sigma(i\omega_n)``
* **Two-Particle Green's Function Data**
    * [solver.g2ph.dat]
    * [solver.g2pp.dat]
    * [solver.h2ph.dat]
    * [solver.h2pp.dat]
    * [solver.v4ph.dat]
    * [solver.v4pp.dat]
* **Susceptibility Data**
    * [solver.sp_t.dat](out_sp_t.md) // Spin-spin correlation function in ``\tau`` space
    * [solver.sp_w.dat](out_sp_w.md) // Spin-spin correlation function in Matsubara frequency space
    * [solver.ch_t.dat](out_ch_t.md) // Orbital-orbital correlation in ``\tau`` space
    * [solver.ch_w.dat](out_ch_w.md) // Orbital-orbital correlation in Matsubara frequency space
* **Physical Observables Data**
    * [solver.hist.dat](out_hist.md) // Histogram for perturbation expansion order
    * [solver.prob.dat](out_prob.md) // Atomic state probability
    * [solver.paux.dat]
    * [solver.szpw.dat]
    * [solver.nmat.dat](out_nmat.md) // Occupation number
    * [solver.kmat.dat](out_kmat.md) // ``\langle k \rangle`` and ``\langle k^2 \rangle``
    * [solver.lrmm.dat](out_lrmm.md) // ``\langle k_L \rangle`` and ``\langle k_R \rangle``
    * [solver.ac_f.dat]
* **Miscellaneous Data**
    * [solver.kernel.dat](out_kern.md) // Frequency-dependent interaction
    * [solver.status.dat](out_stat.md) // Status of the quantum impurity solvers
    * [solver.diag.dat](out_diag.md) // Snapshot of the current configuration for diagrammatic perturbation expansion
