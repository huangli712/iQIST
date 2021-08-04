## Parameters

In this section, we will introduce all of the valid parameters which can be used in the *solver.ctqmc.in* (for **AZALEA**, **GARDENIA**, **NARCISSUS**, **BEGONIA**, **LAVENDER**, **CAMELLIA**, **PANSY**, and **MANJUSHAKA** components) and *solver.hfqmc.in* (for **DAISY** component) files. For more information, the users can also refer to the comments in corresponding *ctqmc\_control.f90* file.

### Classification of the parameters

* ** Parameters for the dynamical mean-field theory engine **
    * [isscf](p_isscf.md)
    * [isbin](p_isbin.md)
    * [niter](p_niter.md)
    * [alpha](p_alpha.md)
* ** Parameters for the strongly correlated models **
    * [issun](p_issun.md)
    * [isspn](p_isspn.md)
    * [isscr](p_isscr.md)
    * [nband](p_nband.md)
    * [nspin](p_nspin.md)
    * [norbs](p_norbs.md)
    * [ncfgs](p_ncfgs.md)
    * [nzero](p_nzero.md)
    * [U](p_u.md)
    * [Uc](p_uc.md)
    * [Uv](p_uv.md)
    * [Jz](p_jz.md)
    * [Js](p_js.md)
    * [Jp](p_jp.md)
    * [lc](p_lc.md)
    * [wc](p_wc.md)
    * [mune](p_mune.md)
    * [beta](p_beta.md)
    * [part](p_part.md)
* ** Parameters for the quantum impurity solvers **
    * [ifast](p_ifast.md)
    * [itrun](p_itrun.md)
    * [mkink](p_mkink.md)
    * [mstep](p_mstep.md)
    * [mfreq](p_mfreq.md)
    * [nsing](p_nsing.md)
    * [ntime](p_ntime.md)
    * [nvect](p_nvect.md)
    * [nleja](p_nleja.md)
    * [npart](p_npart.md)
    * [nflip](p_nflip.md)
    * [ntherm](p_ntherm.md)
    * [nsweep](p_nsweep.md)
    * [nwrite](p_nwrite.md)
    * [nclean](p_nclean.md)
    * [ncarlo](p_ncarlo.md)
    * [nmonte](p_nmonte.md)
* ** Parameters for the Monte Carlo sampling **
    * [isort](p_isort.md)
    * [issus](p_issus.md)
    * [isvrt](p_isvrt.md)
    * [lemax](p_lemax.md)
    * [legrd](p_legrd.md)
    * [chmax](p_chmax.md)
    * [chgrd](p_chgrd.md)
    * [nffrq](p_nffrq.md)
    * [nbfrq](p_nbfrq.md)
    * [nfreq](p_nfreq.md)

---

### Concerning the speed

The following parameters have big influences on the computational efficiency of the quantum impurity solvers.

* [ifast](p_ifast.md)
* [itrun](p_itrun.md)
* [mstep](p_mstep.md)
* [nvect](p_nvect.md)
* [nleja](p_nleja.md)
* [npart](p_npart.md)

---

### Concerning the accuracy

The following parameters will affect the computational accuracy.

* [ifast](p_ifast.md)
* [itrun](p_itrun.md)
* [ntime](p_ntime.md)
* [nvect](p_nvect.md)
* [ntherm](p_ntherm.md)
* [nsweep](p_nsweep.md)
* [nclean](p_nclean.md)
* [ncarlo](p_ncarlo.md)
* [nmonte](p_nmonte.md)
* [isort](p_isort.md)
* [isvrt](p_isvrt.md)
* [lemax](p_lemax.md)
* [legrd](p_legrd.md)
* [chmax](p_chmax.md)
* [chgrd](p_chgrd.md)

---

!!! tip

    1. All of the parameters are case-insensitive.
    2. All of the parameters have default values. You can override them via the *solver.ctqmc.in* and *solver.hfqmc.in* files.
    3. The quantum impurity solvers won't check the correctness, rationality of the parameters. They won't adjust the parameters automatically. For example, if you setup *nband = 2* in the *solver.ctqmc.in*, you have to setup *norbs*, *ncfgs* parameters by yourself at the same time. The quantum impurity solvers won't do that.

**See also**:

* [solver.ctqmc.in](in_ctqmc.md) // File format for *solver.ctqmc.in*
* [solver.hfqmc.in](in_hfqmc.md) // File format for *solver.hfqmc.in*