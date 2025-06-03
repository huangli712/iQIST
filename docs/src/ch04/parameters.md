# Summary

In this section, we will introduce all of the valid parameters which can be used in the *solver.ctqmc.in* (for **NARCISSUS** and **MANJUSHAKA** components) file. For more information, the users can also refer to the comments in corresponding *ctqmc\_control.f90* file.

### Classification of the parameters

* **Parameters for the dynamical mean-field theory engine**
    * [isscf](p_isscf.md)
    * [niter](p_niter.md)
    * [alpha](p_alpha.md)
* **Parameters for the strongly correlated models**
    * [issun]
    * [isspn](p_isspn.md)
    * [isscr]
    * [nband](p_nband.md)
    * [nspin](p_nspin.md)
    * [norbs](p_norbs.md)
    * [ncfgs](p_ncfgs.md)
    * [nzero]
    * [U](p_u.md)
    * [Uc](p_uc.md)
    * [Uv](p_uv.md)
    * [Jz](p_jz.md)
    * [Js](p_js.md)
    * [Jp](p_jp.md)
    * [lc]
    * [wc]
    * [mune](p_mune.md)
    * [beta](p_beta.md)
    * [part](p_part.md)
* **Parameters for the quantum impurity solvers**
    * [ifast]
    * [itrun]
    * [mkink](p_mkink.md)
    * [mfreq](p_mfreq.md)
    * [ntime](p_ntime.md)
    * [npart]
    * [nflip](p_nflip.md)
    * [ntherm](p_ntherm.md)
    * [nsweep](p_nsweep.md)
    * [nwrite](p_nwrite.md)
    * [nclean](p_nclean.md)
    * [ncarlo](p_ncarlo.md)
    * [nmonte](p_nmonte.md)
* **Parameters for the Monte Carlo sampling**
    * [isort]
    * [issus]
    * [isvrt]
    * [lemax]
    * [legrd]
    * [nffrq](p_nffrq.md)
    * [nbfrq](p_nbfrq.md)
    * [nfreq](p_nfreq.md)

---

### Concerning the speed

The following parameters have big influences on the computational efficiency of the quantum impurity solvers.

* [ifast]
* [itrun]
* [npart]

---

### Concerning the accuracy

The following parameters will affect the computational accuracy.

* [ifast]
* [itrun]
* [ntime](p_ntime.md)
* [ntherm](p_ntherm.md)
* [nsweep](p_nsweep.md)
* [nclean](p_nclean.md)
* [ncarlo](p_ncarlo.md)
* [nmonte](p_nmonte.md)
* [isort]
* [isvrt]
* [lemax]
* [legrd]

---

!!! tip

    1. All of the parameters are case-insensitive.
    2. All of the parameters have default values. You can override them via the *solver.ctqmc.in* and *solver.hfqmc.in* files.
    3. The quantum impurity solvers won't check the correctness, rationality of the parameters. They won't adjust the parameters automatically. For example, if you setup *nband = 2* in the *solver.ctqmc.in*, you have to setup *norbs*, *ncfgs* parameters by yourself at the same time. The quantum impurity solvers won't do that.

**See also**:

* [solver.ctqmc.in](in_ctqmc.md) // File format for *solver.ctqmc.in*
