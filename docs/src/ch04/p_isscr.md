### Parameter: isscr

**Definition**

It is a key control flag, which specifies what kind of impurity Hamiltonian model need to be solved. Usually, the quantum impurity solvers in the iQIST software package only support the Hubbard model/Anderson impurity model. But with the help of the **NARCISSUS** component, we have several other choices.

**Type**

Integer

**Default value**

1

**Component**

Only for the **NARCISSUS** component.

**Behavior**

There are five possible values for the *isscr* parameters so far:

* *isscr* = 1, normal Hubbard model/Anderson impurity model (``U`` is static)

* *isscr* = 2, Holstein-Hubbard model

* *isscr* = 3, Dynamic screening effect (``U`` is frequency-dependent), palsmon pole model

* *isscr* = 4, Dynamic screening effect (``U`` is frequency-dependent), ohmic model

* *isscr* =99, Dynamic screening effect (``U`` is frequency-dependent), realistic materials

When *isscr* = 1, *lc* and *wc* are ignored. 

When *isscr* = 2, *lc* means the electron-phonon coupling constant ``\lambda``, and wc phonon vibration frequency ``\omega_{0}``. The dynamical ``U`` is built analytically.

When *isscr* = 3, *lc* and *wc* just mean the control parameters ``\lambda`` and ``\omega^{'}``, respectively. The dynamical ``U`` is built analytically.

When *isscr* = 4, *lc* and *wc* just mean the control parameters ``\alpha`` and ``\omega_{c}``, respectively. The dynamical ``U`` is built analytically.

When *isscr* = 99, *wc* is ignored and *lc* means the shift for interaction matrix and chemical potential. The dynamical ``U`` is read from the *solver.ktau.in* file. Please see [solver.ktau.in](in_ktau.md) for more details.

**Comment**

When you want to use the quantum impurity solvers contained in iQIST software package in extended-DMFT calculations, the *isscr* parameter must be considered.

When *isscr* = 1, 3, 4, or 99, you can use the combination of improved estimator and orthogonal polynomial technology to measure ``G`` and ``\Sigma``. In other words, in such cases, the *isort* parameter can be any values (*isort* ``\in [1,6]``). On the other hand, when *isscr* = 2, the orthogonal polynomial technology is useful as well (``G`` is accurate), but the improved estimator for ``\Sigma`` and vertex function does not work any more! So in this case, you can not setup *isort* to 4, 5, or 6. See [isort](p_isort.md) parameter for more details.

*isscr* = 2 is not compatible with the *p* = 3 bit of the *isvrt* parameter. Then if you want to study the two-particle Green's function and vertex function of the Holstein-Hubbard model, you can only set the *p* = 2 bit of the *isvrt* parameter. Consult [isvrt](p_isvrt.md) parameter for more details.

Please refer to [lc](p_lc.md) and [wc](p_wc.md) parameters for more details.