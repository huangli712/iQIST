# Parameter: isscr

**Definition**

> It is a key control flag, which specifies what kind of impurity Hamiltonian model need to be solved. Usually, the quantum impurity solvers in the iQIST software package only support the Hubbard model/Anderson impurity model. But with the help of the **NARCISSUS** component, we have several other choices.

**Type**

> Integer

**Default value**

> 1

**Component**

> Only for the **NARCISSUS** component.

**Behavior**

> There are five possible values for the *isscr* parameters so far:
>
> * *isscr* = 1, normal Hubbard model/Anderson impurity model (``U`` is static)
>
> * *isscr* = 2, Dynamic screening effect (``U`` is frequency-dependent), palsmon pole model
>
> * *isscr* = 3, Dynamic screening effect (``U`` is frequency-dependent), ohmic model
>
> * *isscr* = 4, Dynamic screening effect (``U`` is frequency-dependent), realistic materials
>
> When *isscr* = 1, *lc* and *wc* are ignored. 
>
> When *isscr* = 2, *lc* and *wc* just mean the control parameters ``\lambda`` and ``\omega^{'}``, respectively. The dynamical ``U`` is built analytically.
>
> When *isscr* = 3, *lc* and *wc* just mean the control parameters ``\alpha`` and ``\omega_{c}``, respectively. The dynamical ``U`` is built analytically.
>
> When *isscr* = 4, *lc* and *wc* are ignored.

**Comment**

> When you want to use the quantum impurity solvers contained in iQIST software package in extended-DMFT calculations, the *isscr* parameter must be considered.
>
> When *isscr* = 1, 2, 3, or 4, you can use the combination of improved estimator and orthogonal polynomial technology to measure ``G`` and ``\Sigma``. In other words, in such cases, the *isort* parameter can be any values (*isort* ``\in [1,3]``). See [isort](p_isort.md) parameter for more details.
>
> Please refer to [lc](p_lc.md) and [wc](p_wc.md) parameters for more details.
