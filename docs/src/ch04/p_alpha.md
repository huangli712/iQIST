# Parameter: alpha

**Definition**

> Mixing parameter ``\alpha`` for dynamical mean-field theory self-consistent engine.

**Type**

> Float, double precision

**Default value**

> 0.7

**Component**

> ALL

**Behavior**

> In the dynamical mean-field theory self-consistent engine, the new hybridization function and self-energy function are mixed using the following equations.
>
> ```math
> \begin{equation}
> \Delta_{\text{new}} \leftarrow \Delta_{\text{old}}(1-\alpha) + \alpha \Delta_{\text{new}}
> \end{equation}
> ```
>
> ```math
> \begin{equation}
> \Sigma_{\text{new}} \leftarrow \Sigma_{\text{old}}(1-\alpha) + \alpha \Sigma_{\text{new}}
> \end{equation}
> ```
>
> Here, ``\Delta_{\text{old}}`` and ``\Sigma_{\text{old}}`` are obtained in the previous iteration, while ``\Delta_{\text{new}}`` and ``\Sigma_{\text{new}}`` are derived in the current iteration.

**Comment**

> If you have trouble in achieving convergence, you can decrease it to 0.5 or even 0.2. If you are doing one-shot calculation (*isscf* = 1), this parameter is useless. It is useful only when *isscf* = 2.
>
> See also [isscf](p_isscf.md) parameter for more details.
