# Parameter: ibasis

**Definition**

> Key control flag, it is used to specify the source for natural basis which is the eigenstate of crystal field splitting plus spin-orbit coupling.

**Type**

> Integer

**Default value**

> 1

**Component**

Only for the **JASMINE** component.

**Behavior**

> Now there are two possible values for the *ibasis* parameter:
>
> * *ibasis* = 1, make natural basis inside of this program. The code will build crystal field splitting ``\Delta_{\alpha,\beta}`` and spin-orbit coupling terms ``\Delta_{\text{SOC}}`` at first, and then use them to generate the transformation matrix ``\mathcal{T}_{\alpha,\beta}`` and on-site impurity level ``E_{\alpha,\beta}``.
>
> * *ibasis* = 2, make natural basis outside of this program. The code will read the on-site impurity level ``E_{\alpha,\beta}`` (in the *atom.emat.in* file) and transformation matrix ``\mathcal{T}_{\alpha,\beta}`` (in the *atom.tmat.in* file) from external files directly. See also [atom.emat.in](in_emat.md) and [atom.tmat.in](in_tmat.md) for more details.

**Comment**

> See also [isoc](p_isoc.md) and [lambda](p_lambda.md) parameters for more details.
