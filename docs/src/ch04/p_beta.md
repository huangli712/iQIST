# Parameter: beta

**Definition**

> Inversion of temperature ``\beta\ (\equiv 1/T)``.

**Type**

> Float, double precision

**Default value**

> 8.0

**Component**

> ALL

**Behavior**

> It is used to determine the system temperature. The larger ``\beta`` is, the slower computation is.

**Comment**

> To convert beta to real temperature, you can use the following formula:
>
> ```math
> T = 11604.505008098 / \beta
> ```
>
> The unit is Kelvin (``K``). For ``\beta = 40.0``, the corresponding temperature (290.11262520245 ``K``) is roughly the room temperature.
