### Parameter: ifast

**Definition**

It is a key control flag, used to choose the efficient algorithms for the calculation of the operator trace. 

**Type**

Integer

**Default value**

1

**Component**

Only for the **MANJUSHAKA** component.

**Behavior**

According to the *ifast* parameter, the quantum impurity solvers will choose suitable algorithm to evaluate the operator trace.

There are three possible values for *ifast* parameter so far:

* *ifast* = 1, the divide-and-conquer algorithm will be used. At this time, you have to adjust *npart* parameter carefully to obtain good computational efficiency. (see [npart](p_npart.md) parameter as well)

* *ifast* = 2, the classic time evolution algorithm will be used. 
> WARNING: 
>
> This algorithm is not implemented in the public version. We are sorry for that.

* *ifast* = 3, the skip listing algorithm will be used. 
> WARNING: 
>
> This algorithm is not implemented in the public version. We are sorry for that.

All in all, now only *ifast* = 1 is valid. Even you set *ifast* = 2 or *ifast* = 3, the **MANJUSHAKA** code will still use the divide-and-conquer algorithm to calculate the operator trace.

**Comment**

In the **MANJUSHAKA** code, the Lazy trace evaluation trick[^1] is always used. As for the divide-and-conquer algorithm and the classic time evolution algorithm, please refer to Emanuel Gull's PhD thesis. 

**Reference**

[^1]P. Sémon, Chuck-Hou Yee, Kristjan Haule, and A.-M. S. Tremblay, *Phys. Rev. B* **90**, 075149 (2014).