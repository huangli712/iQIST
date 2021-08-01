### solver.ktau.in

**Introduction**

The *solver.ktau.in* is used to store the screening function $$K(\tau)$$ and its derivates $$K'(\tau)$$ in imaginary-time space. It is only useful for the **NARCISSUS** component.

> NOTE:

> Only when *isscr* = 99, the *solver.ktau.in* is essential for the **NARCISSUS** component. See [isscr](p_isscr.md) for more details.

**Format**

The format of the *solver.ktau.in* file is as follows:

---

*column 1*: imaginary-time point, $$\tau$$, double precision

*column 2*: screening function, $$K(\tau)$$, double precision

*column 3*: the first derivates of screening function, $$K'(\tau)$$, double precision

---

In principle, $$K(\tau)$$ and $$K'(\tau)$$ should be orbital-dependent. However in the **NARCISSUS** component, for the sake of simplicity, we treat them as orbital-independent vectors.

**Code**

The corresponding Fortran code block for the reading of *solver.ktau.in* file is as follows:

```fortran
open(mytmp, file='solver.ktau.in', form='formatted', status='unknown')
read(mytmp,*) ! skip one line
do i=1,ntime
    read(mytmp,*) rtmp, ktau(i), ptau(i)
enddo ! over i={1,ntime} loop
close(mytmp)
```

For some typical models, the *solver.ktau.in* file can be generated using the toolbox/makescr tool which is included in the **HIBISCUS** component. See [toolbox/makescr](../ch07/scr.md) for more details.