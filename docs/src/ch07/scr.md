### toolbox/makescr

**Introduction**

The *makescr* code is often used to calculate the $$K(\tau)$$ from classic models or from the screening spectral function $$W(\omega)$$ (*scr.frq.dat* file). The results are necessary input for the **NARCISSUS** code. 

Retarded interaction function/screening function $$K(\tau)$$:

$$
K(\tau) = \int^{\infty}_0 \frac{d\omega}{\pi} \frac{\Im W(\omega)}{\omega^2} [\mathcal{B}(\omega,\tau)-\mathcal{B}(\omega,0)]
$$

Shifted interaction $$U_{\text{scr}}$$:

$$
U_{\text{scr}} = U + 2 \int^{\infty}_0 \frac{d\omega}{\pi}\frac{\Im W(\omega)}{\omega}
$$

Shifted chemical potential $$\mu_{\text{scr}}$$:

$$
\mu_{\text{scr}} = \mu +  \int^{\infty}_0 \frac{d\omega}{\pi}\frac{\Im W(\omega)}{\omega}
$$

**Usage**

```
$ ./mscr
```

**Input**

* See the terminal prompt
* *scr.frq.dat* (optional)
 
**Output**

* *scr.tau.dat*

**Comment**

In this code and **NARCISSUS** code, we assume that $$K(\tau)$$ is degenerated for multi-orbital system.

The $$W(\omega)$$ data are often obtained by the cRPA calculations. This feature is unfortunately not included in the iQIST software package.

In order to be compatible with the **NARCISSUS** code, you have to rename the output file (i.e., *scr.tau.dat*) to *solver.ktau.in*.

See also [solver.ktau.in](../ch04/in_ktau.md) for more details.