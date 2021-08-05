### toolbox/makekra

**Introduction**

Once the analytical continuation is finished, we can obtain the spectral function ``A(\omega)`` and the imaginary part of the real-frequency Green's function ``\Im G(\omega)``,

```math
A(\omega) = -\frac{\Im G(\omega)}{\pi}
```

From the well-known Kramers-Kronig transformation, the real part of ``G(\omega)`` can be determined as well:

```math
\Re G(\omega) = -\frac{1}{\pi} \int^{\infty}_{-\infty} d\omega' \frac{\Im G(\omega)}{\omega - \omega'}
```

In the **HIBISCUS** component, we offer a Fortran code to do this job. This is the toolbox/makekra code.

The *makekra* code reads the density of states data, and then calculate the imaginary part of Matsubara green's function. And then using the Kramers-Kronig transformation, we can calculate the real part of the Matsubara green's function easily. So the complete green's function is obtained, which can be used to calculate the self-energy function on real axis by the invert Hilbert transformation.

**Usage**

```sh
$ ./mkra
```

**Input**

* See the terminal prompt
* *mem.dos.dat* (necessary)

**Output**

* *kra.grn.dat*

**Comment**

Now this code is interfaced with **HIBISCUS**/entropy code merely. It can read the *mem.dos.dat* file as input data. While to interface it with the **HIBISCUS**/stoch code is very simple. What you need to do is to rename *sac.imsum.dat* to *mem.dos.dat* file, and then supplement the data for different orbitals.