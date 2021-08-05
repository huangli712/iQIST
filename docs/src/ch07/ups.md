### toolbox/makeups

**Introduction**

The *makeups* code is often used to postprocess the spectral function data to compare with the XAS and UPS experiments. 

The Fermi-Dirac function:

```math
f_{T}(\omega) = \frac{1}{1+e^{\frac{\omega-\mu}{k_BT}}}
```

The theoretical PES:

```math
A_{\text{PES}}(\omega) = A(\omega) f_{T}(\omega)
```

The theoretical PES with Gaussian broadening:

```math
\tilde{A}_{\text{PES}}(\omega) = \int d\omega' \frac{1}{\sqrt{2\pi}\sigma} e^{\frac{(\omega-\omega')^2}{2\sigma^2}}A_{\text{PES}}(\omega')
```

The theoretical XAS:

```math
A_{\text{XAS}}(\omega) = A(\omega)f_{T}(-\omega)
```

**Usage**

```
$ ./mups
```

**Input**:

* See the terminal prompt
* *mem.dos.dat* (necessary)

!!! note

    About the smearing parameter ``\sigma``: 
    The standard deviation smearing parameter was chosen to be in the same range as estimates of experimental resolution (which are around 0.1 for high resolution PES, and approximately 0.2 to 0.4 for XAS. A good test to decide if the broadening is correct is the comparison of the Fermi edge in theory and experiment.

    About the beta parameter ``\beta (\equiv 1/T)``:
    The beta parameter practically plays no role if one uses the Fermi function at the experimental temperature or at the temperature of the CT-QMC calculations.

**Output**

* *ups.pes.dat*
* *ups.xas.dat*

**Comment**

Now this code is interfaced with the **HIBISCUS**/entropy code merely. It can read the *mem.dos.dat* file as input data. While to interface it with the others code, such as **HIBISCUS**/stoch code, is straightforward. What you need to do is just to rename *sac.imsum.dat* to *mem.dos.dat* file, and then supplement the lost orbital-dependent data.