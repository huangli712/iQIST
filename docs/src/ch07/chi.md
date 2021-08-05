### toolbox/makechi

**Introduction**

The *makechi* code is often used to deal with the spin-spin correlation function (in imaginary time). It can output the magnetic susceptibility and effective magnetic moment. So far only the **GARDENIA** and **NARCISSUS** codes can build the spin-spin correlation function in imaginary time space. If the spin-spin correlation function is in frequency space, this code can not be used directly.

Spin-spin correlation function ``\chi_{\text{spin}}(\tau)``:

```math
\chi_{\text{spin}}(\tau) = \langle S_z(0) S_z(\tau)\rangle
```

Magnetic susceptibility ``\chi_{\text{loc}}``:

```math
\chi_{\text{loc}} = \int^{\beta}_0 d\tau \chi_{\text{spin}}(\tau)
```

Effective local magnetic moment ``M_e``:

```math
M_e = \sqrt{T\chi_{\text{loc}}}
```

**Usage**

```sh
$ ./mchi
```

**Input**

* See the terminal prompt
* *solver.schi.dat* (necessary)
 
**Output**

See the terminal output.

**Comment**

N/A