## Summary

The CT-HYB quantum impurity solvers usually work in the imaginary time axis, so that the calculated single-particle and two-particle correlation functions are in the imaginary-time axis [such as ``G(\tau)``] or Matsubara frequency space [such as ``G(i\omega_n)`` and ``\Sigma(i\omega_n)``]. In order to compare them with the experimental results, subtly analytical continuation procedures for them are necessary. After the analytical continuations, the correlation functions are converted into real axis.

In this section, we will illustrate that how to use the scripts/programs contained in the **ACFlow** code to perform analytical continuation for ``G(\tau)`` and ``\Sigma(i\omega_n)``. Through the analytical continuation, we can extract spectral function ``A(\omega)`` from ``G(\tau)``, and ``\Sigma(\omega)`` from ``\Sigma(i\omega_n)``. ``A(\omega)`` can be compared with the photoemission spectra directly, and ``\Sigma(\omega)`` is useful in the calculations of angle-resolved photoemission spectra ``A(\vec{k},\omega)``, optical conductivity ``\sigma(\omega)``, and transport properties, etc.

* [Analytical continuation for imaginary-time Green's function]
* [Analytical continuation for Matsubara self-energy function]
