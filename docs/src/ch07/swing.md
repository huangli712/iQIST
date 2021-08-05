## Analytical continuation for self-energy

### Introduction

The calculated results for the self-energy function on the real axis using Pade approximation strongly depend on the numerical accuracy of the original self-energy data. However, the CT-HYB/DMFT calculations usually yield a Matsubara self-energy function ``\Sigma(i\omega_n)`` with numerical noises. In this case, the Pade approximation does not work well. To overcome this problem, Haule *et al*[^1]. suggested to split the Matsubara self-energy function into a low-frequency part and high-frequency tail. The low-frequency part is fitted by some sort of model functions which depends on whether the system is metallic or insulating, and the high-frequency part is fitted by modified Gaussian polynomials. It was shown that their trick works quite well even when the original self-energy function is noisy, and is superior to the Pade approximation in all cases. Thus, in the **HIBISCUS** component, we also implemented this algorithm. It has broad applications in the DFT + DMFT calculations.

The **HIBISCUS**/swing code is used to continue self-energy analytically from imaginary axis to real axis using K. Haule's strategy. 

[^1]: K. Haule, C.-H. Yee, and K. Kim, Phys. Rev. B 81, 195107 (2010)

!!! note

    This code was written by K. Haule originally. See
    ```
    http://hauleweb.rutgers.edu/downloads/
    ```
    And then we adapted this code such that it can be used with the iQIST software package.

### Usage

```
$ python ./swing_main.py [options]
```

!!! note

    1. In order to run this code properly, you need to ensure *scipy* and *numpy* were correctly installed on your system. This code was tested on *scipy* 0.14.0 and *numpy* 1.7.0 only. So for the older versions of scipy and numpy, we can not guarantee that it can work always.
    2. To run this code, please use the 'make' command in the source directory (i.e., *iqist/src/tools/hibiscus/swing*) to build the dynamical library *swing_fast.so* at first. The *f2py* package must be installed and configured correctly in advance. 

### Input

The options for the **HIBISCUS**/swing code can be as follows:
```
-sig filename
  mandatory, input self-energy on imaginary axis (default:Sig.out)

-nom number
  mandatory, number of matsubara points used (default:100)

-beta float
  mandatory, inverse temperature (default:100.0)

-FL bool
  mandatory, low energy expansion of a Fermi liquid of Mott
  insulator (default:True)

-poles [[y0,A0,B0],[y1,A1,B1],...]
  optional, poles of self-energy determined from spectral function,
  poles will be forced at y0,y1,... this is achieved by contribution
  to chi2+= Al/((x_i-yl)**2+w_i*Bl) where x_i are positions and w_i
  weights of lorentzians to be determined

-b float:[0-1]
  optional, basis functions, b parameter to determin family of basis
  functions (default:0.8)

-Ng number
  optional, number of basis functions (default:12)

-wexp float:[1-2]
  optional, wexp^i is position where basis functions are
  peaked (default:1.5)

-ifit number[4-..]
  optional, low energy fit, number of matsubara points used to fit
  the low energy expansion

-alpha3 float[0-1]
  optional, weight for the normalization in functional to be
  minimized (default:0.01)

-alpha4 float[0-1]
  optional, weight for the low energy expansion in the functional
  to be minimized (default:0.001)

-maxsteps number
  optional, maximum number of function evaluations in minimization
  routine (default:400000)
```

### Output

The possible output files are as follows:

```
sigr.out
  final self-energy function on real axis

sigr_linear.out
  final self-energy function on fine real axis, it is used to
  interface with DFT + DMFT code

siom.nnn
  current function on imaginary axis (nnn counts every 40000 steps)

sres.nnn
  current analytic continuation to real axis (nnn counts every
  40000 steps)

gaus.nnn
  current configuration for gaussians (nnn counts every 40000 steps)
```

### Recipe: how to convert ``\Sigma(i\omega)`` to ``\Sigma(\omega)`` using the **HIBISCUS**/swing code

**Step 1**: 

Perform CT-HYB or HF-QMC calculations, generate a *solver.sgm.dat* file or multiple *solver.sgm.dat.``*``* files.

**Step 2**:

Using the **HIBISCUS**/toolbox/makestd to post-process the *solver.sgm.dat file* or *solver.sgm.dat.``*``* files. The output should be *std.sgm.dat* file.

**Step 3**:

Guess whether the system is metallic or insulating. Edit the corresponding Bash shell scripts (*metal.sh* or *insulator.sh*). See *iqist/working/tools/hibiscus/swing* directory for some examples.

**Step 4**:

Execute the **HIBISCUS**/swing code via the Bash shell script.

**Step 5**:

Validate the real-frequency self-energy function ``\Sigma(\omega)`` in the *sigr.dat* file or *sigr_linear.dat* file. That is what you need.

!!! note

    The **HIBISCUS**/swing code does not support multi-orbital models. So if you want to use it to post-process multi-orbital systems, you have to split the self-energy function at first by youself. Once the analytical continuation is finished, you have to combine the different *sigr.dat* or *sigr_linear.dat* files to a single file. It is a trivial task. 