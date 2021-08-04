## Maximum entropy method

### Introduction

In the Monte Carlo community, the maximum entropy method[^1] is often employed to extract the impurity spectral function $$A(\omega)$$ from the imaginary-time Green's function $$G(\tau)$$. Thus, in the **HIBISCUS** component, we implemented the standard maximum entropy algorithm. In the E-DMFT calculations, sometimes we have to perform analytical continuation for the retarded interaction function $$U(i\nu)$$ to obtain $$U(\nu)$$. So we developed a modified version of the maximum entropy method to enable this calculation.

The **HIBISCUS**/entropy code is often used to perform the analytical continuation to build impurity spectral function from imaginary-time Green's function using the well-known maximum entropy method. In principle, it solves the Laplace transformation

```math
    G(\tau) = \int K(\tau,\omega) A(\omega) d\omega
```

where $$K(\tau,\omega)$$ is the so-called kernel function. Its definition is as follows:

```math
    K(\tau,\omega) = \frac{ \exp{(-\tau\omega)} }{1.0+\exp{(-\beta\omega)}}
```

[^1]: M. Jarrell and J. Gubernatis, Phys. Rep. 269, 133 (1996).

!!! note

    The code is originally written by
    ```
    **Anders W. Sandvik**
    *Akademi University, Finland*
    *email:asandvik@ra.abo.fi*
    ```
    and modified and improved by the iQIST Developer Team using Fortran 90 language.

### Usage

```
$ ./entropy
```

or

```
$ mpiexec -n number_of_cores ./entropy
```

!!! note

    The **HIBISCUS**/entropy code also support the MPI parallelism. So you can apply MPI to improve the computational accuracy of it.

### Input

* *tau.grn.dat* (necessary)
* *entropy.in* (necessary)

The *tau.grn.dat* file contains the $$G(\tau)$$ data. It has to be generated using the **HIBISCUS**/toolbox/maketau code. 

See also [toolbox/maketau](tau.md) for more details.

The *entropy.in* file contains all of the necessary control parameters for the **HIBISCUS**/entropy code. The syntax of it is the same with the *solver.ctqmc.in* file. As for the valid control parameters, please see the following text.

See also [solver.ctqmc.in](../ch04/in_ctqmc.md) for more details.

### Output

* *mem.dos.dat*
* *mem.sum.dat*

The impurity spectral function $$A(\omega)$$ is stored in the *mem.dos.dat* file. In the *mem.sum.dat* file, the sum-rules for the impurity spectral function are examined.

### Parameters

In the following, we will show the original definitions for the control parameters:

```fortran
!!========================================================================
!!>>> integer variables                                                <<<
!!========================================================================

! number of imaginary time slices sampling by continuous time or hirsh-fye
! quantum Monte Carlo quantum impurity solver
     integer, public, save :: ntime = 129

! number of frequency points on half axis, energy range can be expressed by
! [ -wstep * nwmax, wstep * nwmax ]
     integer, public, save :: nwmax = 200

! maximum number of cycles for classic maximum entropy method
     integer, public, save :: niter = 20

! number of smooth runs for classic maximum entropy method
     integer, public, save :: ntune = 20

! number of annealing steps per classic maximum entropy method cycle
     integer, public, save :: nstep = 4000

! number of bands
     integer, public, save :: nband = 1

! number of orbitals
     integer, public, save :: norbs = 2

! the way the default model function is build
! if ntype == 0, gaussian model
! if ntype == 1, flat model
     integer, public, save :: ntype = 1

!!========================================================================
!!>>> real variables                                                   <<<
!!========================================================================

! initial alpha parameter
     real(dp), public, save :: ainit = 1200._dp

! it is the deviation from the average green's function
     real(dp), public, save :: devia = 0.001_dp

! \beta, inversion of real temperature
     real(dp), public, save :: beta  = 10.00_dp

! gauss broadening parameter, used to build the default model
     real(dp), public, save :: sigma = 1.600_dp

! delta frequency, step of real frequency grid
     real(dp), public, save :: wstep = 0.025_dp
```

### Recipe: how to convert $$G(\tau)$$ to $$A(\omega)$$ using the **HIBISCUS**/entropy code

**Step 1**: 

Perform CT-HYB or HF-QMC calculations, generate a *solver.green.dat* file or multiple *solver.green.dat.$$*$$* files.

**Step 2**:

Using the **HIBISCUS**/toolbox/maketau to post-process the *solver.green.dat file* or *solver.green.dat.$$*$$* files. The output should be *tau.grn.dat* file.

**Step 3**:

Edit the *entropy.in* file, setup reasonable control parameters.

**Step 4**:

Execute the **HIBISCUS**/entropy code, the *tau.grn.dat* and *entropy.in* files are necessary inputs.

**Step 5**:

Validate the impurity spectral function $$A(\omega)$$ in the *mem.dos.dat* file. That is what you need.

**Step 6**:

Now you can use the data in *mem.dos.dat* file to generate beautiful figures, or use the other tools to postprocess it again.