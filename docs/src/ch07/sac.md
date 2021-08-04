## Stochastic analytical continuation

### Introduction

An alternative way to extract the $$A(\omega)$$ from $$G(\tau)$$ is using the stochastic analytical continuation[^1]. Unlike the maximum entropy method, the stochastic analytical continuation does not depend on any *a priori* parameters. It has been argued that the stochastic analytical continuation can produce more accurate spectral functions with more subtle structures. 

In the **HIBISCUS** component, we also implemented the stochastic analytical continuation which can be viewed as a useful complementary procedure to the maximum entropy method. Since the stochastic analytical continuation is computationally much heavier than the maximum entropy method, we parallelized it with MPI and OpenMP.

The **HIBISCUS**/stoch code is often used to perform the analytical continuation to build impurity spectral function from imaginary-time Green's function using the modern stochastic analytic continuation method. In principle, it solves the Laplace transformation

$$
    G(\tau) = \int K(\tau,\omega) A(\omega) d\omega
$$

where $$K(\tau,\omega)$$ is the so-called kernel function. Its definition is as follows:

$$
    K(\tau,\omega) = \frac{ \exp{(-\tau\omega)} }{1.0+\exp{(-\beta\omega)}}
$$

[^1]: K. S. D. Beach, arXiv:0403055 [cond-mat]

!!! note

    This code is based originally on Dr. QuanSheng. Wu's code. Now QuanSheng is a postdoc in ETH Zurich, Switzerland. 

### Usage

```
$ ./sac
```

or

```
$ mpiexec -n number_of_cores ./sac
```

> NOTE:

> The **HIBISCUS**/stoch code also support the MPI parallelism. So you can apply MPI to improve the computational accuracy of it.

### Input

* *tau.grn.dat* (necessary)
* *sac.in* (necessary)

The *tau.grn.dat* file contains the $$G(\tau)$$ data. It has to be generated using the **HIBISCUS**/toolbox/maketau code. 

See also [toolbox/maketau](tau.md) for more details.

The *sac.in* file contains all of the necessary control parameters for the **HIBISCUS**/stoch code. The syntax of it is the same with the *solver.ctqmc.in* file. As for the valid control parameters, please see the following text.

See also [solver.ctqmc.in](../ch04/in_ctqmc.md) for more details.

### Output

* *sac.image.dat*
* *sac.imsum.dat*
* *sac.move.dat*
* *sac.swap.dat*

The $$\alpha$$-resolved and $$\alpha$$-averaged impurity spectral functions are stored in the *sac.image.dat* and *sac.imsum.dat* files, respectively. In the *sac.move.dat* and *sac.swap.dat* files, the statistics for the Monte Carlo move/swap operations are recorded.

### Parameters

In the following, we will show the original definitions for the control parameters:

```fortran
!!========================================================================
!!>>> integer variables                                                <<<
!!========================================================================

! number of imaginary time slices sampling by continuous time or hirsh-fye
! quantum Monte Carlo quantum impurity solver
     integer, public, save :: ntime = 1024

! number of frequency points on half axis, energy range can be expressed by
! [ -wstep * nwmax, wstep * nwmax ]
     integer, public, save :: nwmax = 128

! number of slices of x in [0,1]
     integer, public, save :: ngrid = 10001

! number of configurations, dimension for r_{\gamma} and a_{\gamma}
     integer, public, save :: ngamm = 1024

! number of alpha parameters used in parallel tempering
! note: it must be an even number, since we need to exchange configurations
! between different alpha channel
     integer, public, save :: nalph = 10

! maximum number of thermalization steps
     integer, public, save :: nwarm = 4000

! maximum number of quantum Monte Carlo sampling steps
     integer, public, save :: nstep = 4000000

! output period for stochastic analytic continuation code
     integer, public, save :: ndump = 40000

! measurement scheme
! if ltype == 1, normal measurement
! if ltype == 2, using legendre polynomial representation
     integer, public, save :: ltype = 1

! maximum order for legendre polynomial
     integer, public, save :: lemax = 64

! number of mesh points for legendre polynomial in [-1,1] range
     integer, public, save :: legrd = 20001

!!========================================================================
!!>>> real variables                                                   <<<
!!========================================================================

! initial alpha parameter
     real(dp), public, save :: ainit = 1.00_dp

! \alpha_(p+1) / \alpha_p = R, used to build alpha parameter list
     real(dp), public, save :: ratio = 2.00_dp

! \beta, inversion of real temperature
     real(dp), public, save :: beta  = 10.0_dp

! lorentz broadening parameter \eta, used to represent delta function
     real(dp), public, save :: eta1  = 0.02_dp
     real(dp), public, save :: eta2  = 4E-4_dp

! gauss broadening parameter, used to build the default model
     real(dp), public, save :: sigma = 1.00_dp

! frequency step, used to build the frequency mesh
     real(dp), public, save :: wstep = 0.05_dp
```

### Recipe: how to convert $$G(\tau)$$ to $$A(\omega)$$ using the **HIBISCUS**/stoch code

**Step 1**: 

Perform CT-HYB or HF-QMC calculations, generate a *solver.green.dat* file or multiple *solver.green.dat.$$*$$* files.

**Step 2**:

Using the **HIBISCUS**/toolbox/maketau to post-process the *solver.green.dat file* or *solver.green.dat.$$*$$* files. The output should be *tau.grn.dat* file.

**Step 3**:

Edit the *sac.in* file, setup reasonable control parameters.

**Step 4**:

Execute the **HIBISCUS**/sac code, the *tau.grn.dat* and *sac.in* files are necessary inputs.

**Step 5**:

Validate the impurity spectral function $$A(\omega)$$ in the *sac.imsum.dat* file. That is what you need.

> NOTE:

> The **HIBISCUS**/stoch code does not support multi-orbital models. So if you want to use it to post-process multi-orbital systems, you have to split the Green's function at first by yourself. Once the analytical continuation is finished, you have to combine the different *sac.imsum.dat* files to a single file. It is a trivial task. 

**Step 6**:

Now you can use the data in *sac.imsum.dat* file to generate beautiful figures, or use the other tools to postprocess it again.