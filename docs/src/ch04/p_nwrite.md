### Parameter: nwrite

**Definition**

Output period for the quantum impurity solvers.

**Type**

Integer

**Default value**

2000000

**Component**

ALL, except for the **DAISY** component.

**Behavior**

This parameter controls the output frequency of the quantum impurity solvers. The quantum impurity solvers will output the intermediate results to disk files and run-time information to the terminal for every *nwrite* Monte Carlo sampling step. As for the continuous-time quantum Monte Carlo impurity solvers, the so-called intermediate results mean the imaginary-time Green's function $$G(\tau)$$ and histogram for the perturbation expansion series.

Especially, as

$$
\frac{\text{nsweep}}{\text{nwrite}} \geq 100,
$$

the CT-HYB quantum impurity solvers will save the current diagram configurations to the *solver.diag.dat* file every *nwrite* Monte Carlo sampling step. Then, you can use the *u_animator.py* to make an animation movie.

**Comment**

When the data binning mode is activated (*isbin = 2*), then the *nwrite* and *nsweep* parameters rise tenfold. 

See also [isbin](p_isbin.md) and [nsweep](p_nsweep.md) parameters, [solver.diag.dat](out_diag.md) file, [script/u_animator.py](../ch07/animator.md) tool for more details.