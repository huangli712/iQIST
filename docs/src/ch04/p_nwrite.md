# Parameter: nwrite

**Definition**

> Output period for the quantum impurity solvers.

**Type**

> Integer

**Default value**

> 2000000

**Component**

> ALL

**Behavior**

> This parameter controls the output frequency of the quantum impurity solvers. The quantum impurity solvers will output the intermediate results to disk files and run-time information to the terminal for every *nwrite* Monte Carlo sampling step. As for the continuous-time quantum Monte Carlo impurity solvers, the so-called intermediate results mean the histogram for the perturbation expansion series.
>
> Especially, the CT-HYB quantum impurity solvers will save the current diagram configurations to the *solver.diag.dat* file every *nwrite* Monte Carlo sampling step. Then, you can use the *u_movie.py* script to make an animation movie.

**Comment**

> See also the [nsweep](p_nsweep.md) parameter, [solver.diag.dat](out_diag.md) file, [script/u_movie.py](../ch06/movie.md) tool for more details.
