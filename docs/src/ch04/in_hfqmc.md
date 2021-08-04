### solver.hfqmc.in

**Introduction**

The only configuration file for the HF-QMC impurity solver in the iQIST software package is the *solver.hfqmc.in*. Just like the CT-QMC impurity solvers, since all of the input parameters have default values, the HF-QMC impurity solver can run without any input files. But if you want to use it to solve a specific problem, a well-prepared *solver.hfqmc.in* file is necessary.

**Format**

The *solver.hfqmc.in* file shares the same format with the *solver.ctqmc.in* file. See [solver.ctqmc.in](in_ctqmc.md) for more details.

!!! warning

    The quantum impurity solver will not check whether the settings in the *solver.hfqmc.in* file is reasonable and correct. It is the user's responsibility.

**Code**

N/A