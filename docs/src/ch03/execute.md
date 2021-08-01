## Execute codes

There are also several ways to execute the quantum impurity solvers in the iQIST software package.

### Parallelism mode

```
$ mpiexec -n num_of_cores iqist/build/ctqmc
```

### Sequential mode

```
$ ./iqist/build/ctqmc
```

### Fortran library mode

```
$ mpiexec -n num_of_cores your_fortran_program
```

or

```
$ ./your_fortran_program
```

### Python module mode

```
mpiexec -n num_of_cores your_python_script.py
```

or

```
$ ./your_python_script.py
```

> NOTE:

> 1. Be sure all of the input files are ready and in correct position. Sometimes the quantum impurity solvers can run without any inputs, but the results may be completely wrong.

> 2. Be sure the message passing interface works properly. If not, please solve it immediately.

> 3. If you are using a queue system to manage the computational jobs, please edit the submit script by yourself. We could not provide further help.

> 4. Sometimes, you want to use the multi-thread technology (OpenMP) to accelerate the quantum impurity solvers in the share memory architecture, it is wise to setup the number of threads in advance. For example:
```
$ export OMP_NUM_THREADS=number_of_cpu_cores_per_node
```

> 5. As for the Python script, the mpi4py package can be used to implement parallelism.

> 6. Even the iQIST is compiled as a paralleled program, it still can run in a sequential mode.