# Execute The Codes

There are also several ways to execute the quantum impurity solvers in the iQIST software package.

**Parallelism Mode**

```shell
$ mpiexec -n num_of_cores iqist/src/ctseg/ctqmc
```

```shell
$ mpiexec -n num_of_cores iqist/src/cthyb/ctqmc
```

**Sequential Mode**

```shell
$ ./iqist/src/ctseg/ctqmc
```

```shell
$ ./iqist/src/cthyb/ctqmc
```

!!! tip

    1. Be sure all of the input files are ready and in correct position. Sometimes the quantum impurity solvers can run without any inputs, but the results may be completely wrong.
    2. Be sure the message passing interface works properly. If not, please solve it immediately.
    3. If you are using a queue system to manage the computational jobs, please edit the submit script by yourself. We could not provide further help.
    4. Sometimes, you want to use the multi-thread technology (OpenMP) to accelerate the quantum impurity solvers in the share memory architecture, it is wise to setup the number of threads in advance. For example:
    ```
    $ export OMP_NUM_THREADS=number_of_cpu_cores_per_node
    ```

!!! note

    Even the iQIST is compiled as a paralleled program, it still can run in a sequential mode.
