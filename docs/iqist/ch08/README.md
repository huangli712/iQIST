# Application programming interfaces

The users can not only execute the components of the iQIST software package directly, but also invoke them in their own programs. To achieve this, we provide simple application programming interfaces (APIs) for most of the components in the iQIST software package in the Fortran and Python languages. With these well-defined APIs, one can easily setup, start, and stop the CT-HYB/HF-QMC quantum impurity solvers. For example, one can use the following Python script fragment to start the CT-HYB impurity solver:

```python
import mpi4py
import pyiqist as iqist
...
iqist.init_ctqmc(myid = 0, num_procs = 10)
iqist.exec_ctqmc(iter = 20)
iqist.stop_ctqmc()
```

When the computations are finished, one can also collect and analyze the calculated results with Python scripts. Using these APIs, the users enjoy more freedom to design and implement very complex computational procedures and to adapt them to their own requirements.

* [Fortran APIs](fortran.md)
* [Python APIs](python.md)

Besides the APIs, we also provide some Python modules to facilitate the development of your own computational codes/scripts.

* [script/u_atomic.py](../ch07/atomic.md) // For the **JASMINE** component.
* [script/u_ctqmc.py](../ch07/ctqmc.md) // For the CT-HYB impurity solver components.
* [script/u_hfqmc.py](../ch07/hfqmc.md) // For the **DAISY** component.
* [script/u_reader.py](../ch07/reader.md) // Read the output data.
* [script/u_writer.py](../ch07/writer.md) // Prepare the input data.