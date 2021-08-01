### Examples

In these section, we will use some typical examples to show you how to apply the Python binding to control the atomic eigenvalue problem solver, CT-HYB and HF-QMC quantum impurity solvers.

**Case 1**: **JASMINE** component (atomic eigenvalue problem solver)

In the following, we will use an example to show you how to use API to control the **JASMINE** code. When you want to run your Python code, you have to ensure that *pyjasmine.so* is in correct **PATH**. Or else the Python will complain that it can not find iQIST.

(1) Import *pyjasmine*

```python
import pyjasmine
```

You have to ensure that the *pyjasmine* package is in the *sys.path*. For example, you can use the following code to modify *sys.path*

```python
sys.path.append('../../src/tools/jasmine/')
```

(2) Configure the atomic eigenvalue problem solver

You have to setup the key parameters for the atomic eigenvalue problem solver, and write them down to the *atom.config.in* file. Now you can do that manually. On the other hand, we provide a powerful Python module to facilitate this work (see *iqist/src/tools/hibiscus/script/u_atomic.py*).

(3) Initialize the atomic eigenvalue problem solver

```python
pyjasmine.cat_init_atomic() # there is no parameter for cat_init_atomic()
```

(4) Start the atomic eigenvalue problem solver

```python
pyjasmine.cat_exec_atomic()
```

(5) Close the atomic eigenvalue problem solver

```python
pyjasmine.cat_stop_atomic()
```

(6) Access the computational results

You have to write your own Python codes to access the results.

---

**Case 2**: **CT-HYB** quantum impurity solvers

In the following, we will use the **AZALEA** code as an example to show how to use API to control it. When you want to run your Python code, you have to ensure that *pyiqist.so* is in correct *PATH*. Or else the Python will complain that it can not find iQIST.

(1) Import MPI support

Here we can use the *pyalps.mpi* package or *mpi4py* package to provide the MPI support, such as,

```python
import pyalps.mpi
```

or

```python
from mpi4py import MPI
```

We recommend to use the *mpi4py* package since it is included in the *scipy* package already. The above code will also start the MPI running environment implicitly.

(2) Import *pyiqist*

```python
import pyiqist
```

You have to ensure that the *pyiqist* package is in the *sys.path*. For example, you can use the following code to modify *sys.path*

```python
sys.path.append('../../src/ctqmc/azalea/')
```

(3) Configure the CT-HYB quantum impurity solver

You have to setup the parameters for the CT-HYB quantum impurity solver, and write them down to the *solver.ctqmc.in* file. Now you can do that manually. On the other hand, we provide a powerful Python module to facilitate this work (see *iqist/src/tools/hibiscus/script/u_ctqmc.py*).

(4) Initialize the CT-HYB quantum impurity solver

```python
pyiqist.cat_init_ctqmc(my_id, num_procs)
```

Here *my_id* means the rank for current process, and *num_procs* means number of processes. If you are using the *mpi4py* package to provide MPI support, then you can use the following code to init the CT-HYB quantum impurity solver.

```python
comm = MPI.COMM_WORLD
pyiqist.cat_init_ctqmc(comm.rank, comm.size)
```

(5) setup *hybf*, *symm*, *eimp*, and *ktau*

For examples:

```python
mfreq = 8193
norbs = 2
size_t = mfreq * norbs * norbs
hybf = numpy.zeros(size_t, dtype = numpy.complex, order = 'F')
...
pyiqist.cat_set_hybf(size_t, hybf)
```

> NOTE: 

> 1. We strongly recommend to select the Fortran rule to manage the *numpy* array memory.

> 2. This step is optional, because the CT-HYB quantum impurity solver will provide default values for *hybf*, *symm*, *eimp*, and *ktau* or read them from external disk files.

(6) Start the CT-HYB quantum impurity solver

```python
pyiqist.cat_exec_ctqmc(i)
```

Here $$i$$ is the current iteration number.

(7) Retrieve the calculated results

For examples:

```python
size_t = norbs
nmat = pyiqist.cat_get_nmat(size_t)
print nmat

size_t = mfreq * norbs * norbs
grnf = pyiqist.cat_get_grnf(size_t)
grnf = numpy.reshape(grnf, (mfreq, norbs, norbs), order = 'F')
```

> NOTE: 

> You have to pay attention to the array order when you try to use numpy.reshape() to convert a 1-D array to a 3-D array.

(8) Close the CT-HYB quantum impurity solver

```python
pyiqist.cat_stop_ctqmc()
```

---

**Case 3**: **HF-QMC** quantum impurity solver

In the following, we will use the **DAISY** code as an example to show how to use API to control it. When you want to run your Python code, you have to ensure that *pydaisy.so* is in correct *PATH*. Or else the Python will complain that it can not find iQIST.

(1) Import MPI support

Here we can use the *pyalps.mpi* package or *mpi4py* package to provide the MPI support, such as,

```python
import pyalps.mpi
```

or

```python
from mpi4py import MPI
```

We recommend to use the *mpi4py* package since it is included in the *scipy* package already. The above code will also start the MPI running environment implicitly.

(2) Import *pydaisy*

```python
import pydaisy
```

You have to ensure that the *pydaisy* package is in the *sys.path*. For example, you can use the following code to modify *sys.path*

```python
sys.path.append('../../src/hfqmc/daisy/')
```

(3) Configure the HF-QMC quantum impurity solver

You have to setup the parameters for the HF-QMC quantum impurity solver, and write them down to the *solver.hfqmc.in* file. Now you can do that manually. On the other hand, we provide a powerful Python module to facilitate this work (see *iqist/src/tools/hibiscus/script/u_hfqmc.py*).

(4) Initialize the HF-QMC quantum impurity solver

```python
pydaisy.cat_init_hfqmc(my_id, num_procs)
```

Here *my_id* means the rank for current process, and *num_procs* means number of processes. If you are using the *mpi4py* package to provide *MPI* support, then you can use the following code to init the HF-QMC quantum impurity solver.

```python
comm = MPI.COMM_WORLD
pydaisy.cat_init_hfqmc(comm.rank, comm.size)
```

(5) Setup *wssf*, *symm*, *eimp*, and *ktau*

For examples:

```python
mfreq = 8193
norbs = 2
size_t = mfreq * norbs
wssf = numpy.zeros(size_t, dtype = numpy.complex, order = 'F')
...
pydaisy.cat_set_wssf(size_t, wssf)
```

> NOTE:

> 1. We strongly recommend to select the Fortran rule to manage the *numpy* array memory.

> 2. This step is optional, because the HF-QMC quantum impurity solver will provide default values for *wssf*, *symm*, *eimp*, and *ktau* or read them from external disk files.

(6) Start the HF-QMC quantum impurity solver

```python
pydaisy.cat_exec_hfqmc(i)
```

Here $$i$$ is the current iteration number.

(7) Retrieve the calculated results

For examples:

```python
size_t = norbs
nmat = pydaisy.cat_get_nmat(size_t)
print nmat

size_t = mfreq * norbs
grnf = pydaisy.cat_get_grnf(size_t)
grnf = numpy.reshape(grnf, (mfreq, norbs), order = 'F')
```

> NOTE:

> You have to pay attention to the array order when you try to use ```numpy.reshape()``` to convert a 1-D array to a 2-D array.

(8) Close the HF-QMC quantum impurity solver

```python
pydaisy.cat_stop_hfqmc()
```