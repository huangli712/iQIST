# script/u_ctqmc.py

**Introduction**

The purpose of this script is to generate essential input file (*solver.ctqmc.in*) for the quantum impurity solver components. Note that you can not use it to control these codes directly.

**Type**

Python module

**APIs**

```python
class p_ctqmc_solver(object):
    """ This class can be used to generate the config file for the quantum
        impurity solver components.
    """

    def __init__(self, solver):
        """ define the class variables
        """

    def setp(self, **kwargs):
        """ setup the parameters using a series of key-value pairs
        """

    def check(self):
        """ check the correctness of input parameters
        """

    def write(self):
        """ write the parameters to the config file: solver.ctqmc.in
        """
```

**Examples**

```python
# import this module
from u_ctqmc import *

# create an instance
p = p_ctqmc_solver('manjushaka')

# setup the parameters
p.setp(isscf = 2, isort = 1, nsweep = 10000000)
p.setp(mune = 2.0, nmaxi = 10)
p.setp()
p.setp(isscf = 1)

# verify the parameters
p.check()

# generate the solver.ctqmc.in file
p.write()

# destroy the instance
del p
```

!!! note

    You can not execute *u_ctqmc.py* in the terminal or Python environment directly, like this:
    ```sh
    $ python u_ctqmc.py
    ```

**Comment**

N/A
