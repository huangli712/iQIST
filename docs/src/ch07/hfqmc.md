### script/u_hfqmc.py

**Introduction**

The purpose of this script is to generate essential input file (*solver.hfqmc.in*) for the **DAISY** code. Note that you can not use it to control the **DAISY** code directly.

**Type**

Python module

**APIs**

```python
class p_hfqmc_solver(object):
    """ This class can be used to generate the config file for the quantum
        impurity solver daisy code.
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
        """ write the parameters to the config file: solver.hfqmc.in
        """
```

**Examples**

```python
# import this module
from u_hfqmc import *

# create an instance
p = p_hfqmc_solver('daisy')

# setup the parameters
p.setp(isscf = 2, nsweep = 10000000)
p.setp(mune = 2.0, mstep = 10)
p.setp()
p.setp(isscf = 1)

# verify the parameters
p.check()

# generate the solver.hfqmc.in file
p.write()

# destroy the instance
del p
```

!!! note

    You can not execute *u_hfqmc.py* in the terminal or Python environment directly, like this:
    ```
    $ python u_hfqmc.py
    ```

**Comment**

N/A