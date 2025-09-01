# u_atomic.py

**Introduction**

The purpose of this script is to generate essential input file (*solver.atomic.in*) for the **JASMINE** code. Note that you can not use it to control the **JASMINE** code directly.

**Type**

Python module

**APIs**

```python
class p_atomic_solver(object):
    """ This class can be used to generate the config file for the jasmine.
    """

    def __init__(self):
        """ define the class variables
        """

    def setp(self, **kwargs):
        """ setup the parameters using a series of key-value pairs
        """

    def check(self):
        """ check the correctness of input parameters
        """

    def write(self):
        """ write the parameters to the config file: solver.atomic.in
        """
```

**Examples**

```python
# import this module
from u_atomic import *

# create an instance
p = p_atomic_solver()

# setup the parameters
p.setp(ibasis = 2, Uv = 2.0)
p.setp(icu = 30) # invalid parameter
p.setp(icu = 1)
p.setp()

# verify the parameters
p.check()

# generate the solver.atomic.in file
p.write()

# destroy the instance
del p
```

!!! note

    You can not execute *u_atomic.py* in the terminal or Python environment directly, like this:
    ```sh
    $ python u_atomic.py
    ```

**Comment**

N/A
