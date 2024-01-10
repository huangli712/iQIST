# u_writer.py

**Introduction**

The purpose of this script is provide an easy-to-use interface to write/dump necessary input files for the quantum impurity solver components.

**Type**

Python module

**APIs**

```python
class iqistWriter(object):
    """ This class provide a few static methods which are used to write
        the necessary input data for the ctqmc impurity solvers and hfqmc
        impurity solver.

        Why do we need this class? Because sometimes it is not convenient
        to call the Python API for iQIST directly. Using this class, we
        can ensure the input file format is correct.
    """

    @staticmethod
    def out_hyb(norbs, mfreq, rmesh, hybf, fileName = None):
        """ try to write the hybridization function to the solver.hyb.in
            file, only suitable for the ctqmc impurity solver
        """

    @staticmethod
    def out_wss(norbs, mfreq, rmesh, wssf, fileName = None):
        """ try to write the bath weiss's function to the solver.wss.in
            file, only suitable for the hfqmc impurity solver
        """

    @staticmethod
    def out_eimp(norbs, symm, eimp, fileName = None):
        """ try to write the impurity levels and symmetry vector to the
            solver.eimp.in file
        """

    @staticmethod
    def out_umat(norbs, umat, fileName = None):
        """ try to write the Coulomb matrix to the solver.umat.in file,
            only suitable for the ctqmc impurity solver
        """

    @staticmethod
    def out_ktau(ntime, tmesh, ktau, ptau, fileName = None):
        """ try to write the screening function K(\tau) and its first
            derivates to the solver.ktau.in file, only suitable for the
            ctqmc impurity solver (narcissus)
        """
```

**Examples**

```python
# import this module
from u_writer import *

# setup parameters
norbs = 2
ntime = 1024
mfreq = 8193

# build rmesh, hybf, wssf, symm, eimp, umat, tmesh, ktau and ptau
...

# write the data
iqistWriter.out_hyb(norbs, mfreq, rmesh, hybf)
iqistWriter.out_wss(norbs, mfreq, rmesh, wssf)
iqistWriter.out_eimp(norbs, symm, eimp)
iqistWriter.out_umat(norbs, umat)
iqistWriter.out_ktau(ntime, tmesh, ktau, ptau)
```

!!! note

    You can not execute *u_writer.py* in the terminal or Python environment directly, like this:
    ```
    $ python u_writer.py
    ```

**Comment**

N/A
