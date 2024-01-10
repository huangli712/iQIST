# u_reader.py

**Introduction**

The purpose of this script is provide an easy-to-use interface to read in and analyze the output data of the quantum impurity solver components.

**Type**

Python module

**APIs**

```python
class iqistReader(object):
    """ This class provide a few static methods which are used to extract
        the data from the ouput files of ctqmc impurity solvers and hfqmc
        impurity solver.
    """

    @staticmethod
    def get_green(norbs, ntime, fileName = None):
        """ try to read the solver.green.dat or solver.green.bin.nnn file
            to return the imaginary time Green's function G(\tau) data
        """

    @staticmethod
    def get_grn(norbs, mfreq, fileName = None):
        """ try to read the solver.grn.dat file to return the matsubara
            Green's function G(i\omega) data
        """

    @staticmethod
    def get_weiss(norbs, ntime, fileName = None):
        """ try to read the solver.weiss.dat file to return the imaginary
            time Weiss's function \mathcal{G}(\tau) data
        """

    @staticmethod
    def get_wss(norbs, mfreq, fileName = None):
        """ try to read the solver.wss.dat file to return the matsubara
            Weiss's function \mathcal{G}(i\omega) data
        """

    @staticmethod
    def get_hybri(norbs, ntime, fileName = None):
        """ try to read the solver.hybri.dat file to return the imaginary
            time hybridization function \Delta(\tau) data
        """

    @staticmethod
    def get_hyb(norbs, mfreq, fileName = None):
        """ try to read the solver.hyb.dat file to return the matsubara
            hybridization function \Delta(i\omega) data
        """

    @staticmethod
    def get_sgm(norbs, mfreq, fileName = None):
        """ try to read the solver.sgm.dat file to return the matsubara
            self-energy function \Sigma(i\omega) data
        """

    @staticmethod
    def get_hub(norbs, mfreq, fileName = None):
        """ try to read the solver.hub.dat file to return the matsubara
            Hubbard-I self-energy function \Sigma_{hub}(i\omega) data and
            Green's function data
        """

    @staticmethod
    def get_hist(mkink, fileName = None):
        """ try to read the solver.hist.dat file to return the histogram
            data for diagrammatic perturbation expansion
        """

    @staticmethod
    def get_prob(ncfgs, nsect = 0, fileName = None):
        """ try to read the solver.prob.dat file to return the atomic
            state probability P_{\Gamma} data
        """

    @staticmethod
    def get_nmat(norbs, fileName = None):
        """ try to read the solver.nmat.dat file to return the occupation
            number < N_i > and double occupation number < N_i N_j > data
        """

    @staticmethod
    def get_kmat(norbs, fileName = None):
        """ try to read the solver.kmat.dat file to return the required
            perturbation order data: < k > and < k^2 >
        """

    @staticmethod
    def get_lmat(norbs, fileName = None):
        """ try to read the solver.lmat.dat file to return the fidelity
            susceptibility data: < k_l >, < k_r >, and < k_l k_r >
        """

    @staticmethod
    def get_schi(nband, ntime, fileName = None):
        """ try to read the solver.schi.dat file to return the spin-spin
            correlation function < S_z(0) S_z(\tau) > data
        """

    @staticmethod
    def get_sfom(nband, nbfrq, fileName = None):
        """ try to read the solver.sfom.dat file to return the spin-spin
            correlation function data
        """

    @staticmethod
    def get_ochi(norbs, ntime, fileName = None):
        """ try to read the solver.ochi.dat file to return the orbital-
            orbital correlation function < N_i(0) N_j(\tau) > data
        """

    @staticmethod
    def get_ofom(norbs, nbfrq, fileName = None):
        """ try to read the solver.ofom.dat file to return the orbital-
            orbital correlation function data
        """

    @staticmethod
    def get_twop(norbs, nffrq, nbfrq, fileName = None):
        """ try to read the solver.twop.dat file to return the two-particle
            Green's function data
        """

    @staticmethod
    def get_vrtx(norbs, nffrq, nbfrq, fileName = None):
        """ try to read the solver.vrtx.dat file to return the two-particle
            Green's function data
        """

    @staticmethod
    def get_pair(norbs, nffrq, nbfrq, fileName = None):
        """ try to read the solver.pair.dat file to return the pair
            susceptibility data
        """

    @staticmethod
    def get_kernel(ntime, fileName = None):
        """ try to read the solver.kernel.dat file to return the screening
            function K(\tau) and its first derivates
        """
```

**Examples**

```python
# import this module
from u_reader import *

# setup parameters
norbs = 2
ntime = 1024
mfreq = 8193

# read the data
(tmesh, gtau) = iqistReader.get_green(norbs, ntime)
(tmesh, gbin) = iqistReader.get_green(norbs, ntime, "solver.green.bin.10")
(rmesh, grnf) = iqistReader.get_grn(norbs, mfreq)
```

!!! note

    You can not execute *u_reader.py* in the terminal or Python environment directly, like this:
    ```sh
    $ python u_reader.py
    ```

**Comment**

N/A
