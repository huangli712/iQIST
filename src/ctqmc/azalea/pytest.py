import pyiqist
import sys
from mpi4py import MPI

pyiqist.cat_init_ctqmc(0,1)
print pyiqist.cat_get_nmat(2)
print pyiqist.cat_get_nnmat(4)
print pyiqist.cat_get_grnf(8193*4)
#sys.exit(-1)
pyiqist.cat_exec_ctqmc(10)
print pyiqist.cat_get_nmat(2)
print pyiqist.cat_get_nnmat(4)
print pyiqist.cat_get_grnf(8193*4)
pyiqist.cat_stop_ctqmc()
