#!/usr/bin/env python

import pyiqist
import sys
from mpi4py import MPI

comm = MPI.COMM_WORLD

pyiqist.cat_init_ctqmc(comm.rank,comm.size)
pyiqist.cat_exec_ctqmc(10)
if comm.rank == 0:
    print pyiqist.cat_get_nmat(2)
    print pyiqist.cat_get_nnmat(4)
    print pyiqist.cat_get_grnf(8193*4)
comm.Barrier()
pyiqist.cat_stop_ctqmc()
