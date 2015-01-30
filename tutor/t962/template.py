#!/usr/bin/env python

import sys
import numpy

# import mpi support
from mpi4py import MPI

# modify sys.path
sys.path.append('../../src/tools/hibiscus/script/')

# import the writer for ctqmc configuration file
from u_ctqmc import *

# import iqist software package
from pyiqist import api as ctqmc

comm = MPI.COMM_WORLD

# check the status of ctqmc impurity solver
if ctqmc.solver_id() == 101:
    if comm.rank == 0 : 
        print "Hello world! This is the AZALEA code"
else:
    if comm.rank == 0 : 
        print "Where is the AZALEA code"
    sys.exit(-1)
if ctqmc.solver_status() != 1 :
    print "I am sorry. This ctqmc impurity solver is not ready."
    sys.exit(-1)
comm.Barrier()

# prepare the input file
if comm.rank == 0:
    # create an instance
    p = p_ctqmc_solver('azalea')

    # setup the parameters
    p.setp(isscf = 2, isbin = 1, niter = 20, U = 4.0, Uc = 4.0, Uv = 4.0, mune = 2.0, beta = 10.0)

    # verify the parameters
    p.check()

    # generate the solver.ctqmc.in file
    p.write()

    # destroy the instance
    del p
comm.Barrier()

ctqmc.init_ctqmc(comm.rank, comm.size)
a = ctqmc.get_grnf(8193 * 2 * 2)
print type(a[0])
print a
sys.exit(-1)

for i in range(1):
    ctqmc.exec_ctqmc(i+1)
ctqmc.stop_ctqmc()
comm.Barrier()
