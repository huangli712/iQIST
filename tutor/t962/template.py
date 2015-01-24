#!/usr/bin/env python

import sys
import numpy
import pyalps.mpi

from pyiqist import api as ctqmc
from u_ctqmc import *

# check the status of ctqmc impurity solver
if ctqmc.solver_id() == 101:
    if pyalps.mpi.rank == 0 : 
        print "Hello world! This is the AZALEA code"
else:
    if pyalps.mpi.rank == 0 : 
        print "Where is the AZALEA code"
    sys.exit(-1)

if ctqmc.solver_status() != 1 :
    print "I am sorry. This ctqmc impurity solver is not ready."
    sys.exit(-1)
pyalps.mpi.world.barrier()

# prepare the input file
if pyalps.mpi.rank == 0:
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
pyalps.mpi.world.barrier()

ctqmc.init_ctqmc(pyalps.mpi.rank, pyalps.mpi.size)

mfreq = 8193
norbs = 2
size_t = mfreq * norbs * norbs
hybf = numpy.zeros((mfreq*norbs*norbs), dtype = numpy.complex)
print size_t
#ctqmc.set_hybf(size_t, hybf)
uumat = numpy.zeros(4, dtype = numpy.float)
ctqmc.set_uumat(4, uumat)
sys.exit(-1)

for i in range(1):
    ctqmc.exec_ctqmc(i+1)

ctqmc.stop_ctqmc()
pyalps.mpi.world.barrier()
