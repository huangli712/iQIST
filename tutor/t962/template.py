#!/usr/bin/env python

import sys
import numpy
import pyalps.mpi

#from pyiqist import api as ctqmc

# modify sys.path
sys.path.append('../../src/tools/hibiscus/script/')

# import the writer for ctqmc configuration file
from u_ctqmc import *

# check the status of ctqmc impurity solver
#if ctqmc.solver_id() == 101:
#    if pyalps.mpi.rank == 0 : 
#        print "Hello world! This is the AZALEA code"
#else:
#    if pyalps.mpi.rank == 0 : 
#        print "Where is the AZALEA code"
#    sys.exit(-1)
#
#if ctqmc.solver_status() != 1 :
#    print "I am sorry. This ctqmc impurity solver is not ready."
#    sys.exit(-1)
#pyalps.mpi.world.barrier()

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

sys.exit(-1)


ctqmc.init_ctqmc(pyalps.mpi.rank, pyalps.mpi.size)

for i in range(2):
    ctqmc.exec_ctqmc(i+1)
print 'hdh'
#ctqmc.stop_ctqmc()
#pyalps.mpi.world.barrier()
