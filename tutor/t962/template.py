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

#for i in range(1):
#    ctqmc.exec_ctqmc(i+1)
#ctqmc.stop_ctqmc()
#comm.Barrier()

norbs = 2
mfreq = 8193

#size_t = norbs * norbs
#nnmat = ctqmc.get_nnmat(size_t)
#print type(nnmat[0])
#print nnmat

#size_t = norbs
#nmat = ctqmc.get_nmat(size_t)
#print type(nmat[0])
#print nmat

#size_t = mfreq * norbs * norbs
#sigf = ctqmc.get_sigf(size_t)
#print type(sigf[0])
#print sigf

#size_t = mfreq * norbs * norbs
#grnf = ctqmc.get_grnf(size_t)
#print type(grnf[0])
#print grnf

#size_t = norbs * norbs
#uumat = numpy.zeros((norbs,norbs), dtype = numpy.float)
#uumat = uumat + 4.0
#uumat = numpy.reshape(uumat, (size_t))
#print uumat
#ctqmc.set_uumat(size_t, uumat)

size_t = norbs
eimp = numpy.zeros((norbs), dtype = numpy.float)
eimp = eimp - 1.0
print eimp
ctqmc.set_eimp(size_t, eimp)
sys.exit(-1)
