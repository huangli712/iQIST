#!/usr/bin/env python

import sys
import numpy

# import mpi support
from mpi4py import MPI

# modify sys.path
sys.path.append('../../src/tools/hibiscus/')

# import the writer for ctqmc configuration file
from u_ctqmc import *

# modify sys.path
sys.path.append('../../src/ctqmc/narcissus/')

# import iqist software package
from pyiqist import *

# get mpi communicator
comm = MPI.COMM_WORLD

# check the status of ctqmc impurity solver
if cat_solver_id() == 101:
    if comm.rank == 0 :
        print "Hello world! This is the NARCISSUS code."
else:
    if comm.rank == 0 :
        print "Where is the NARCISSUS code?"
    sys.exit(-1)
if cat_solver_status() != 1 :
    print "I am sorry. This ctqmc impurity solver is not ready."
    sys.exit(-1)

# mpi barrier
comm.Barrier()

# prepare the input file
if comm.rank == 0:
    # create an instance
    p = p_ctqmc_solver('narcissus')

    # setup the parameters
    p.setp(isscf = 1, isbin = 1, U = 4.0, Uc = 4.0, Uv = 4.0, mune = 2.0, beta = 10.0)

    # verify the parameters
    p.check()

    # generate the solver.ctqmc.in file
    p.write()

    # destroy the instance
    del p

# mpi barrier
comm.Barrier()

# setup parameters
mfreq = 8193 # number of matsubara frequency points
norbs = 2    # number of orbitals
niter = 20   # number of iterations
size_t = mfreq * norbs * norbs

# allocate memory
hybf = numpy.zeros(size_t, dtype = numpy.complex)
grnf = numpy.zeros(size_t, dtype = numpy.complex)
grnf_s = numpy.zeros(size_t, dtype = numpy.complex)

# init ctqmc impurity solver
cat_init_ctqmc(comm.rank, comm.size)

# try to implement the DMFT self-consistent loop
for i in range(niter):
    cat_exec_ctqmc(i+1)
    grnf = cat_get_grnf(size_t)
    hybf = 0.25 * grnf
    cat_set_hybf(size_t, hybf)
    print 'MAX_ERROR:', (numpy.absolute(grnf - grnf_s)).max()
    grnf_s = (grnf + grnf_s)/2.0

# stop ctqmc impurity solver
cat_stop_ctqmc()

# mpi barrier
comm.Barrier()

# deallocate memory
del hybf, grnf, grnf_s
