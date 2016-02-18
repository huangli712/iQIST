#!/usr/bin/env python

import sys
import numpy

# import mpi support
from mpi4py import MPI

# modify sys.path
sys.path.append('../../src/tools/hibiscus/script/')

# import the writer for hfqmc configuration file
from u_hfqmc import *

# modify sys.path
sys.path.append('../../src/hfqmc/daisy/')

# import daisy software package
from pydaisy import *

def do_dmft_loop(mfreq, norbs, grnf):
    """ implement the DMFT self-consistent condition for bethe lattice
        t = 0.5, beta = 10.0, half-filling case
    """
    size_t = mfreq * norbs
    rmesh = numpy.zeros(mfreq, dtype = numpy.float)
    for i in range(mfreq):
        rmesh[i] = (2 * i + 1.0) * numpy.pi / 10.0
    grnf_t = numpy.reshape(grnf, (mfreq, norbs), order = 'F')
    wssf_t = grnf_t * 0.0
    for i in range(norbs):
        for j in range(mfreq):
            wssf_t[j,i] = 1.0 / (1j * rmesh[j] - 0.25 * grnf_t[j,i])
    return numpy.reshape(wssf_t, size_t, order = 'F')

# get mpi communicator
comm = MPI.COMM_WORLD

# check the status of hfqmc impurity solver
if cat_solver_id() == 901:
    if comm.rank == 0 :
        print "Hello world! This is the DAISY code."
else:
    if comm.rank == 0 :
        print "Where is the DAISY code?"
    sys.exit(-1)
if cat_solver_status() != 1 :
    print "I am sorry. This hfqmc impurity solver is not ready."
    sys.exit(-1)

# mpi barrier
comm.Barrier()

# prepare the input file
if comm.rank == 0:
    # create an instance
    p = p_hfqmc_solver('daisy')

    # setup the parameters
    p.setp(isscf = 1, isbin = 1, Uc = 4.0, mune = 2.0, beta = 10.0)

    # verify the parameters
    p.check()

    # generate the solver.hfqmc.in file
    p.write()

    # destroy the instance
    del p

# mpi barrier
comm.Barrier()

# setup parameters
mfreq = 8193 # number of matsubara frequency points
norbs = 2    # number of orbitals
niter = 20   # number of iterations
size_t = mfreq * norbs

# allocate memory
wssf = numpy.zeros(size_t, dtype = numpy.complex, order = 'F')
grnf = numpy.zeros(size_t, dtype = numpy.complex, order = 'F')
grnf_s = numpy.zeros(size_t, dtype = numpy.complex, order = 'F')

# init hfqmc impurity solver
cat_init_hfqmc(comm.rank, comm.size)

# try to implement the DMFT self-consistent loop
for i in range(niter):
    cat_exec_hfqmc(i+1)
    grnf = cat_get_grnf(size_t)
    wssf = do_dmft_loop(mfreq, norbs, grnf)
    cat_set_wssf(size_t, wssf)
    print 'MAX_ERROR:', (numpy.absolute(grnf - grnf_s)).max()
    grnf_s = (grnf + grnf_s)/2.0

# stop hfqmc impurity solver
cat_stop_hfqmc()

# mpi barrier
comm.Barrier()

# deallocate memory
del wssf, grnf, grnf_s
