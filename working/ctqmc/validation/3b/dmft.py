#!/usr/bin/env python

import sys
import numpy

# import mpi support
from mpi4py import MPI

# modify sys.path
sys.path.append('../../../../src/tools/hibiscus/')

# import the writer for ctqmc configuration file
from u_ctqmc import *

# modify sys.path
sys.path.append('../../../../src/ctqmc/manjushaka/')

# import iqist software package
import pyiqist as ctqmc

def do_dmft_loop(mfreq, norbs, grnf):
    """ implement the DMFT self-consistent condition for bethe lattice
        t1 = 1.0, t2 = 2.0, half-filling case
    """
    size_t = mfreq * norbs * norbs
    grnf_t = numpy.reshape(grnf, (mfreq, norbs, norbs), order = 'F')
    hybf_t = grnf_t * 0.0
    part_1 = 1.0
    part_2 = 2.0
    hybf_t[:,0,0] = grnf_t[:,0,0] * part_1 * part_1 # band 1, spin up
    hybf_t[:,1,1] = grnf_t[:,1,1] * part_2 * part_2 # band 2, spin up
    hybf_t[:,2,2] = grnf_t[:,2,2] * part_1 * part_1 # band 1, spin dn
    hybf_t[:,3,3] = grnf_t[:,3,3] * part_2 * part_2 # band 2, spin dn
    return numpy.reshape(hybf_t, size_t, order = 'F')

# get mpi communicator
comm = MPI.COMM_WORLD

# check the status of ctqmc impurity solver
if ctqmc.cat_solver_id() == 201:
    if comm.rank == 0 :
        print "Hello world! This is the MANJUSHAKA code."
else:
    if comm.rank == 0 :
        print "Where is the MANJUSHAKA code?"
    sys.exit(-1)
if ctqmc.cat_solver_status() != 1 :
    print "I am sorry. This ctqmc impurity solver is not ready."
    sys.exit(-1)

# mpi barrier
comm.Barrier()

# prepare the input file
if comm.rank == 0:
    # create an instance
    p = p_ctqmc_solver('manjushaka')

    # setup the parameters
    p.setp(isscf = 1, issun = 1, isbin = 1)
    p.setp(nband = 2, norbs = 4, ncfgs = 16)
    p.setp(mune = 5.25, part = 1.0, beta = 50.0)
    p.setp(nsweep = 200000000)

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
norbs = 4    # number of orbitals
niter = 20   # number of iterations
size_t = mfreq * norbs * norbs

# allocate memory
hybf = numpy.zeros(size_t, dtype = numpy.complex)
grnf = numpy.zeros(size_t, dtype = numpy.complex)
grnf_s = numpy.zeros(size_t, dtype = numpy.complex)

# init ctqmc impurity solver
ctqmc.cat_init_ctqmc(comm.rank, comm.size)

# try to implement the DMFT self-consistent loop
for i in range(niter):
    ctqmc.cat_exec_ctqmc(i+1)
    grnf = ctqmc.cat_get_grnf(size_t)
    hybf = do_dmft_loop(mfreq, norbs, grnf)
    ctqmc.cat_set_hybf(size_t, hybf)
    print 'MAX_ERROR:', (numpy.absolute(grnf - grnf_s)).max()
    grnf_s = (grnf + grnf_s)/2.0

# stop ctqmc impurity solver
ctqmc.cat_stop_ctqmc()

# mpi barrier
comm.Barrier()

# deallocate memory
del hybf, grnf, grnf_s
