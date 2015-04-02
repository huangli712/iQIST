#!/usr/bin/env python

import sys
import numpy

# import mpi support
from mpi4py import MPI

# modify sys.path
sys.path.append('../../../../src/tools/hibiscus/script/')

# import the writer for ctqmc configuration file
from u_ctqmc import *

# modify sys.path
sys.path.append('../../../../src/api/')

# import iqist software package
from pyiqist import api as ctqmc

def do_dmft_loop(mfreq, norbs, grnf):
    """ implement the DMFT self-consistent condition for bethe lattice
        t1 = 1.0, t2 = 2.0, half-filling case
    """
    size_t = mfreq * norbs * norbs
    grnf_t = numpy.reshape(grnf, (mfreq, norbs, norbs), order = 'F')
    part = 0.25
    hybf_t = grnf_t * part * part
    return numpy.reshape(hybf_t, size_t, order = 'F')

# get mpi communicator
comm = MPI.COMM_WORLD

# check the status of ctqmc impurity solver
if ctqmc.solver_id() == 102:
    if comm.rank == 0 :
        print "Hello world! This is the GARDENIA code."
else:
    if comm.rank == 0 :
        print "Where is the GARDENIA code?"
    sys.exit(-1)
if ctqmc.solver_status() != 1 :
    print "I am sorry. This ctqmc impurity solver is not ready."
    sys.exit(-1)

# mpi barrier
comm.Barrier()

# setup parameters
mfreq = 8193 # number of matsubara frequency points
norbs = 6    # number of orbitals
niter = 20   # number of iterations
mune  = 1.0  # initial chemical potential
occup = 1.6  # required occupation number
size_t = mfreq * norbs * norbs

# allocate memory
hybf = numpy.zeros(size_t, dtype = numpy.complex)
grnf = numpy.zeros(size_t, dtype = numpy.complex)
grnf_s = numpy.zeros(size_t, dtype = numpy.complex)

# mpi barrier
comm.Barrier()

# try to implement the DMFT self-consistent loop
for i in range(niter):

    # prepare the input file
    if comm.rank == 0:
        # create an instance
        p = p_ctqmc_solver('gardenia')

        # setup the parameters
        p.setp(isscf = 1, isbin = 1, issus = 2)
        p.setp(nband = 3, norbs = 6, ncfgs = 64)
        p.setp(Uc = 1.0, Uv = 0.5, Jz = 0.25)
        p.setp(mune = mune, part = 0.25, beta = 200.0)
        p.setp(nsweep = 20000000, nmonte = 50, ncarlo = 50)

        # verify the parameters
        p.check()

        # generate the solver.ctqmc.in file
        p.write()

        # destroy the instance
        del p

    # mpi barrier
    comm.Barrier()

    # init ctqmc impurity solver
    ctqmc.init_ctqmc(comm.rank, comm.size)

    # setup hybridization function
    # for the first iteration, we just use the default value
    if i > 0: ctqmc.set_hybf(size_t, hybf)

    # exec ctqmc impurity solver
    ctqmc.exec_ctqmc(i+1)

    # get green's function
    grnf = ctqmc.get_grnf(size_t)

    # get new hybridization function
    hybf = do_dmft_loop(mfreq, norbs, grnf)

    # get occupation number
    nmat = ctqmc.get_nmat(norbs)

    # calculate new chemical potential
    mune = mune + 0.3 * float( occup - sum(nmat) )

    # convergence analysis
    grnf_s = (grnf + grnf_s) / 2.0
    if comm.rank == 0:
        print 'MAX_ERROR:', (numpy.absolute(grnf - grnf_s)).max()
        print 'curr_mune:', mune, 'curr_occu:', sum(nmat), 'need_occu:', occup

    # stop ctqmc impurity solver
    ctqmc.stop_ctqmc()

# mpi barrier
comm.Barrier()

# deallocate memory
del hybf, grnf, grnf_s
