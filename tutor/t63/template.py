#!/usr/bin/env python

import sys

# import mpi support
from mpi4py import MPI

# modify sys.path
sys.path.append('../../src/tools/hibiscus/')

# import the writer for atomic configuration file
from u_atomic import *

# modify sys.path
sys.path.append('../../src/tools/jasmine/')

# import jasmine software package
from pyjasmine import *

# get mpi communicator
comm = MPI.COMM_WORLD

# mpi barrier
comm.Barrier()

# prepare the input file
if comm.rank == 0:
    # create an instance
    p = p_atomic_solver()

    # setup the parameters
    p.setp(nband = 2, ictqmc = 3, icu = 1, Uc = 4.0, Uv = 2.0, Jz = 1.0, Js = 1.0, Jp = 1.0)

    # verify the parameters
    p.check()

    # generate the atom.config.in file
    p.write()

    # destroy the instance
    del p

# mpi barrier
comm.Barrier()

# only the master node can do this job
if comm.rank == 0:
    # init the atomic eigenvalue problem solver
    cat_init_atomic()

    # execute the atomic eigenvalue problem solver
    cat_exec_atomic()

    # stop the atomic eigenvalue problem solver
    cat_stop_atomic()

# mpi barrier
comm.Barrier()
