#!/usr/bin/env python

import sys

# import mpi support
from mpi4py import MPI

# modify sys.path
sys.path.append('../../src/tools/hibiscus/script/')

# import the writer for atomic configuration file
from u_atomic import *

# modify sys.path
sys.path.append('../../src/tools/jasmine/')

# import jasmine software package
from pyjasmine import japi as atomic

# get mpi communicator
comm = MPI.COMM_WORLD

# mpi barrier
comm.Barrier()

# prepare the input file
if comm.rank == 0:
    # create an instance
    p = p_atomic_solver()

    # setup the parameters
    p.setp(ibasis = 2, Uv = 2.0, icu = 2)

    # verify the parameters
    p.check()

    # generate the atom.config.in file
    p.write()

    # destroy the instance
    del p

# mpi barrier
comm.Barrier()

atomic.init_atomic()
atomic.exec_atomic()
atomic.stop_atomic()

# mpi barrier
comm.Barrier()
