#!/usr/bin/env python

import sys
import math
import numpy

sys.path.append("/home/huangl/iqist/iqist2015/src/tools/hibiscus/script")

from u_writer import *
from u_ctqmc import *

nband = 1
norbs = 2
ntime = 1024
mfreq = 1024
beta = 20.0 # have to CHECK IT
mune = 0.0

tmesh = numpy.zeros(ntime, dtype = numpy.float)
rmesh = numpy.zeros(mfreq, dtype = numpy.float)
ktau = numpy.zeros(ntime, dtype = numpy.float)
ptau = numpy.zeros(ntime, dtype = numpy.float)
hybf = numpy.zeros((mfreq, norbs,norbs), dtype = numpy.complex)

# prepare tmesh and rmesh
for i in range(ntime):
    tmesh[i] = float(i) * beta / ( ntime - 1 )

for i in range(mfreq):
    rmesh[i] = ( 2.0 * i + 1.0 ) * math.pi / beta

# prepare hybridization function
f = open('c_hybf.dat', 'r')
for i in range(norbs):
    for j in range(mfreq):
        spl = f.readline().split()
        hybf[j,i,i] = float( spl[3] ) + float( spl[4] ) * 1j
    f.readline()
    f.readline()
f.close()

# prepare K(\tau)
f = open('K.dat', 'r')
for i in range(ntime):
    spl = f.readline().split()
    ktau[i] = float( spl[1] )
    ptau[i] = float( spl[2] )
f.close()

# get chemical potential
f = open('mune.dat', 'r')
spl = f.readline().split()
mune = float( spl[0] ) + ptau[0]

# prepare solver.ctqmc.in
# create an instance
p = p_ctqmc_solver('narcissus')

# setup the parameters
p.setp(isscf = 1, isbin = 2, isvrt = 32, isscr = 99)
p.setp(nband = 1, norbs = norbs, ncfgs = 4)
p.setp(ntime = ntime, mfreq = mfreq, nbfrq = 1, nffrq = 128)
p.setp(U = 3.652818090610, Uc = 3.652818090610, Uv = 3.652818090610, Jz = 0.0)
p.setp(lc = float(2.0 * ptau[0]))
p.setp(beta = beta, mune = float(mune))

# verify the parameters
p.check()

# generate the solver.ctqmc.in file
p.write()

# destroy the instance
del p

# prepare solver.hyb.in and solver.ktau.in
iqistWriter.out_hyb(norbs, mfreq, rmesh, hybf)
iqistWriter.out_ktau(ntime, tmesh, ktau, ptau)
