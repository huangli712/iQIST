#!/usr/bin/env python

##
##
## Introduction
## ============
##
## It is a python script. The purpose of this script is to split the data
## file into several pieces. We usually use it to deal with the following
## files:
##     solver.green.dat
##     solver.sgm.dat
##     solver.sp_w.dat
##     solver.sp_t.dat
##     solver.ch_w.dat
##     solver.ch_t.dat
## so that they can be postprocessed by the other analytical continuation
## codes, such as SpM or OmegaMaxEnt.
##
## Usage
## =====
##
## ./u_cut.py [-h] [-b NBLOCK] [-l NLINE] fn
##
## Here fn is the original data file. NBLOCK means number of blocks, and
## NLINE means number of lines for each block. The output files should be
## fn.*, * means the index of data block. In other words, each data block
## in the original data file should be converted into a single file. The
## following shows a few concrete examples:
##
## Solit the imaginary-time green's function
##     ./u_cut -b 2 -l 1024 solver.green.dat
##
## Split the matsubara self-energy function
##     ./u_cut -b 2 -l 8193 solver.sgm.dat
##
#fn = 'solver.hyb.in'
#epsilon = [-1.0, 0.0, 1.0]
#V = 0.5
#beta = 10.0
#mfreq = 8193
#norbs = 2
##
## Author
## ======
##
## This python script is designed, created, and maintained by
##
## Li Huang // email: lihuang.dmft@gmail.com
##
## History
## =======
##
## 03/28/2015 by li huang (created)
## 07/25/2017 by li huang (last modified)
##
##

import numpy
import argparse
import sys

# define command-line arguments
parser = argparse.ArgumentParser(description='Try to generate a valid hybridization function.')
parser.add_argument('fn', help='the file name')
parser.add_argument('-e', dest = 'epsilon', type = str, help = 'energy levels for bath')
parser.add_argument('-v', dest = 'V', type = float, help = 'hybridization strength')
parser.add_argument('-b', dest = 'beta', type = float, help = 'inverse temperature')
parser.add_argument('-m', dest = 'mfreq', type = int, help = 'number of matsubara frequency points')
parser.add_argument('-n', dest = 'norbs', type = int, help = 'number of orbitals')
args = parser.parse_args()

# extract command-line arguments
fn = str(args.fn)
epsilon = [float(item) for item in args.epsilon.split(',')]
V = args.V
beta = args.beta
mfreq = args.mfreq
norbs = args.norbs

# allocate memory
rmesh = numpy.zeros(mfreq, dtype = numpy.float)
hybf = numpy.zeros((mfreq,norbs,norbs), dtype = numpy.complex)

# build matsubara frequency mesh
for i in range(mfreq):
    rmesh[i] = (2*i + 1) * numpy.pi / beta

# build hybridization function
for i in range(mfreq):
    for e in epsilon:
        for b in range(norbs):
            hybf[i,b,b] = hybf[i,b,b] + V*V / ( 1j * rmesh[i] - e )

# write the hybridization function
with open(fn, 'w') as f:
    for i in range(norbs):
        for j in range(mfreq):
            print >> f, '%6d %16.8f %16.8f %16.8f %16.8f %16.8f' % \
            ( i+1, rmesh[j], hybf[j,i,i].real, hybf[j,i,i].imag, 0.0, 0.0 )
        print >> f
        print >> f
