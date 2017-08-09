#!/usr/bin/env python

##
##
## Introduction
## ============
##
## It is a python script. The purpose of this script is to generate valid
## hybridization function by using the following formula:
##
##     \Delta(i\omega_n) = \sum_{\alpha} \frac{V^2}{i\omega_n - \epsilon_{\alpha}}
##
## Here V is the hybridization strength, \epsilon_{\alpha} is the energy
## level of bath, \omega_n is the matsubara frequency.
##
##
## Usage
## =====
##
## u_hyb.py [-h] [-e EPSILON] [-v V] [-b BETA] [-m MFREQ] [-n NORBS] fn
##
## Here fn is the data filename, NORBS the number of orbitals, MFREQ the
## number of matsubara frequency points, BETA the inverse temperature, V
## the hybridization strength, EPSILON the energy level of bath. Note that
## there are no default values for these arguments. And the form of the
## EPSILON argument looks like the list structure in Python language. The
## following shows a few concrete examples:
##
## ./u_hyb.py -e [-1.0,1.0] -v 1.0 -b 100.0 -m 8193 -n 2 solver.hyb.in
## ./u_hyb.py -e [-1.0,1.0] -v 0.5 -b 40.0 -m 8193 -n 4 solver.hyb.in
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
## 08/09/2017 by li huang (last modified)
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
epsilon = args.epsilon.rstrip(']').lstrip('[')
epsilon = [float(item) for item in epsilon.split(',')]
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

# build hybridization function (only diagonal part)
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
