#!/usr/bin/env python

import numpy
import math
import argparse

# define command-line arguments
parser = argparse.ArgumentParser(description='Try to generate a valid screening function.')
parser.add_argument('fn', help='the file name')
parser.add_argument('-m', dest = 'model', type = int, help = 'model for screening function (1: ppm, 2: om)')
parser.add_argument('-t', dest = 'ntime', type = int, help = 'number of time slices')
parser.add_argument('-b', dest = 'beta', type = float, help = 'inverse temperature')
parser.add_argument('-l', dest = 'lc', type = float, help = 'experienced parameter for screening function lc')
parser.add_argument('-w', dest = 'wc', type = float, help = 'experienced parameter for screening function wc')
args = parser.parse_args()

# extract command-line arguments
fn = str(args.fn)
model = args.model
ntime = args.ntime
beta = args.beta
lc = args.lc
wc = args.wc

# allocate memory
ktau = numpy.zeros(ntime, dtype = numpy.float)
ptau = numpy.zeros(ntime, dtype = numpy.float)

# build linear imaginary time mesh
kmesh = numpy.linspace(0.0, beta, ntime)

if model == 1: # plasmon pole model
    for i in range(ntime):
        ktau[i] = (lc / wc)**2 / math.sinh(beta * wc / 2.0)
        ktau[i] = ktau[i] * ( math.cosh(beta * wc / 2.0) - math.cosh(beta * wc / 2.0 - kmesh[i] * wc) )
        ptau[i] = (lc / wc)**2 / sinh(beta * wc / 2.0)
        ptau[i] = ptau[i] * sinh(beta * wc / 2.0 - kmsh(i) * wc) * wc

if model == 2: # ohmic model
    #do i=1,ntime
    #    ktau(i) = lc * log(one + beta * wc * sin(pi * kmsh(i) / beta) / pi)
    #    ptau(i) = lc * wc * cos(pi * kmsh(i) / beta) / (one + beta * wc * sin(pi * kmsh(i) / beta) / pi)
    #enddo ! over i={1,ntime} loop
    pass
