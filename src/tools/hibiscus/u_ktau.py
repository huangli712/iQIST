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

# try to build screening function
if model == 1: # plasmon pole model
    p = beta * wc / 2.0
    q = (lc / wc)**2
    for i in range(ntime):
        ktau[i] = q / math.sinh(p) * ( math.cosh(p) - math.cosh(p - kmesh[i] * wc) )
        ptau[i] = q / math.sinh(p) * math.sinh(p - kmesh[i] * wc) * wc

# try to build screening function
if model == 2: # ohmic model
    q = beta * wc / math.pi
    for i in range(ntime):
        p = math.pi * kmesh[i] / beta
        ktau[i] = lc * math.log(1.0 + q * math.sin(p))
        ptau[i] = lc * wc * math.cos(p) / (1.0 + q * math.sin(p))

# write the screening function
with open(fn, 'w') as f:
    print >> f, '# u shift: %16.8f mu shift: %16.8f' % ( 2.0*ptau[0], ptau[0] )
    for i in range(ntime):
        print >> f, '%16.8f %16.8f %16.8f' % ( kmesh[i], ktau[i], ptau[i] )
