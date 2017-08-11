#!/usr/bin/env python

import numpy
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
kmesh = numpy.zeros(ntime, dtype = numpy.float)
ktau = numpy.zeros(ntime, dtype = numpy.float)
ptau = numpy.zeros(ntime, dtype = numpy.float)
