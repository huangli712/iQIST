#!/usr/bin/env python

import numpy
import argparse

# define command-line arguments
parser = argparse.ArgumentParser(description='Try to generate a valid hybridization function.')
parser.add_argument('fn', help='the file name')
parser.add_argument('-e', dest = 'epsilon', type = str, help = 'energy levels for bath')
parser.add_argument('-v', dest = 'V', type = float, help = 'hybridization strength')
parser.add_argument('-b', dest = 'beta', type = float, help = 'inverse temperature')
parser.add_argument('-m', dest = 'mfreq', type = int, help = 'number of matsubara frequency points')
parser.add_argument('-n', dest = 'norbs', type = int, help = 'number of orbitals')
args = parser.parse_args()
