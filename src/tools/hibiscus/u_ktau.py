#!/usr/bin/env python

import numpy
import argparse

# define command-line arguments
parser = argparse.ArgumentParser(description='Try to generate a valid screening function.')
parser.add_argument('fn', help='the file name')
parser.add_argument('-m', dest = 'model', type = int, help = 'hybridization strength')
parser.add_argument('-t', dest = 'ntime', type = int, help = 'hybridization strength')
parser.add_argument('-b', dest = 'beta', type = float, help = 'hybridization strength')
parser.add_argument('-l', dest = 'lc', type = float, help = 'hybridization strength')
parser.add_argument('-w', dest = 'wc', type = float, help = 'hybridization strength')
args = parser.parse_args()
