#!/usr/bin/env python

import argparse
import sys

# define command-line arguments
parser = argparse.ArgumentParser(description='Split a data file into several files.')
parser.add_argument('fn', help='the file name')
parser.add_argument('-b', dest = 'nblock', type = int, help = 'number of blocks')
parser.add_argument('-l', dest = 'nline', type = int, help = 'number of lines for each data block')
args = parser.parse_args()

# extract command-line arguments
nblock =  args.nblock
nline = args.nline
fn = str(args.fn)

# try to split the given file into nblock files. each file should have
# nline lines
with open(fn, 'r') as f:
    for ib in range(nblock):
        with open(fn + '.' + str(ib + 1), 'w') as fb:
            for il in range(nline):
                print >> fb, f.readline().rstrip()
        f.readline() # there should be two empty lines after each block 
        f.readline()
