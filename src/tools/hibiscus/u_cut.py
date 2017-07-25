#!/usr/bin/env python

import argparse
import sys

parser = argparse.ArgumentParser(description='Process some integers.')
parser.add_argument('-b', dest = 'nblock', type = int, help = 'number of blocks')
parser.add_argument('-l', dest = 'nline', type = int, help = 'number of lines for each data block')
parser.add_argument('filename', help='the file name')
args = parser.parse_args()

nblock =  args.nblock
nline = args.nline
filename = str(args.filename)

print nblock, nline, filename

f = open(filename, 'r')

for ib in range(nblock):
    fb = open(filename + '.' + str(ib + 1), 'w')

    for il in range(nline):
        read_data = f.readline()
        print >> fb, read_data

    fb.close()
    f.readline()
    f.readline()

f.close()
