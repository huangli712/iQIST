#!/usr/bin/env python

##
##
## Introduction
## ============
##
## It is a python script. The purpose of this script is to generate the
## animation movie using the data contained in the solver.diag.dat.
##
## Usage
## =====
##
## edit the configuration parameter carefully, and then execute
##
## ./u_movie.py movie.mp4
##
## Here movie.mp4 is the output file. We can use the VLC to play it. If you
## don't supply any filename, the default output should be diag.mp4.
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
## 06/05/2017 by li huang (last modified)
##
##

import sys
import argparse

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
