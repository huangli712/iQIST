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
