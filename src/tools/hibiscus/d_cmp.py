#!/usr/bin/env python

##
##
## Introduction
## ============
##
## It is a python script. The purpose of this script is to compare the
## output files between two cases when testing solvers.
##
## This script should be used by the developer only.
##
## Usage
## =====
##
## ./d_cmp.py  directory_of_case_A directory_of_case_B
##
## Author
## ======
##
## This python script is designed, created, implemented, and maintained by
##
## Yilin Wang // email: qhwyl2006@126.com
##
## History
## =======
##
## 01/12/2015 by yilin wang (created)
## 08/17/2015 by li huang (last modified)
##
##

import sys
import os
import numpy as np

from u_reader import *

def get_control(filename):
    """ read some key control parameters from input file solver.ctqmc.in,
        which will be used when reading output files.
    """
    f = open(filename,'r')
    lines = [elm.strip() for elm in f.readlines()]
    f.close()

    for elm in lines:
        if elm.startswith('norbs'):
            norbs = int( elm.split(':')[1] )

        if elm.startswith('ntime'):
            ntime = int( elm.split(':')[1] )

        if elm.startswith('mfreq'):
            mfreq = int( elm.split(':')[1] )

        if elm.startswith('mkink'):
            mkink = int( elm.split(':')[1] )

    return (norbs, ntime, mfreq, mkink)

def print_result(case, func, diff, eps):
    """ print the results of comparing
    """
    if diff < eps:
        print case, func, "diff: %10.6f  eps: %10.6f  PASSED" % (diff, eps)
    else:
        print case, func, "diff: %10.6f  eps: %10.6f  FAILED" % (diff, eps)
 
if __name__ == '__main__':

    # set the tolerance for the average relative error
    hist_eps = 1E-2
    gtau_eps = 1E-2
    grnf_eps = 1E-2
    sigf_eps = 1E-2

    # just compare a few points of low frequency for grnf and sigf
    NGF = 50
    NSF = 10

    # a is the case to be checked, b is the base case 
    dir_a, dir_b = sys.argv[1], sys.argv[2]
    print "compare ", dir_a, " with ", dir_b

    # loop for each case
    cases = [elm for elm in os.listdir(dir_a) if os.path.isdir(os.path.join(dir_a,elm))]
    cases.sort()
    for elm in cases:
        dir_case_a = os.path.join(dir_a,elm) 
        dir_case_b = os.path.join(dir_b,elm) 
        print "compare case: ", elm

        # get the key control parameters for this case
        norbs, ntime, mfreq, mkink = get_control(os.path.join(dir_case_a,'solver.ctqmc.in'))

        # compare hist
        hist_a = iqistReader.get_hist(mkink, os.path.join(dir_case_a,'solver.hist.dat'))
        hist_b = iqistReader.get_hist(mkink, os.path.join(dir_case_b,'solver.hist.dat'))
        diff = 2.0 * np.sum(np.abs(hist_a - hist_b)) / np.sum(np.abs(hist_a + hist_b))
        print_result(elm, 'hist', diff, hist_eps)
 
        # compare green
        _, gtau_a = iqistReader.get_green(norbs, ntime, os.path.join(dir_case_a,'solver.green.dat'))
        _, gtau_b = iqistReader.get_green(norbs, ntime, os.path.join(dir_case_b,'solver.green.dat'))
        diff = 2.0 * np.sum(np.abs(gtau_a - gtau_b)) / np.sum(np.abs(gtau_a + gtau_b))
        print_result(elm, 'gtau', diff, gtau_eps)

        # compare grnf
        _, grnf_a = iqistReader.get_grn(norbs, mfreq, os.path.join(dir_case_a,'solver.grn.dat'))
        _, grnf_b = iqistReader.get_grn(norbs, mfreq, os.path.join(dir_case_b,'solver.grn.dat'))
        diff = 2.0 * np.sum(np.abs(grnf_a[0:NGF,:,:] - grnf_b[0:NGF,:,:])) / np.sum(np.abs(grnf_a[0:NGF,:,:] + grnf_b[0:NGF,:,:])) 
        print_result(elm, 'grnf', diff, grnf_eps)

        # compare sigf
        _, sigf_a = iqistReader.get_sgm(norbs, mfreq, os.path.join(dir_case_a,'solver.sgm.dat'))
        _, sigf_b = iqistReader.get_sgm(norbs, mfreq, os.path.join(dir_case_b,'solver.sgm.dat'))
        diff = 2.0 * np.sum(np.abs(sigf_a[0:NSF,:,:] - sigf_b[0:NSF,:,:])) / np.sum(np.abs(sigf_a[0:NSF,:,:] + sigf_b[0:NSF,:,:]))
        print_result(elm, 'sigf', diff, sigf_eps)
