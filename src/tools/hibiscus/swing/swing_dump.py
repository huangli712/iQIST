#!/usr/bin/env python

##
##
## Introduction
## ============
##
## It is a python script. The purpose of this script is provide the basic
## writing functions for the swing code. Now it implements the following
## python functions/classes:
##
##     def swing_print_header
##     def swing_print_footer
##     def swing_dump_hist
##     def swing_dump_keys
##     def swing_dump_sigr
##     def swing_dump_gaus
##     def swing_dump_siom
##     def swing_dump_sres
##
## Usage
## =====
##
## Sorry, it can not be invoked manually.
##
## Author
## ======
##
## This python script is designed, created, implemented, and maintained by
##
## Li Huang // email: huangli712@gmail.com
##
## History
## =======
##
## 12/20/2014 by li huang
##
##

from scipy import *
from scipy import interpolate

def swing_print_header():
    """ to display the header for the hibiscus/swing code to the screen
    """
    print '  HIBISCUS/swing'
    print '  >>> A Stochastic Analytic Continuation Code for Self-Energy Data'
    print
    print '  Version: 2014.10.11T'
    print '  Develop: by li huang (at IOP/CAS & SPCLab/CAEP & UNIFR)'
    print '  Support: huangli712@gmail.com'
    print '  License: GNU General Public License version 3'
    print

def swing_print_footer(tot_time):
    """ to display the footer for the hibiscus/swing code to the screen
    """
    print '  HIBISCUS/swing >>> total time spent:', tot_time, 's'
    print '  HIBISCUS/swing >>> I am tired and want to go to bed. Bye!'
    print '  HIBISCUS/swing >>> happy ending'

def swing_dump_hist(argv):
    """ to record the command history in hist.dat file
    """
    fs = open('hist.dat','a')
    for cmd in argv:
        print >> fs, cmd,
    print >> fs
    fs.close()

def swing_dump_keys(params):
    """ to display the parameter lists to the screen
    """
    print '  HIBISCUS/swing >>> parameters list:'
    for var in params.keys():
        print '%s %-8s %s %-6s %s' % ('   ', var, ':', params[var][0], params[var][1])
    print

def swing_dump_sigr(om, vary, fixed, gweigh, gwfix, rfunc, sinfty):
    """ dump the final self-energy function in real axis to sigr.out file
        and sigr_linear.out file
    """
    # zsum is used to store the self-energy function on real axis
    zsum = []
    for im in range(len(om)):
        zsum.append(0.0)

    # calculate self-energy function on real logarithm axis
    fs = open('sigr.out','w')
    for im in range(len(om)):
        for i in range(len(gweigh)):
            zsum[im] += rfunc[vary[i]-1,im] * gweigh[i]
        for i in range(len(fixed)):
            zsum[im] += rfunc[fixed[i]-1,im] * gwfix[i]
        print >> fs, '%16.8f %16.8f %16.8f' % (om[im], zsum[im].real+sinfty, zsum[im].imag)
    fs.close()

    # build real symmetry linear axis
    step = 0.04
    nfrq = 400
    new_om = []
    for im in range(2*nfrq+1):
        new_om.append(step*(im - nfrq))

    # interpolate self-energy function from om mesh to new_om mesh
    # the results are dumped into sigr_linear.out
    spl_re = interpolate.splrep(om, real(zsum), k=3, s=0.0)
    spl_im = interpolate.splrep(om, imag(zsum), k=3, s=0.0)
    fs = open('sigr_linear.out','w')
    for im in range(2*nfrq+1):
        sig_re = interpolate.splev(new_om[im], spl_re, der=0) + sinfty
        sig_im = interpolate.splev(new_om[im], spl_im, der=0)
        print >> fs, '%6d %16.8f %16.8f %16.8f' % (im, new_om[im], sig_re, sig_im)
    fs.close()

def swing_dump_gaus(it, gpos, gweigh):
    """ dump the modified gaussian function and their weights to gaus.nnn
    """
    fs = open('gaus.'+str(it),'w')
    for i in range(len(gweigh)):
        print >> fs, '%4d %16.8f %16.8f' % (i, gpos[i], gweigh[i])
    fs.close()

def swing_dump_siom(it, iom, vary, fixed, gweigh, gwfix, ifunr, ifuni):
    """ dump the fitted self-energy function in matsubara axis to siom.nnn
    """
    fs = open('siom.'+str(it),'w')
    for im in range(len(iom)):
        gc = 0j
        for i in range(len(gweigh)):
            gc += (ifunr[vary[i]-1,im]  + 1j * ifuni[vary[i]-1,im]) * gweigh[i]
        for i in range(len(fixed)):
            gc += (ifunr[fixed[i]-1,im] + 1j * ifuni[fixed[i]-1,im]) * gwfix[i]
        print >> fs, '%16.8f %16.8f %16.8f' % (iom[im], gc.real, gc.imag)
    fs.close()

def swing_dump_sres(it, om, vary, fixed, gweigh, gwfix, rfunc):
    """ dump the fitted self-energy function in real axis to sres.nnn
    """
    fs = open('sres.'+str(it),'w')
    for im in range(len(om)):
        zsum = 0.0
        for i in range(len(gweigh)):
            zsum += rfunc[vary[i]-1,im] * gweigh[i]
        for i in range(len(fixed)):
            zsum += rfunc[fixed[i]-1,im] * gwfix[i]
        print >> fs, '%16.8f %16.8f %16.8f' % (om[im], zsum.real, zsum.imag)
    fs.close()
