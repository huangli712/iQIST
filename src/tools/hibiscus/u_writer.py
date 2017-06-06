#!/usr/bin/env python

##
##
## Introduction
## ============
##
## It is a python script. The purpose of it is to provide an easy-to-use
## interface to write/dump necessary input files for various quantum
## impurity solver components.
##
## Usage
## =====
##
## see the document string
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
## 02/05/2015 by li huang (created)
## 06/06/2017 by li huang (last modified)
##
##

import os
import sys
import numpy

class iqistWriter(object):
    """ This class provide a few static methods which are used to write
        the necessary input data for various ctqmc impurity solvers.

        Why do we need this class? 

        With this class, we can ensure the input file format is correct.

        typical usage:
        # import this module
        from u_writer import *

        # setup parameters
        norbs = 2
        ntime = 1024
        mfreq = 8193

        # build rmesh, hybf, wssf, symm, eimp, umat, tmesh, ktau and ptau
        ...

        # write the data
        iqistWriter.out_hyb(norbs, mfreq, rmesh, hybf)
        iqistWriter.out_eimp(norbs, symm, eimp)
        iqistWriter.out_umat(norbs, umat)
        iqistWriter.out_ktau(ntime, tmesh, ktau, ptau)
    """

    @staticmethod
    def out_hyb(norbs, mfreq, rmesh, hybf, fileName = None):
        """ try to write the hybridization function to the solver.hyb.in
            file, only suitable for the ctqmc impurity solver
        """
        if fileName is None:
            f = open("solver.hyb.in","w")
        else:
            f = open(fileName,"w")

        for i in range(norbs):
            for j in range(mfreq):
                print >> f, '%6d %16.8f %16.8f %16.8f %16.8f %16.8f' % \
                ( i+1, rmesh[j], hybf[j,i,i].real, hybf[j,i,i].imag, 0.0, 0.0 )
            print >> f
            print >> f

        f.close()

    @staticmethod
    def out_eimp(norbs, symm, eimp, fileName = None):
        """ try to write the impurity levels and symmetry vector to the
            solver.eimp.in file, only suitable for the ctqmc impurity solver
        """
        if fileName is None:
            f = open("solver.eimp.in","w")
        else:
            f = open(fileName,"w")

        for i in range(norbs):
            print >> f, '%6d %16.8f %6d' % ( i+1, eimp[i], symm[i] )

        f.close()

    @staticmethod
    def out_umat(norbs, umat, fileName = None):
        """ try to write the Coulomb matrix to the solver.umat.in file,
            only suitable for the ctqmc impurity solver
        """
        if fileName is None:
            f = open("solver.umat.in","w")
        else:
            f = open(fileName,"w")

        for i in range(norbs):
            for j in range(norbs):
                print >> f, '%6d %6d %16.8f' % ( i+1, j+1, umat[i,j] )

        f.close()

    @staticmethod
    def out_ktau(ntime, tmesh, ktau, ptau, fileName = None):
        """ try to write the screening function K(\tau) and its first
            derivates to the solver.ktau.in file, only suitable for the
            ctqmc impurity solver (narcissus)
        """
        if fileName is None:
            f = open("solver.ktau.in","w")
        else:
            f = open(fileName,"w")

        print >> f, '%s %16.8f %s %16.8f' % ( '# u shift:', 2.0 * ptau[0], '  mu shift:', ptau[0] )
        for i in range(ntime):
            print >> f, '%16.8f %16.8f %16.8f' % ( tmesh[i], ktau[i], ptau[i] )

        f.close()
