#!/usr/bin/env python

##
##
## Introduction
## ============
##
## It is a python script. The purpose of this script is provide an easy-
## to-use interface to write/dump necessary input files for the quantum
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
## This python script is designed, created, implemented, and maintained by
##
## Li Huang // email: huangli712@gmail.com
##
## History
## =======
##
## 02/05/2015 by li huang
##
##

import os
import sys
import numpy

class iqistWriter(object):
    """ This class provide a few static methods which are used to write
        the necessary input data for the ctqmc impurity solvers and hfqmc
        impurity solver.

        Why do we need this class? Because sometimes it is not convenient
        to call the Python API for iQIST directly. Using this class, we
        can ensure the input file format is correct.

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
        iqistWriter.out_wss(norbs, mfreq, rmesh, wssf)
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

        nband = norbs / 2
        for i in range(nband):
            for j in range(mfreq):
                print >> f, '%6d %16.8f %16.8f %16.8f %16.8f %16.8f' % \
                ( i+1, rmesh[j], hybf[j,i,i].real, hybf[j,i,i].imag,   \
                                       hybf[j,i+nband,i+nband].real,   \
                                       hybf[j,i+nband,i+nband].imag )
            print >> f
            print >> f

        f.close()

    @staticmethod
    def out_wss(norbs, mfreq, rmesh, wssf, fileName = None):
        """ try to write the bath weiss's function to the solver.wss.in
            file, only suitable for the hfqmc impurity solver
        """
        if fileName is None:
            f = open("solver.wss.in","w")
        else:
            f = open(fileName,"w")

        nband = norbs / 2
        for i in range(nband):
            for j in range(mfreq):
                print >> f, '%6d %16.8f %16.8f %16.8f %16.8f %16.8f' % \
                    ( i+1, rmesh[j], wssf[j,i].real, wssf[j,i].imag,   \
                         wssf[j,i+nband].real, wssf[j,i+nband].imag )
            print >> f
            print >> f

        f.close()

    @staticmethod
    def out_eimp(norbs, symm, eimp, fileName = None):
        """ try to write the impurity levels and symmetry vector to the
            solver.eimp.in file
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

        for i in range(ntime):
            print >> f, '%16.8f %16.8f %16.8f' % ( tmesh[i], ktau[i], ptau[i] )

        f.close()

        return (tmesh, ktau, ptau)
