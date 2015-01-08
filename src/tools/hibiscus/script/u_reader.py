#!/usr/bin/env python

##
##
## Introduction
## ============
##
## It is a python script. The purpose of this script is provide an easy-
## to-use interface to read in and analyze the output data of the quantum
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
## 01/08/2015 by li huang
##
##

import os
import sys
import numpy

class iqistReader(object):
    """
    """

    @staticmethod
    def get_green(norbs, ntime, fileName = None):
        """ try to read the solver.green.dat or solver.green.bin.nnn file
            to return the imaginary time Green's function G(\tau) data
        """
        if fileName is None:
            f = open("solver.green.dat","r")
        else:
            f = open(fileName,"r")

        nband = norbs / 2
        tmesh = numpy.zeros((ntime), dtype = numpy.float)
        gtau = numpy.zeros((ntime,norbs,norbs), dtype = numpy.float)
        for i in range(nband):
            for j in range(ntime):
                spl = f.readline().split()
                tmesh[j] = float( spl[2] )
                gtau[j,i,i] = float( spl[3] )
                gtau[j,i+nband,i+nband] = float( spl[4] )
            f.readline() # skip two blank lines
            f.readline()

        f.close()

        return (tmesh, gtau)

    @staticmethod
    def get_grn(norbs, mfreq, fileName = None):
        """ try to read the solver.grn.dat file to return the matsubara
            Green's function G(i\omega) data
        """
        if fileName is None:
            f = open("solver.grn.dat","r")
        else:
            f = open(fileName,"r")

        nband = norbs / 2
        rmesh = numpy.zeros((mfreq), dtype = numpy.float)
        grnf = numpy.zeros((mfreq,norbs,norbs), dtype = numpy.complex)
        for i in range(nband):
            for j in range(mfreq):
                spl = f.readline().split()
                rmesh[j] = float( spl[1] )
                grnf[j,i,i] = float( spl[2] ) + 1j * float( spl[3] )
                grnf[j,i+nband,i+nband] = float( spl[4] ) + 1j * float( spl[5] )
            f.readline() # skip two blank lines
            f.readline()

        f.close()

        return (rmesh, grnf)

    @staticmethod
    def get_weiss(norbs, ntime, fileName = None):
        """ try to read the solver.weiss.dat file to return the imaginary
            time Weiss's function \mathcal{G}(\tau) data
        """
        if fileName is None:
            f = open("solver.weiss.dat","r")
        else:
            f = open(fileName,"r")

        nband = norbs / 2
        tmesh = numpy.zeros((ntime), dtype = numpy.float)
        wtau = numpy.zeros((ntime,norbs,norbs), dtype = numpy.float)
        for i in range(nband):
            for j in range(ntime):
                spl = f.readline().split()
                tmesh[j] = float( spl[2] )
                wtau[j,i,i] = float( spl[3] )
                wtau[j,i+nband,i+nband] = float( spl[4] )
            f.readline() # skip two blank lines
            f.readline()

        f.close()

        return (tmesh, wtau)

    @staticmethod
    def get_wss(norbs, mfreq, fileName = None):
        """ try to read the solver.wss.dat file to return the matsubara
            Weiss's function \mathcal{G}(i\omega) data
        """
        if fileName is None:
            f = open("solver.wss.dat","r")
        else:
            f = open(fileName,"r")

        nband = norbs / 2
        rmesh = numpy.zeros((mfreq), dtype = numpy.float)
        wssf = numpy.zeros((mfreq,norbs,norbs), dtype = numpy.complex)
        for i in range(nband):
            for j in range(mfreq):
                spl = f.readline().split()
                rmesh[j] = float( spl[1] )
                wssf[j,i,i] = float( spl[2] ) + 1j * float( spl[3] )
                wssf[j,i+nband,i+nband] = float( spl[4] ) + 1j * float( spl[5] )
            f.readline() # skip two blank lines
            f.readline()

        f.close()

        return (rmesh, wssf)

    @staticmethod
    def get_hybri(norbs, ntime, fileName = None):
        """ try to read the solver.hybri.dat file to return the imaginary
            time hybridization function \Delta(\tau) data
        """
        if fileName is None:
            f = open("solver.hybri.dat","r")
        else:
            f = open(fileName,"r")

        nband = norbs / 2
        tmesh = numpy.zeros((ntime), dtype = numpy.float)
        htau = numpy.zeros((ntime,norbs,norbs), dtype = numpy.float)
        for i in range(nband):
            for j in range(ntime):
                spl = f.readline().split()
                tmesh[j] = float( spl[2] )
                htau[j,i,i] = float( spl[3] )
                htau[j,i+nband,i+nband] = float( spl[4] )
            f.readline() # skip two blank lines
            f.readline()

        f.close()

        return (tmesh, htau)

    @staticmethod
    def get_hyb(norbs, mfreq, fileName = None):
        """ try to read the solver.hyb.dat file to return the matsubara
            hybridization function \Delta(i\omega) data
        """
        if fileName is None:
            f = open("solver.hyb.dat","r")
        else:
            f = open(fileName,"r")

        nband = norbs / 2
        rmesh = numpy.zeros((mfreq), dtype = numpy.float)
        hybf = numpy.zeros((mfreq,norbs,norbs), dtype = numpy.complex)
        for i in range(nband):
            for j in range(mfreq):
                spl = f.readline().split()
                rmesh[j] = float( spl[1] )
                hybf[j,i,i] = float( spl[2] ) + 1j * float( spl[3] )
                hybf[j,i+nband,i+nband] = float( spl[4] ) + 1j * float( spl[5] )
            f.readline() # skip two blank lines
            f.readline()

        f.close()

        return (rmesh, hybf)

    @staticmethod
    def get_sgm(norbs, mfreq, fileName = None):
        """ try to read the solver.sgm.dat file to return the matsubara
            self-energy function \Sigma(i\omega) data
        """
        if fileName is None:
            f = open("solver.sgm.dat","r")
        else:
            f = open(fileName,"r")

        nband = norbs / 2
        rmesh = numpy.zeros((mfreq), dtype = numpy.float)
        sig2 = numpy.zeros((mfreq,norbs,norbs), dtype = numpy.complex)
        for i in range(nband):
            for j in range(mfreq):
                spl = f.readline().split()
                rmesh[j] = float( spl[1] )
                sig2[j,i,i] = float( spl[2] ) + 1j * float( spl[3] )
                sig2[j,i+nband,i+nband] = float( spl[4] ) + 1j * float( spl[5] )
            f.readline() # skip two blank lines
            f.readline()

        f.close()

        return (rmesh, sig2)

    @staticmethod
    def get_hub(norbs, mfreq, fileName = None):
        """ try to read the solver.hub.dat file to return the matsubara
            Hubbard-I self-energy function \Sigma_{hub}(i\omega) data
        """
        pass

    @staticmethod
    def get_hist():
        """ try to read the solver.hist.dat file to return the histogram
            data for diagrammatic perturbation expansion
        """
        pass

    @staticmethod
    def get_prob():
        """ try to read the solver.prob.dat file to return the atomic
            state probability P_{\Gamma} data
        """
        pass

    @staticmethod
    def get_nmat():
        """ try to read the solver.nmat.dat file to return the occupation
            number <N_i> and double occupation number <N_i N_j> data
        """
        pass

    @staticmethod
    def get_schi():
        """ try to read the solver.schi.dat file to return the spin-spin
            correlation function <S_z(0) S_z(\tau)> data
        """
        pass

    @staticmethod
    def get_ochi():
        """ try to read the solver.ochi.dat file to return the orbital-
            orbital correlation function <N_i(0) N_j(\tau)> data
        """
        pass

    @staticmethod
    def get_twop():
        """ try to read the solver.twop.dat file to return the two-particle
            Green's function data
        """
        pass

    @staticmethod
    def get_vrtx():
        """ try to read the solver.vrtx.dat file to return the two-particle
            Green's function data
        """
        pass

    @staticmethod
    def get_pair():
        """ try to read the solver.pair.dat file to return the pair
            susceptibility data
        """
        pass

    @staticmethod
    def get_kernel():
        """ try to read the solver.kernel.dat file to return the screening
            function K(\tau)
        """
        pass

if __name__ == '__main__':
    print "hehe"

    norbs = 2
    ntime = 1024
    mfreq = 8193
    #(tmesh, gtau) = iqistReader.get_green(norbs, ntime, "solver.green.bin.10")
    #print gtau[:,1,1]
    #(rmesh, grnf) = iqistReader.get_grn(norbs, mfreq)
    #print grnf[:,0,0]
    #(rmesh, hybf) = iqistReader.get_hyb(norbs, mfreq, "solver.grn.dat")
    #print hybf[:,1,1]
    (rmesh, sig2) = iqistReader.get_sgm(norbs, mfreq)
    print sig2[:,1,1]
