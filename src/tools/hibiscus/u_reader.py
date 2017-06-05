#!/usr/bin/env python

##
##
## Introduction
## ============
##
## It is a python script. The purpose of it is to provide an easy-to-use
## interface to read in and analyze the output data of various quantum
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
## 08/15/2015 by li huang (created)
## 06/05/2017 by li huang (last modified)
##
##

import os
import sys
import numpy

class iqistReader(object):
    """ This class provide a few static methods which are used to extract
        the data from the ouput files of various ctqmc impurity solvers.

        typical usage:
        # import this module
        from u_reader import *

        # setup parameters
        norbs = 2
        ntime = 1024
        mfreq = 8193

        # read the data
        (tmesh, gtau) = iqistReader.get_gtau(norbs, ntime)
        (rmesh, grnf) = iqistReader.get_grnf(norbs, mfreq)
    """

    @staticmethod
    def get_hist(mkink, fileName = None):
        """ try to read the solver.hist.dat file to return the histogram
            data for diagrammatic perturbation expansion
        """
        if fileName is None:
            f = open("solver.hist.dat","r")
        else:
            f = open(fileName,"r")

        hist = numpy.zeros((mkink), dtype = numpy.float)
        f.readline() # skip one comment line
        for i in range(mkink):
            spl = f.readline().split()
            hist[i] = float( spl[2] )

        f.close()

        return hist

    @staticmethod
    def get_prob(ncfgs, nsect = 0, fileName = None):
        """ try to read the solver.prob.dat file to return the atomic
            state probability P_{\Gamma} data
        """
        if fileName is None:
            f = open("solver.prob.dat","r")
        else:
            f = open(fileName,"r")

        prob = numpy.zeros((ncfgs), dtype = numpy.float)
        f.readline() # skip one comment line
        # read atomic state probability (prob)
        for i in range(ncfgs):
            spl = f.readline().split()
            prob[i] = float( spl[1] )
        if nsect > 0:
            sprob = numpy.zeros((nsect), dtype = numpy.float)
            f.readline() # skip one comment line
            # read sector probability (sprob)
            for j in range(nsect):
                spl = f.readline().split()
                sprob[j] = float( spl[2] )

        f.close()

        if nsect > 0:
            return (prob, sprob)
        else:
            return prob

    @staticmethod
    def get_nmat(norbs, fileName = None):
        """ try to read the solver.nmat.dat file to return the occupation
            number < N_i > and double occupation number < N_i N_j > data
        """
        if fileName is None:
            f = open("solver.nmat.dat","r")
        else:
            f = open(fileName,"r")

        nimp = numpy.zeros((norbs), dtype = numpy.float)
        nmat = numpy.zeros((norbs,norbs), dtype = numpy.float)
        f.readline() # skip one comment line
        # read nmat
        for i in range(norbs):
            spl = f.readline().split()
            nimp[i] = float( spl[1] )
        f.readline() # skip four lines
        f.readline()
        f.readline()
        f.readline()
        # read nnmat
        for i in range(norbs):
            for j in range(norbs):
                spl = f.readline().split()
                nmat[i,j] = float( spl[2] )

        f.close()

        return (nimp, nmat)

    @staticmethod
    def get_gtau(norbs, ntime, fileName = None):
        """ try to read the solver.green.dat file to return the imaginary
            time Green's function G(\tau) data
        """
        if fileName is None:
            f = open("solver.green.dat","r")
        else:
            f = open(fileName,"r")

        tmesh = numpy.zeros((ntime), dtype = numpy.float)
        gtau = numpy.zeros((ntime,norbs,norbs), dtype = numpy.float)
        for i in range(norbs):
            for j in range(ntime):
                spl = f.readline().split()
                tmesh[j] = float( spl[2] )
                gtau[j,i,i] = float( spl[3] )
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

        rmesh = numpy.zeros((mfreq), dtype = numpy.float)
        grnf = numpy.zeros((mfreq,norbs,norbs), dtype = numpy.complex)
        for i in range(norbs):
            for j in range(mfreq):
                spl = f.readline().split()
                rmesh[j] = float( spl[1] )
                grnf[j,i,i] = float( spl[2] ) + 1j * float( spl[3] )
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

        tmesh = numpy.zeros((ntime), dtype = numpy.float)
        wtau = numpy.zeros((ntime,norbs,norbs), dtype = numpy.float)
        for i in range(norbs):
            for j in range(ntime):
                spl = f.readline().split()
                tmesh[j] = float( spl[2] )
                wtau[j,i,i] = float( spl[3] )
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

        rmesh = numpy.zeros((mfreq), dtype = numpy.float)
        wssf = numpy.zeros((mfreq,norbs,norbs), dtype = numpy.complex)
        for i in range(norbs):
            for j in range(mfreq):
                spl = f.readline().split()
                rmesh[j] = float( spl[1] )
                wssf[j,i,i] = float( spl[2] ) + 1j * float( spl[3] )
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

        tmesh = numpy.zeros((ntime), dtype = numpy.float)
        htau = numpy.zeros((ntime,norbs,norbs), dtype = numpy.float)
        for i in range(norbs):
            for j in range(ntime):
                spl = f.readline().split()
                tmesh[j] = float( spl[2] )
                htau[j,i,i] = float( spl[3] )
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

        rmesh = numpy.zeros((mfreq), dtype = numpy.float)
        hybf = numpy.zeros((mfreq,norbs,norbs), dtype = numpy.complex)
        for i in range(norbs):
            for j in range(mfreq):
                spl = f.readline().split()
                rmesh[j] = float( spl[1] )
                hybf[j,i,i] = float( spl[2] ) + 1j * float( spl[3] )
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

        rmesh = numpy.zeros((mfreq), dtype = numpy.float)
        sig2 = numpy.zeros((mfreq,norbs,norbs), dtype = numpy.complex)
        for i in range(norbs):
            for j in range(mfreq):
                spl = f.readline().split()
                rmesh[j] = float( spl[1] )
                sig2[j,i,i] = float( spl[2] ) + 1j * float( spl[3] )
            f.readline() # skip two blank lines
            f.readline()

        f.close()

        return (rmesh, sig2)

    @staticmethod
    def get_kmat(norbs, fileName = None):
        """ try to read the solver.kmat.dat file to return the required
            perturbation order data: < k > and < k^2 >
        """
        if fileName is None:
            f = open("solver.kmat.dat","r")
        else:
            f = open(fileName,"r")

        kmat = numpy.zeros((norbs), dtype = numpy.float)
        kkmat = numpy.zeros((norbs,norbs), dtype = numpy.float)
        f.readline() # skip one comment line
        # read kmat
        for i in range(norbs):
            spl = f.readline().split()
            kmat[i] = float( spl[1] )
        f.readline() # skip two lines
        f.readline()
        # read kkmat
        for i in range(norbs):
            for j in range(norbs):
                spl = f.readline().split()
                kkmat[i,j] = float( spl[2] )

        f.close()

        return (kmat, kkmat)

    @staticmethod
    def get_lmat(norbs, fileName = None):
        """ try to read the solver.lmat.dat file to return the fidelity
            susceptibility data: < k_l >, < k_r >, and < k_l k_r >
        """
        if fileName is None:
            f = open("solver.lmat.dat","r")
        else:
            f = open(fileName,"r")

        lmat = numpy.zeros((norbs), dtype = numpy.float)
        rmat = numpy.zeros((norbs), dtype = numpy.float)
        lrmat = numpy.zeros((norbs,norbs), dtype = numpy.float)
        f.readline() # skip one comment line
        # read lmat and rmat
        for i in range(norbs):
            spl = f.readline().split()
            lmat[i] = float( spl[1] )
            rmat[i] = float( spl[2] )
        f.readline() # skip three lines
        f.readline()
        f.readline()
        # read lrmat
        for i in range(norbs):
            for j in range(norbs):
                spl = f.readline().split()
                lrmat[i,j] = float( spl[2] )

        f.close()

        return (lmat, rmat, lrmat)

    @staticmethod
    def get_schi(nband, ntime, fileName = None):
        """ try to read the solver.schi.dat file to return the spin-spin
            correlation function < S_z(0) S_z(\tau) > data
        """
        if fileName is None:
            f = open("solver.schi.dat","r")
        else:
            f = open(fileName,"r")

        tmesh = numpy.zeros((ntime), dtype = numpy.float)
        schi = numpy.zeros((ntime), dtype = numpy.float)
        sschi = numpy.zeros((ntime,nband), dtype = numpy.float)
        # read sschi
        for i in range(nband):
            f.readline() # skip one comment line
            for j in range(ntime):
                spl = f.readline().split()
                sschi[j,i] = float( spl[1] )
            f.readline() # skip two blank lines
            f.readline()
        f.readline() # skip one comment line
        # read schi
        for i in range(ntime):
            spl = f.readline().split()
            tmesh[i] = float( spl[0] )
            schi[i] = float( spl[1] )

        f.close()

        return (tmesh, schi, sschi)

    @staticmethod
    def get_sfom(nband, nbfrq, fileName = None):
        """ try to read the solver.sfom.dat file to return the spin-spin
            correlation function data
        """
        if fileName is None:
            f = open("solver.sfom.dat","r")
        else:
            f = open(fileName,"r")

        bmesh = numpy.zeros((nbfrq), dtype = numpy.float)
        ssfom = numpy.zeros((nbfrq,nband), dtype = numpy.float)
        # read ssfom
        for i in range(nband):
            f.readline() # skip one comment line
            for j in range(nbfrq):
                spl = f.readline().split()
                bmesh[j] = float( spl[0] )
                ssfom[j,i] = float( spl[1] )
            f.readline() # skip two blank lines
            f.readline()

        f.close()

        return (bmesh, ssfom)

    @staticmethod
    def get_ochi(norbs, ntime, fileName = None):
        """ try to read the solver.ochi.dat file to return the orbital-
            orbital correlation function < N_i(0) N_j(\tau) > data
        """
        if fileName is None:
            f = open("solver.ochi.dat","r")
        else:
            f = open(fileName,"r")

        tmesh = numpy.zeros((ntime), dtype = numpy.float)
        ochi = numpy.zeros((ntime), dtype = numpy.float)
        oochi = numpy.zeros((ntime,norbs,norbs), dtype = numpy.float)
        # read oochi
        for i in range(norbs):
            for j in range(norbs):
                f.readline() # skip one comment line
                for k in range(ntime):
                    spl = f.readline().split()
                    oochi[k,j,i] = float( spl[1] )
                f.readline() # skip two blank lines
                f.readline()
        f.readline() # skip one comment line
        # read ochi
        for i in range(ntime):
            spl = f.readline().split()
            tmesh[i] = float( spl[0] )
            ochi[i] = float( spl[1] )

        f.close()

        return (tmesh, ochi, oochi)

    @staticmethod
    def get_ofom(norbs, nbfrq, fileName = None):
        """ try to read the solver.ofom.dat file to return the orbital-
            orbital correlation function data
        """
        if fileName is None:
            f = open("solver.ofom.dat","r")
        else:
            f = open(fileName,"r")

        bmesh = numpy.zeros((nbfrq), dtype = numpy.float)
        oofom = numpy.zeros((nbfrq,norbs,norbs), dtype = numpy.float)
        # read oofom
        for i in range(norbs):
            for j in range(norbs):
                f.readline() # skip one comment line
                for k in range(nbfrq):
                    spl = f.readline().split()
                    bmesh[k] = float( spl[0] )
                    oofom[k,j,i] = float( spl[1] )
                f.readline() # skip two blank lines
                f.readline()

        f.close()

        return (bmesh, oofom)

    @staticmethod
    def get_twop(norbs, nffrq, nbfrq, fileName = None):
        """ try to read the solver.twop.dat file to return the two-particle
            Green's function data
        """
        if fileName is None:
            f = open("solver.twop.dat","r")
        else:
            f = open(fileName,"r")

        g2 = numpy.zeros((nffrq,nffrq,nbfrq,norbs,norbs), dtype = numpy.complex)
        f2 = numpy.zeros((nffrq,nffrq,nbfrq,norbs,norbs), dtype = numpy.complex)
        for m in range(norbs):
            for n in range(m+1):
                for k in range(nbfrq):
                    f.readline() # skip three comment lines
                    f.readline()
                    f.readline()
                    for j in range(nffrq):
                        for i in range(nffrq):
                            spl = f.readline().split()
                            g2[i,j,k,n,m] = float( spl[2] ) + 1j * float( spl[3] )
                            g2[i,j,k,m,n] = float( spl[2] ) + 1j * float( spl[3] )
                            f2[i,j,k,n,m] = float( spl[8] ) + 1j * float( spl[9] )
                            f2[i,j,k,m,n] = float( spl[8] ) + 1j * float( spl[9] )
                    f.readline() # skip two blank lines
                    f.readline()

        f.close()

        return (g2, f2)

    @staticmethod
    def get_vrtx(norbs, nffrq, nbfrq, fileName = None):
        """ try to read the solver.vrtx.dat file to return the two-particle
            Green's function data
        """
        if fileName is None:
            f = open("solver.vrtx.dat","r")
        else:
            f = open(fileName,"r")

        g2 = numpy.zeros((nffrq,nffrq,nbfrq,norbs,norbs), dtype = numpy.complex)
        f2 = numpy.zeros((nffrq,nffrq,nbfrq,norbs,norbs), dtype = numpy.complex)
        for m in range(norbs):
            for n in range(m+1):
                for k in range(nbfrq):
                    f.readline() # skip three comment lines
                    f.readline()
                    f.readline()
                    for j in range(nffrq):
                        for i in range(nffrq):
                            spl = f.readline().split()
                            g2[i,j,k,n,m] = float( spl[2] ) + 1j * float( spl[3] )
                            g2[i,j,k,m,n] = float( spl[2] ) + 1j * float( spl[3] )
                            f2[i,j,k,n,m] = float( spl[8] ) + 1j * float( spl[9] )
                            f2[i,j,k,m,n] = float( spl[8] ) + 1j * float( spl[9] )
                    f.readline() # skip two blank lines
                    f.readline()

        f.close()

        return (g2, f2)

    @staticmethod
    def get_pair(norbs, nffrq, nbfrq, fileName = None):
        """ try to read the solver.pair.dat file to return the pair
            susceptibility data
        """
        if fileName is None:
            f = open("solver.pair.dat","r")
        else:
            f = open(fileName,"r")

        p2 = numpy.zeros((nffrq,nffrq,nbfrq,norbs,norbs), dtype = numpy.complex)
        for m in range(norbs):
            for n in range(m+1):
                for k in range(nbfrq):
                    f.readline() # skip three comment lines
                    f.readline()
                    f.readline()
                    for j in range(nffrq):
                        for i in range(nffrq):
                            spl = f.readline().split()
                            p2[i,j,k,n,m] = float( spl[2] ) + 1j * float( spl[3] )
                            p2[i,j,k,m,n] = float( spl[2] ) + 1j * float( spl[3] )
                    f.readline() # skip two blank lines
                    f.readline()

        f.close()

        return p2

    @staticmethod
    def get_kernel(ntime, fileName = None):
        """ try to read the solver.kernel.dat file to return the screening
            function K(\tau) and its first derivates
        """
        if fileName is None:
            f = open("solver.kernel.dat","r")
        else:
            f = open(fileName,"r")

        tmesh = numpy.zeros((ntime), dtype = numpy.float)
        ktau = numpy.zeros((ntime), dtype = numpy.float)
        ptau = numpy.zeros((ntime), dtype = numpy.float)
        ksed = numpy.zeros((ntime), dtype = numpy.float)
        psed = numpy.zeros((ntime), dtype = numpy.float)
        for i in range(ntime):
            spl = f.readline().split()
            tmesh[i] = float( spl[1] )
            ktau[i] = float( spl[2] )
            ptau[i] = float( spl[3] )
            ksed[i] = float( spl[4] )
            psed[i] = float( spl[5] )

        f.close()

        return (tmesh, ktau, ptau, ksed, psed)
