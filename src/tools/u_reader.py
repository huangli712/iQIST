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
## See the document string
##
## Author
## ======
##
## This python script is designed, created, and maintained by
##
## Li Huang // email: huangli@caep.cn
##
## History
## =======
##
## 08/15/2015 by li huang (created)
## 01/10/2024 by li huang (last modified)
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
    def get_ac_f(ntime, fileName = None):
        """ try to read the solver.ac_f.dat file to return the auto-
            correlation function data which can be used to estimate
            reasonable nmonte parameter
        """
        if fileName is None:
            f = open("solver.ac_f.dat","r")
        else:
            f = open(fileName,"r")

        ac_f = numpy.zeros((ntime), dtype = numpy.float)
        for i in range(ntime):
            spl = f.readline().split()
            ac_f[i] = float( spl[1] )

        f.close()

        return ac_f

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
    def get_paux(fileName = None):
        """ try to read the solver.paux.dat file to return the auxiliary
            physical observables
        """
        if fileName is None:
            f = open("solver.paux.dat","r")
        else:
            f = open(fileName,"r")

        paux = numpy.zeros((9), dtype = numpy.float)
        for i in range(9):
            spl = f.readline().split()
            paux[i] = float( spl[1] )
        f.close()

        return paux

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
        # read nimp
        for i in range(norbs):
            spl = f.readline().split()
            nimp[i] = float( spl[1] )
        f.readline() # skip four lines
        f.readline()
        f.readline()
        f.readline()
        # read nmat
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
    def get_ftau(norbs, ntime, fileName = None):
        """ try to read the solver.fcorr.dat file to return the imaginary
            time auxiliary correlation function F(\tau) data
        """
        if fileName is None:
            f = open("solver.fcorr.dat","r")
        else:
            f = open(fileName,"r")

        tmesh = numpy.zeros((ntime), dtype = numpy.float)
        ftau = numpy.zeros((ntime,norbs,norbs), dtype = numpy.float)
        for i in range(norbs):
            for j in range(ntime):
                spl = f.readline().split()
                tmesh[j] = float( spl[2] )
                ftau[j,i,i] = float( spl[3] )
            f.readline() # skip two blank lines
            f.readline()

        f.close()

        return (tmesh, ftau)

    @staticmethod
    def get_htau(norbs, ntime, fileName = None):
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
    def get_wtau(norbs, ntime, fileName = None):
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
    def get_ktau(ntime, fileName = None):
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

    @staticmethod
    def get_grnf(norbs, mfreq, fileName = None):
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
    def get_frnf(norbs, mfreq, fileName = None):
        """ try to read the solver.frn.dat file to return the matsubara
            auxiliary correlation function F(i\omega) data
        """
        if fileName is None:
            f = open("solver.frn.dat","r")
        else:
            f = open(fileName,"r")

        rmesh = numpy.zeros((mfreq), dtype = numpy.float)
        frnf = numpy.zeros((mfreq,norbs,norbs), dtype = numpy.complex)
        for i in range(norbs):
            for j in range(mfreq):
                spl = f.readline().split()
                rmesh[j] = float( spl[1] )
                frnf[j,i,i] = float( spl[2] ) + 1j * float( spl[3] )
            f.readline() # skip two blank lines
            f.readline()

        f.close()

        return (rmesh, frnf)

    @staticmethod
    def get_hybf(norbs, mfreq, fileName = None):
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
    def get_wssf(norbs, mfreq, fileName = None):
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
    def get_sig2(norbs, mfreq, fileName = None):
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

        knop = numpy.zeros((norbs), dtype = numpy.float)
        kmat = numpy.zeros((norbs,norbs), dtype = numpy.float)
        f.readline() # skip one comment line
        # read knop
        for i in range(norbs):
            spl = f.readline().split()
            knop[i] = float( spl[1] )
        f.readline() # skip two lines
        f.readline()
        # read kmat
        for i in range(norbs):
            for j in range(norbs):
                spl = f.readline().split()
                kmat[i,j] = float( spl[2] )

        f.close()

        return (knop, kmat)

    @staticmethod
    def get_lrmm(norbs, fileName = None):
        """ try to read the solver.lrmm.dat file to return the fidelity
            susceptibility data: < k_l >, < k_r >, and < k_l k_r >
        """
        if fileName is None:
            f = open("solver.lrmm.dat","r")
        else:
            f = open(fileName,"r")

        lnop = numpy.zeros((norbs), dtype = numpy.float)
        rnop = numpy.zeros((norbs), dtype = numpy.float)
        lrmm = numpy.zeros((norbs,norbs), dtype = numpy.float)
        f.readline() # skip one comment line
        # read lnop and rnop
        for i in range(norbs):
            spl = f.readline().split()
            lnop[i] = float( spl[1] )
            rnop[i] = float( spl[2] )
        f.readline() # skip three lines
        f.readline()
        f.readline()
        # read lrmm
        for i in range(norbs):
            for j in range(norbs):
                spl = f.readline().split()
                lrmm[i,j] = float( spl[2] )

        f.close()

        return (lnop, rnop, lrmm)

    @staticmethod
    def get_szpw(nband, norbs, fileName = None):
        """ try to read the solver.szpw.dat file to return the powers of
            local magnetization
        """
        if fileName is None:
            f = open("solver.szpw.dat","r")
        else:
            f = open(fileName,"r")

        szpw = numpy.zeros((4,norbs), dtype = numpy.float)
        for i in range(nband+1):
            f.readline() # skip one comment line
            for j in range(4):
                spl = f.readline().split()
                szpw[j,i] = float( spl[1] )
            f.readline() # skip two blank lines
            f.readline()

        f.close()

        return szpw

    @staticmethod
    def get_sp_t(nband, ntime, fileName = None):
        """ try to read the solver.sp_t.dat file to return the spin-spin
            correlation function < S_z(0) S_z(\tau) > data
        """
        if fileName is None:
            f = open("solver.sp_t.dat","r")
        else:
            f = open(fileName,"r")

        tmesh = numpy.zeros((ntime), dtype = numpy.float)
        schi = numpy.zeros((ntime), dtype = numpy.float)
        sp_t = numpy.zeros((ntime,nband), dtype = numpy.float)
        # read sp_t
        for i in range(nband):
            f.readline() # skip one comment line
            for j in range(ntime):
                spl = f.readline().split()
                sp_t[j,i] = float( spl[1] )
            f.readline() # skip two blank lines
            f.readline()
        f.readline() # skip one comment line
        # read schi
        for i in range(ntime):
            spl = f.readline().split()
            tmesh[i] = float( spl[0] )
            schi[i] = float( spl[1] )

        f.close()

        return (tmesh, schi, sp_t)

    @staticmethod
    def get_sp_w(nband, nbfrq, fileName = None):
        """ try to read the solver.sp_w.dat file to return the spin-spin
            correlation function data
        """
        if fileName is None:
            f = open("solver.sp_w.dat","r")
        else:
            f = open(fileName,"r")

        bmesh = numpy.zeros((nbfrq), dtype = numpy.float)
        sp_w = numpy.zeros((nbfrq,nband), dtype = numpy.float)
        # read sp_w
        for i in range(nband):
            f.readline() # skip one comment line
            for j in range(nbfrq):
                spl = f.readline().split()
                bmesh[j] = float( spl[0] )
                sp_w[j,i] = float( spl[1] )
            f.readline() # skip two blank lines
            f.readline()

        f.close()

        return (bmesh, sp_w)

    @staticmethod
    def get_ch_t(norbs, ntime, fileName = None):
        """ try to read the solver.ch_t.dat file to return the orbital-
            orbital correlation function < N_i(0) N_j(\tau) > data
        """
        if fileName is None:
            f = open("solver.ch_t.dat","r")
        else:
            f = open(fileName,"r")

        tmesh = numpy.zeros((ntime), dtype = numpy.float)
        cchi = numpy.zeros((ntime), dtype = numpy.float)
        ch_t = numpy.zeros((ntime,norbs,norbs), dtype = numpy.float)
        # read ch_t
        for i in range(norbs):
            for j in range(norbs):
                f.readline() # skip one comment line
                for k in range(ntime):
                    spl = f.readline().split()
                    ch_t[k,j,i] = float( spl[1] )
                f.readline() # skip two blank lines
                f.readline()
        f.readline() # skip one comment line
        # read cchi
        for i in range(ntime):
            spl = f.readline().split()
            tmesh[i] = float( spl[0] )
            cchi[i] = float( spl[1] )

        f.close()

        return (tmesh, cchi, ch_t)

    @staticmethod
    def get_ch_w(norbs, nbfrq, fileName = None):
        """ try to read the solver.ch_w.dat file to return the orbital-
            orbital correlation function data
        """
        if fileName is None:
            f = open("solver.ch_w.dat","r")
        else:
            f = open(fileName,"r")

        bmesh = numpy.zeros((nbfrq), dtype = numpy.float)
        ch_w = numpy.zeros((nbfrq,norbs,norbs), dtype = numpy.float)
        # read ch_w
        for i in range(norbs):
            for j in range(norbs):
                f.readline() # skip one comment line
                for k in range(nbfrq):
                    spl = f.readline().split()
                    bmesh[k] = float( spl[0] )
                    ch_w[k,j,i] = float( spl[1] )
                f.readline() # skip two blank lines
                f.readline()

        f.close()

        return (bmesh, ch_w)

    @staticmethod
    def get_g2ph(norbs, nffrq, nbfrq, fileName = None):
        """ try to read the solver.g2ph.dat file to return the two-particle
            Green's function data
        """
        if fileName is None:
            f = open("solver.g2ph.dat","r")
        else:
            f = open(fileName,"r")

        g2ph = numpy.zeros((nffrq,nffrq,nbfrq,norbs,norbs), dtype = numpy.complex)
        for m in range(norbs):
            for n in range(m+1):
                for k in range(nbfrq):
                    f.readline() # skip three comment lines
                    f.readline()
                    f.readline()
                    for j in range(nffrq):
                        for i in range(nffrq):
                            spl = f.readline().split()
                            g2ph[i,j,k,n,m] = float( spl[2] ) + 1j * float( spl[3] )
                            g2ph[i,j,k,m,n] = float( spl[2] ) + 1j * float( spl[3] )
                    f.readline() # skip two blank lines
                    f.readline()

        f.close()

        return g2ph

    @staticmethod
    def get_h2ph(norbs, nffrq, nbfrq, fileName = None):
        """ try to read the solver.h2ph.dat file to return the two-particle
            Green's function data
        """
        if fileName is None:
            f = open("solver.h2ph.dat","r")
        else:
            f = open(fileName,"r")

        h2ph = numpy.zeros((nffrq,nffrq,nbfrq,norbs,norbs), dtype = numpy.complex)
        for m in range(norbs):
            for n in range(m+1):
                for k in range(nbfrq):
                    f.readline() # skip three comment lines
                    f.readline()
                    f.readline()
                    for j in range(nffrq):
                        for i in range(nffrq):
                            spl = f.readline().split()
                            h2ph[i,j,k,n,m] = float( spl[2] ) + 1j * float( spl[3] )
                            h2ph[i,j,k,m,n] = float( spl[2] ) + 1j * float( spl[3] )
                    f.readline() # skip two blank lines
                    f.readline()

        f.close()

        return h2ph

    @staticmethod
    def get_v4ph(norbs, nffrq, nbfrq, fileName = None):
        """ try to read the solver.v4ph.dat file to return the two-particle
            vertex function data
        """
        if fileName is None:
            f = open("solver.v4ph.dat","r")
        else:
            f = open(fileName,"r")

        chic = numpy.zeros((nffrq,nffrq,nbfrq,norbs,norbs), dtype = numpy.complex)
        chig = numpy.zeros((nffrq,nffrq,nbfrq,norbs,norbs), dtype = numpy.complex)
        for m in range(norbs):
            for n in range(m+1):
                for k in range(nbfrq):
                    f.readline() # skip three comment lines
                    f.readline()
                    f.readline()
                    for j in range(nffrq):
                        for i in range(nffrq):
                            spl = f.readline().split()
                            chic[i,j,k,n,m] = float( spl[2] ) + 1j * float( spl[3] )
                            chic[i,j,k,m,n] = float( spl[2] ) + 1j * float( spl[3] )
                            chig[i,j,k,n,m] = float( spl[4] ) + 1j * float( spl[5] )
                            chig[i,j,k,m,n] = float( spl[4] ) + 1j * float( spl[5] )
                    f.readline() # skip two blank lines
                    f.readline()

        f.close()

        return (chic, chig)
