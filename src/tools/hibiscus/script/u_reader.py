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
## 01/07/2015 by li huang
##
##

import os
import sys
import numpy

class iqistReader(object):
    """
    """

    @staticmethod
    def get_green():
        """ try to read the solver.green.dat or solver.green.bin.nnn file
            to return the imaginary time Green's function G(\tau) data
        """
        pass

    @staticmethod
    def get_grn():
        """ try to read the solver.grn.dat file to return the matsubara
            Green's function G(i\omega) data
        """
        pass

    @staticmethod
    def get_hybri():
        """ try to read the solver.hybri.dat file to return the imaginary
            time hybridization function \Delta(\tau) data
        """
        pass

    @staticmethod
    def get_hyb():
        """ try to read the solver.hyb.dat file to return the matsubara
            hybridization function \Delta(i\omega) data
        """
        pass

    @staticmethod
    def get_weiss():
        """ try to read the solver.weiss.dat file to return the imaginary
            time Weiss's function \mathcal{G}(\tau) data
        """
        pass

    @staticmethod
    def get_wss():
        """ try to read the solver.wss.dat file to return the matsubara
            Weiss's function \mathcal{G}(i\omega) data
        """
        pass

    @staticmethod
    def get_sgm():
        """ try to read the solver.sgm.dat file to return the matsubara
            self-energy function \Sigma(i\omega) data
        """
        pass

    @staticmethod
    def get_hub():
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
