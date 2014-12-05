#!/usr/bin/env python

##
##
## Introduction
## ============
##
## It is a python script. The purpose of this script is generate essential
## input file (solver.hfqmc.in) for the daisy code. Note that you can not
## use it to control the daisy code.
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
## 12/05/2014 by li huang
##
##

import sys

class p_hfqmc_solver(object):
    """ This class can be used to generate the config file for the quantum
        impurity solver daisy code.

        typical usage:
            # import this module
            from p_hfqmc import *

            # create an instance
            p = p_hfqmc_solver()

            # setup the parameters
            p.setp(isscf = 2, nsweep = 10000000)
            p.setp(mune = 2.0, mstep = 10)
            p.setp()
            p.setp(isscf = 1)

            # verify the parameters
            p.check()

            # generate the solver.hfqmc.in file
            p.write()

            # destroy the instance
            del p
    """

    def __init__(self, solver):
        """ define the class variables
        """
        # __p_cmp: the official parameter dict
        # here the default values are just used to verify the data type of
        # user's input data
        self.__p_cmp = {
            'isscf'  : 2       ,
            'issun'  : 2       ,
            'isspn'  : 1       ,
            'isbin'  : 2       ,
            'nband'  : 1       ,
            'nspin'  : 2       ,
            'norbs'  : 2       ,
            'niter'  : 20      ,
            'mstep'  : 16      ,
            'mfreq'  : 8193    ,
            'nsing'  : 1       ,
            'ntime'  : 64      ,
            'ntherm' : 100     ,
            'nsweep' : 240000  ,
            'nclean' : 100     ,
            'ncarlo' : 10      ,
            'Uc'     : 4.00    ,
            'Jz'     : 0.00    ,
            'mune'   : 2.00    ,
            'beta'   : 8.00    ,
            'part'   : 0.50    ,
            'alpha'  : 0.70    ,
        }

        # __p_inp: the user-input parameter dict
        self.__p_inp = {}

    def setp(self, **kwargs):
        """ setup the parameters using a series of key-value pairs
        """
        if kwargs is not None:
            for key, value in kwargs.iteritems():
                self.__p_inp[key] = value

    def check(self):
        """ check the correctness of input parameters
        """
        for key in self.__p_inp.iterkeys():
            # check whether the key is valid
            if key not in self.__p_cmp:
                sys.exit('FATAL ERROR: wrong key ' + key)
            # check the data type of key's value
            if type( self.__p_inp[key] ) is not type( self.__p_cmp[key] ):
                sys.exit('FATAL ERROR: wrong value ' + key + ' = ' + str(self.__p_inp[key]))

    def write(self):
        """ write the parameters to the config file: solver.hfqmc.in
        """
        f = open('solver.hfqmc.in','w')
        for key in self.__p_inp.iterkeys():
            empty = ( 8 - len(key) ) * ' ' + ': '
            f.write(key + empty + str(self.__p_inp[key]) + '\n')
        f.close()
