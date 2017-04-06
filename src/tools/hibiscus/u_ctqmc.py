#!/usr/bin/env python

##
##
## Introduction
## ============
##
## It is a python script. The purpose of it is to generate essential input
## file (solver.ctqmc.in) for various quantum impurity solver components.
## Note that you can not use it to control these codes directly.
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
## Li Huang // email: lihuang.dmft@gmail.com
##
## History
## =======
##
## 03/28/2015 by li huang (created)
## 04/06/2017 by li huang (last modified)
##
##

import sys

class p_ctqmc_solver(object):
    """ This class can be used to generate the config file for the quantum
        impurity solver components.

        typical usage:
            # import this module
            from u_ctqmc import *

            # create an instance
            p = p_ctqmc_solver('manjushaka')

            # setup the parameters
            p.setp(isscf = 2, isort = 1, nsweep = 10000000)
            p.setp(mune = 2.0)
            p.setp()
            p.setp(isscf = 1)

            # verify the parameters
            p.check()

            # generate the solver.ctqmc.in file
            p.write()

            # destroy the instance
            del p
    """

    def __init__(self, solver):
        """ define the class variables
        """
        # __p_cmp_solver: the official parameter dict for generic solver
        self.__p_cmp_solver = {
            'isscf'  : 2       ,
            'issun'  : 2       ,
            'isspn'  : 1       ,
            'isbin'  : 2       ,
            'nband'  : 1       ,
            'nspin'  : 2       ,
            'norbs'  : 2       ,
            'ncfgs'  : 4       ,
            'niter'  : 20      ,
            'mkink'  : 1024    ,
            'mfreq'  : 8193    ,
            'nfreq'  : 128     ,
            'ntime'  : 1024    ,
            'nflip'  : 20000   ,
            'ntherm' : 200000  ,
            'nsweep' : 20000000,
            'nwrite' : 2000000 ,
            'nclean' : 100000  ,
            'nmonte' : 10      ,
            'ncarlo' : 10      ,
            'U'      : 4.00    ,
            'Uc'     : 4.00    ,
            'Uv'     : 4.00    ,
            'Jz'     : 0.00    ,
            'Js'     : 0.00    ,
            'Jp'     : 0.00    ,
            'mune'   : 2.00    ,
            'beta'   : 8.00    ,
            'part'   : 0.50    ,
            'alpha'  : 0.70    ,
        }

        # __p_cmp_narcissus: the official parameter dict for narcissus
        self.__p_cmp_narcissus = self.__p_cmp_solver.copy()
        self.__p_cmp_narcissus['isort'] = 1
        self.__p_cmp_narcissus['issus'] = 1
        self.__p_cmp_narcissus['isvrt'] = 1
        self.__p_cmp_narcissus['isscr'] = 1
        self.__p_cmp_narcissus['lemax'] = 32
        self.__p_cmp_narcissus['legrd'] = 20001
        self.__p_cmp_narcissus['chmax'] = 32
        self.__p_cmp_narcissus['chgrd'] = 20001
        self.__p_cmp_narcissus['nffrq'] = 32
        self.__p_cmp_narcissus['nbfrq'] = 8

        # __p_cmp_manjushaka: the official parameter dict for manjushaka
        self.__p_cmp_manjushaka = self.__p_cmp_solver.copy()
        self.__p_cmp_manjushaka['isort'] = 1
        self.__p_cmp_manjushaka['issus'] = 1
        self.__p_cmp_manjushaka['isvrt'] = 1
        self.__p_cmp_manjushaka['ifast'] = 1
        self.__p_cmp_manjushaka['itrun'] = 1
        self.__p_cmp_manjushaka['lemax'] = 32
        self.__p_cmp_manjushaka['legrd'] = 20001
        self.__p_cmp_manjushaka['chmax'] = 32
        self.__p_cmp_manjushaka['chgrd'] = 20001
        self.__p_cmp_manjushaka['nffrq'] = 32
        self.__p_cmp_manjushaka['nbfrq'] = 8
        self.__p_cmp_manjushaka['npart'] = 4

        # __p_cmp: the official parameter dict
        self.__p_cmp = {}

        # __p_inp: the user-input parameter dict
        self.__p_inp = {}

        # config __p_cmp according to the selected solver
        for case in switch( solver.lower() ):
            if case ('narcissus'):
                self.__p_cmp = self.__p_cmp_narcissus.copy()
                break

            if case ('manjushaka'):
                self.__p_cmp = self.__p_cmp_manjushaka.copy()
                break

            if case ():
                sys.exit('FATAL ERROR: unrecognize quantum impurity solver')

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
        """ write the parameters to the config file: solver.ctqmc.in
        """
        f = open('solver.ctqmc.in','w')
        for key in self.__p_inp.iterkeys():
            empty = ( 8 - len(key) ) * ' ' + ': '
            f.write(key + empty + str(self.__p_inp[key]) + '\n')
        f.close()

class switch(object):
    """ This class provides the functionality we want. You only need to
        look at this if you want to know how this works. It only needs to
        be defined once, no need to muck around with its internals.
    """
    def __init__(self, value):
        """ class constructor for the switch
        """
        self.__value = value
        self.__fall = False

    def __iter__(self):
        """ return the match method once, then stop
        """
        yield self.match
        raise StopIteration

    def match(self, *args):
        """ indicate whether or not to enter a case suite
        """
        if self.__fall or not args:
            return True
        elif self.__value in args:
            self.__fall = True
            return True
        else:
            return False
