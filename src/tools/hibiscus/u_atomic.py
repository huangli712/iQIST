#!/usr/bin/env python

##
##
## Introduction
## ============
##
## It is a python script. The purpose of it is to generate essential input
## file (atom.config.in) for the jasmine code. Note that you can not use
## it to control the jasmine code directly.
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
## 12/02/2014 by li huang (created)
## 04/06/2017 by li huang (last modified)
##
##

import sys

class p_atomic_solver(object):
    """ This class can be used to generate the config file for the jasmine.

        typical usage:
            # import this module
            from u_atomic import *

            # create an instance
            p = p_atomic_solver()

            # setup the parameters
            p.setp(ibasis = 2, Uv = 2.0)
            p.setp(icu = 30) # invalid parameter
            p.setp(icu = 1)
            p.setp()

            # verify the parameters
            p.check()

            # generate the atom.config.in file
            p.write()

            # destroy the instance
            del p
    """

    def __init__(self):
        """ define the class variables
        """
        # __p_cmp: the official parameter dict
        # here the default values are just used to verify the data type of
        # user's input data
        self.__p_cmp = {
            'ibasis' : 1   ,
            'ictqmc' : 1   ,
            'icu'    : 1   ,
            'icf'    : 0   ,
            'isoc'   : 0   ,
            'nband'  : 1   ,
            'nspin'  : 2   ,
            'norbs'  : 2   ,
            'ncfgs'  : 4   ,
            'nmini'  : 0   ,
            'nmaxi'  : 2   ,
            'Uc'     : 2.00,
            'Uv'     : 2.00,
            'Jz'     : 0.00,
            'Js'     : 0.00,
            'Jp'     : 0.00,
            'Ud'     : 2.00,
            'Jh'     : 0.00,
            'mune'   : 0.00,
            'lambda' : 0.00,
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
        """ write the parameters to the config file: atom.config.in
        """
        f = open('atom.config.in','w')
        for key in self.__p_inp.iterkeys():
            empty = ( 8 - len(key) ) * ' ' + ': '
            f.write(key + empty + str(self.__p_inp[key]) + '\n')
        f.close()
