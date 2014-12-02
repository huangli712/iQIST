#!/usr/bin/env python

##
##
## Introduction
## ============
##
## It is a python script. The purpose of this script is generate essential
## input file (atom.config.in) for the jasmine code. Note that you can not
## use it to control the jasmine code.
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
## 12/02/2014 by li huang
##
##

import sys

class p_atomic_solver(object):
    """
    """

    def __init__(self):
        """
        """
        self._p_cmp = {
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

        self._p_inp = {}

    def setp(self, **kwargs):
        """
        """
        if kwargs is not None:
            for key, value in kwargs.iteritems():
                self._p_inp[key] = value

    def check(self):
        """
        """
        for key in self._p_inp.iterkeys():
            if key not in self._p_cmp:
                sys.exit('FATAL ERROR: wrong key ' + key)
            if type( self._p_inp[key] ) is not type( self._p_cmp[key] ):
                sys.exit('FATAL ERROR: wrong value ' + key + ' = ' + str(self._p_inp[key]))

    def write(self):
        """
        """
        f = open('atom.config.in','w')
        for key in self._p_inp.iterkeys():
            empty = ( 8 - len(key) ) * ' ' + ': '
            f.write(key + empty + str(self._p_inp[key]) + '\n')
        f.close()
