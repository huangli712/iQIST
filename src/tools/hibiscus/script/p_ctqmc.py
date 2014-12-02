#!/usr/bin/env python

##
##
## Introduction
## ============
##
## Usage
## =====
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

class p_ctqmc_solver(object):
    """
    """

    def __init__(self, solver):
        """
        """
        pass

    def setp(self, aa):
        """
        """
        pass

    def check(self):
        """
        """
        pass

    def write(self):
        """
        write the parameters to the config file: solver.ctqmc.in
        """
        f = open('solver.ctqmc.in','w')
        for key in self.__p_inp.iterkeys():
            empty = ( 8 - len(key) ) * ' ' + ': '
            f.write(key + empty + str(self.__p_inp[key]) + '\n')
        f.close()
