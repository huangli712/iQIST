#!/usr/bin/env python

##
##
## Introduction
## ============
##
## It is a python script. The purpose of this script is generate essential
## input file (solver.ctqmc.in) for the quantum impurity solver components.
## Note that you can not use it to control these codes.
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

class p_ctqmc_solver(object):
    """
    """

    def __init__(self, solver):
        """ define the class variables
        """
        self.__p_cmp_azalea = {}
        self.__p_cmp_gardenia = {}
        self.__p_cmp_narcissus = {}
        self.__p_cmp_begonia = {}
        self.__p_cmp_lavender = {}
        self.__p_cmp_pansy = {}
        self.__p_cmp_manjushaka = {}

        # __p_cmp: the official parameter dict
        self.__p_cmp = {}

        # __p_inp: the user-input parameter dict
        self.__p_inp = {}

        for case in switch( solver.lower() ):
            if case ('azalea'):
                self.__p_cmp = self.__p_cmp_azalea
                break

            if case ('gardenia'):
                self.__p_cmp = self.__p_cmp_gardenia
                break

            if case ('narcissus'):
                self.__p_cmp = self.__p_cmp_narcissus
                break

            if case ('begonia'):
                self.__p_cmp = self.__p_cmp_begonia
                break

            if case ('lavender'):
                self.__p_cmp = self.__p_cmp_lavender
                break

            if case ('pansy'):
                self.__p_cmp = self.__p_cmp_pansy
                break

            if case ('manjushaka'):
                self.__p_cmp = self.__p_cmp_manjushaka
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
