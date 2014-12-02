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

    def setp(self, **kwargs):
        """
        setup the parameters using a series of key-value pairs
        """
        if kwargs is not None:
            for key, value in kwargs.iteritems():
                self.__p_inp[key] = value

    def check(self):
        """
        check the correctness of input parameters
        """
        for key in self.__p_inp.iterkeys():
            # check whether the key is valid
            if key not in self.__p_cmp:
                sys.exit('FATAL ERROR: wrong key ' + key)
            # check the data type of key's value
            if type( self.__p_inp[key] ) is not type( self.__p_cmp[key] ):
                sys.exit('FATAL ERROR: wrong value ' + key + ' = ' + str(self.__p_inp[key]))

    def write(self):
        """
        write the parameters to the config file: solver.ctqmc.in
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
