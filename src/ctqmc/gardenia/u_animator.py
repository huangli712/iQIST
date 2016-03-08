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
## Li Huang // email: lihuang.dmft@gmail.com
##
## History
## =======
##
## 03/28/2015 by li huang (created)
## 08/17/2015 by li huang (last modified)
##
##

import sys
import numpy

num_iter = 20
num_diag = 200
num_band = 1
num_orbs = 2
max_pair = 100

time_s = numpy.zeros((max_pair,num_orbs,num_diag,num_iter), dtype = numpy.float)

f = open('solver.diag.dat', 'r')
f.close()
