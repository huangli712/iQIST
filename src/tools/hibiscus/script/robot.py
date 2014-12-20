#!/usr/bin/env python

##
##
## Introduction
## ============
##
## It is a python script. The purpose of this script is to perform some
## automatic tests for the ctqmc quantum impurity solvers.
##
## This script should be used by the developer only.
##
## Usage
## =====
##
## ./robot.py job_file
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
## 12/20/2014 by li huang
##
##

import os
import sys

# setup the directory we want to scan
argu = sys.argv[1:]
if ( len(argu) > 0 ):
    scan_dir = argu[0]
else:
    scan_dir = '.'
