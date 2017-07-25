#!/usr/bin/env python

##
##
## Introduction
## ============
##
## It is a python script. The purpose of this script is to check whether
## the lines in a given file are ended with blanks.
##
## This script should be used by the developer only.
##
## Usage
## =====
##
## ./d_check.py file_name
##
## Author
## ======
##
## This python script is designed, created, and maintained by
##
## Li Huang // email: lihuang.dmft@gmail.com
##
## History
## =======
##
## 11/13/2014 by li huang (created)
## 07/26/2017 by li huang (last modified)
##
##

import sys

# access the command-line argument
argu = sys.argv[1:]

with open(argu[0], 'r') as f:
    i = 0 # line counter
    while True:
        line = f.readline()
        if len(line) == 0:
            break
        i = i + 1
        if len(line) != len(line.rstrip()) + 1:
            print 'line number:', i
            print '--->', line.rstrip()
