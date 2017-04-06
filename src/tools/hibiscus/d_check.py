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
## This python script is designed, created, implemented, and maintained by
##
## Li Huang // email: lihuang.dmft@gmail.com
##
## History
## =======
##
## 11/13/2014 by li huang (created)
## 08/17/2015 by li huang (last modified)
##
##

import sys

argu = sys.argv[1:]
f = file(argu[0])

i = 0 # line counter
while True:
    line = f.readline()
    if len(line) == 0:
        break
    i = i + 1
    if len(line) != len(line.rstrip()) + 1:
        print 'check space:', i, line,

f.close()
