#!/usr/bin/env python

##
##
## Introduction
## ============
##
## It is a python script. The purpose of this script is to clean the dull
## .DS_Store files for the Mac OS X filesystem.
##
## This script should be used by the developer only.
##
## Usage
## =====
##
## ./d_clean.py directory_name
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
## 08/17/2015 by li huang (last modified)
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

# scan the directory and delete all of the found .DS_Store files directly
for path, dirs, files in os.walk(scan_dir):
    for f in files:
        if f.endswith('DS_Store'):
            print path + '/' + f
            os.remove(path + '/' + f)
