#!/bin/bash

##
##
## Introduction
## ============
##
## It is a shell script. The purpose of this script is to remove the
## trailing whitespaces in the given file.
##
## This script should be used by the developer only.
##
## Usage
## =====
##
## ./d_trim.sh file_name
##
## For macOS system, the grammar for sed is (we don't generate backup)
##     sed -i '' ...
##
## However, for Linux-based system, the grammar for sed is
##     sed -i ...
##
## Author
## ======
##
## This shell script is designed, created, and maintained by
##
## Li Huang // email: lihuang.dmft@gmail.com
##
## History
## =======
##
## 11/13/2014 by li huang (created)
## 06/05/2017 by li huang (last modified)
##
##

sed -i '' -e's/[[:space:]]*$//' $1
