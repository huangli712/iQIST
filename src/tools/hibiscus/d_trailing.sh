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
## ./d_trailing.sh file_name
##
## Author
## ======
##
## This shell script is designed, created, implemented, and maintained by
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

sed -i '' -e's/[[:space:]]*$//' $1
