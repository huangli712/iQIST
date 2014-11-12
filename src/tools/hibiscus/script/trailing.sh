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
## ./trailing.sh file_name
##
## Author
## ======
##
## This shell script is designed, created, implemented, and maintained by
##
## Li Huang // email: huangli712@gmail.com
##
## History
## =======
##
## 11/13/2014 by li huang
##
##
sed -i '' -e's/[[:space:]]*$//' $1
