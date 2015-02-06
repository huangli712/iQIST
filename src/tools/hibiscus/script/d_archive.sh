#!/bin/bash

##
##
## Introduction
## ============
##
## It is a shell script. The purpose of this script is to create archive
## of files (in the working directory) from the current repo branch.
##
## This script should be used by the developer only.
##
## Usage
## =====
##
## ./d_archive.sh
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
## 02/06/2015 by li huang
##
##

git archive -o latest.tar.gz HEAD
