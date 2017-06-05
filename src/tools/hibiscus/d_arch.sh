#!/bin/bash

##
##
## Introduction
## ============
##
## It is a shell script. The purpose of this script is to create archive
## of files (in the working directory) from the current repo branch.
##
## The name of the output archive should like this:
##     iqist_43e2cbb_1441276643.tar.gz
## Here 43e2cbb is the abbreviated commit hash, and 1441276643 is the
## UNIX timestamp when this commit was committed.
##
## This script should be used by the developer only.
##
## Usage
## =====
##
## ./d_arch.sh
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
## 02/06/2015 by li huang (created)
## 06/05/2017 by li huang (last modified)
##
##

short_hash_tag=`git show -s --format='%h_%at' HEAD`
archive_name=iqist_$short_hash_tag.tar.gz
git archive -o $archive_name HEAD
