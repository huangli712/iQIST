#!/bin/bash

##
##
## Introduction
## ============
##
## It is a bash shell script. The purpose of this script is to make symbolic
## links for the executable codes in the bin directory (current directory). 
##
## Usage
## =====
##
## ./setup.sh
##
## Author
## ======
##
## This shell script is designed, created, implemented, and maintained
## by Li Huang (email: huangli712@gmail.com).
##
## Department of Physics, University of Fribourg, Fribourg CH-1700, Switzerland
##
## History
## =======
##
## 09/17/2014 by li huang
##
##

# define my ln function
function mln {
    name=$(echo $2 | tr '[:lower:]' '[:upper:]')
    if [ -e "$1" ]
    then
        echo "[$name]: found"
        ln -fs $1 $2.x
        echo "[$name]: setup OK"
    fi
}

# loop over the components
for solver in azalea gardenia narcissus begonia lavender pansy manjushaka
do
    dir=$(echo ../src/ctqmc/$solver/ctqmc)
    mln $dir $solver
done
