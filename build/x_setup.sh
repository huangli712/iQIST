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
## ./x_setup.sh
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
## 11/10/2014 by li huang (created)
## 03/30/2017 by li huang (last modified)
##
##

# define my own ln function
function make_link {
    name=$(echo $2 | tr '[:lower:]' '[:upper:]')
    if [ -e "$1" ]
    then
        echo "[$name]: found"
        ln -fs $1 $2.x
        echo "[$name]: setup OK"
    fi
}

# loop over the ctqmc components
for component in narcissus manjushaka
do
    dir=$(echo ../src/ctqmc/$component/ctqmc)
    make_link $dir $component
done

# loop over the jasmine components
for component in jasmine
do
    dir=$(echo ../src/tools/jasmine/atomic)
    make_link $dir $component
done

# loop over the hibiscus components
for component in mchi mdos mscr
do
    dir=$(echo ../src/tools/hibiscus/toolbox/$component)
    make_link $dir $component
done
