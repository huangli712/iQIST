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
## This shell script is designed, created, implemented, and maintained by
##
## Li Huang // email: huangli712@gmail.com
##
## History
## =======
##
## 11/10/2014 by li huang
##
##

# define my own ln function
function mln {
    name=$(echo $2 | tr '[:lower:]' '[:upper:]')
    if [ -e "$1" ]
    then
        echo "[$name]: found"
        ln -fs $1 $2.x
        echo "[$name]: setup OK"
    fi
}

# loop over the ct-qmc components
for component in azalea gardenia narcissus begonia lavender pansy manjushaka
do
    dir=$(echo ../src/ctqmc/$component/ctqmc)
    mln $dir $component
done

# loop over the hf-qmc components
for component in daisy
do
    dir=$(echo ../src/hfqmc/$component/hfqmc)
    mln $dir $component
done

# loop over the jasmine components
for component in jasmine
do
    dir=$(echo ../src/tools/jasmine/atomic)
    mln $dir $component
done

# loop over the hibiscus components
for component in entropy
do
    dir=$(echo ../src/tools/hibiscus/entropy1/$component)
    mln $dir $component
done
for component in maxent
do
    dir=$(echo ../src/tools/hibiscus/entropy2/$component)
    mln $dir $component
done
for component in sac
do
    dir=$(echo ../src/tools/hibiscus/stoch/$component)
    mln $dir $component
done
for component in mchi mdos mkra mscr msig mstd mtau mups
do
    dir=$(echo ../src/tools/hibiscus/toolbox/$component)
    mln $dir $component
done
