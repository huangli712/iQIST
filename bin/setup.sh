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
## 08/09/2014 by li huang
##
##

# azalea code
if [ -e "../src/ctqmc/azalea/ctqmc" ]
then
    echo "[AZALEA]: found"
    ln -fs ../src/ctqmc/azalea/ctqmc azalea.x
    echo "[AZALEA]: setup OK"
fi

# gardenia code
if [ -e "../src/ctqmc/gardenia/ctqmc" ]
then
    echo "[GARDENIA]: found"
    ln -fs ../src/ctqmc/gardenia/ctqmc gardenia.x
    echo "[GARDENIA]: setup OK"
fi

# narcissus code
if [ -e "../src/ctqmc/narcissus/ctqmc" ]
then
    echo "[NARCISSUS]: found"
    ln -fs ../src/ctqmc/narcissus/ctqmc narcissus.x
    echo "[NARCISSUS]: setup OK"
fi

# begonia code
if [ -e "../src/ctqmc/begonia/ctqmc" ]
then
    echo "[BEGONIA]: found"
    ln -fs ../src/ctqmc/begonia/ctqmc begonia.x
    echo "[BEGONIA]: setup OK"
fi

# lavender code
if [ -e "../src/ctqmc/lavender/ctqmc" ]
then
    echo "[LAVENDER]: found"
    ln -fs ../src/ctqmc/lavender/ctqmc lavender.x
    echo "[LAVENDER]: setup OK"
fi

# pansy code
if [ -e "../src/ctqmc/pansy/ctqmc" ]
then
    echo "[PANSY]: found"
    ln -fs ../src/ctqmc/pansy/ctqmc pansy.x
    echo "[PANSY]: setup OK"
fi

# manjushaka code
if [ -e "../src/ctqmc/manjushaka/ctqmc" ]
then
    echo "[MANJUSHAKA]: found"
    ln -fs ../src/ctqmc/manjushaka/ctqmc manjushaka.x
    echo "[MANJUSHAKA]: setup OK"
fi
