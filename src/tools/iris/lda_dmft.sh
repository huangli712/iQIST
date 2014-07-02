#!/bin/bash

#-------------------------------------------------------------------------
# lda_dmft script is used to boost the lda + dmft self-consistent loop.
# it is written in Bash shell.
#
# author  : li huang (email:huangli712@yahoo.com.cn)
# version : 0.0.8
# status  : unstable
#-------------------------------------------------------------------------

# define color escape sequence for terminal system
#-------------------------------------------------------------------------
BLACK='\033[30;1m'
RED='\033[31;1m'
GREEN='\033[32;1m'
YELLOW='\033[33;1m'
BLUE='\033[34;1m'
MAGENTA='\033[35;1m'
CYAN='\033[36;1m'
WHITE='\033[37;1m'
RESET='\033[0m'
#-------------------------------------------------------------------------

# define root directory
#-------------------------------------------------------------------------
ROOT=`pwd`
#-------------------------------------------------------------------------

# define mpi running envirnoment
#-------------------------------------------------------------------------
# mpi running driver
MPID=mpiexec

# number of processors used by lda + dmft engine
NUMP=4

# number of processors used by quantum impurity solver
NUMQ=8
#-------------------------------------------------------------------------

# define exectuable code
#-------------------------------------------------------------------------
# lda + dmft engine based on maximally localized wannier function
WANN=~/Work/sunset-forge/wann

# continuous time quantum monte carlo quantum impurity solver
CTQMC=~/Work/azalea-forge/ctqmc
#-------------------------------------------------------------------------

# define control flag for lda + dmft self-consistent loop
#-------------------------------------------------------------------------
# maximum allowable lda + dmft iteration number
NITER=20

# once the iteration number reach CPIMP, quantum impurity solver should
# provide impurity occupation number for lda + dmft engine
CPIMP=999
#-------------------------------------------------------------------------

# start lda + dmft self-consistent loop
for ((i=1; i<=$NITER; ++i))
do

# print iteration information
    echo -ne $RED ">>> LDA + DMFT ITER <<<" $RESET $i; echo

# call lda + dmft engine
    echo -ne $GREEN "< STAGE 1 > LDA + DMFT Engine: WANN" $RESET; echo
    cd $ROOT/dmft
    echo -ne $GREEN "current directory:" $RESET $YELLOW `pwd` $RESET; echo
    $MPID -n $NUMP $WANN < /dev/null
    cd ..
    echo -ne $GREEN "current directory:" $RESET $YELLOW `pwd` $RESET; echo

# prepare necessary input data for quantum impurity solver
    cp $ROOT/dmft/dmft.hybf.dat $ROOT/imp/solver.hyb.in
    cp $ROOT/dmft/dmft.eimp.dat $ROOT/imp/solver.eimp.in

# special treatment for self-energy function
    if [ -e $ROOT/dmft/dmft.smat.in ]
    then
        cp $ROOT/dmft/dmft.smat.in $ROOT/dmft/dmft.smat.it
    fi

# print two blank lines
    echo
    echo

# call quantum impurity solver
    echo -ne $GREEN "< STAGE 2 > Quantum Impurity Solver: CTQMC" $RESET; echo
    cd $ROOT/imp
    echo -ne $GREEN "current directory:" $RESET $YELLOW `pwd` $RESET; echo
    $MPID -n $NUMQ $CTQMC < /dev/null
    cd ..
    echo -ne $GREEN "current directory:" $RESET $YELLOW `pwd` $RESET; echo

# prepare necessary input data for lda + dmft engine
    cp $ROOT/imp/solver.hyb.dat $ROOT/dmft/dmft.hmat.in
    cp $ROOT/imp/solver.grn.dat $ROOT/dmft/dmft.gmat.in
    cp $ROOT/imp/solver.sgm.dat $ROOT/dmft/dmft.smat.in

# special treatment for impurity occupation number
    if [ $i -ge $CPIMP ]
    then
        cp $ROOT/imp/solver.nmat.dat $ROOT/dmft/dmft.zimp.in
    fi

# print two blank lines
    echo
    echo

done
#<<< finish lda + dmft self-consistent loop
