#!/bin/bash

#-------------------------------------------------------------------------
# ctqmc_tester script is used to boost the testing farm.
# it is written in Bash shell.
#
# author  : li huang (email:huangli712@yahoo.com.cn)
# version : 0.0.5
# status  : unstable
#-------------------------------------------------------------------------

# define mpi running envirnoment
#-------------------------------------------------------------------------
# mpi running driver
mpid=mpiexec

# number of processors used by quantum impurity solver
# please choose appropriate number
nump=8
#-------------------------------------------------------------------------

# define exectuable code
#-------------------------------------------------------------------------
# continuous time quantum monte carlo quantum impurity solver
# please setup it carefully
ctqmc=~/Work/azalea-forge/ctqmc
#-------------------------------------------------------------------------

# define testing farm
#-------------------------------------------------------------------------
# testing case list
# please check the necessary input file before start formal testing
wlst='1* 2* 3* 4* 5* 6* 7*'

# log file
clog=ctqmc.log
#-------------------------------------------------------------------------

for directory in $wlst
do

# enter correct directory
   cd $directory; pwd

# cleaning current directory
   rm -f *.log *.out *.dat *.bin* .DS_Store

# power on the quantum impurity solver
   $mpid -n $nump $ctqmc < /dev/null > ctqmc.out

# print successful information
   echo 'job' $directory 'is done' >> ../$clog

# goto the top directory
   cd ..; pwd

done
