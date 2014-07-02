#!/bin/bash

#-------------------------------------------------------------------------
# ctqmc_cleaner script is used to clean the testing farm.
# it is written in Bash shell.
#
# author  : li huang (email:huangli712@yahoo.com.cn)
# version : 0.0.3
# status  : unstable
#-------------------------------------------------------------------------

# define testing farm
#-------------------------------------------------------------------------
# testing case list
wlst='1* 2* 3* 4* 5* 6* 7*'

# log file
clog=clean.log
#-------------------------------------------------------------------------

for directory in $wlst
do

# enter correct directory
   cd $directory; pwd

# cleaning current directory
   rm -f *.log *.out *.dat *.bin* .DS_Store

# print successful information
   echo 'job' $directory 'is cleaned' >> ../$clog

# goto the top directory
   cd ..; pwd

done
