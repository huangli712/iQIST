#!/bin/bash

#-------------------------------------------------------------------------
# ctqmc_ploter script is used to compare and verify those results from
# testing farm quickly, it is written in Bash shell.
#
# author  : li huang (email:huangli712@yahoo.com.cn)
# version : 0.0.3
# status  : unstable
#-------------------------------------------------------------------------

# source directory
proj_a=azalea

# destination directory
proj_b=begonia

# compare directory
comp_d=1

for i in ../$proj_a/$comp_d*
do

# get rid of the ../path/ substring
    list=${i/..\/*\//}; echo $list

# dump the plot script for gnuplot
# for impurity green's function in matsubara frequency
    cat > grnf.gp << !
        set terminal postscript enhanced color solid
        set output "test1.eps"
        plot [0:10]'../$proj_a/$list/solver.grn.dat' u 2:4 w lp,'../$proj_b/$list/solver.grn.dat' u 2:4 w lp
!

# for self-energy function in matsubara frequency
    cat > sigf.gp << !
        set terminal postscript enhanced color solid
        set output "test2.eps"
        plot [0:10]'../$proj_a/$list/solver.sgm.dat' u 2:4 w lp,'../$proj_b/$list/solver.sgm.dat' u 2:4 w lp
!

# for impurity green's function in imaginary time
    cat > gtau.gp << !
        set terminal postscript enhanced color solid
        set output "test3.eps"
        plot '../$proj_a/$list/solver.green.dat' u 2:4 w lp,'../$proj_b/$list/solver.green.dat' u 2:4 w lp
!

# draw it
   gnuplot grnf.gp
   gnuplot sigf.gp
   gnuplot gtau.gp

# copy the figure
   mv test1.eps $list.grnf.eps
   mv test2.eps $list.sigf.eps
   mv test3.eps $list.gtau.eps

done
