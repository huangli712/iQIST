#!/usr/bin/env python

##
##
## Introduction
## ============
##
## It is a python script. The purpose of this script is to generate the
## animation movie using the data contained in the solver.diag.dat.
##
## Usage
## =====
##
## edit the configuration parameter carefully, and then execute
##
## ./u_movie.py movie.mp4
##
## Here movie.mp4 is the output file. We can use the VLC to play it. If you
## don't supply any filename, the default output should be diag.mp4.
##
## Author
## ======
##
## This python script is designed, created, and maintained by
##
## Li Huang // email: lihuang.dmft@gmail.com
##
## History
## =======
##
## 03/28/2015 by li huang (created)
## 06/05/2017 by li huang (last modified)
##
##

import sys
import numpy
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.animation as ani

# configuration parameter
# you have to fix them carefully before you try to execute this script
zero = 0.0      # zero constant
beta = 8.00     # inverse temperature (it = beta)
#               #
sel_iter = 0    # which animation movie is produced (0 <= it <= num_iter - 1)
num_iter = 20   # number of animation movies (it = niter)
num_diag = 200  # number of frames per animation (it = nsweep / nwrite) 
num_orbs = 2    # number of orbitals (it = norbs)
max_pair = 128  # maximum number of operator pairs per orbital

# setup colormap
# we use colors to distinguish the operators of various orbitals
cmap = ['red',         # 01
        'sienna',      # 02
        'tan',         # 03
        'orange',      # 04
        'darkkhaki',   # 05
        'greenyellow', # 06
        'green',       # 07
        'teal',        # 08
        'skyblue',     # 09
        'slategrey',   # 10
        'blue',        # 11
        'purple',      # 12
        'magenta',     # 13
        'hotpink']     # 14

# creator operators (x points)
time_s = numpy.zeros((max_pair,num_orbs,num_diag,num_iter), dtype = numpy.float)

# destroy operators (x points)
time_e = numpy.zeros((max_pair,num_orbs,num_diag,num_iter), dtype = numpy.float)

# y points (its elements are all zero)
time_z = numpy.zeros((max_pair,num_orbs,num_diag,num_iter), dtype = numpy.float)

# numboer of operator pairs
rank_t = numpy.zeros((num_orbs,num_diag,num_iter), dtype = numpy.int)

# read the data
f = open('solver.diag.dat', 'r')
for iter in range(num_iter):
    for diag in range(num_diag):
        line = f.readline()
        line = f.readline()
        for orbs in range(num_orbs):
            line = f.readline().split()
            rank_t[orbs,diag,iter] = int(line[4])
            if rank_t[orbs,diag,iter] >= max_pair: 
                print 'ERR: max_pair is too small'
                sys.exit(-1)
            for pair in range( rank_t[orbs,diag,iter] ):
                line = f.readline().split()
                time_s[pair,orbs,diag,iter] = float(line[1]) 
                time_e[pair,orbs,diag,iter] = float(line[2])
        line = f.readline()
        line = f.readline()
f.close()

# create the figure
fig = plt.figure()

# setup static imagine objects
# sometimes you have to adjust the below parameters to obtain a better visualization
plt.gca().set_aspect(2.0)
plt.gca().set_yticks([])
plt.axis([zero, beta, -0.2, 0.2])
plt.text(zero-0.1, 0.3, r'0', fontsize = 18)
plt.text(beta-0.1, 0.3, r'$\beta$', fontsize = 18)
plt.plot([zero,zero], [-0.2,0.2], color = 'black', lw = 5)
plt.plot([beta,beta], [-0.2,0.2], color = 'black', lw = 5)
plt.plot([zero,beta], [-0.0,0.0], '--', color = 'black', lw = 1)
for orbs in range(num_orbs):
    plt.plot([orbs*beta/num_orbs], [-0.8], 'o', ms = 8, mfc = cmap[orbs], mec = cmap[orbs], mew = 2, clip_on = False)
    plt.plot([orbs*beta/num_orbs], [-1.0], 'o', ms = 8, mfc = 'white',    mec = cmap[orbs], mew = 2, clip_on = False)

# build the dynamic imagine objects
lines = []
for orbs in range(num_orbs):
    c_l, = plt.plot([], [], 'o', ms = 8, mfc = cmap[orbs], mec = cmap[orbs], mew = 2, alpha = 0.8)
    lines.append(c_l)
    d_l, = plt.plot([], [], 'o', ms = 8, mfc = 'white',    mec = cmap[orbs], mew = 2, alpha = 0.8)
    lines.append(d_l)
time_text = plt.text(zero, 0.6, '', fontsize = 18)
lines.append(time_text)

# select which movie should be produced
iter = sel_iter

# draw the figure
def draw_fig(diag):
    print 'diag:', diag
    for orbs in range(num_orbs):
        pair = rank_t[orbs,diag,iter]
        print 'orbs:', orbs, 'pair:', pair
        lines[2*orbs+0].set_data(time_s[0:pair,orbs,diag,iter], time_z[0:pair,orbs,diag,iter]+0.1)
        lines[2*orbs+1].set_data(time_e[0:pair,orbs,diag,iter], time_z[0:pair,orbs,diag,iter]-0.1)
    lines[2*num_orbs].set_text('snapshot: ' + str(diag) + '/' + str(num_diag))
    print ''
    return lines

# produce the animation
movie = ani.FuncAnimation(fig, draw_fig, num_diag, interval=500, blit=True)

# setup the filename we want to use
argu = sys.argv[1:]
if ( len(argu) > 0 ):
    mov_name = argu[0]
else:
    mov_name = 'diag.mp4'

# save the movie file
movie.save(mov_name)
