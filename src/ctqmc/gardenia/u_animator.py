#!/usr/bin/env python

##
##
## Introduction
## ============
##
## It is a python script. The purpose of this script is generate essential
## input file (solver.ctqmc.in) for the quantum impurity solver components.
## Note that you can not use it to control these codes.
##
## Usage
## =====
##
## see the document string
##
## Author
## ======
##
## This python script is designed, created, implemented, and maintained by
##
## Li Huang // email: lihuang.dmft@gmail.com
##
## History
## =======
##
## 03/28/2015 by li huang (created)
## 08/17/2015 by li huang (last modified)
##
##

import sys
import numpy
import matplotlib
matplotlib.use('TKAgg')
#matplotlib.use("pdf") # setup backend
import matplotlib.pyplot as plt
import matplotlib.animation as ani

sel_iter = 1

num_iter = 20
num_diag = 200
num_band = 1
num_orbs = 2
max_pair = 100
zero = 0.0
beta = 8.00
cmap = ['red','green','purple']

time_s = numpy.zeros((max_pair,num_orbs,num_diag,num_iter), dtype = numpy.float)
time_e = numpy.zeros((max_pair,num_orbs,num_diag,num_iter), dtype = numpy.float)
time_z = numpy.zeros((max_pair,num_orbs,num_diag,num_iter), dtype = numpy.float)
rank_t = numpy.zeros((num_orbs,num_diag,num_iter), dtype = numpy.int)

f = open('solver.diag.dat', 'r')
for iter in range(num_iter):
    for diag in range(num_diag):
        line = f.readline()
        line = f.readline()
        for orbs in range(num_orbs):
            line = f.readline().split()
            rank_t[orbs,diag,iter] = int(line[4])
            for pair in range( rank_t[orbs,diag,iter] ):
                line = f.readline().split()
                time_s[pair,orbs,diag,iter] = float(line[1]) 
                time_e[pair,orbs,diag,iter] = float(line[2])
        line = f.readline()
        line = f.readline()
        #print iter, num_iter
f.close()

fig = plt.figure()
time_text = plt.text(zero, 0.5, '', fontsize = 18)
plt.text(zero-0.1, -0.5, r'0', fontsize = 18)
plt.text(beta-0.1, -0.5, r'$\beta$', fontsize = 18)
plt.gca().set_aspect(2.0)
plt.gca().set_yticks([])
plt.plot([zero,zero],[-0.2,0.2], color = 'black', lw = 5)
plt.plot([beta,beta],[-0.2,0.2], color = 'black', lw = 5)
plt.plot([zero,beta],[-0.0,0.0], '--', color = 'black', lw = 1)

iter = sel_iter
def update_figure(diag):
    print 'diag:', diag
    for orbs in range(num_orbs):
        pair = rank_t[orbs,diag,iter]
        lines[2*orbs+0].set_data(time_s[0:pair,orbs,diag,iter], time_z[0:pair,orbs,diag,iter])
        lines[2*orbs+1].set_data(time_e[0:pair,orbs,diag,iter], time_z[0:pair,orbs,diag,iter])
        print 'orbs:', orbs, 'pair:', pair
        print 'c--->', time_s[0:pair,orbs,diag,iter]
        print 'd--->', time_e[0:pair,orbs,diag,iter]
    lines[2*num_orbs].set_text('snapshot: ' + str(diag))
    print ''
    return lines

lines = []
for orbs in range(num_orbs):
    c_l, = plt.plot([], [], 'o', ms = 8, mfc = cmap[orbs], mec = cmap[orbs], mew = 2, alpha = 0.8)
    lines.append(c_l)
    d_l, = plt.plot([], [], 'o', ms = 8, mfc = 'white',    mec = cmap[orbs], mew = 2, alpha = 0.8)
    lines.append(d_l)
lines.append(time_text)

line_ani = ani.FuncAnimation(fig, update_figure, num_diag, interval=500, repeat = False, blit=True)
line_ani.save('im.mp4')
