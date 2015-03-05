#!/usr/bin/env python

import sys
import numpy
import matplotlib.pyplot as plt

sys.path.append('../..//tools/hibiscus/script')
from u_reader import *

#norbs = 2
#nbfrq = 4
#nffrq = 48
#beta = 10.0

#g2, f2 = iqistReader.get_twop(norbs, nffrq, nbfrq)
#data = f2[:,28,2,0,0]
#x = numpy.zeros(nffrq, dtype = numpy.float)
#for i in range(nffrq):
#    x[i] = (2*float(i+1) - nffrq - 1) * numpy.pi / beta
#    print i,x[i]
#
#lines = plt.plot(x, data.real, alpha = 0.5, clip_on = True)
#plt.xlim(-20,20)
#plt.show()

#norbs = 4
#nband = 2
#ntime = 1024
#beta = 50.0
#schi, sschi = iqistReader.get_schi(nband, ntime)
#
#schi = schi * nband / 4.0
#x = numpy.zeros(ntime, dtype = numpy.float)
#for i in range(ntime):
#    x[i] = i * beta / float(ntime - 1)
#lines = plt.plot(x, schi, alpha = 0.5, clip_on = True)
#plt.xlim(0,25)
#plt.ylim(0,0.9)
#plt.show()

norbs = 4
nbfrq = 2
nffrq = 48
beta = 50.0

g2, f2 = iqistReader.get_twop(norbs, nffrq, nbfrq)
data_1_1 = f2[:,23,1,0,0]
data_2_2 = f2[:,23,1,1,1]
data_3_3 = f2[:,23,1,2,2]
data_4_4 = f2[:,23,1,3,3]
data_1_3 = f2[:,23,1,0,2]
data_1_2 = f2[:,23,1,0,1]
data_1_4 = f2[:,23,1,0,3]
x = numpy.zeros(nffrq, dtype = numpy.float)
for i in range(nffrq):
    x[i] = (2*float(i+1) - nffrq - 1) * numpy.pi / beta
    print i,x[i]

lines = plt.plot(x, (data_1_1 + data_1_3).real, x, (data_1_1 - data_1_3).real, alpha = 0.5, clip_on = True)
#lines = plt.plot(x, (data_1_1 + data_1_3).imag, x, (data_1_1 - data_1_3).imag, alpha = 0.5, clip_on = True)
#lines = plt.plot(x, (data_1_2 + data_1_4).real, x, (data_1_2 - data_1_4).real, alpha = 0.5, clip_on = True)
#lines = plt.plot(x, (data_1_2 + data_1_4).imag, x, (data_1_2 - data_1_4).imag, alpha = 0.5, clip_on = True)
plt.xlim(-3,3)
plt.ylim(-15,45)
plt.ylim(-40,40)
plt.show()
