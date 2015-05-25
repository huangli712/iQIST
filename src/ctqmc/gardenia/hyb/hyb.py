#!/usr/bin/env python

import sys
import numpy
import math
import cmath

nffrq = 8193   # number of matsubara frequencies
nfreq = 2000   # number of mesh points
norbs = 2      # number of orbitals
width = 4.0    # bandwidth
beta  = 100.   # inversion temperature
step  = 0.0    # energy step
V     = 0.3    # hybridization strength 

# real axis
wmesh = numpy.zeros(2*nfreq + 1, numpy.float)

# semicircle density of states
dos_bethe = numpy.zeros(2*nfreq + 1, numpy.float)

for i in range(2*nfreq + 1):
    wmesh[i] = (i - nfreq)* width / (2.0 * nfreq)
    dos_bethe[i] = 2.0 / 4.0 / math.pi * math.sqrt(4.0 - wmesh[i] * wmesh[i])

dos_bethe[0] = 0.0
dos_bethe[2*nfreq] = 0.0

# matsubara frequency mesh
fmesh = numpy.zeros(nffrq, numpy.float)

# hybridization function for semicircular density of states
hyb_bethe = numpy.zeros(nffrq, numpy.complex)

step = wmesh[1] - wmesh[0]
scale = step * V * V
for i in range(nffrq):
    print 'index:', i
    fmesh[i] = ( 2.0 * i + 1.0 ) * math.pi / beta
    for j in range(2*nfreq + 1):
        hyb_bethe[i] = hyb_bethe[i] + dos_bethe[j] * scale / ( 1j * fmesh[i] - wmesh[j] )

# dump the data
f_hyb = open('solver.hyb.in', 'w')
for orb in range(norbs/2):
    for i in range(nffrq):
        f_hyb.write(" %6u %16.12s %16.12E %16.12E %16.12E %16.12E\n" % 
                     (orb + 1, fmesh[i], hyb_bethe[i].real, hyb_bethe[i].imag, hyb_bethe[i].real, hyb_bethe[i].imag))
    f_hyb.write("\n")
    f_hyb.write("\n")
f_hyb.close()
