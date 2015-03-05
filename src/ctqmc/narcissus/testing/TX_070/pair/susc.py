#!/usr/bin/env python

import sys
import math
import numpy
import numpy.linalg
sys.path.append("/home/huangl/iqist/iqist2015/src/tools/hibiscus/script/")
from u_reader import *

def cal_gk(nband, norbs, nkpts, mfreq, mune, fmesh, hk, slat):
    """ calculate G_{lat}(i\omega) by the dyson equation
        mune  : chemical potential
        fmesh : matsubara frequency
        hk    : lda hamiltonian
    """
    # reset glat
    glat = numpy.zeros((mfreq, norbs, nkpts), dtype = numpy.complex)

    # unit matrix
    I = numpy.eye(nband)

    # dummy matrix, used to store self-energy
    saux = numpy.zeros(nband, dtype = numpy.complex)

    # loop over frequencies and k-points, it is parallelized
    for f in range(mfreq):
        fM = (1j * fmesh[f] + mune) * I
        for k in range(nkpts):
            # calculate \Sigma - \Sigma_{dc}, only consider the spin up part
            for s in range(nband):
                saux[s] = slat[f,s,s]
            # calculate G^-1 = i\omega + \mu - H - (\Sigma - \Sigma_{dc})
            A =  fM - hk[:,:,k] - numpy.diag(saux)
            # note: you can calculate inverse of A directly, B = numpy.linalg.inv(A)
            # but what I do is more efficient.
            B = numpy.linalg.solve(A, I)
            # setup the final results
            for s in range(nband):
                glat[f,s,k]   = B[s,s] # spin up part
                glat[f,s+nband,k] = B[s,s] # spin down part

    # release memory
    del saux

    return glat

def get_hk(nband, norbs, nkpts, Z):
    """ build the lda hamiltonian and corresponding impurity level
        nband : number of bands
        norbs : number of orbitals
        nkpts : number of k-points
        Z     : renormalization factor for hamiltonian
    """
    # lda hamiltonian
    hk = numpy.zeros((nband, nband, nkpts), dtype = numpy.complex)

    # open data file : in.hk
    f = open("../in.hk", "r")

    # check the dimension size of hamiltonian
    spl = f.readline().split()
    assert ( nkpts == float(spl[0]) )
    assert ( nband == float(spl[1]) )

    # read the data
    # loop over k-points
    for ikpt in range(nkpts):
        # check the index of k-point
        spl = f.readline().split()
        assert ( ikpt + 1 == int(spl[0]) )
        # loop over band index
        for m in range(nband):
            for n in range(nband):
                line = f.readline()
                spl = map(float, line.split())
                hk[m,n,ikpt] = spl[0] + 1j * spl[1]

    # close data file
    f.close()

    # we just assume it is a hermitian
    # overwrite the input matrix with its hermitean part
    for ikpt in range(nkpts):
        for m in range(nband):
            for n in range(m,nband):
                hk[m,n,ikpt] = ( hk[m,n,ikpt] + hk[n,m,ikpt].conjugate() ) * 0.5
                hk[n,m,ikpt] = hk[m,n,ikpt].conjugate()

    # renormalize the hamiltonian
    # note: the Z is equal to 1.0 at default
    hk = Z * hk

    return hk

# setup parameters
norbs = 2
nband = 1
nbfrq = 1
nffrq = 128
nfreq = 1024
nkpts = 216
mfreq = 1024
beta = 20.0 # have to be fixed
mune = 0.303235740115 # have to be fixed
bose = 0.629867541391

# read my solver.pair.dat
p2 = iqistReader.get_pair(norbs, nffrq, nbfrq)

# read G(i\omega)
rmesh, grnf = iqistReader.get_grn(norbs, mfreq)

# read self-energy function
rmesh, sigf = iqistReader.get_sgm(norbs, mfreq)

# read Hk
hk = get_hk(nband, norbs, nkpts, bose)

# build G_k
gk = cal_gk(nband, norbs, nkpts, mfreq, mune, rmesh, hk, sigf)

# \chi_{loc}
chi_loc = p2[:,:,0,1,0] / beta

# \chi^{0}_{loc}
chi_loc_z = numpy.zeros(nfreq, dtype = numpy.complex)
for i in range(nfreq):
    if i < nfreq/2:
        chi_loc_z[i] = grnf[nfreq/2-1-i,0,0].conjugate() * grnf[nfreq/2-1-i,0,0]
    else:
        chi_loc_z[i] = grnf[i-nfreq/2,0,0] * grnf[i-nfreq/2,0,0].conjugate()

# \chi^{0}_{q}
chi_q_z = numpy.zeros(nfreq, dtype = numpy.complex)
for i in range(nfreq):
    kaux = 0.0
    if i < nfreq/2:
        for k in range(nkpts):
            kaux = kaux + gk[nfreq/2-1-i,0,k].conjugate() * gk[nfreq/2-1-i,0,k]
    else:
        for k in range(nkpts):
            kaux = kaux + gk[i-nfreq/2,0,k] * gk[i-nfreq/2,0,k].conjugate()
    chi_q_z[i] = kaux / float(nkpts)

# calculate final results: \chi_q
gamma = 1.0 / chi_q_z - 1.0 / chi_loc_z
start = nfreq / 2 - nffrq / 2
end = nfreq / 2 + nffrq / 2 - 1
A = numpy.linalg.inv(chi_loc)
for i in range(nffrq):
    A[i,i] = A[i,i] + gamma[start + i]
A = numpy.linalg.inv(A)
gamma = chi_q_z
tail = numpy.sum(gamma[0:start]) + numpy.sum(gamma[end+1:nfreq])
print 'chi^-1: ', (beta / (numpy.sum(A) + tail)).real
