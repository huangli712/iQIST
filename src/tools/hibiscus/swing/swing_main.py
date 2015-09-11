#!/usr/bin/env python

##
##
## HIBISCUS/swing @ iQIST
##
## version: v2015.01.06T
## status : WARNING: IN TESTING STAGE, USE IT IN YOUR RISK
##
## Introduction
## ============
##
## This code is used to continue self-energy analytically from imaginary
## axis to real axis. This algorithm is invented by K. Haule, please refer
## to original paper: Phys. Rev. B 81, 195107 (2010).
##
## This code was written by K. Haule originally. And then we adapted this
## code such that it can be used with the iQIST software package.
##
## Usage
## =====
##
## python swing_main.py [options]
##
## The options can be as follows:
##    -sig filename
##      mandatory, input self-energy on imaginary axis (default:Sig.out)
##
##    -nom number
##      mandatory, number of matsubara points used (default:100)
##
##    -beta float
##      mandatory, inverse temperature (default:100.0)
##
##    -FL bool
##      mandatory, low energy expansion of a Fermi liquid of Mott
##      insulator (default:True)
##
##    -poles [[y0,A0,B0],[y1,A1,B1],...]
##      optional, poles of self-energy determined from spectral function,
##      poles will be forced at y0,y1,... this is achieved by contribution
##      to chi2+= Al/((x_i-yl)**2+w_i*Bl) where x_i are positions and w_i
##      weights of lorentzians to be determined
##
##    -b float:[0-1]
##      optional, basis functions, b parameter to determin family of basis
##      functions (default:0.8)
##
##    -Ng number
##      optional, number of basis functions (default:12)
##
##    -wexp float:[1-2]
##      optional, wexp^i is position where basis functions are
##      peaked (default:1.5)
##
##    -ifit number[4-..]
##      optional, low energy fit, number of matsubara points used to fit
##      the low energy expansion
##
##    -alpha3 float[0-1]
##      optional, weight for the normalization in functional to be
##      minimized (default:0.01)
##
##    -alpha4 float[0-1]
##      optional, weight for the low energy expansion in the functional
##      to be minimized (default:0.001)
##
##    -maxsteps number
##      optional, maximum number of function evaluations in minimization
##      routine (default:400000)
##
## The possible output files are as follows:
##    sigr.out
##      final self-energy function on real axis
##
##    sigr_linear.out
##      final self-energy function on fine real axis, it is used to
##      interface with LDA + DMFT code
##
##    siom.nnn
##      current function on imaginary axis (nnn counts every 40000 steps)
##
##    sres.nnn
##      current analytic continuation to real axis (nnn counts every
##      40000 steps)
##
##    gaus.nnn
##      current configuration for gaussians (nnn counts every 40000 steps)
##
## Comment: 
##    In order to run this code properly, you need to ensure scipy and
##    numpy were correctly installed on your system. This code was tested
##    on scipy 0.14.0 and numpy 1.7.0 only. So for the older versions of
##    scipy and numpy, we can not guarantee that it can work always.
##
##    To run this code, please use 'make' command to build the dynamical
##    library swing_fast.so at first.
##
## Author
## ======
##
## The original code was developed by K. Haule
## see http://hauleweb.rutgers.edu/downloads/
##
## This python script is designed, created, implemented, and maintained by
##
## Li Huang // email: lihuang.dmft@gmail.com
##
## History
## =======
##
## 01/06/2015 by li huang (created)
## 08/17/2015 by li huang (last modified)
##
##

import sys
import time

from scipy import *
from scipy import optimize

from  swing_dump import *
from  swing_chi2 import *
from  swing_mesh import *
from swing_model import *

if __name__ == '__main__':

    # print the running header at first
    swing_print_header()

    # record start time
    start_time = time.clock()

    # parameters in this module
    params = {
        'FL'      : [True,      '# low energy expansion of a Fermi liquid of Mott insulator'],
        'sig'     : ['Sig.out', '# input self-energy on imaginary axis'],
        'nom'     : [100,       '# number of matsubara points used'],
        'lcut'    : [0.0,       '# the lowest frequency lorentzian position'],
        'beta'    : [100.,      '# inverse temperature'],
        'b'       : [0.8,       '# b parameter to determin family of basis functions'],
        'Ng'      : [12,        '# number of basis functions'],
        'wexp'    : [1.5,       '# wexp^i is position where basis functions are peaked'],
        'ifit'    : [4,         '# number of matsubara points used to fit the low energy expansion'],
        'alpha3'  : [0.01,      '# weight for the normalization in functional to be minimized'],
        'alpha4'  : [0.001,     '# weight for the low energy expansion in the functional to be minimized'],
        'p0'      : [1,         '# how fast should the low energy function fall off'],
        'poles'   : [[],        '# additional poles'],
        'vunity'  : [0.05,      '# value of low-energy function at unity'],
        'L0'      : [15,        '# high energy cut-off for plotting'],
        'N0'      : [400,       '# number of frequency points for plotting'],
        'a0'      : [5e-3,      '# frequency mesh start'],
        'b0'      : [0.5,       '# frequency mesh parameter'],
        'rps'     : [1,         '# the lowest lorentzian is usualy at pi*T. this factor can bring it closer to zero.'],
        'maxsteps': [400000,    '# maximum number of function evaluations in minimization routine'],
    }

    # from command line
    arguments = sys.argv[1:]

    # if no arguments is given, print the document string
    if (len(arguments)==0 or arguments[0]=="-h" or arguments[0]=="--help"):
        print __doc__
        sys.exit(0)

    # dump the historic commands to hist.dat
    swing_dump_hist(sys.argv)

    # analyze the arguments, and update the default params
    for i in range(len(arguments)/2):
        if (arguments[2*i][1:] in params.keys()):
            # get variable name
            var = arguments[2*i][1:]
            # new value is assigned to string parameters
            if (type(params[var][0]) is str):
                params[var][0] = arguments[2*i+1]
            # new value is assigned to non-string parameters
            else:
                params[var][0] = eval(arguments[2*i+1])

    # we go over all parameters and evaluate them
    # let the parameters to be a variables in current namespace
    for var in params.keys():
        val = params[var][0]
        if (type(val) is str):
            val = '"'+str(val)+'"'
        else:
            val = str(val)
        # variable is set to its current value
        exec(var + '=' + val)

    # dump all the parameters to the screen
    swing_dump_keys(params)

    # reading qmc self-energy data
    f = open(sig)
    iom0 = []
    sqmc0 = []
    for line in f:
        spl = map(float,line.split())
        if not spl: break
        # matsubara frequency
        iom0.append(spl[0])
        # spl[1]: real part
        # spl[2]: imag part
        # spl[3]: standard deviation for real part
        # spl[4]: standard deviation for imag part
        sqmc0.append([spl[1], spl[2], spl[3], spl[4]])
    f.close()

    # subtracting Sigma_{infty}
    sinfty = sqmc0[-1][0]
    for i in range(len(sqmc0)):
        sqmc0[i][0] -= sinfty
    sqmc0 = array(sqmc0)

    # calculate integral of self-energy function
    Nstart = int(size(iom0) * 0.1)
    if Nstart < nom:
        Nstart = nom+1
    Nend = int(size(iom0) * 0.7)
    htail = [-sqmc0[i][1] * iom0[i] for i in range(Nstart,Nend)]
    intg = sum(htail) / len(htail)

    # setup temperature
    T = 1.0 / beta

    # position of gaussians
    gpos=[]
    for i in range(Ng):
        x0 = -rps*pi*T*wexp**(Ng-1-i)
        if abs(x0) > lcut:
            gpos.append(x0)
    for i in range(Ng):
        x0 = rps*pi*T*wexp**i
        if abs(x0) > lcut:
            gpos.append(x0)

    # creat real-frequency mesh, it is not a linear mesh, but a symmetric
    # logarithm and tan mesh
    (om, dh) = swing_make_mesh(N0, a0, b0, L0)

    # create matsubara mesh
    iom = [(2*n+1)*pi/beta for n in range(nom)]
    if (abs(iom[0]-iom0[0]) < T):
        sqmc = array(sqmc0[0:len(iom),  :], order='F')
    else:
        sqmc = array(sqmc0[1:len(iom)+1,:], order='F')

    # fitting low energy self-energy function on imaginary axis
    # the broad function at low frequency is prepared on real axis
    if (FL == True):
        lwe = lowEnergyFermiLiquid(ifit, iom0, sqmc0, p0, vunity)
    else:
        lwe = lowEnergyInsulator(ifit, iom0, sqmc0, 2*T)

    expandr = [lwe.expan_r[0],  lwe.expan_i[1], -lwe.expan_r[2]]
    expandi = [lwe.expan_i[0], -lwe.expan_r[1], -lwe.expan_i[2]]
    expand = array(expandr + expandi, order='F')

    # the low energy basis function
    A0 = lwe.Func(om)

    # the low energy function on matsubara axis
    (F0r, F0i, weight0) = lwe.Matsubara(iom)

    # real parts
    (Gw0, ders0) = lwe.RealParts(om)

    # build all gaussians on high energy
    hge = highEnergy(b, gpos)

    # basis functions
    Arest = hge.Func(om)

    # rest of the functions on matsubara axis
    (Frest_r, Frest_i) = hge.Matsubara(iom)

    # real parts
    (Gw_rest, ders_rest) = hge.RealParts(om)

    # assigning weights to each basis functions
    gweigh = []
    for i in range(len(gpos)):
        gweigh.append( abs(gpos[i])**1 )
    gweigh = array(gweigh)
    sm = sum(gweigh)
    # normalizing to start with the correct weight
    gweigh *= ( intg - weight0 ) / sm

    # putting them rogether
    ifunr = array(Frest_r   + [F0r], order='F')
    ifuni = array(Frest_i   + [F0i], order='F')
    rfun  = array(Arest     + [A0])
    rfunc = array(Gw_rest   + [Gw0])
    ders  = array(ders_rest + [ders0], order='F')

    # preparing index arrays for minimization
    gwfix= array([weight0], order='F')
    abounds = [(1e-12,None)] * len(gpos)
    alphas = [alpha3,alpha4]

    vary=array(range(1,len(gpos) + 1), order='F')
    fixed=array([len(gpos) + 1], order='F')

    # minimization routine is called, the well known L-BFGS-B algorithm
    # is used, which is included in scipy.optimize package
    #
    # please pay attention to the chi2 in the output log. it should be
    # less than 1. If it is large than 1 when the code is terminated,
    # please increase maxsteps and decrease factr, and run it again.
    #
    # approx_grad: bool, whether to approximate the gradient numerically
    #
    # bounds: list, (min, max) pairs for each element in x, defining the
    # bounds on that parameter. Use None for one of min or max when there
    # is no bound in that direction.
    #
    # maxfun: maximum number of function evaluations.
    (gweigh, fmin, dinf) = optimize.fmin_l_bfgs_b(swing_cchi, gweigh, 
       approx_grad=True, factr=1000, bounds=abounds, maxfun=maxsteps,
       args=(vary, gwfix, fixed, sqmc, ifunr, ifuni, iom, intg, om, 
             rfun, rfunc, expand, ders, alphas, gpos, poles))

    # dump the final self-energy function on real axis, please refer to
    # sigr.out and sigr_linear.out file
    swing_dump_sigr(om, vary, fixed, gweigh, gwfix, rfunc, sinfty)

    # record the end time
    end_time = time.clock()

    # print the final footer for this code
    swing_print_footer(end_time - start_time)
