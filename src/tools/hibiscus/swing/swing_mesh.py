#!/usr/bin/env python

##
##
## Introduction
## ============
##
## It is a python script. The purpose of this script is build logarithm
## and tan mesh for the swing code. Now it implements the following
## python functions/classes:
##
##     def swing_ltan_mesh
##     def swing_make_mesh
##
## Usage
## =====
##
## Sorry, it can not be invoked manually.
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
## 12/20/2014 by li huang (created)
## 08/17/2015 by li huang (last modified)
##
##

from scipy import *

def swing_ltan_mesh(N, x0, x1, x2):
    """ build logarithm and tan mesh
    """
    eta = log(x1 / x0) / (x2 / x1 - 1)
    N1 = int((1 + eta * N) / (1 + eta) + 0.5)
    if (N1 > N - 2): N1 = N - 2
    N2 = N - N1
    xt = x2 / x1
    dwt = N2 * (log(x1) - log(x0)) / (N1 - 1)
    ut = 1e-5

    a = arctan(tan(ut) / xt)
    b = dwt * sin(a) * cos(a)
    w = x1 / tan(a)

    om=[]
    # build logarithm mesh
    for i in range(N1):
        om.append( exp( log(x0) + i * ( log(x1) - log(x0) ) / ( N1 - 1 ) ) )

    # build tan mesh
    for i in range(N2):
        om.append( w * tan( a + (i + 1) * b / N2 ) )

    return om

def swing_make_mesh(N, x0, x1, x2):
    """ for building mesh which is logarithmic at small frequency and tan
        at large frequency. mesh is symmetric around zero frequency
    """
    N2 = N / 2
    om1 = swing_ltan_mesh(N2, x0, x1, x2)

    # negative half axis
    om=[]
    for i in range(N2):
        om.append(-om1[N2-i-1])

    # positive half axis
    for i in range(N2):
        om.append(om1[i])

    # build \delta h
    dh=[]
    dh.append( 0.5 * ( om[1] - om[0] ) )
    for i in range(1,len(om)-1):
        dh.append( 0.5 * ( om[i+1] - om[i-1] ) )
    dh.append( 0.5 * ( om[-1] - om[-2] ) )

    return (array(om), array(dh))
