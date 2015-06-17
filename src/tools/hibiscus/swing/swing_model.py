#!/usr/bin/env python

##
##
## Introduction
## ============
##
## It is a python script. The purpose of this script is define the basis
## functions for different models (Fermi liquid case and Mott insulator
## case). Now it implements the following python functions/classes:
##
##     def fparab
##
##     class lowEnergyFermiLiquid
##     ---> def Fun
##     ---> def Func
##     ---> def Matsubara
##     ---> def RealParts
##
##     class lowEnergyInsulator
##     ---> def Fun
##     ---> def Funr
##     ---> def Func
##     ---> def Matsubara
##     ---> def RealParts
##
##     class highEnergy
##     ---> def F0
##     ---> def Func
##     ---> def Matsubara
##     ---> def RealParts
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
## 12/20/2014 by li huang
##
##

from scipy import *
from scipy import optimize
from scipy import integrate
from scipy import interpolate

from swing_fast import *
from swing_dump import *
from swing_mesh import *

def fparab(par, x, data):
    """ the fitting function for self-energy, please refer to eq.(113)
    """
    chi2 = 0
    for i in range(len(x)):
        chi2 += ( (par[0] + x[i]*par[1] + x[i]**2*par[2]) - data[i] )**2
    return chi2

class lowEnergyFermiLiquid(object):
    """ low energy model for Fermi liquid system
    """
    def __init__(self, ifit, iom0, isig, p0, valueAtUnity=0.05):
        """ find out the fitting parameters
        """
        # fetch matsubara axis
        xfit = iom0[0:ifit]

        # imag part for self-energy function
        yfit_i = array(isig)[0:ifit,1]

        # real part for self-energy function
        yfit_r = array(isig)[0:ifit,0]

        # fit the real part and imag part respectively
        self.expan_i = optimize.fmin_powell(fparab, [0,0,0], args=(xfit,yfit_i), disp=0)
        self.expan_r = optimize.fmin_powell(fparab, [0,0,0], args=(xfit,yfit_r), disp=0)
        self.expan = self.expan_r.tolist() + self.expan_i.tolist()

        # fitting parameters
        self.a0 = abs(self.expan_i[0])
        self.a1 = self.expan_r[1]
        self.a2 = abs(self.expan_i[2])

        # parabola would otherwise becomes negative at some point
        if abs(self.a1) > sqrt(2*self.a0*self.a2):
            self.a1 = sqrt(2*self.a0*self.a2)*sign(self.a1)
        p0n = ((self.a0 + self.a1 + 0.5*self.a2) / valueAtUnity - 1) * (2/self.a2)**4

        # function should be smaller than valueAtUnity at x=1.0
        if p0 > p0n: self.p0 = p0
        else: self.p0 = p0n

        # normalization for the function
        self.pf = sqrt(2.)*self.a2**2*self.p0**(3/4.) / (2*pi*(2+self.a0*self.a2*sqrt(self.p0)))

    def Fun(self, x):
        """ the polynomial function, please refer to eq.(115)
        """
        return self.pf*(self.a0 + self.a1*x + 0.5*self.a2*x**2) / (1 + self.p0*(self.a2*x/2)**4)

    def Func(self, om):
        """ the polynomial function, please refer to eq.(115)
        """
        return array([self.Fun(x) for x in om])
        
    def Matsubara(self, iom):
        """ evaluate self-energy at matsubara axis
        """
        gm = 1/self.a2
        (x0, dh0) = swing_make_mesh(500, 1e-5*gm, 300*gm, gm)
        F0 = array([self.Fun(x) for x in x0])
        F00 = self.Fun(0.0)
        weigh0 = abs(self.expan_i[0]) / (pi*F00)

        datai=zeros(len(x0),dtype=float)
        datar=zeros(len(x0),dtype=float)
        F0i=[]
        F0r=[]
        for n in range(len(iom)):
            omn = iom[n]
            if (omn<0.3): subtract=1
            else: subtract=0
            for i in range(len(F0)):
                datai[i] = (F0[i]-F00*subtract)/(omn**2+x0[i]**2)
                datar[i] = F0[i]*x0[i]/(omn**2+x0[i]**2)
            wi = -(omn*integrate.simps(datai, x0) + F00*pi*subtract)
            wr = -integrate.simps(datar, x0)
            F0i.append(wi)
            F0r.append(wr)
            
        return (array(F0r), array(F0i), weigh0)
    
    def RealParts(self, om):
        """ evaluate self-energy at real axis
        """
        gm = 1/self.expan_i[2]
        (x0, dh0) = swing_make_mesh(500, 1e-5*gm, 300*gm, gm)
        F0 = array([self.Fun(x) for x in x0])
        F0j = array([self.Fun(x) for x in om])

        Frc = kramskron(om, F0j, F0, x0, dh0, 0.0)
        
        spl = interpolate.splrep(om, real(Frc), k=3, s=0.0)
        dersr = [interpolate.splev(0.0, spl, der=0),
                interpolate.splev(0.0, spl, der=1),
                interpolate.splev(0.0, spl, der=2)]
        spl = interpolate.splrep(om, imag(Frc), k=3, s=0.0)
        dersi = [interpolate.splev(0.0, spl, der=0),
                interpolate.splev(0.0, spl, der=1),
                interpolate.splev(0.0, spl, der=2)]
        
        ders = array(dersr + dersi)
        return (Frc, ders)

class lowEnergyInsulator(object):
    """ low energy model for mott insulator system
    """
    def __init__(self, ifit, iom0, isig, width):
        """ find out the fitting parameters
        """
        # fetch matsubara axis
        xfit = iom0[0:ifit]

        # imag part for self-energy function
        yfit_i = array(isig)[0:ifit,1]

        # real part for self-energy function
        yfit_r = array(isig)[0:ifit,0]

        # fit the real part and imag part respectively
        self.expan_i = optimize.fmin_powell(fparab, [0,0,0], args=(xfit,yfit_i), disp=0)
        self.expan_r = optimize.fmin_powell(fparab, [0,0,0], args=(xfit,yfit_r), disp=0)
        self.expan = self.expan_r.tolist() + self.expan_i.tolist()

        # fitting parameters
        self.a0 = abs(self.expan_i[0])
        self.width = width

    def Fun(self, x):
        """ the polynomial function, please refer to eq.(115)
        """
        return (self.width/pi)/(x**2 + self.width**2)

    def Funr(self, x):
        return 1/(x + self.width*1j)

    def Func(self, om):
        """ the polynomial function, please refer to eq.(115)
        """
        return array([self.Fun(x) for x in om])

    def Matsubara(self, iom):
        """ evaluate self-energy at matsubara axis
        """
        F00 = self.Fun(0.0)
        weigh0 = abs(self.expan_i[0])/(pi*F00)

        F0i=[]
        for n in range(len(iom)): F0i.append(-1/(iom[n]+self.width))

        F0r = zeros(len(iom), dtype=float)
        return (F0r, array(F0i), weigh0)
    
    def RealParts(self, om):
        """ evaluate self-energy at real axis
        """
        Frc = array([self.Funr(x) for x in om])
        dersr = [0.0, 1/self.width**2, 0.0]
        dersi = [-1/self.width, 0.0, 2/self.width**3]
        ders = array(dersr + dersi)
        return (Frc, ders)

class highEnergy(object):
    """ high energy model for self-energy function
    """
    def __init__(self, b_, pos_):
        self.b = b_
        self.pos = pos_

    def F0(self, om, En):
        """ define modified gaussians
        """
        if (En*om > 0):
            return 1./(self.b*abs(En)*sqrt(pi))*exp(-self.b**2/4.-(log(om/En)/self.b)**2)    
        else:
            return 0.0

    def Func(self, om):
        """ calculate all the gaussians
        """
        F=[]
        for i in range(len(self.pos)):
            En = self.pos[i]
            F.append(array([self.F0(x, En) for x in om]))
        return F
        
    def Matsubara(self, iom):
        """ calculate functions on imaginary axis
        """
        Funcr=[]
        Funci=[]
        for i in range(len(self.pos)):
            En = self.pos[i]

            (x0, dh0) = swing_make_mesh(300, 1e-2*abs(En), 2*abs(En), 7*abs(En))

            wb = []
            for x in x0:
                wb.append(self.F0(x+En, En))
            wb = array(wb)
            
            (gr, gi) = matsum(En, iom, x0, dh0, wb)
            
            Funcr.append(array(gr))
            Funci.append(array(gi))
            
        return (Funcr, Funci)

    def RealParts(self, om):
        """ calculate functions on real axis
        """
        Func=[]
        derivs=[]
        for i in range(len(self.pos)):
            En = self.pos[i]

            (x0, dh0) = swing_make_mesh(300, 1e-2*abs(En), 2*abs(En), 7*abs(En))

            wb = array([self.F0(x+En, En) for x in x0])
            
            F0j = array([self.F0(x, En) for x in om])

            Frc = kramskron(om, F0j, wb, x0, dh0, En)

            spl = interpolate.splrep(om, real(Frc), k=3, s=0.0)            
            dersr = [interpolate.splev(0.0, spl, der=0),
                     interpolate.splev(0.0, spl, der=1),
                     interpolate.splev(0.0, spl, der=2)]

            spl = interpolate.splrep(om, imag(Frc), k=3, s=0.0)
            dersi = [interpolate.splev(0.0, spl, der=0),
                     interpolate.splev(0.0, spl, der=1),
                     interpolate.splev(0.0, spl, der=2)]

            Func.append(Frc)
            ders = array(dersr + dersi)
            derivs.append(array(ders))
        return (Func, derivs)
