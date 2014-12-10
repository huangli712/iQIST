#!/usr/bin/env python

#
#
# call the alps impurity solver to solve an impurity model
#
#

# python standard library
import os
import sys
import time
import math
import shutil
import ctypes

# python 3rd-party library
import numpy
import scipy.weave
import scipy.interpolate

# python binding for ALPS
import pyalps.mpi
import pyalps.cthyb

class switch(object):
    """ This class provides the functionality we want. You only need to look
        at this if you want to know how this works. It only needs to be defined
        once, no need to muck around with its internals.
    """
    def __init__(self, value):
        """ class constructor for the switch
        """
        self.__value = value
        self.__fall = False

    def __iter__(self):
        """ return the match method once, then stop
        """
        yield self.match
        raise StopIteration

    def match(self, *args):
        """ indicate whether or not to enter a case suite
        """
        if self.__fall or not args:
            return True
        elif self.__value in args:
            self.__fall = True
            return True
        else:
            return False

class Symm(object):
    """ class Symm is used to symmetrize the input data over imaginary time
        or spin. it can be also used to get rid of the real part or imaginary
        part of the input data.
    """
    @staticmethod
    def symm_over_time(time_data):
        """ perform symmetrization over imaginary time
            time_data : imaginary time data
        """
        # get number of elements
        n = len(time_data)
        for i in range(n):
            raux = 0.5 * ( time_data[i] + time_data[n-1-i] )
            time_data[i] = raux
            time_data[n-1-i] = raux

    @staticmethod
    def symm_over_spin(time_data_up, time_data_dn):
        """ perform symmetrization over spin
            time_data_up : imaginary time data for spin up
            time_data_dn : imaginary time data for spin dn
        """
        raux = 0.5 * ( time_data_up + time_data_dn )
        time_data_up = raux.copy()
        time_data_dn = raux.copy()
        del raux

    @staticmethod
    def symm_over_band(time_data):
        """ perform symmetrization over band/orbital
            time_data : imaginary time data. it should be a numpy array
        """
        norbs = len(time_data)
        mean = numpy.mean(time_data, axis = 0)
        for i in range(norbs):
            time_data[i] = mean.copy()
        del mean

    @staticmethod
    def symm_over_real(freq_data):
        """ get rid of the imaginary part of the input data
            freq_data : matsubara frequency data
        """
        # get number of elements
        n = len(freq_data)
        for i in range(n):
            freq_data[i] = freq_data[i].real + 0j

    @staticmethod
    def symm_over_imag(freq_data):
        """ get rid of the real part of the input data
            freq_data : matsubara frequency data
        """
        # get number of elements
        n = len(freq_data)
        for i in range(n):
            freq_data[i] = 0.0 + 1j * freq_data[i].imag

    @staticmethod
    def symm_smth_data(data, naver):
        """ used to smooth 1-d float data
            data : 1-d array
            naver : used to define the smooth region
        """
        # determine the length of 1-d data
        n = len(data)

        # create a temporary array
        data_s = numpy.zeros_like(data, dtype = numpy.float)

        for i in range(n):
            # bracketing the smooth region [i1, i2]
            i1 = i - naver
            if i1 < 0: i1 = 0

            i2 = i + naver
            if i2 > n - 1: i2 = n - 1

            # accumulate the data
            count = 0
            for j in range(i1, i2):
                count = count + 1
                data_s[i] = data_s[i] + data[j]

            # calculate the average
            data_s[i] = data_s[i] / float(count)

        # copy smoothed data to original data
        for i in range(n):
            data[i] = data_s[i]

        # release memory
        del data_s

class Mesh(object):
    """ class Mesh is used to define fermionic and bosonic meshes, imaginary
        time mesh is also defined in this class.
        note : the k-point mesh is defined in the Lattice class.
    """
    def __init__(self, p):
        """ class constructor for the Mesh
            p : param dictionary
            beta  : inversion temperature
            ntime : number of imaginary time points
            nffrq : number of frequencies for fermionic mesh
            nbfrq : number of frequencies for bosonic mesh
        """
        self.__beta  = p["beta"][0]
        self.__ntime = p["ntime"][0]
        self.__nffrq = p["nffrq"][0]
        self.__nbfrq = p["nbfrq"][0]

    def create_mesh(self, mesh_type):
        """ construct imaginary time, fermionic or bosonic mesh according
            to the mesh_type variable.
            mesh_type : which type of mesh should be constructed
        """
        if   mesh_type == 1: # construct imaginary time mesh
            tmesh = numpy.zeros(self.__ntime, dtype = numpy.float)
            for i in range(len(tmesh)):
                tmesh[i] = float(i) * self.__beta / ( self.__ntime - 1 )
            return tmesh

        elif mesh_type == 2: # construct fermionic mesh
            fmesh = numpy.zeros(self.__nffrq, dtype = numpy.float)
            for i in range(len(fmesh)):
                fmesh[i] = ( 2.0 * i + 1.0 ) * math.pi / self.__beta
            return fmesh
        else:                # construct bosonic mesh
            bmesh = numpy.zeros(self.__nbfrq, dtype = numpy.float)
            for i in range(len(bmesh)):
                bmesh[i] = ( 2.0 * i + 0.0 ) * math.pi / self.__beta
            return bmesh

    @staticmethod
    def create_time_mesh(beta, ntime):
        """ construct the imaginary time mesh
            beta  : inversion temperature
            ntime : number of imaginary time points
        """
        tmesh = numpy.zeros(ntime, dtype = numpy.float)
        for i in range(len(tmesh)):
            tmesh[i] = float(i) * beta / ( ntime - 1)
        return tmesh

    @staticmethod
    def create_fermion_mesh(beta, nffrq):
        """ construct the fermionic mesh
            beta  : inversion temperature
            nffrq : number of frequencies for fermionic mesh
        """
        fmesh = numpy.zeros(nffrq, dtype = numpy.float)
        for i in range(len(fmesh)):
            fmesh[i] = ( 2.0 * i + 1.0 ) * math.pi / beta
        return fmesh

    @staticmethod
    def create_boson_mesh(beta, nbfrq):
        """ create the bosonic mesh
            beta  : inversion temperature
            nbfrq : number of frequencies for bosonic mesh
        """
        bmesh = numpy.zeros(nbfrq, dtype = numpy.float)
        for i in range(len(bmesh)):
            bmesh[i] = ( 2.0 * i + 0.0 ) * math.pi / beta
        return bmesh

class Dump(object):
    """ class Dump is used to dump the key variables to the disk file.
    """
    @staticmethod
    def dump_time_data3(file_name, time_mesh, time_data):
        """ dump the imaginary time data to the disk file, multi-orbital version
            file_name : file name
            time_mesh : imaginary time mesh
            time_data : array in imaginary time mesh. it is in 2d form
        """
        norbs = len(time_data) # get number of orbitals
        f = open(file_name, "w")
        for orb in range(norbs):
            for i in range(len(time_mesh)):
                f.write( "%6u %6u %16.8s %16.8E\n" %
                    (orb + 1, i + 1, time_mesh[i], time_data[orb][i]) )
            f.write("\n") # dump two blank lines to separate the data for different orbitals
            f.write("\n")
        f.close()

    @staticmethod
    def dump_freq_data3(file_name, freq_mesh, freq_data):
        """ dump the matsubara frequency data to the disk file, multi-orbital version
            file_name : file name
            freq_mesh : matsubara frequency mesh
            freq_data : array in matsubara frequency mesh. it is in 2d form
        """
        norbs = len(freq_data) # get number of orbitals
        f = open(file_name, "w")
        for orb in range(norbs):
            for i in range(len(freq_mesh)):
                f.write( "%6u %6u %16.8s %16.8E %16.8E\n" %
                    (orb + 1, i + 1, freq_mesh[i], freq_data[orb][i].real, freq_data[orb][i].imag) )
            f.write("\n") # dump two blank lines to separate the data for different orbitals
            f.write("\n")
        f.close()

    @staticmethod
    def dump_htau_data(file_name, time_mesh, htau_data):
        """ dump the hybridization function \Delta(\tau) to the disk file,
            multi-orbital version, this function is designed for the ctqmc
            impurity solver contained in ALPS package
            file_name : file name
            time_mesh : imaginary time mesh
            htau_data : array in imaginary time mesh
        """
        norbs = len(htau_data) # get number of orbitals
        f = open(file_name, "w")
        for i in range(len(time_mesh)):
            f.write( "%6u " % (i) )
            for orb in range(norbs):
                f.write( "%16.8E " % (htau_data[orb][i]) )
            f.write("\n") # start a new line
        f.close()

    @staticmethod
    def dump_ktau_data(file_name, time_mesh, ktau_data, ptau_data):
        """ dump the K(\tau) and K'(\tau) to the disk file, this function is
            designed for the ctqmc impurity solver contained in ALPS package
            file_name : file name
            time_mesh : imaginary time mesh
            ktau_data : array in imaginary time mesh
            ptau_data : array in imaginary time mesh
        """
        # note: for ALPS ctqmc impurity solver, the K(\tau) and K'(\tau)
        # are equal for all orbitals. So, here we only output the data
        # for band 1. If you want to see all of the orbital-dependent
        # data for K(\tau) and K'(\tau), please check c_ktau.dat.* and
        # c_ptau.dat.* files.
        nband = len(ktau_data) # get number of bands
        f = open(file_name, "w")
        for i in range(len(time_mesh)):
            f.write( "%6u " % (i) )
            for orb in range(1):
                f.write( "%16.8E %16.8E " % (ktau_data[orb][i], ptau_data[orb][i]) )
            f.write("\n") # start a new line
        f.close()

    @staticmethod
    def dump_kdep_data(file_name, nfreq, nkx, nky, kdep_data):
        """ dump k-dependent data in matsubara frequency axis to the disk file
            file_name : file name
            nfreq     : number of frequency points
            nkx       : number of k-points
            nky       : number of k-points
            kdep_data : array in matsubara frequency axis
        """
        is_complex = isinstance(kdep_data[0,0,0], complex)
        fk = open(file_name, "w")
        for f in range(nfreq):
            if is_complex: # for complex data
                # write real part
                for ikx in range(nkx):
                    for iky in range(nky):
                        fk.write( "%16.8E" % ( kdep_data[f,ikx,iky].real ) )
                    fk.write("\n")
                fk.write("\n") # write two blank lines
                fk.write("\n")

                # write imaginary part
                for ikx in range(nkx):
                    for iky in range(nky):
                        fk.write( "%16.8E" % ( kdep_data[f,ikx,iky].imag ) )
                    fk.write("\n")
                fk.write("\n") # write two blank lines
                fk.write("\n")
            else: # for real data
                for ikx in range(nkx):
                    for iky in range(nky):
                        fk.write( "%16.8E" % ( kdep_data[f,ikx,iky] ) )
                    fk.write("\n")
                fk.write("\n") # write two blank lines
                fk.write("\n")
        fk.close()

class Fourier(object):
    """ class Fourier is used to perform forward and backward fourier
        transformations between imaginary time and matsubara frequency
        spaces.
    """
    def __init__(self, p, ftype):
        """ class constructor for the Fourier
            p : param dictionary
            ftype : fourier type
            beta  : inversion temperature
            ntime : number of imaginary time points
            nfreq : number of frequencies for fermionic/bosonic mesh
        """
        self.__beta  = p["beta"][0]
        self.__ntime = p["ntime"][0]

        if ftype == 1: # for fermionic system
            # set nfreq to number of fermionic frequency
            self.__nfreq = p["nffrq"][0]
        else:          # for bosonic system
            # set nfreq to number of bosonic frequency
            self.__nfreq = p["nbfrq"][0]

    def __tails(self, freq_mesh, freq_data):
        """ calculate high frequency tails using K. Haule's trick
            freq_mesh : matsubara frequency grid
            freq_data : function on matsubara frequency space
        """
        Sn = 0.0
        Sx = 0.0
        Sy = 0.0

        Sxx = 0.0
        Sxy = 0.0

        for j in range(self.__nfreq - 128, self.__nfreq):
            Sn = Sn + 1.0
            Sx = Sx + 1.0 / freq_mesh[j]**2
            Sy = Sy + freq_data[j].imag * freq_mesh[j]
            Sxx = Sxx + 1.0 / freq_mesh[j]**4
            Sxy = Sxy + freq_data[j].imag * freq_mesh[j] / freq_mesh[j]**2

        return (Sx * Sxy - Sxx * Sy) / (Sn * Sxx - Sx * Sx)

    def freq_to_time_F(self, time_mesh, time_data, freq_mesh, freq_data, fast = True):
        """ backward fourier transformation from matsubara frequency to imaginary time,
            this function is only suitable for the fermionic system
            time_mesh : imaginary time mesh
            time_data : function on imaginary time axis
            freq_mesh : matsubara frequency mesh
            freq_data : function on matsubara frequency axis
            fast      : whether scipy.weave is used to accelerate the code
        """
        # calculate high frequency tails need to be subtracted
        tail = self.__tails(freq_mesh, freq_data)

        # perform backward fourier transformation
        if fast: # scipy weave version
            # build weave arguments
            ntime = self.__ntime
            nfreq = self.__nfreq
            beta  = self.__beta
            ftail = float(tail) # attention: tail is in numpy.float64 type, while ftail is in float type

            # define c++ code
            code = """
            #include <complex>

            double raux;
            for ( int i = 0; i < ntime; i++ ) {
                raux = 0.0;
                for ( int j = 0; j < nfreq; j++ ) {
                    raux = raux + cos( freq_mesh(j) * time_mesh(i) ) *   real(freq_data(j));
                    raux = raux + sin( freq_mesh(j) * time_mesh(i) ) * ( imag(freq_data(j)) + ftail / freq_mesh(j));
                }
                time_data(i) = 2.0 * raux / beta - 0.5 * ftail;
            }
            """

            # inline c++ code using weave
            scipy.weave.inline(code, ['ntime', 'nfreq', 'beta', 'ftail', 'freq_mesh', 'freq_data', 'time_mesh', 'time_data'],
                               type_converters=scipy.weave.converters.blitz)
        else: # pure python + numpy version
            for i in range(self.__ntime):
                raux = 0.0
                for j in range(self.__nfreq):
                    raux = raux + math.cos( freq_mesh[j] * time_mesh[i] ) *   freq_data[j].real
                    raux = raux + math.sin( freq_mesh[j] * time_mesh[i] ) * ( freq_data[j].imag + tail / freq_mesh[j])
                time_data[i] = 2.0 * raux / self.__beta - 0.5 * tail

        # corrections for the boundary point
        raux = freq_data[self.__nfreq-1].real * freq_mesh[self.__nfreq-1] / math.pi
        time_data[0] = time_data[0] + raux
        time_data[self.__ntime-1] = time_data[self.__ntime-1] - raux

    def time_to_freq_F(self, time_mesh, time_data, freq_mesh, freq_data, fast = True):
        """ forward fourier transformation from imaginary time to matsubara frequency,
            this function is only suitable for the fermionic system
            time_mesh : imaginary time mesh
            time_data : function on imaginary time axis
            freq_mesh : matsubara frequency mesh
            freq_data : function on matsubara frequency axis
            fast      : whether scipy.weave is used to accelerate the code
        """
        # perform fourier transformation
        if fast : # scipy weave version
            # build weave arguments
            ntime = self.__ntime
            nfreq = self.__nfreq

            # define c++ code
            code = """
            #include <complex>

            std::complex<double> ci(0.0, 1.0);
            for ( int i = 0; i < nfreq; i++ ) {
                double sre = 0.0;
                double sim = 0.0;
                for ( int j = 0; j < ntime - 1; j++ ) {
                    double c0 = cos( time_mesh(j)   * freq_mesh(i) );
                    double c1 = cos( time_mesh(j+1) * freq_mesh(i) );
                    double s0 = sin( time_mesh(j)   * freq_mesh(i) );
                    double s1 = sin( time_mesh(j+1) * freq_mesh(i) );
                    double g0 = time_data(j);
                    double g1 = time_data(j+1);
                    double dg = ( g1 - g0 ) / ( time_mesh(j+1) - time_mesh(j) );
                    sim = sim + ( c0 * g0 - c1 * g1 + dg * (s1 - s0) / freq_mesh(i) ) / freq_mesh(i);
                    sre = sre + ( s1 * g1 - s0 * g0 + dg * (c1 - c0) / freq_mesh(i) ) / freq_mesh(i);
                }
                freq_data(i) = sre + ci*sim;
            }
            """

            # inline c++ code using weave
            scipy.weave.inline(code, ['ntime', 'nfreq', 'freq_mesh', 'freq_data', 'time_mesh', 'time_data'],
                               type_converters=scipy.weave.converters.blitz)
        else : # pure python + numpy version
            for i in range(self.__nfreq):
                sre = 0.0
                sim = 0.0
                for j in range(self.__ntime-1):
                    c0 = math.cos( time_mesh[j]   * freq_mesh[i] )
                    c1 = math.cos( time_mesh[j+1] * freq_mesh[i] )
                    s0 = math.sin( time_mesh[j]   * freq_mesh[i] )
                    s1 = math.sin( time_mesh[j+1] * freq_mesh[i] )
                    g0 = time_data[j]
                    g1 = time_data[j+1]
                    dg = ( g1 - g0 ) / ( time_mesh[j+1] - time_mesh[j] )
                    sim = sim + ( c0 * g0 - c1 * g1 + dg * (s1 - s0) / freq_mesh[i] ) / freq_mesh[i]
                    sre = sre + ( s1 * g1 - s0 * g0 + dg * (c1 - c0) / freq_mesh[i] ) / freq_mesh[i]
                freq_data[i] = sre + 1j*sim

    def freq_to_time_B(self, time_mesh, time_data, freq_mesh, freq_data, fast = True):
        """ backward fourier transformation from matsubara frequency to imaginary time,
            this function is only suitable for the bosonic function
            time_mesh : imaginary time mesh
            time_data : function on imaginary time axis
            freq_mesh : matsubara frequency mesh
            freq_data : function on matsubara frequency axis
            fast      : whether scipy.weave is used to accelerate the code
        """
        # perform backward fourier transformation
        if fast: # scipy weave version
            # build weave arguments
            ntime = self.__ntime
            nfreq = self.__nfreq
            beta  = self.__beta

            # define c++ code
            code = """
            #include <complex>

            std::complex<double> ci(0.0, 1.0);
            std::complex<double> caux;
            for ( int t = 0; t < ntime; t++ ) {
                caux = freq_data(0);
                for ( int f = 1; f < nfreq; f++ ) {
                    caux = caux + 2.0 * freq_data(f) * exp(-ci * freq_mesh(f) * time_mesh(t) );
                }
                time_data(t) = real(caux) / beta;
            }
            """

            # inline c++ code using weave
            scipy.weave.inline(code, ['ntime', 'nfreq', 'beta', 'freq_mesh', 'freq_data', 'time_mesh', 'time_data'],
                               type_converters=scipy.weave.converters.blitz)
        else: # pure python + numpy version
            for t in range(self.__ntime):
                caux = freq_data[0]
                for f in range(1,self.__nfreq):
                    caux = caux + 2.0 * freq_data[f] * numpy.exp(-1j * freq_mesh[f] * time_mesh[t])
                time_data[t] = caux.real / self.__beta

    def time_to_freq_B(self, time_mesh, time_data, freq_mesh, freq_data, fast = True):
        """ forward fourier transformation from imaginary time to matsubara frequency,
            this function is only suitable for the bosonic function
            time_mesh : imaginary time mesh
            time_data : function on imaginary time axis
            freq_mesh : matsubara frequency mesh
            freq_data : function on matsubara frequency axis
            fast      : whether scipy.weave is used to accelerate the code
        """
        fit = scipy.interpolate.InterpolatedUnivariateSpline(time_mesh, time_data)
        ntime_dense = 4 * self.__ntime # used to build a dense imaginary time mesh
        # denser time mesh
        time_mesh_dense = Mesh.create_time_mesh(self.__beta, ntime_dense)
        # denser time data
        time_data_dense = fit(time_mesh_dense)
        for im in range(self.__nfreq):
            faux = time_data_dense * numpy.exp(1j * freq_mesh[im] * time_mesh_dense)
            # calculate \int^{\beta}_{0} f(\tau) e^{i \omega \tau} d\tau
            # now faux = f(\tau) e^{i \omega \tau}
            freq_data[im] = numpy.trapz(faux, time_mesh_dense)

    def time_to_freq_C(self, time_mesh, time_data, freq_mesh, freq_data):
        """ forward fourier transformation from imaginary time to matsubara frequency,
            this function is only suitable for the bosonic function \chi
            time_mesh : imaginary time mesh
            time_data : function on imaginary time axis
            freq_mesh : matsubara frequency mesh
            freq_data : function on matsubara frequency axis
        """
	# spline interpolation to evaluate the high-frequency tail
        fit = scipy.interpolate.InterpolatedUnivariateSpline(time_mesh, time_data)
        deriv_0 = fit.derivatives(0.0)
        deriv_beta = fit.derivatives(self.__beta)
        c1 = deriv_beta[0] - deriv_0[0]
        c2 = -(deriv_beta[1] - deriv_0[1])
        c3 = deriv_beta[2] - deriv_0[2]
        for im in range(1,self.__nfreq): # im = 0 will cause invalid arithmetic opeartion
            freq_data[im] = c1 / (1j*freq_mesh[im])    \
                          + c2 / (1j*freq_mesh[im])**2 \
                          + c3 / (1j*freq_mesh[im])**3

	# contribution from the rest part
        ntime_dense = 2 * self.__ntime # used to build a dense imaginary time mesh
        time_mesh_dense = Mesh.create_time_mesh(self.__beta, ntime_dense)
        time_data_dense = fit(time_mesh_dense)
        for i in range(ntime_dense):
            time_data_dense[i] = time_data_dense[i] \
                               - c1 * ( time_mesh_dense[i] - self.__beta/2 ) / self.__beta \
                               + c2 * ( time_mesh_dense[i] - self.__beta/2 )**2 / (2.0*self.__beta) - c2 * self.__beta / 24.0 \
                               - c3 * ( time_mesh_dense[i] - self.__beta/2 )**3 / (6.0*self.__beta)

        cutoff = min(self.__ntime, self.__nfreq)
        for im in range(cutoff):
            faux = time_data_dense * numpy.exp(1j * freq_mesh[im] * time_mesh_dense)
            # calculate \int^{\beta}_{0} f(\tau) e^{i \omega \tau} d\tau
            # now faux = f(\tau) e^{i \omega \tau}
            freq_data[im] += numpy.trapz(faux, time_mesh_dense)

class AlpsSolver(object):
    """ class AlpsSolver is used to drive the hybridization expansion impurity
        solver contained in the ALPS package.
    """
    def __init__(self, params):
        """ class constructor for the AlpsSolver
            params : user-supply input parameters
        """
        # setup default parameters for the quantum impurity solver
        # note : MAX_TIME IS TOO SMALL. IT SHOULD BE MODIFIED BEFORE CALCULATION.
        self.__params = {
            'N_MEAS'                    : 100        ,
            'SWEEPS'                    : 200000000  ,
            'THERMALIZATION'            : 1000       ,
            'MAX_TIME'                  : 240        ,
            'N_ORBITALS'                : 2          ,
            'N_TAU'                     : 1023       ,
            'N_nn'                      : 256        ,
            'N_MATSUBARA'               : 512        ,
            'N_LEGENDRE'                : 48         ,
            'N_HISTOGRAM_ORDERS'        : 128        ,
            'MU'                        : 4.0        ,
            'U'                         : 8.0        ,
            'J'                         : 0.0        ,
            'BETA'                      : 100.0      ,
            'MEASURE_time'              : 1          ,
            'MEASURE_freq'              : 1          ,
            'MEASURE_nn'                : 1          ,
            'MEASURE_nnt'               : 16         ,
            'MEASURE_legendre'          : 1          ,
            'TEXT_OUTPUT'               : 1          ,
            #'SPINFLIP'                  : 1          ,
            'DELTA'                     : "D.dat"    ,
            #'RET_INT_K'                 : "K.dat"    ,
            'SEED'                      : int( time.time() ) + ( pyalps.mpi.rank + 3 ) * 1981
        }

        # setup additional keys and values
        for key in params.iterkeys():
            self.__params[key] = params[key]

    def setup(self):
        """ prepare necessary inputs and perform check to ensure everything is OK
        """
        # only master node can do it
        if pyalps.mpi.rank != 0:
            return

        # check whether the delta.dat is available
        if not self.__params.has_key("DELTA"):
            sys.exit("FATAL ERROR : please provide the hybridization function")

        if not os.path.isfile(self.__params["DELTA"]):
            try:
                raise IOError(self.__params["DELTA"])
            except IOError as e:
                print ('IOError : ' + str(e) + " does not exist!")
                sys.exit("FATAL ERROR : please check it and then run this code again")

        return
        # check whether the K.dat is available
        if not self.__params.has_key("RET_INT_K"):
            sys.exit("FATAL ERROR : please provide the K(\tau) and K'(\tau) data")

        if not os.path.isfile(self.__params["RET_INT_K"]):
            try:
                raise IOError(self.__params["RET_INT_K"])
            except IOError as e:
                print ('IOError : ' + str(e) + " does not exist!")
                sys.exit("FATAL ERROR : please check it and then run this code again")

    def start(self, curr_iter, params = {}):
        """ invoke the quantum impurity solver
            curr_iter : current iteration number
            params    : user-supply input parameters
        """
        # update the parameters dynamically if necessary
        for key in params.iterkeys():
            self.__params[key] = params[key]

        # dump the parameters to external file
        if pyalps.mpi.rank == 0:
            f = open("ctqmc.param." + str(curr_iter), "w")
            for key in self.__params.iterkeys():
                print >> f, key, ":", self.__params[key]
            f.close()

        # mpi barrier
        pyalps.mpi.world.barrier()
        if pyalps.mpi.rank == 0 : print ("... Quantum Impurity Solver BEGIN ...")

        # just do it
        pyalps.cthyb.solve(self.__params)

        # mpi barrier
        if pyalps.mpi.rank == 0 : print ("... Quantum Impurity Solver END ...")
        pyalps.mpi.world.barrier()

    def final(self, curr_iter):
        """ finialize the quantum impurity solver
            curr_iter : current iteration number
        """
        if pyalps.mpi.rank != 0:
            return

        # rename the output data of the quantum impurity solver
        # only the master node can do it
        try:
            shutil.move("Gt.dat", "t_Gt.dat." + str(curr_iter))
            shutil.move("Gw.dat", "t_Gw.dat." + str(curr_iter))
            shutil.move("Ft.dat", "t_Ft.dat." + str(curr_iter))
            shutil.move("Fw.dat", "t_Fw.dat." + str(curr_iter))
            shutil.move("Sw.dat", "t_Sw.dat." + str(curr_iter))

            shutil.move("nnt.dat", "t_nnt.dat." + str(curr_iter))

            if self.__params["MEASURE_legendre"] == 1: shutil.move("Gtl.dat", "t_Gtl.dat." + str(curr_iter))
            if self.__params["MEASURE_legendre"] == 1: shutil.move("Gwl.dat", "t_Gwl.dat." + str(curr_iter))
            if self.__params["MEASURE_legendre"] == 1: shutil.move("Ftl.dat", "t_Ftl.dat." + str(curr_iter))
            if self.__params["MEASURE_legendre"] == 1: shutil.move("Fwl.dat", "t_Fwl.dat." + str(curr_iter))
            if self.__params["MEASURE_legendre"] == 1: shutil.move("Swl.dat", "t_Swl.dat." + str(curr_iter))
        except Exception:
            print ("NON-FATAL ERROR : shutil move error")
            pass

        # deal with the other files
        try:
            shutil.move("observables.dat", "observables.dat." + str(curr_iter))
            shutil.move("orders.dat", "orders.dat." + str(curr_iter))
            shutil.move("simulation.dat", "simulation.dat." + str(curr_iter))
        except Exception:
            print ("NON-FATAL ERROR : shutil move error")
            pass

    def extract(self, obs_type, p = 0):
        """ get the simulated results when the quantum impurity solver is finished
            obs_type : used to specify the object we want to extract
            p : param dictionary
        """
        if pyalps.mpi.rank != 0:
            return

        # open results.out.h5, which is in HDF5 format
        # only the master node can do it
        ar = pyalps.ngs.h5ar('results.out.h5','r')

        norbs = p["norbs"][0] # extract number of orbitals
        nband = p["nband"][0] # extract number of bands
        nffrq = p["nffrq"][0] # extract number of matsubara frequencies
        for case in switch(obs_type):
            if case(1): # extract the green's function in matsubara frequency space
                obs_data = numpy.zeros((norbs, nffrq), numpy.complex)
                for orb in range(norbs):
                    obs_data[orb] = ar['G_omega/'+str(orb)+'/mean/value']
                return obs_data
                break

            if case(3): # extract the self-energy function in matsubara frequency space
                obs_data = numpy.zeros((norbs, nffrq), numpy.complex)
                for orb in range(norbs):
                    obs_data[orb] = ar['S_omega/'+str(orb)+'/mean/value']
                return obs_data
                break

            if case(4): # extract the green's function in matsubara frequency space
                obs_data = numpy.zeros((norbs, nffrq), numpy.complex)
                for orb in range(norbs):
                    obs_data[orb] = ar['G_l_omega/'+str(orb)+'/mean/value']
                return obs_data
                break

            if case(6): # extract the self-energy function in matsubara frequency space
                obs_data = numpy.zeros((norbs, nffrq), numpy.complex)
                for orb in range(norbs):
                    obs_data[orb] = ar['S_l_omega/'+str(orb)+'/mean/value']
                return obs_data
                break

            if case(99): # extract the charge-charge correlation function
                # get orbital occupation
                occ_i = numpy.zeros(norbs, dtype = numpy.float)
                for i in range(norbs):
                    occ_i[i] = ar['simulation/results/density_' + str(i) + '/mean/value']
                # get charge-charge correlation function
                nnt_i_i = numpy.zeros((norbs, self.__params["N_nn"] + 1), dtype = numpy.float)
                for i in range(norbs):
                    nnt_i_i[i] = ar['simulation/results/nnt_' + str(i) + '_' + str(i) + '/mean/value']
                nnt_i_j = numpy.zeros((nband, self.__params["N_nn"] + 1), dtype = numpy.float)
                for i in range(0,norbs,2):
                    nnt_i_j[i/2] = ar['simulation/results/nnt_' + str(i + 1) + '_' + str(i) + '/mean/value']
                # sum them into nnt_tot
                nnt_tot = numpy.zeros((nband, self.__params["N_nn"] + 1), dtype = numpy.float)
                for i in range(0,norbs,2):
                    nnt_tot[i/2] = nnt_i_i[i] + 2*nnt_i_j[i/2] + nnt_i_i[i+1] - (occ_i[i] + occ_i[i+1])**2
                # symmetrize nnt_tot according to particle-hole symmetry
                if p["sy_ph"][0]:
                    for i in range(nband):
                        Symm.symm_over_time(nnt_tot[i])
                if p["sy_pm"][0]:
                    Symm.symm_over_band(nnt_tot)
                # allocate memory for \chi(i\omega)
                nnw_tot = numpy.zeros((nband, p["nbfrq"][0]), numpy.complex)
                # create a backup dict object
                paux = p.copy()
                paux["ntime"] = [self.__params["N_nn"] + 1, 'a dummy ntime']
                # calculate nnw_tot from nnt_tot by using fourier transformation
                fourier = Fourier(paux, 2)
                tmesh = Mesh.create_time_mesh(paux["beta"][0], paux["ntime"][0])
                bmesh = Mesh.create_boson_mesh(paux["beta"][0], paux["nbfrq"][0])
                for i in range(nband):
                    fourier.time_to_freq_C(tmesh, nnt_tot[i], bmesh, nnw_tot[i])
                del fourier, tmesh, bmesh, paux, nnt_tot, nnt_i_i, nnt_i_j, occ_i
                return nnw_tot
                break

            if case():  # default, could also just omit condition or 'if True'
                sys.exit("FATAL ERROR : this obs_type parameter is not supported yet.")

        # erase the HDF5 object
        del ar

if __name__ == '__main__':

    p = {
        'U'     : [4.00  , "Onsite Coulomb interaction"             ],
        'J'     : [1.00  , "Onsite Hund exchange interaction"       ],
        'mune'  : [11.761, "Chemical potential"                     ],
        'beta'  : [100.0 , "Inverse temperature"                    ],
        'norbs' : [10    , "Number of orbitals"                     ],
        'nband' : [5     , "Number of bands"                        ],
        'nspin' : [2     , "Number of spin projections"             ],
        'ntime' : [1024  , "Number of imaginary time points"        ],
        'nffrq' : [256   , "Number of fermionic frequencies"        ],
        'nbfrq' : [256   , "Number of bosonic frequencies"          ],
    }

    s_para = {
        'U'           : p["U"    ][0],     # Coulomb interaction
        'J'           : p["J"    ][0],     # Hund exchange interaction
        'N_ORBITALS'  : p["norbs"][0],     # number of orbitals
        'BETA'        : p["beta" ][0],     # inversion temperature
        'MU'          : p["mune" ][0],     # chemical potential
        'N_MATSUBARA' : p["nffrq"][0],     # number of matsubara frequencies
        'N_TAU'       : p["ntime"][0] - 1, # number of imaginary time points
    }

    hybf = numpy.zeros((p['norbs'][0],p['nffrq'][0]), dtype = numpy.complex)
    sloc = numpy.zeros((p['norbs'][0],p['nffrq'][0]), dtype = numpy.complex)
    htau = numpy.zeros((p['norbs'][0],p['ntime'][0]), dtype = numpy.float)

    # create imaginary time mesh
    tmesh = Mesh.create_time_mesh(p["beta"][0], p["ntime"][0])

    # create fermionic mesh
    fmesh = Mesh.create_fermion_mesh(p["beta"][0], p["nffrq"][0])

    # create bosonic mesh
    bmesh = Mesh.create_boson_mesh(p["beta"][0], p["nbfrq"][0])

    if pyalps.mpi.rank == 0:
        f = open('solver.hyb.in', 'r')
        for i in range(p['nffrq'][0]):
            spl = f.readline().split()
            hybf[0,i] = float(spl[2]) + 1j * float(spl[3])
            for j in range(1,p['norbs'][0]):
                hybf[j,i] = hybf[0,i]
        f.close()

    hybf = pyalps.mpi.broadcast(pyalps.mpi.world, hybf, 0)
    pyalps.mpi.world.barrier()

    fourier = Fourier(p, 1)
    for s in range(p['norbs'][0]):
        fourier.freq_to_time_F(tmesh, htau[s,:], fmesh, hybf[s,:])
    del fourier

    if pyalps.mpi.rank == 0:
        Dump.dump_htau_data("D.dat", tmesh, htau)
    
    pyalps.mpi.world.barrier()

    solver = AlpsSolver(s_para)
    solver.setup()
    _iter_ = 1
    solver.start(_iter_, s_para)
    solver.final(_iter_)
    gloc = solver.extract(4, p)
    gloc = pyalps.mpi.broadcast(pyalps.mpi.world, gloc, 0)
    pyalps.mpi.world.barrier()

    Symm.symm_over_band(gloc)
    mune = p['mune'][0]
    for s in range(p['norbs'][0]):
        sloc[s,:] = 1j*fmesh + mune - hybf[s,:] - 1.0 / gloc[s,:]
    Symm.symm_over_band(sloc)

    if pyalps.mpi.rank == 0:
        Dump.dump_freq_data3("s_hybf.dat." + str(_iter_), fmesh, hybf)
        Dump.dump_freq_data3("s_gloc.dat." + str(_iter_), fmesh, gloc)
        Dump.dump_freq_data3("s_sloc.dat." + str(_iter_), fmesh, sloc)
