#!/usr/bin/env python

# python standard library
import os
import sys
import time
import math
import shutil
import ctypes

# python 3rd-party library
import numpy
import numpy.linalg
import scipy.weave
import scipy.interpolate

# python binding for ALPS
import pyalps.mpi
import pyalps.cthyb

def lda_edmft_config():
    """ setup the config parameters for the lda + edmft code.
        note: don't modify this function. if you want to change
        the config parameters, please edit the in.param file
    """
    # setup the parameters
    # param : a dict object, used to record all the necessary parameters
    param = {
        'dmft'  : [1     , "dmft (1), dmft + u(w) (2), edmft (3)"   ],
        'sedc'  : [1     , "zero (1), FLL         (2), KAV   (3)"   ],
        'semu'  : [1     , "impurity mune (1), lattice mune (2)"    ],
        'U'     : [2.20  , "Onsite Coulomb interaction"             ],
        'J'     : [0.00  , "Onsite Hund exchange interaction"       ],
        'Z'     : [1.00  , "Renormalization factor for Hamiltonian" ],
        'beta'  : [100.  , "Inverse temperature"                    ],
        'mune'  : [1.10  , "Chemical potential"                     ],
        'occup' : [1.00  , "Impurity occupation number"             ],
        'alpha' : [0.70  , "Mixing parameter"                       ],
        'niter' : [40    , "Maximum iteration number"               ],
        'norbs' : [2     , "Number of orbitals"                     ],
        'nband' : [1     , "Number of bands"                        ],
        'nkpts' : [1000  , "Number of k-points along x axis"        ],
        'ntime' : [1024  , "Number of imaginary time points"        ],
        'nffrq' : [1024  , "Number of fermionic frequencies"        ],
        'nbfrq' : [1024  , "Number of bosonic frequencies"          ],
        'start' : [True  , "Whether start from scratch"             ],
        'sy_ph' : [True  , "Whether enforce symmetrization (PH)"    ],
        'sy_pm' : [True  , "Whether enforce symmetrization (PM)"    ],
    }

    # read config parameters from in.param file
    # only master node can do it
    if ( pyalps.mpi.rank == 0 ) :
        # check whether file in.param exists
        if os.path.isfile("in.param"):
            in_param = {}
            # read data, pay attention to their data type
            f_in = open("in.param", "r")
            for line in f_in:
                A = line.split()
                if   A[2] == 'float':
                    in_param[A[0]] = float(A[4])
                elif A[2] == 'int'  :
                    in_param[A[0]] = int(A[4])
                elif A[2] == 'bool' :
                    in_param[A[0]] = eval(A[4])
                else :
                    sys.exit("FATAL ERROR : wrong data type")
            # close the file handler
            f_in.close()

            # now the config parameters are stored in in_param dict,
            # we need to transfer them to param dict
            for key in in_param.iterkeys():
                param[key][0] = in_param[key]
        # if in.param does not exist, we have to stop the code
        else:
            sys.exit("FATAL ERROR : in.param does not exist")

    # broadcast the param dict to all children processes
    param = pyalps.mpi.broadcast(pyalps.mpi.world, param, 0)
    pyalps.mpi.world.barrier()

    # dump the parameters into the stdout
    if pyalps.mpi.rank == 0:
        print (">>> Listing parameters:")
        for k in sorted(param.keys()):
            print ("    " + param[k][1].ljust(42) + " : " + k.rjust(5) + " : " + str(param[k][0]).rjust(5))
        print ("")

    return param

def lda_edmft_sconfig(p):
    """ extract some key parameters which are useful for the impurity solver
        p : param dictionary
    """
    # note: if U_MATRIX is specified, then the U and J parameters are
    # skipped. But we still need them to determine the double counting term.
    # on the other hand, if MU_VECTOR is specified, then the MU parameter
    # will be skipped as well.
    s_p = {
        'U'           : p["U"    ][0],     # Coulomb interaction
        'J'           : p["J"    ][0],     # Hund exchange interaction
        'MU'          : p["mune" ][0],     # chemical potential
        'BETA'        : p["beta" ][0],     # inversion temperature
        'N_ORBITALS'  : p["norbs"][0],     # number of orbitals
        'N_MATSUBARA' : p["nffrq"][0],     # number of matsubara frequencies
        'N_TAU'       : p["ntime"][0] - 1, # number of imaginary time points
        'DELTA'       : "D.dat"      ,     # file which is used to store the hybridization function
        'U_MATRIX'    : "umat.dat"   ,     # file which is used to store the u-matrix
        'MU_VECTOR'   : "mune.dat"   ,     # file which is used to store the crystal field splitting
    }

    # determine whether we are dealing with frequency-dependent interaction
    if p["dmft"][0] > 1:
        s_p['RET_INT_K'] = "K.dat"

    return s_p

def lda_edmft_header():
    """ print the header to standard terminal
    """
    # only the master node can do it
    if pyalps.mpi.rank != 0:
        return

    print ("LE :: LDA + EDMFT")
    print ("Ver.2014/06/17D (accelerated by mpi + scipy.weave)")
    print ("Developed by Li HUANG (li.huang@unifr.ch)")
    print ("Department of Physics, Fribourg University, Switzerland")
    print ("")

    print (">>> Starting in " + time.asctime())
    print (">>> Using %i processors, current rank is %i" % (pyalps.mpi.size, pyalps.mpi.rank))
    print ("")

def lda_edmft_footer(dmft_type):
    """ print the footer to standard terminal
        dmft_type : the calculating scheme, it may be (1) lda + dmft / (2) lda + dmft + u(w) / (3) lda + edmft
    """
    # only the master node can do it
    if pyalps.mpi.rank != 0:
        return

    for case in switch( dmft_type ):
        if case(3):
            print ("-"*60)
            print (">>> LDA + EDMFT ITER : Finished")
            print ("-"*60)
            break

        if case(2):
            print ("-"*60)
            print (">>> LDA + DMFT + U(W) ITER : Finished")
            print ("-"*60)
            break

        if case(1):
            print ("-"*60)
            print (">>> LDA + DMFT ITER : Finished")
            print ("-"*60)
            break
    print (">>> Ending in " + time.asctime())

def lda_edmft_iter(dmft_type, niter, curr_iter):
    """ print the iteration information to standard terminal
        dmft_type : the calculating scheme, it may be (1) lda + dmft / (2) lda + dmft + u(w) / (3) lda + edmft
        niter     : number of iterations
        curr_iter : current iteration number
    """
    # only the master node can do it
    if pyalps.mpi.rank != 0:
        return

    for case in switch( dmft_type ):
        if case(3):
            print ("-"*60)
            print (">>> LDA + EDMFT ITER : current => %i  maximum => %i" % (curr_iter, niter))
            print ("-"*60)
            print ("")
            break

        if case(2):
            print ("-"*60)
            print (">>> LDA + DMFT + U(W) ITER : current => %i  maximum => %i" % (curr_iter, niter))
            print ("-"*60)
            print ("")
            break

        if case(1):
            print ("-"*60)
            print (">>> LDA + DMFT ITER : current => %i  maximum => %i"  % (curr_iter, niter))
            print ("-"*60)
            print ("")
            break

def lda_edmft_fdirac(beta, omega):
    """ to calculate the Fermi-Dirac distribution function
        beta  : inverse temperature
        omega : complex frequency
    """
    om = omega.real # extract the real frequency
    value = 0.0
    if   beta * om >= +600.0 :
        value = 0.0
    elif beta * om <= -600.0 :
        value = 1.0
    else :
        value = 1.0 / ( 1.0 + math.exp( beta * om ) )
    return value

def lda_edmft_occupy(norbs, nkpts, nwcut, beta, maux, fmesh, eigs, einf):
    """ to determine the electron occupation for a given chemical potential
        norbs : number of orbitals
        nkpts : number of k-points
        nwcut : number of matsubara frequency points
        beta  : inverse temperature 
        maux  : pseudo chemical potential
        fmesh : mesh of matsubara frequency
        eigs  : eigenvalues of H(k) + \Sigma(i\omega)
        einf  : eigenvalues of H(k) + \Sigma(i\omega -> \infty)
    """
    # build complex mesh
    cmesh = 1j * fmesh + maux

    # calculate tail correction to occupation numbers
    focc = numpy.zeros((norbs, nwcut), dtype = numpy.complex)
    for ib in range(norbs):
        for iw in range(nwcut):
            focc[ib,iw] = focc[ib,iw] + numpy.sum( 1.0 / ( cmesh[iw] - eigs[ib,iw,:] ) - 1.0 / ( cmesh[iw] - einf[ib,:] ) )

    # evaulate zocc
    zocc = numpy.zeros(norbs, dtype = numpy.complex)
    for ib in range(norbs):
        zocc[ib] = numpy.sum( focc[ib,:] ) * ( 2.0 / beta )

    # calculate occupation numbers by fermi-dirac distribution
    for ik in range(nkpts):
        for ib in range(norbs):
            zocc[ib] = zocc[ib] + lda_edmft_fdirac( beta, einf[ib,ik] - maux )

    # release the memory and return what we need
    zsum = numpy.sum(zocc.real) / float(nkpts)
    del zocc, focc
    return zsum

def lda_edmft_slevel(norbs, nband, nkpts, nwish, beta, fmesh, hk, slat, sedc):
    """ to determine the chemical potential for a given occupation number
        norbs : number of orbitals
        nband : number of bands
        nkpts : number of k-points
        nwish : expected occupation number
        beta  : inverse temperature 
        fmesh : mesh of matsubara frequency
        hk    : LDA hamiltonian
        slat  : self-energy function without double counting term
        sedc  : double counting term
    """
    nwcut = 256 # number of matsubara frequency points considered

    # a auxiliary array, used to build self-energy matrix
    saux = numpy.zeros(nband, dtype = numpy.complex)

    # eigenvalues of H(k) + \Sigma(i\omega)
    eigs = numpy.zeros((norbs, nwcut, nkpts), dtype = numpy.complex)

    # eigenvalues of H(k) + \Sigma(i\omega -> \infty)
    einf = numpy.zeros((norbs, nkpts), dtype = numpy.complex)

    # evaluate eigs
    for ik in range(nkpts):
        for iw in range(nwcut):
            # get the spin up part of self-energy
            for ib in range(nband):
                saux[ib] = slat[2*ib,iw,ik] - sedc[2*ib]
            # diagonalize the hamiltonian
            ek = numpy.linalg.eigvalsh( hk[:,:,ik] + numpy.diag(saux) )
            # setup the eigenvalues
            for ib in range(nband):
                eigs[2*ib,iw,ik]   = ek[ib] # spin up
                eigs[2*ib+1,iw,ik] = ek[ib] # spin down
    pyalps.mpi.world.barrier()

    # evaluate einf
    for ik in range(nkpts):
        # get the spin up part of self-energy
        for ib in range(nband):
            saux[ib] = slat[2*ib,-1,ik] - sedc[2*ib]
        # diagonalize the hamiltonian
        ek = numpy.linalg.eigvalsh( hk[:,:,ik] + numpy.diag(saux) )
        # setup the eigenvalues
        for ib in range(nband):
            einf[2*ib,ik]   = ek[ib] # spin up
            einf[2*ib+1,ik] = ek[ib] # spin down
    pyalps.mpi.world.barrier()

    # search the fermi level using the bisection algorithm
    lmin = -50.0 # left boundary for the fermi level
    rmax = +50.0 # right boundary for the fermi level
    maux = ( lmin + rmax ) * 0.5 # temporary fermi level
    naux = 0.0   # temporary occupation number
    step = 0     # iteration counter
    while abs( naux - nwish) > 0.000001 :
        # calculate the chemical potential
        naux = lda_edmft_occupy(norbs, nkpts, nwcut, beta, maux, fmesh, eigs, einf)
        # increase the counter
        step = step + 1
        # print the useful information
        if step > 100 :
            sys.exit("FATAL ERROR : can not locate the correct chemical potential")
        else :
            if pyalps.mpi.rank == 0 :
                print ("    step : %3i  mune : %8.7s  Ncurr: %8.7s  Nwish : %8.7s" % (step, maux, naux, nwish))
        # adjust the boundaries
        # if fermi level is too large, so we cut down the right boundary
        # if fermi level is too small, so we move up the left boundary
        if naux > nwish :
            rmax = maux
        else :
            lmin = maux
        # adjust the fermi level
        maux = ( lmin + rmax ) * 0.5
    pyalps.mpi.world.barrier()

    # release the memory
    del saux, eigs, einf

    # return the global chemical potential
    return maux

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
    def symm_keep_real(freq_data):
        """ get rid of the imaginary part of the input data
            freq_data : matsubara frequency data
        """
        # get number of elements
        n = len(freq_data)
        for i in range(n):
            freq_data[i] = freq_data[i].real + 0j

    @staticmethod
    def symm_keep_imag(freq_data):
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
            p     : param dictionary
            beta  : inversion temperature
            ntime : number of imaginary time points
            nffrq : number of frequencies for fermionic mesh
            nbfrq : number of frequencies for bosonic mesh
        """
        self.__beta  = p["beta" ][0]
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
    def dump_time_data(file_name, time_mesh, time_data):
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
    def dump_freq_data(file_name, freq_mesh, freq_data):
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
    def dump_umat_data(file_name, umat, shift = 0.0):
        """ dump Coulomb interaction matrix to the disk file
            file_name : file name
            umat      : Coulomb interaction matrix
            shift     : shift for the Coulomb interaction matrix
        """
        norbs = len(umat)
        f = open(file_name, "w")
        for orb1 in range(norbs):
            for orb2 in range(norbs):
                if orb1 != orb2: # non-diagonal elements
                    f.write(" %16.8s " % ( umat[orb1,orb2] - shift ))
                else:            # diagonal elements
                    f.write(" %16.8s " % ( umat[orb1,orb2] ))
            f.write("\n") # start a new line
        f.close()

    @staticmethod
    def dump_mune_data(file_name, mune, eimp, sedc):
        """ dump the orbital chemical potential to the disk file
            file_name : file name
            mune      : the chemical potential
            eimp      : impurity level
            sedc      : double counting term
        """
        mu_vector = mune - eimp + sedc
        norbs = len(mu_vector)
        f = open(file_name, "w")
        for orb in range(norbs):
            f.write(" %16.8s " % mu_vector[orb].real)
        f.write("\n") # start a new line
        f.close()

class Fourier(object):
    """ class Fourier is used to perform forward and backward fourier
        transformations between imaginary time and matsubara frequency
        spaces.
    """
    def __init__(self, p, ftype):
        """ class constructor for the Fourier
            p     : param dictionary
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

class Lattice(object):
    """ class Lattice is used to define a LDA Hamiltonian.
    """
    @staticmethod
    def create_hk(nband, norbs, nkpts, Z):
        """ build the lda hamiltonian and corresponding impurity level
            nband : number of bands
            norbs : number of orbitals
            nkpts : number of k-points
            Z     : renormalization factor for hamiltonian
        """
        # impurity level
        eimp = numpy.zeros(norbs, dtype = numpy.complex)

        # lda hamiltonian
        hk = numpy.zeros((nband, nband, nkpts), dtype = numpy.complex)

        # read hamiltonian data
        if pyalps.mpi.rank == 0:

            # open data file : in.hk
            f = open("in.hk", "r")

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

            # build the impurity level
            for ikpt in range(nkpts):
                for i in range(nband):
                    eimp[2*i] = eimp[2*i] + hk[i,i,ikpt] / float(nkpts) # spin up part
                    eimp[2*i+1] = eimp[2*i+1] + hk[i,i,ikpt] / float(nkpts) # spin down part

        # broadcast the hk and eimp data from master node to children nodes
        eimp = pyalps.mpi.broadcast(pyalps.mpi.world, eimp, 0)
        hk = pyalps.mpi.broadcast(pyalps.mpi.world, hk, 0)
        pyalps.mpi.world.barrier()

        return eimp, hk

    @staticmethod
    def create_vk(nband, nkpts):
        """ build the k-dependent general bare interaction
            nband : number of bands
            nkpts : number of k-points
        """
        vk = numpy.zeros((nband,nkpts), dtype = numpy.complex)

        # read vk data
        if pyalps.mpi.rank == 0:
            # open data file : in.vk
            f = open("in.vk", "r")

            # loop over k-points
            for k in range(nkpts):
                f.readline()
                for s in range(nband):
                    spl = f.readline().split()
                    vk[s,k] = float( spl[0] ) + float( spl[1] ) * 1j

            # close data file
            f.close()
        
        # broadcast the vk data from master node to children nodes
        vk = pyalps.mpi.broadcast(pyalps.mpi.world, vk, 0)
        pyalps.mpi.world.barrier()

        return vk

    @staticmethod
    def create_umat(norbs):
        """ build the onsite Coulomb interaction matrix
            nband : number of orbitals
        """
        umat = numpy.zeros((norbs,norbs), dtype = numpy.float)

        # read umat data
        if pyalps.mpi.rank == 0:
            # open data file : in.umat
            f = open("in.umat", "r")

            # loop over norbs
            for orb1 in range(norbs):
                spl = f.readline().split()
                for orb2 in range(norbs):
                    umat[orb1,orb2] = float( spl[orb2] )

            # close data file
            f.close()

        # broadcast the umat data from master node to children nodes
        umat = pyalps.mpi.broadcast(pyalps.mpi.world, umat, 0)
        pyalps.mpi.world.barrier()

        return umat

class Mixer(object):
    """ class Mixer is used to mix the old and new vectors to produce a
        newer one.
    """
    def __init__(self, p):
        """ class constructor for the Mixer
            p : param dictionary
            alpha : mixing parameter
        """
        self.__alpha = p["alpha"][0]

    def linear_mixing(self, vec1, vec2):
        """ perform the simplest linear mixing
            vec1 : old vector on input
            vec2 : new vector on input
        """
        return vec1 * ( 1 - alpha ) + vec2 * alpha

    @staticmethod
    def linear_mixing(alpha, vec1, vec2):
        """ perform the simplest linear mixing
            alpha : mixing parameter
            vec1  : old vector on input
            vec2  : new vector on input
        """
        return vec1 * ( 1 - alpha ) + vec2 * alpha

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
            'MAX_TIME'                  : 1200       ,
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
            'SEED'                      : int( time.time() ) + ( pyalps.mpi.rank + 3 ) * 1981,
        }

        # setup additional keys and values
        for key in params.iterkeys():
            self.__params[key] = params[key]

    def setup(self, mode):
        """ prepare necessary inputs and perform check to ensure everything is OK
            mode : control the running mode of impurity solver, if mode > 1, the
                   Coulomb interaction is frequency-dependent, we need to check
                   whether the file K.dat exists 
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

        # check whether we are dealing with frequency-dependent interaction
        if mode == 1:
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
        if pyalps.mpi.rank == 0 : print ("... Quantum Impurity Solver START ...")

        # just do it
        pyalps.cthyb.solve(self.__params)

        # mpi barrier
        if pyalps.mpi.rank == 0 : print ("... Quantum Impurity Solver STOP ...")
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
            # we can not move results.out.h5, because we still need it in extract()
            shutil.copy("results.out.h5", "results.out.h5." + str(curr_iter))
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

            if case(98): # extract the charge-charge correlation function
                # get orbital occupation
                occ_i = numpy.zeros(norbs, dtype = numpy.float)
                for i in range(norbs):
                    occ_i[i] = ar['simulation/results/density_' + str(i) + '/mean/value']
                # get charge-charge correlation function
                nnt_i_j = numpy.zeros((norbs, norbs, self.__params["N_nn"] + 1), dtype = numpy.float)
                for j in range(norbs):
                    for i in range(j,norbs):
                        nnt_i_j[i,j] = ar['simulation/results/nnt_' + str(i) + '_' + str(j) + '/mean/value']
                        nnt_i_j[j,i] = ar['simulation/results/nnt_' + str(i) + '_' + str(j) + '/mean/value']
                # sum them into nnt_tot
                nnt_tot = numpy.zeros(self.__params["N_nn"] + 1, dtype = numpy.float)
                for i in range(norbs):
                    for j in range(norbs):
                        nnt_tot = nnt_tot + nnt_i_j[i,j]
                nnt_tot = nnt_tot - ( numpy.sum(occ_i) )**2
                # symmetrize nnt_tot according to the symmetry
                Symm.symm_over_time(nnt_tot)
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
                    fourier.time_to_freq_C(tmesh, nnt_tot, bmesh, nnw_tot[i])
                del fourier, tmesh, bmesh, paux, nnt_tot, nnt_i_j, occ_i
                return nnw_tot
                break

            if case(100): # extract the occupation number
                # get orbital occupation
                occ_i = numpy.zeros(norbs, dtype = numpy.float)
                for i in range(norbs):
                    occ_i[i] = ar['simulation/results/density_' + str(i) + '/mean/value']
                # get total occupation number
                obs_data = numpy.sum(occ_i)
                del occ_i
                return obs_data
                break

            if case():  # default, could also just omit condition or 'if True'
                sys.exit("FATAL ERROR : this obs_type parameter is not supported yet.")

        # erase the HDF5 object
        del ar

class LocVar(object):
    """ class LocVar is used to store and manipulate the local objects.
    """
    def __init__(self, p):
        """ class constructor for the LocVar
            p     : param dictionary
            nband : number of bands
            norbs : number of orbitals
            nffrq : number of fermionic frequencies
            nbfrq : number of bosonic frequencies
            beta  : inverse temperature
        """
        self.__nband = p["nband"][0]
        self.__norbs = p["norbs"][0]
        self.__nffrq = p["nffrq"][0]
        self.__nbfrq = p["nbfrq"][0]
        self.__beta  = p["beta"] [0]

        # allocate memory for G(i\omega), \Sigma(i\omega), \mathcal{G}(i\omega)
        self.__gloc = numpy.zeros((self.__norbs, self.__nffrq), dtype = numpy.complex)
        self.__sloc = numpy.zeros((self.__norbs, self.__nffrq), dtype = numpy.complex)
        self.__bloc = numpy.zeros((self.__norbs, self.__nffrq), dtype = numpy.complex)

        # allocate memory for double counting term
        self.__sedc = numpy.zeros((self.__norbs              ), dtype = numpy.complex)

        # allocate memory for W(i\omega), \Pi(i\omega), \mathcal{U}(i\omega)
        self.__wloc = numpy.zeros((self.__nband, self.__nbfrq), dtype = numpy.complex)
        self.__ploc = numpy.zeros((self.__nband, self.__nbfrq), dtype = numpy.complex)
        self.__uloc = numpy.zeros((self.__nband, self.__nbfrq), dtype = numpy.complex)

        # allocate memory for \chi_{loc}, which is used to calculate wloc
        self.__cloc = numpy.zeros((self.__nband, self.__nbfrq), dtype = numpy.complex)

        # allocate memory for G(i\omega), \Sigma(i\omega), \Pi(i\omega), used as backup
        self.__gloc_save = numpy.zeros((self.__norbs, self.__nffrq), dtype = numpy.complex)
        self.__sloc_save = numpy.zeros((self.__norbs, self.__nffrq), dtype = numpy.complex)
        self.__ploc_save = numpy.zeros((self.__nband, self.__nbfrq), dtype = numpy.complex)

    @property
    def gloc(self):
        """ getter for __gloc property
        """
        return self.__gloc

    @property
    def sloc(self):
        """ getter for __sloc property
        """
        return self.__sloc

    @property
    def bloc(self):
        """ getter for __bloc property
        """
        return self.__bloc

    @property
    def sedc(self):
        """ getter for __sedc property
        """
        return self.__sedc

    @property
    def wloc(self):
        """ getter for __wloc property
        """
        return self.__wloc

    @property
    def ploc(self):
        """ getter for __ploc property
        """
        return self.__ploc

    @property
    def uloc(self):
        """ getter for __uloc property
        """
        return self.__uloc

    @property
    def cloc(self):
        """ getter for __cloc property
        """
        return self.__cloc

    @gloc.setter
    def gloc(self, gloc):
        """ setter for __gloc property
        """
        self.__gloc = gloc

    @sloc.setter
    def sloc(self, sloc):
        """ setter for __sloc property
        """
        self.__sloc = sloc

    @bloc.setter
    def bloc(self, bloc):
        """ setter for __bloc property
        """
        self.__bloc = bloc

    @sedc.setter
    def sedc(self, sedc):
        """ setter for __sedc property
        """
        self.__sedc = sedc

    @wloc.setter
    def wloc(self, wloc):
        """ setter for __wloc property
        """
        self.__wloc = wloc

    @ploc.setter
    def ploc(self, ploc):
        """ setter for __ploc property
        """
        self.__ploc = ploc

    @uloc.setter
    def uloc(self, uloc):
        """ setter for __uloc property
        """
        self.__uloc = uloc

    @cloc.setter
    def cloc(self, cloc):
        """ setter for __cloc property
        """
        self.__cloc = cloc

    def reload_sloc(self):
        """ reload local self-energy function from in.sloc file
        """
        # read data from in.sloc file, only master node can do it
        if pyalps.mpi.rank == 0 :
            print ("    reloading sloc data from in.sloc ...")
            if os.path.isfile("in.sloc"):
                f_sig = open("in.sloc", 'r')
                for orb in range(self.__norbs):
                    for f in range(self.__nffrq):
                        line = f_sig.readline()
                        spl = map(float, line.split())
                        self.__sloc[orb][f] = spl[3] + 1j * spl[4] # build complex numbers
                    f_sig.readline() # skip two empty lines
                    f_sig.readline()
                f_sig.close()
            else:
                sys.exit("FATAL ERROR : in.sloc does not exist")

        # broadcast the self.__sloc data to children processes
        self.__sloc = pyalps.mpi.broadcast(pyalps.mpi.world, self.__sloc, 0)
        pyalps.mpi.world.barrier()

        # copy sloc to sloc_save
        self.__sloc_save = self.__sloc.copy()

    def reload_ploc(self):
        """ reload local bosonic self-energy function from in.ploc file
        """
        # read data from in.ploc file, only master node can do it
        if pyalps.mpi.rank == 0 :
            print ("    reloading ploc data from in.ploc ...")
            if os.path.isfile("in.ploc"):
                f_prod = open("in.ploc", 'r')
                for orb in range(self.__nband):
                    for b in range(self.__nbfrq):
                        line = f_prod.readline()
                        spl = map(float, line.split())
                        self.__ploc[orb][b] = spl[3] + 1j * spl[4] # build a complex number
                    f_prod.readline() # skip two empty lines
                    f_prod.readline()
                f_prod.close()
            else:
                sys.exit("FATAL ERROR : in.ploc does not exist")

        # broadcast the self.__ploc data to children processes
        self.__ploc = pyalps.mpi.broadcast(pyalps.mpi.world, self.__ploc, 0)
        pyalps.mpi.world.barrier()

        # copy ploc to ploc_save
        self.__ploc_save = self.__ploc.copy()

    def mix_gloc(self, curr_iter, alpha):
        """ check the convergence of gloc, and then mix it with gloc_save
            curr_iter : current iteration number
            alpha : mixing parameter
        """
        error = (numpy.absolute(self.__gloc - self.__gloc_save)).max()
        # for the first iteration, we do not mix them
        if curr_iter > 1 :
            self.__gloc  = Mixer.linear_mixing(alpha, self.__gloc_save, self.__gloc)
        self.__gloc_save = self.__gloc.copy()
        return error

    def mix_sloc(self, curr_iter, alpha):
        """ check the convergence of sloc, and then mix it with sloc_save
            curr_iter : current iteration number
            alpha : mixing parameter
        """
        error = (numpy.absolute(self.__sloc - self.__sloc_save)).max()
        # for the first iteration, we do not mix them
        if curr_iter > 1 :
            self.__sloc  = Mixer.linear_mixing(alpha, self.__sloc_save, self.__sloc)
        self.__sloc_save = self.__sloc.copy()
        return error

    def mix_ploc(self, curr_iter, alpha):
        """ check the convergence of ploc, and then mix it with ploc_save
            curr_iter : current iteration number
            alpha : mixing parameter
        """
        error = (numpy.absolute(self.__ploc - self.__ploc_save)).max()
        # for the first iteration, we do not mix them
        if curr_iter > 1 :
            self.__ploc  = Mixer.linear_mixing(alpha, self.__ploc_save, self.__ploc)
        self.__ploc_save = self.__ploc.copy()
        return error

    def cal_ksum(self, v_lat):
        """ core subroutine used to calculate the k-summing
            note : this function is used internally
            v_lat : k-dependent function
        """
        return numpy.mean(v_lat)

    def cal_gloc_by_ksum(self, glat):
        """ calculate G_{loc} by G_{lat}
            glat  : lattice green's function
        """
        for s in range(self.__norbs):
            for f in range(self.__nffrq):
                self.__gloc[s,f] = self.cal_ksum( glat[s,f,:] )

    def cal_wloc_by_ksum(self, wlat):
        """ calculate W_{loc} by W_{lat}
            wlat  : lattice screening function
        """
        for s in range(self.__nband):
            for b in range(self.__nbfrq):
                self.__wloc[s,b] = self.cal_ksum( wlat[s,b,:] )

    def cal_sloc_by_ksum(self, slat):
        """ calculate \Sigma_{loc} by \Sigma_{lat}
            slat  : lattice self-energy function
        """
        for s in range(self.__norbs):
            for f in range(self.__nffrq):
                self.__sloc[s,f] = self.cal_ksum( slat[s,f,:] )

    def cal_ploc_by_ksum(self, plat):
        """ calculate \Pi_{loc} by \Pi_{lat}
            plat  : lattice bosonic self-energy function
        """
        for s in range(self.__nband):
            for b in range(self.__nbfrq):
                self.__ploc[s,b] = self.cal_ksum( plat[s,b,:] )

    def cal_bloc_by_dyson(self):
        """ calculate \mathcal{G}_{loc} by using the dyson equation
        """
        self.__bloc = 1.0 / ( 1.0 / self.__gloc + self.__sloc )

    def cal_uloc_by_dyson(self):
        """ calculate \mathcal{U}_{loc} by using the dyson equation
        """
        self.__uloc = 1.0 / ( 1.0 / self.__wloc + self.__ploc )

    def cal_sloc_by_dyson(self):
        """ calculate \Sigma_{loc} by using the dyson equation
        """
        self.__sloc = 1.0 / self.__bloc - 1.0 / self.__gloc

    def cal_ploc_by_dyson(self):
        """ calculate \Pi_{loc} by using the dyson equation
        """
        self.__ploc = 1.0 / self.__uloc - 1.0 / self.__wloc

    def cal_wloc_by_dyson(self):
        """ calculate W_{loc} by using the dyson equation
        """
        self.__wloc = self.__uloc - self.__uloc * self.__cloc * self.__uloc

    def cal_ploc_by_direct(self):
        """ calculate \Pi_{loc} via \chi_{loc} and \mathcal{U}_{loc}
            note: it is equivalent to cal_ploc_by_dyson()
        """
        self.__ploc = self.__cloc / ( self.__uloc * self.__cloc - 1.0 )

    def cal_sedc_by_fll(self, U, J, occupy):
        """ calculate the double counting term using FLL scheme
            U : Coulomb interaction
            J : Hund's exchange interaction
            occupy : occupation number
        """
        for s in range(self.__norbs):
            self.__sedc[s] = U * ( occupy - 0.5 ) - J * ( occupy / 2.0 - 0.5 )

    def cal_sedc_by_kav(self, U, J, occupy):
        """ calculate the double counting term using KAV scheme
            U : Coulomb interaction
            J : Hund's exchange interaction
            occupy : occupation number
        """
        ubar = U + float( self.__nband - 1 ) * ( 2.0 * U - 5.0 * J )
        for s in range(self.__norbs):
            self.__sedc[s] = ubar * occupy / float(self.__norbs)

    def fix_sloc_by_dyson(self, fmesh, mune, eimp, hybf):
        """ fix the real part of self-energy function using the Dyson equation
            fmesh : frequency mesh
            mune  : chemical potential
            eimp  : impurity level
            hybf  : hybridization function
        """
        # calculate auxiliary self-energy function using dyson equation
        saux = numpy.zeros((self.__norbs, self.__nffrq), dtype = numpy.complex)
        for s in range(self.__norbs):
            for f in range(self.__nffrq):
                saux[s,f] = fmesh[f]*1j + mune - eimp[s] + self.__sedc[s] - 1.0 / self.__gloc[s,f] - hybf[s,f]

        # combine saux and self.__sloc
        # for the imaginary part: use self.__sloc's contribution
        # for the real part
        # [0,half], use saux's contribution
        # [half,half+half/2], intermediate zone
        # [half+half/2,nffrq], use self.__sloc's contribution
        half = self.__nffrq / 2
        f_start = half - 1
        f_end = half + half / 2
        for s in range(self.__norbs):
            for f in range(half):
                self.__sloc[s,f] = saux[s,f].real + self.__sloc[s,f].imag * 1j
            s_start = saux[s,f_start].real
            s_end = self.__sloc[s,f_end].real
            f_begin = fmesh[f_start]
            for f in range(f_start, f_end+1):
                real_p = s_start + (s_end - s_start) * ( (fmesh[f] - f_begin) / (fmesh[f_end] - f_begin) )**2
                imag_p = self.__sloc[s,f].imag
                self.__sloc[s,f] = real_p + imag_p * 1j
       
        # smooth the real part of self-energy function 
        for s in range(self.__norbs):
            stmp = self.__sloc[s,half:self.__nffrq].real
            for smth in range(16):
                Symm.symm_smth_data(stmp, 8)
            for f in range(half,self.__nffrq):
                self.__sloc[s,f] = stmp[f-half] + self.__sloc[s,f].imag * 1j

        # release memory
        del saux

    def sym_gloc(self, sy_pm, sy_ph):
        """ symmetrize the gloc over spin, we also remove its real part
            sy_pm : spin polarization or non spin polarization
            sy_ph : particle-hole symmetry or not
        """
        if sy_pm: Symm.symm_over_band(self.__gloc)
        if sy_ph: Symm.symm_keep_imag(self.__gloc)

    def sym_sloc(self, sy_pm):
        """ symmetrize the sloc over spin
            sy_pm : spin polarization or non spin polarization
        """
        if sy_pm: Symm.symm_over_band(self.__sloc)

    def sym_bloc(self, sy_pm):
        """ symmetrize the bloc over spin
            sy_pm : spin polarization or non spin polarization
        """
        if sy_pm: Symm.symm_over_band(self.__bloc)

    def sym_wloc(self, sy_pm):
        """ ensure \W_{loc} to satisfy some basis properties
            sy_pm : spin polarization or non spin polarization
        """
        # W_{loc} should be real and not less than 0
        #for orb in range(self.__nband):
        #    for b in range(self.__nbfrq):
        #        if self.__wloc[orb,b].real < 0.0 :
        #            self.__wloc[orb,b] = 0.0 + 0.0j
        #            #sys.exit("FATAL ERROR : wrong wloc data")
        Symm.symm_keep_real(self.__wloc)
        if sy_pm: Symm.symm_over_band(self.__wloc)

    def sym_ploc(self, sy_pm):
        """ ensure \Pi_{loc} to satisfy some basis properties
            sy_pm : spin polarization or non spin polarization
        """
        # \Pi_{loc} should be real and not larger than 0
        for orb in range(self.__nband):
            for b in range(self.__nbfrq):
                if self.__ploc[orb,b].real > 0.0 :
                    self.__ploc[orb,b] = 0.0 + 0.0j
                    #sys.exit("FATAL ERROR : wrong ploc data")
        for orb in range(self.__nband):
            bmesh = Mesh.create_boson_mesh(self.__beta, self.__nbfrq)
            fit = scipy.interpolate.InterpolatedUnivariateSpline(bmesh[6:20], self.__ploc[orb,6:20].real, k = 2)
            self.__ploc[orb,0:6] = fit( bmesh[0:6] ) + 1j * self.__ploc[orb,0:6].imag
        Symm.symm_keep_real(self.__ploc)
        if sy_pm: Symm.symm_over_band(self.__ploc)

    def sym_uloc(self, sy_pm):
        """ ensure \U_{loc} to satisfy some basis properties
            sy_pm : spin polarization or non spin polarization
        """
        Symm.symm_keep_real(self.__uloc)
        if sy_pm: Symm.symm_over_band(self.__uloc)

    def sym_cloc(self, sy_pm):
        """ ensure \chi_{loc} to satisfy some basis properties
            sy_pm : spin polarization or non spin polarization
        """
        # \chi_{loc} should be real
        #for orb in range(self.__nband):
        #    scale = self.__cloc[orb,0] * self.__uloc[orb,0]
        #    if scale > 1.0:
        #        if pyalps.mpi.rank == 0:
        #            print ("    sym_cloc: rescale cloc %i -> %16.8s" % (orb, scale.real))
        #        self.__cloc[orb] = self.__cloc[orb] / ( scale + 0.01 )
        Symm.symm_keep_real(self.__cloc)
        if sy_pm: Symm.symm_over_band(self.__cloc)

class LatVar(object):
    """ class LatVar is used to store and manipulate the lattice K-dependent objects.
    """
    def __init__(self, p):
        """ class constructor for the LocVar
            p     : param dictionary
            nband : number of bands
            norbs : number of orbitals
            nffrq : number of fermionic frequencies
            nbfrq : number of bosonic frequencies
            nkpts : number of k-points in x axis
        """
        self.__nband = p["nband"][0]
        self.__norbs = p["norbs"][0]
        self.__nffrq = p["nffrq"][0]
        self.__nbfrq = p["nbfrq"][0]
        self.__nkpts = p["nkpts"][0]

        # allocate memory for G_{lat}(i\omega), \Sigma_{lat}(i\omega)
        self.__glat = numpy.zeros((self.__norbs, self.__nffrq, self.__nkpts), dtype = numpy.complex)
        self.__slat = numpy.zeros((self.__norbs, self.__nffrq, self.__nkpts), dtype = numpy.complex)

        # allocate memory for W_{lat}(i\omega), \Pi_{lat}(i\omega)
        self.__wlat = numpy.zeros((self.__nband, self.__nbfrq, self.__nkpts), dtype = numpy.complex)
        self.__plat = numpy.zeros((self.__nband, self.__nbfrq, self.__nkpts), dtype = numpy.complex)

    @property
    def glat(self):
        """ getter for __glat property
        """
        return self.__glat

    @property
    def slat(self):
        """ getter for __slat property
        """
        return self.__slat

    @property
    def wlat(self):
        """ getter for __wlat property
        """
        return self.__wlat

    @property
    def plat(self):
        """ getter for __plat property
        """
        return self.__plat

    @glat.setter
    def glat(self, glat):
        """ setter for __glat property
        """
        self.__glat = glat

    @slat.setter
    def slat(self, slat):
        """ setter for __slat property
        """
        self.__slat = slat

    @wlat.setter
    def wlat(self, wlat):
        """ setter for __wlat property
        """
        self.__wlat = wlat

    @plat.setter
    def plat(self, plat):
        """ setter for __plat property
        """
        self.__plat = plat

    def sloc_to_slat(self, sloc):
        """ update slat with sloc, ignore the k-dependence
            sloc : k-independent self-energy function
        """
        for ikpt in range(self.__nkpts):
            self.__slat[:,:,ikpt] = sloc.copy()

    def ploc_to_plat(self, ploc):
        """ update plat with ploc, ignore the k-dependence
            ploc : k-independent polarization function
        """
        for ikpt in range(self.__nkpts):
            self.__plat[:,:,ikpt] = ploc.copy()

    def cal_glat_by_dyson(self, mune, fmesh, hk, sedc):
        """ calculate G_{lat}(i\omega) by the dyson equation
            mune  : chemical potential
            fmesh : matsubara frequency
            hk    : lda hamiltonian
            sedc  : double counting term
        """
        # reset self.__glat
        self.__glat.fill(0j)

        # unit matrix
        I = numpy.eye(self.__nband)

        # dummy matrix, used to store self-energy
        saux = numpy.zeros(self.__nband, dtype = numpy.complex)

        # loop over frequencies and k-points, it is parallelized
        for f in range(pyalps.mpi.rank, self.__nffrq, pyalps.mpi.size):
            fM = (1j * fmesh[f] + mune) * I
            for k in range(self.__nkpts):
                # calculate \Sigma - \Sigma_{dc}, only consider the spin up part
                for s in range(self.__nband):
                    saux[s] = self.__slat[2*s,f,k] - sedc[2*s]
                # calculate G^-1 = i\omega + \mu - H - (\Sigma - \Sigma_{dc})
                A =  fM - hk[:,:,k] - numpy.diag(saux)
                # note: you can calculate inverse of A directly, B = numpy.linalg.inv(A)
                # but what I do is more efficient.
                B = numpy.linalg.solve(A, I)
                # setup the final results
                for s in range(self.__nband):
                    self.__glat[2*s,f,k]   = B[s,s] # spin up part
                    self.__glat[2*s+1,f,k] = B[s,s] # spin down part

        # reduce the values from all children nodes
        self.__glat = pyalps.mpi.all_reduce(pyalps.mpi.world, self.__glat, lambda x,y : x + y)
        pyalps.mpi.world.barrier()

        # release memory
        del saux

    def cal_wlat_by_dyson(self, vk):
        """ calculate W_{lat}(i\omega) by the dyson equation
            vk   : k-dependent interaction
        """
        # reset self.__wlat
        self.__wlat.fill(0j)

        # loop over frequencies and orbitals, it is parallelized
        for s in range(self.__nband):
            for b in range(pyalps.mpi.rank, self.__nbfrq, pyalps.mpi.size):
                self.__wlat[s,b,:] = 1.0 / ( 1.0 / vk[s,:] - self.__plat[s,b,:] )
        
        # reduce the values from all children nodes
        self.__wlat = pyalps.mpi.all_reduce(pyalps.mpi.world, self.__wlat, lambda x,y : x + y)
        pyalps.mpi.world.barrier()

class EDMFT(object):
    """ class EDMFT is used to defined the self-consistent equation.
    """
    def __init__(self, p):
        """ class constructor for the EDMFT
            p     : param dictionary
            beta  : inversion temperature
            nband : number of bands
            norbs : number of orbitals
            ntime : number of imaginary time points
            nffrq : number of fermionic frequencies
            nbfrq : number of bosonic frequencies
        """
        self.__beta  = p["beta" ][0]

        self.__nband = p["nband"][0]
        self.__norbs = p["norbs"][0]
        self.__ntime = p["ntime"][0]
        self.__nffrq = p["nffrq"][0]
        self.__nbfrq = p["nbfrq"][0]

        # allocate memory for \Delta(i\omega) and \Delta(\tau)
        self.__hybf = numpy.zeros((self.__norbs, self.__nffrq), dtype = numpy.complex)
        self.__htau = numpy.zeros((self.__norbs, self.__ntime), dtype = numpy.float)

        # allocate memory for K(\tau), K'(\tau)
        self.__ktau = numpy.zeros((self.__nband, self.__ntime), dtype = numpy.float)
        self.__ptau = numpy.zeros((self.__nband, self.__ntime), dtype = numpy.float)

    @property
    def hybf(self):
        """ getter for __hybf property
        """
        return self.__hybf

    @property
    def htau(self):
        """ getter for __htau property
        """
        return self.__htau

    @property
    def ktau(self):
        """ getter for __ktau property
        """
        return self.__ktau

    @property
    def ptau(self):
        """ getter for __ptau property
        """
        return self.__ptau

    @hybf.setter
    def hybf(self, hybf):
        """ setter for __hybf property
        """
        self.__hybf = hybf

    @htau.setter
    def htau(self, htau):
        """ setter for __htau property
        """
        self.__htau = htau

    @ktau.setter
    def ktau(self, ktau):
        """ setter for __ktau property
        """
        self.__ktau = ktau

    @ptau.setter
    def ptau(self, ptau):
        """ setter for __ptau property
        """
        self.__ptau = ptau

    def cal_hybf_by_edmft(self, mune, fmesh, eimp, sedc, sloc, gloc):
        """ calculate \Delta(i\omega) by using the self-consistent equation
            mune  : chemical potential
            fmesh : fermionic mesh
            eimp  : impurity level
            sedc  : double counting term
            sloc  : \Sigma(i\omega)
            gloc  : G(i\omega)
        """
        for s in range(self.__norbs):
            self.__hybf[s,:] = 1j*fmesh + mune - eimp[s] - sloc[s,:] + sedc[s] - 1.0 / gloc[s,:]
            # suppress the real part of hybridization function
            #for f in range(self.__nffrq):
            #    self.__hybf[s,f] = self.__hybf[s,f] - self.__hybf[s, self.__nffrq - 1].real

    def cal_htau_by_fourier(self, p, tmesh, fmesh):
        """ fourier \Delta(i\omega) to \Delta(\tau)
            p : param dictionary
            tmesh : imaginary time mesh
            fmesh : fermionic mesh
        """
        fourier = Fourier(p, 1) # just for fermionic system
        for s in range(self.__norbs):
            fourier.freq_to_time_F(tmesh, self.__htau[s,:], fmesh, self.__hybf[s,:])
        del fourier

        # check its casuality
        for s in range(self.__norbs):
            for t in range(self.__ntime):
                self.__htau[s,t] = min(self.__htau[s,t], -0.00001)

        # symmetrize over imaginary time
        if p["sy_ph"][0]:
            for s in range(self.__norbs):
                Symm.symm_over_time(self.__htau[s])

        # symmetrize over (spin-)orbital index
        if p["sy_pm"][0]: Symm.symm_over_band(self.__htau)

    def cal_ktau_by_fourier(self, U, uloc, tmesh, bmesh):
        """ calculate K(i\omega), and then fourier K(i\omega) to K(\tau)
            U     : static Coulomb interaction
            uloc  : \mathcal{U}(i\omega)
            tmesh : imaginary time mesh
            bmesh : bosonic mesh
        """
        D = uloc - U
        for s in range(self.__nband):
            for t in range(self.__ntime):
                rtmp1 = 0.0
                rtmp2 = 0.0
                for b in range(1,self.__nbfrq):
                    texp = numpy.exp(1j * bmesh[b] * tmesh[t])
                    rtmp1 = rtmp1 + ( D[s,b] * ( 1.0 - texp ) / bmesh[b]**2 ).real
                    rtmp2 = rtmp2 + ( -1j * texp * D[s,b] / bmesh[b] ).real
                self.__ktau[s,t] = 2.0 * rtmp1 / self.__beta
                self.__ptau[s,t] = 2.0 * rtmp2 / self.__beta
                # add the D(0) tail correction term for K(\tau)
                tail = 0.5 * D[s,0].real / self.__beta * tmesh[t] * ( self.__beta - tmesh[t] )
                self.__ktau[s,t] = self.__ktau[s,t] - tail
                # add the D(0) tail correction term for K'(\tau)
                tail = D[s,0].real / self.__beta * ( tmesh[t] - self.__beta/2 )
                self.__ptau[s,t] = self.__ptau[s,t] + tail
        del D

        # enforce positive
        for s in range(self.__nband):
            for t in range(1,self.__ntime-1):
                if self.__ktau[s,t] < 1e-6 :
                    self.__ktau[s,t] = 1e-6
                    self.__ptau[s,t] = 0.0
            self.__ktau[s,0] = 0.0
            self.__ktau[s,self.__ntime - 1] = 0.0

        # symmetrize it?
        Symm.symm_over_band(self.__ktau)
        Symm.symm_over_band(self.__ptau)

    def reload_ktau(self):
        """ read in K(\tau) and K'(\tau) data from external file
        """
        if pyalps.mpi.rank == 0:
            f = open("in.K", "r")
            for t in range(self.__ntime):
                spl = f.readline().split()
                for s in range(self.__nband):
                    self.__ktau[s,t] = float( spl[1] )
                    self.__ptau[s,t] = float( spl[2] )
            f.close()

        # broadcast the ktau and ptau to all children processes
        self.__ktau = pyalps.mpi.broadcast(pyalps.mpi.world, self.__ktau, 0)
        self.__ptau = pyalps.mpi.broadcast(pyalps.mpi.world, self.__ptau, 0)
        pyalps.mpi.world.barrier()

        # enforce positive
        for s in range(self.__nband):
            for t in range(1,self.__ntime-1):
                if self.__ktau[s,t] < 1e-6 :
                    self.__ktau[s,t] = 1e-6
                    self.__ptau[s,t] = 0.0
            self.__ktau[s,0] = 0.0
            self.__ktau[s,self.__ntime - 1] = 0.0

        # symmetrize it?
        Symm.symm_over_band(self.__ktau)
        Symm.symm_over_band(self.__ptau)

# let's go
if __name__ == '__main__':

    #---------------------------------------------------------------------

    # print the header for the lda + edmft code
    lda_edmft_header()

    # setup config parameters
    p = lda_edmft_config()

    # extract some key parameters which are useful for the impurity solver
    s_p = lda_edmft_sconfig(p)

    #---------------------------------------------------------------------

    if pyalps.mpi.rank == 0:
        print (">>> Create class instance (lat, loc, dmft, solver) ...")
        t_start = time.clock()

    # create lattice instance
    lat = LatVar(p)

    # create local instance
    loc = LocVar(p)

    # create edmft instance
    dmft = EDMFT(p)

    # create quantum impurity solver instance
    solver = AlpsSolver(s_p)

    if pyalps.mpi.rank == 0:
        t_end = time.clock()
        print ("    status: OK    time: " + str(t_end - t_start) + "s")
        print ("")

    #---------------------------------------------------------------------

    if pyalps.mpi.rank == 0:
        print (">>> Build mesh (tmesh, fmesh, bmesh) ...")
        t_start = time.clock()

    # create imaginary time mesh
    tmesh = Mesh.create_time_mesh(p["beta"][0], p["ntime"][0])

    # create fermionic mesh
    fmesh = Mesh.create_fermion_mesh(p["beta"][0], p["nffrq"][0])

    # create bosonic mesh
    bmesh = Mesh.create_boson_mesh(p["beta"][0], p["nbfrq"][0])

    if pyalps.mpi.rank == 0:
        t_end = time.clock()
        print ("    status: OK    time: " + str(t_end - t_start) + "s")
        print ("")

    #---------------------------------------------------------------------

    if pyalps.mpi.rank == 0:
        print (">>> Build hamiltonian (hk, vk, eimp, umat) ...")
        t_start = time.clock()

    # create hamiltonian for the model
    eimp, hk = Lattice.create_hk(p["nband"][0], p["norbs"][0], p["nkpts"][0], p["Z"][0])

    # create general interaction for the model
    if p["dmft"][0] == 3 : vk = Lattice.create_vk(p["nband"][0], p["nkpts"][0])

    # create onsite Coulomb interaction for the model
    umat = Lattice.create_umat(p["norbs"][0])

    if pyalps.mpi.rank == 0:
        t_end = time.clock()
        print ("    status: OK    time: " + str(t_end - t_start) + "s")
        print ("")

    #---------------------------------------------------------------------

    if pyalps.mpi.rank == 0:
        print (">>> Reload iteration variables (sloc, ploc, slat, plat) ...")
        t_start = time.clock()

    # if we do not need to start from scratch, just try to initialize loc.sloc
    # and loc.ploc from existing data
    # use loc.sloc to initialize lat.slat, and loc.ploc to initialize lat.plat
    if p["start"][0] == False :
        loc.reload_sloc()
        if p["dmft"][0] == 3 : loc.reload_ploc()

    if p["start"][0] == False :
        lat.sloc_to_slat(loc.sloc)
        if p["dmft"][0] == 3 : lat.ploc_to_plat(loc.ploc)

    if pyalps.mpi.rank == 0:
        t_end = time.clock()
        print ("    status: OK    time: " + str(t_end - t_start) + "s")
        print ("")

    #---------------------------------------------------------------------

    g_conv = False     # convergence flag for gloc, impurity green's function
    s_conv = False     # convergence flag for sloc, self-energy function
    p_conv = False     # convergence flag for ploc, bosonic self-energy function
    conv_value = 0.002 # convergence condition
    for _iter_ in range (1, p["niter"][0] + 1):

        # print the iteration information
        lda_edmft_iter(p["dmft"][0], p["niter"][0], _iter_)

    #---------------------------------------------------------------------

        if pyalps.mpi.rank == 0:
            print (">>> Adjust the double counting term (sedc) ...")
            t_start = time.clock()

        # extract current occupation number
        if _iter_ == 1: # for the first iteration
            curr_occu = p["occup"][0]
        else: # when _iter_ > 1
            curr_occu = solver.extract(100, p)
            curr_occu = pyalps.mpi.broadcast(pyalps.mpi.world, curr_occu, 0)
            pyalps.mpi.world.barrier()

        # calculate the double counting term
        for case in switch( p["sedc"][0] ):
            if case(5): # fixed occupation for KAV scheme
                loc.cal_sedc_by_kav(p["U"][0], p["J"][0], p["occup"][0])
                break

            if case(4): # fixed occupation for FLL scheme
                loc.cal_sedc_by_fll(p["U"][0], p["J"][0], p["occup"][0])
                break

            if case(3): # standard KAV scheme
                loc.cal_sedc_by_kav(p["U"][0], p["J"][0], curr_occu)
                break

            if case(2): # standard FLL scheme
                loc.cal_sedc_by_fll(p["U"][0], p["J"][0], curr_occu)
                break

            if case(1): # zero double counting term
                break

        # special consideration for the first iteration
        if _iter_ == 1 and p["start"][0]:
            loc.sedc = loc.sedc * 0.0

        # print the double counting term
        if pyalps.mpi.rank == 0:
            print ("    curr_sedc: %s  mode_sedc: %i" % (loc.sedc[0].real, p["sedc"][0]) )

        if pyalps.mpi.rank == 0:
            t_end = time.clock()
            print ("    status: OK    time: " + str(t_end - t_start) + "s")
            print ("")

    #---------------------------------------------------------------------

        if pyalps.mpi.rank == 0:
            print (">>> Adjust the chemical potential (mune) ...")
            t_start = time.clock()

        # determine the chemical potential using the different algorithms
        for case in switch ( p["semu"][0] ):
            if case(2): # adjust chemical potential for the lattice problem
                p["mune"][0] = lda_edmft_slevel(p["norbs"][0], p["nband"][0], p["nkpts"][0], p["occup"][0], p["beta"][0], fmesh, hk, lat.slat, loc.sedc)
                break

            if case(1): # adjust chemical potential for the impurity model
                p["mune"][0] = p["mune"][0] + ( p["occup"][0] - curr_occu )
                break

        # print the chemical potential
        if pyalps.mpi.rank == 0:
            print ("    curr_occu: %s  curr_mune: %s" % (curr_occu, p["mune"][0]) )

        if pyalps.mpi.rank == 0:
            t_end = time.clock()
            print ("    status: OK    time: " + str(t_end - t_start) + "s")
            print ("")

    #---------------------------------------------------------------------

        if pyalps.mpi.rank == 0:
            print (">>> Build lattice variables (glat, wlat) ...")
            t_start = time.clock()

        # build glat from scratch
        lat.cal_glat_by_dyson(p["mune"][0], fmesh, hk, loc.sedc)

        # build wlat from scratch
        if p["dmft"][0] == 3 : lat.cal_wlat_by_dyson(vk)

        if pyalps.mpi.rank == 0:
            t_end = time.clock()
            print ("    status: OK    time: " + str(t_end - t_start) + "s")
            print ("")

    #---------------------------------------------------------------------

        if pyalps.mpi.rank == 0:
            print (">>> Build local variables (gloc, wloc, sloc, ploc, bloc, uloc) ...")
            t_start = time.clock()

        # build gloc from glat
        loc.cal_gloc_by_ksum(lat.glat)

        # build sloc from slat
        loc.cal_sloc_by_ksum(lat.slat)

        # build bloc from gloc and sloc
        loc.cal_bloc_by_dyson()

        # build wloc from wlat
        if p["dmft"][0] == 3 : loc.cal_wloc_by_ksum(lat.wlat)

        # build ploc from plat
        if p["dmft"][0] == 3 : loc.cal_ploc_by_ksum(lat.plat)

        # build uloc from wloc and ploc
        if p["dmft"][0] == 3 : loc.cal_uloc_by_dyson()

        # dump the data, only the master node can do it
        if pyalps.mpi.rank == 0:
            Dump.dump_freq_data("c_gloc.dat." + str(_iter_), fmesh, loc.gloc)
            Dump.dump_freq_data("c_sloc.dat." + str(_iter_), fmesh, loc.sloc)
            Dump.dump_freq_data("c_bloc.dat." + str(_iter_), fmesh, loc.bloc)

            if p["dmft"][0] == 3 : Dump.dump_freq_data("c_wloc.dat." + str(_iter_), bmesh, loc.wloc)
            if p["dmft"][0] == 3 : Dump.dump_freq_data("c_ploc.dat." + str(_iter_), bmesh, loc.ploc)
            if p["dmft"][0] == 3 : Dump.dump_freq_data("c_uloc.dat." + str(_iter_), bmesh, loc.uloc)

        if pyalps.mpi.rank == 0:
            t_end = time.clock()
            print ("    status: OK    time: " + str(t_end - t_start) + "s")
            print ("")

    #---------------------------------------------------------------------

        if pyalps.mpi.rank == 0:
            print (">>> Build dmft variables (hybf, htau, ktau, ptau) ...")
            t_start = time.clock()

        # build ktau and ptau from fourier transformation
        for case in switch ( p["dmft"][0] ):
            if case(1):
                break;

            if case(2):
                dmft.reload_ktau()
                break

            if case(3):
                dmft.cal_ktau_by_fourier(p["U"][0], loc.uloc, tmesh, bmesh)
                break

        # build hybf from dmft self-consistent equation
        dmft.cal_hybf_by_edmft(p["mune"][0], fmesh, eimp, loc.sedc, loc.sloc, loc.gloc)

        # build htau from fourier transformation
        dmft.cal_htau_by_fourier(p, tmesh, fmesh)

        # dump the htau, ktau and ptau data to disk file, used by ctqmc impurity solver
        if pyalps.mpi.rank == 0:
            Dump.dump_htau_data("D.dat", tmesh, dmft.htau)
            if p["dmft"][0] > 1 : Dump.dump_ktau_data("K.dat", tmesh, dmft.ktau, dmft.ptau)

        # dump the data, only master node can do it
        if pyalps.mpi.rank == 0:
            Dump.dump_freq_data("c_hybf.dat." + str(_iter_), fmesh, dmft.hybf)
            Dump.dump_time_data("c_htau.dat." + str(_iter_), tmesh, dmft.htau)
            if p["dmft"][0] > 1 : Dump.dump_time_data("c_ktau.dat." + str(_iter_), tmesh, dmft.ktau)
            if p["dmft"][0] > 1 : Dump.dump_time_data("c_ptau.dat." + str(_iter_), tmesh, dmft.ptau)

        if pyalps.mpi.rank == 0:
            t_end = time.clock()
            print ("    status: OK    time: " + str(t_end - t_start) + "s")
            print ("")

    #---------------------------------------------------------------------

        if pyalps.mpi.rank == 0:
            print (">>> Launch quantum impurity solver (alps::cthyb) ...")
            t_start = time.clock()

        # check the status of quantum impurity solver
        solver.setup(mode = p["dmft"][0])

        # setup the orbital chemical potential and shift for Coulomb interaction
        s_p["MU"] = p["mune"][0]
        if p["dmft"][0] == 3:
            shift = 2.0 * dmft.ptau[0,0]
        else:
            shift = 0.0

        # dump the orbital chemical potential to the disk file
        pyalps.mpi.world.barrier()
        if pyalps.mpi.rank == 0:
            Dump.dump_umat_data(s_p["U_MATRIX"], umat, shift)
            Dump.dump_mune_data(s_p["MU_VECTOR"], s_p["MU"], eimp, loc.sedc)
        pyalps.mpi.world.barrier()

        # boot the quantum impurity solver
        pyalps.mpi.world.barrier()
        solver.start(_iter_, s_p)
        pyalps.mpi.world.barrier()

        # finalize the quantum impurity solver
        solver.final(_iter_)

        if pyalps.mpi.rank == 0:
            t_end = time.clock()
            print ("    status: OK    time: " + str(t_end - t_start) + "s")
            print ("")

    #---------------------------------------------------------------------

        if pyalps.mpi.rank == 0:
            print (">>> Rebuild local variables (cloc, wloc, ploc, gloc, sloc) ...")
            t_start = time.clock()

        # extract the measurement results: local green's function
        # only master node can do it
        loc.gloc = solver.extract(4, p)
        loc.gloc = pyalps.mpi.broadcast(pyalps.mpi.world, loc.gloc, 0)
        pyalps.mpi.world.barrier()
        loc.sym_gloc(p["sy_pm"][0], p["sy_ph"][0])

        # extract the measurement results: self-energy function
        # only master node can do it
        loc.sloc = solver.extract(6, p)
        loc.sloc = pyalps.mpi.broadcast(pyalps.mpi.world, loc.sloc, 0)
        loc.fix_sloc_by_dyson(fmesh, p["mune"][0], eimp, dmft.hybf)
        pyalps.mpi.world.barrier()
        loc.sym_sloc(p["sy_pm"][0])

        # extract the measurement results: charge susceptibility
        # only master node can do it
        loc.cloc = solver.extract(98, p)
        loc.cloc = pyalps.mpi.broadcast(pyalps.mpi.world, loc.cloc, 0)
        pyalps.mpi.world.barrier()
        loc.sym_cloc(p["sy_pm"][0])

        if p["dmft"][0] == 3 :
            # update wloc by using dyson equation (need uloc and cloc)
            loc.cal_wloc_by_dyson()
            loc.sym_wloc(p["sy_pm"][0])

            # update ploc by using dyson equation (need uloc and cloc)
            #loc.cal_ploc_by_dyson()
            loc.cal_ploc_by_direct()
            loc.sym_ploc(p["sy_pm"][0])

        # dump the data, only master node can do it
        if pyalps.mpi.rank == 0:
            Dump.dump_freq_data("s_gloc.dat." + str(_iter_), fmesh, loc.gloc)
            Dump.dump_freq_data("s_sloc.dat." + str(_iter_), fmesh, loc.sloc)
            Dump.dump_freq_data("s_cloc.dat." + str(_iter_), bmesh, loc.cloc)
            if p["dmft"][0] == 3 : Dump.dump_freq_data("s_wloc.dat." + str(_iter_), bmesh, loc.wloc)
            if p["dmft"][0] == 3 : Dump.dump_freq_data("s_ploc.dat." + str(_iter_), bmesh, loc.ploc)

        if pyalps.mpi.rank == 0:
            t_end = time.clock()
            print ("    status: OK    time: " + str(t_end - t_start) + "s")
            print ("")

    #---------------------------------------------------------------------

        if pyalps.mpi.rank == 0:
            print (">>> Check and Mix local variables (sloc, ploc) ...")
            t_start = time.clock()

        # mix gloc and check its convergence
        g_err = loc.mix_gloc(_iter_, p["alpha"][0])
        g_conv = False
        if g_err < conv_value : g_conv = True

        # mix sloc and check its convergence
        s_err = loc.mix_sloc(_iter_, p["alpha"][0])
        s_conv = False
        if s_err < conv_value : s_conv = True

        # mix ploc and check its convergence
        if p["dmft"][0] == 3 :
            p_err = loc.mix_ploc(_iter_, p["alpha"][0])
        else :
            p_err = 100.0
        p_conv = False
        if p_err < conv_value : p_conv = True

        # dump the convergence results
        if pyalps.mpi.rank == 0:
            print ("    curr_iter: %i   g_err: %s   s_err: %s   p_err: %s" % (_iter_, g_err, s_err, p_err) )

        # dump the data, only master node can do it
        if pyalps.mpi.rank == 0:
            Dump.dump_freq_data("m_gloc.dat." + str(_iter_), fmesh, loc.gloc)
            Dump.dump_freq_data("m_sloc.dat." + str(_iter_), fmesh, loc.sloc)
            if p["dmft"][0] == 3 : Dump.dump_freq_data("m_ploc.dat." + str(_iter_), bmesh, loc.ploc)

        if pyalps.mpi.rank == 0:
            t_end = time.clock()
            print ("    status: OK    time: " + str(t_end - t_start) + "s")
            print ("")

    #---------------------------------------------------------------------

        if pyalps.mpi.rank == 0:
            print (">>> Apply self-consistent equation (slat, plat) ...")
            t_start = time.clock()

        # update \Sigma(k, i\omega) with \Sigma_{loc}(i\omega)
        lat.sloc_to_slat(loc.sloc)

        # update \Pi(k, i\omega) with \Pi_{loc}(i\omega)
        if p["dmft"][0] == 3 : lat.ploc_to_plat(loc.ploc)

        if pyalps.mpi.rank == 0:
            t_end = time.clock()
            print ("    status: OK    time: " + str(t_end - t_start) + "s")
            print ("")

    #---------------------------------------------------------------------

        # check whether the convergence is obtained
        # note : 10 is minimum iteration loop
        if ( g_conv or s_conv or p_conv ) and _iter_ > 25 :
            if pyalps.mpi.rank == 0:
                print(">>> Congratulations! The convergence is obtained.")
                print ("")
            break

    #---------------------------------------------------------------------

    lda_edmft_footer(p["dmft"][0])
