import numpy

epsilon = [-1.0, 0.0, 1.0]
V = 0.5
mfreq = 8193
beta = 10.0
norbs = 2

def out_hyb(norbs, mfreq, rmesh, hybf, fileName = None):
    """ try to write the hybridization function to the solver.hyb.in
        file, only suitable for the ctqmc impurity solver
    """
    if fileName is None:
        f = open("solver.hyb.in","w")
    else:
        f = open(fileName,"w")

    for i in range(norbs):
        for j in range(mfreq):
            print >> f, '%6d %16.8f %16.8f %16.8f %16.8f %16.8f' % \
            ( i+1, rmesh[j], hybf[j,i,i].real, hybf[j,i,i].imag, 0.0, 0.0 )
        print >> f
        print >> f

    f.close()

rmesh = numpy.zeros(mfreq, dtype = numpy.float)
hybf = numpy.zeros((mfreq,norbs,norbs), dtype = numpy.complex)

for i in range(mfreq):
    rmesh[i] = (2*i + 1) * numpy.pi / beta

for i in range(mfreq):
    for e in epsilon:
        for b in range(norbs):
            hybf[i,b,b] = hybf[i,b,b] + V*V / ( 1j * rmesh[i] - e )
out_hyb(norbs, mfreq, rmesh, hybf)
