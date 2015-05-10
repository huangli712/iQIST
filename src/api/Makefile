##
##
## Introduction
## ============
##
## Makefile for iqist/src/api directory.
##
## Usage
## =====
##
## 'make' or 'make clean'
##
## Author
## ======
##
## This building system is designed, created, implemented, and maintained by
##
## Li Huang // email: huangli712@gmail.com
##
## History
## =======
##
## 02/28/2015 by li huang
##
##

.SUFFIXES: .f90

include ../build/make.sys

mod = ctqmc_api.o hfqmc_api.o atomic_api.o

default: all

all: build_mod

build_mod: $(mod)

f2py:
	cpp -E -DMPY ctqmc_api.f90 > ctqmc_api_py.f90
	cpp -E -DMPY hfqmc_api.f90 > hfqmc_api_py.f90
	cpp -E -DMPY atomic_api.f90 > atomic_api_py.f90

libctqmc:
	cp ../ctqmc/azalea/libctqmc.a .

libhfqmc:
	cp ../hfqmc/daisy/libhfqmc.a .

libatomic:
	cp ../tools/jasmine/libatomic.a .

libMM:
	cp ../common/libMM.a .

pyiqist: f2py libctqmc libMM
	$(F2PY) -c --f90exec=$(F2PYC) -m pyiqist -L. -lctqmc -lMM --link-lapack_opt ctqmc_api_py.f90
	rm -fr *.dSYM

pydaisy: f2py libhfqmc libMM
	f2py -c --fcompiler=intelem hfqmc_api_py.f90 -m pydaisy -L. -lhfqmc -lMM $(MPIL) $(LIBS)

pyjasmine: f2py libatomic libMM
	f2py -c --fcompiler=intelem atomic_api_py.f90 -m pyjasmine -L. -latomic -lMM $(LIBS)

.f90.o:
	$(F90) $(FFLAGS) $*.f90

clean:
	rm -f *.mod
	rm -f *.o
	rm -f *.so
	rm -f libctqmc.a
	rm -f libhfqmc.a
	rm -f libatomic.a
	rm -f libMM.a
	rm -f ctqmc_api_py.f90
	rm -f hfqmc_api_py.f90
	rm -f atomic_api_py.f90

clean-dat:
	rm -f *.dat
	rm -f *.bin.*
	rm -f *.out

clean-all: clean clean-dat
