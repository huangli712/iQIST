##
##
## Introduction
## ============
##
## Makefile for iqist/src/ctqmc/api directory.
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
## 09/21/2014 by li huang
##
##

.SUFFIXES: .f90

include ../../build/make.sys

mod = ctqmc_api.o

default: all

all: build_mod

build_mod: $(mod)

f2py:
	fpp -P -DF2PY ctqmc_api.f90 > py_api.f90

build-azalea-py: f2py
	cp ../azalea/libctqmc.a .
	cp ../../common/libMM.a .
	f2py -c --fcompiler=intelem py_api.f90 -m pyiqist -L. -lctqmc -lMM -L/usr/local/openmpi/lib -lmpi_mpifh $(LIBS) 

.f90.o:
	$(F90) $(FFLAGS) $*.f90

clean:
	rm -f *.mod
	rm -f *.o
	rm -f *.so
	rm -f libctqmc.a
	rm -f libMM.a
	rm -f py_api.f90

clean-dat:
	rm -f *.dat
	rm -f *.bin.*
	rm -f *.out

clean-all: clean clean-dat
