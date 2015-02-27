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
## 02/28/2015 by li huang
##
##

.SUFFIXES: .f90

include ../build/make.sys

mod = ctqmc_api.o

default: all

all: build_mod

build_mod: $(mod)

f2py:
	fpp -P -DMPY ctqmc_api.f90 > ctqmc_api_py.f90

libctqmc:
	cp ../azalea/libctqmc.a .

libMM:
	cp ../../common/libMM.a .

pyiqist: f2py libctqmc libMM
	f2py -c --fcompiler=intelem ctqmc_api_py.f90 -m pyiqist -L. -lctqmc -lMM $(MPIL) $(LIBS)

.f90.o:
	$(F90) $(FFLAGS) $*.f90

clean:
	rm -f *.mod
	rm -f *.o
	rm -f *.so
	rm -f libctqmc.a
	rm -f libMM.a
	rm -f ctqmc_api_py.f90

clean-dat:
	rm -f *.dat
	rm -f *.bin.*
	rm -f *.out

clean-all: clean clean-dat
