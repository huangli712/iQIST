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

.f90.o:
	$(F90) $(FFLAGS) $*.f90

clean:
	rm -f *.mod
	rm -f *.o
	rm -f libMM.a

clean-dat:
	rm -f *.dat
	rm -f *.bin.*
	rm -f *.out

clean-all: clean clean-dat
