##
##
## Introduction
## ============
##
## Makefile for iqist/src/ctqmc/azalea directory.
##
## Usage
## =====
##
## 'make', 'make lib', 'make clean', 'make clean-dat' or 'make clean-all'
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

modc = hfqmc_control.o hfqmc_context.o
dmft = hfqmc_dmft.o
core = hfqmc_solver.o
lev1 = hfqmc_update.o
lev2 = hfqmc_record.o hfqmc_stream.o hfqmc_util.o
lev3 = hfqmc_dump.o hfqmc_print.o
main = hfqmc_main.o
MAPI = hfqmc_api.o libMM.a

objects = $(modc) $(dmft) $(core) $(lev1) $(lev2) $(lev3) $(main) $(MAPI)

default: all

all: exe

api: hfqmc_api.o

cssl:
	cp ../../common/libMM.a .

csml:
	cp ../../common/constants.mod .
	cp ../../common/mmpi.mod .
	cp ../../common/spring.mod .
	cp ../../common/parser.mod .
	cp ../../common/linkedlist.mod .

exe: api cssl csml hfqmc-exe
hfqmc-exe: $(objects)
	$(LINKER) $(objects) -o hfqmc $(LFLAGS) $(LIBS)

lib: api cssl csml hfqmc-lib
hfqmc-lib: $(objects)
	$(ARCHIVER) libhfqmc.a $(objects)

f2py:
	fpp -P -DMPY hfqmc_api.f90 > hfqmc_api_py.f90

pydaisy: f2py lib
	f2py -c --fcompiler=intelem hfqmc_api_py.f90 -m pydaisy -L. -lhfqmc -lMM $(MPIL) $(LIBS)

.f90.o:
	$(F90) $(FFLAGS) $*.f90

clean:
	rm -f *.mod
	rm -f *.o
	rm -f *.so
	rm -f hfqmc
	rm -f libhfqmc.a
	rm -f libMM.a
	rm -f hfqmc_api_py.f90

clean-dat:
	rm -f *.dat
	rm -f *.bin.*
	rm -f *.out

clean-all: clean clean-dat
