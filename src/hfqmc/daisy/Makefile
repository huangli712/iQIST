##
##
## Introduction
## ============
##
## Makefile for iqist/src/hfqmc/daisy directory.
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
## Li Huang // email: lihuang.dmft@gmail.com
##
## History
## =======
##
## 02/28/2015 by li huang (created)
## 08/17/2015 by li huang (last modified)
##
##

.SUFFIXES: .f90

include ../../../build/make.sys

modc = hfqmc_control.o hfqmc_context.o
dmft = hfqmc_dmft.o
core = hfqmc_solver.o
lev1 = hfqmc_update.o
lev2 = hfqmc_record.o hfqmc_stream.o hfqmc_util.o
lev3 = hfqmc_dump.o hfqmc_print.o
main = hfqmc_main.o
open = hfqmc_open.o

CAPI = hfqmc_api.o
CLIB = libMM.a

objects = $(modc) $(dmft) $(core) $(lev1) $(lev2) $(lev3)
exe_obj = $(objects) $(main) $(CLIB)
lib_obj = $(objects) $(open) $(CAPI)

default: all

all: exe

api:
	cp ../../capi/dapi.mod .
	cp ../../capi/hfqmc_api.o .

cssl:
	cp ../../base/libMM.a .

csml:
	cp ../../base/constants.mod .
	cp ../../base/mmpi.mod .
	cp ../../base/spring.mod .
	cp ../../base/parser.mod .
	cp ../../base/linkedlist.mod .

exe: api cssl csml hfqmc-exe
hfqmc-exe: $(exe_obj)
	$(LINKER) $(exe_obj) -o hfqmc $(LFLAGS) $(LIBS)

lib: api cssl csml hfqmc-lib
hfqmc-lib: $(lib_obj)
	$(ARCHIVER) libhfqmc.a $(lib_obj)

pylib: lib
	cpp -E -DPYAPI hfqmc_open.f90 > hfqmc_open_py.f90
	$(F2PY) -c --f90exec=$(F2PYC) -m pydaisy -L. -lhfqmc -lMM $(F2PYL) hfqmc_open_py.f90
	rm -fr *.dSYM hfqmc_open_py.f90

.f90.o:
	$(F90) $(FFLAGS) $*.f90

clean:
	rm -f *.mod
	rm -f *.o
	rm -f *.so
	rm -f hfqmc
	rm -f libhfqmc.a
	rm -f libMM.a

clean-dat:
	rm -f *.dat
	rm -f *.out

clean-all: clean clean-dat
