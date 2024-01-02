##
##
## Introduction
## ============
##
## Makefile for iqist/src/atomic directory.
##
## Usage
## =====
##
## 'make', 'make clean', 'make clean-dat' or 'make clean-all'
##
## Author
## ======
##
## This building system is designed, created, and maintained by
##
## Li Huang // email: huangli@caep.cn
##
## History
## =======
##
## 02/28/2015 by li huang (created)
## 01/02/2024 by li huang (last modified)
##
##

.SUFFIXES: .f90

include ../../build/make.inc

modc = atomic_control.o atomic_context.o
core = atomic_driver.o
lev1 = atomic_full.o atomic_sector.o
lev2 = atomic_stream.o atomic_util.o
lev3 = atomic_dump.o atomic_print.o
main = atomic_main.o
open = atomic_open.o

CAPI = i_ctqmc.o
CLIB = libflink.a

objects = $(modc) $(core) $(lev1) $(lev2) $(lev3)
exe_obj = $(objects) $(main) $(CLIB)
lib_obj = $(objects) $(open) $(CAPI)

default: all

all: exe

api:
	#cp ../../capi/capi.mod .
	#cp ../../capi/i_ctqmc.o .

cssl:
	#cp ../../base/libMM.a .
	cp $(FLINK)/libflink.a .

csml:
	cp $(FLINK)/constants.mod .
	cp $(FLINK)/parser.mod .
	cp $(FLINK)/linkedlist.mod .
	#cp $(FLINK)/version.mod .

exe: api cssl csml atomic-exe
atomic-exe: $(exe_obj)
	$(LINKER) $(exe_obj) -o atomic $(LFLAGS) $(LIBS)

lib: api cssl csml atomic-lib
atomic-lib: $(lib_obj)
	$(ARCHIVER) libatomic.a $(lib_obj)

pylib: lib
	cpp -E -DPYAPI atomic_open.f90 > atomic_open_py.f90
	$(F2PY) -c --f90exec=$(F2PYC) -m pyjasmine -L. -latomic -lMM $(F2PYL) atomic_open_py.f90
	rm -fr *.dSYM atomic_open_py.f90

.f90.o:
	$(F90) $(FFLAGS) $*.f90

clean:
	rm -f *.mod
	rm -f *.o
	rm -f *.so
	rm -f atomic
	rm -f libatomic.a
	rm -f libMM.a

clean-dat:
	rm -f *.cix
	rm -f *.dat
	rm -f *.out
	rm -f solver.umat.in

clean-all: clean clean-dat
