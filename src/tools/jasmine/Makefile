##
##
## Introduction
## ============
##
## Makefile for iqist/src/tools/jasmine directory.
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

modc = atomic_control.o atomic_context.o
core = atomic_driver.o
lev1 = atomic_full.o atomic_sector.o
lev2 = atomic_stream.o atomic_util.o
lev3 = atomic_dump.o atomic_print.o
main = atomic_main.o
MAPI = atomic_api.o libMM.a

objects = $(modc) $(core) $(lev1) $(lev2) $(lev3) $(main) $(MAPI)

default: all

all: exe

api: atomic_api.o

cssl:
	cp ../../common/libMM.a .

csml:
	cp ../../common/constants.mod .
	cp ../../common/parser.mod .
	cp ../../common/linkedlist.mod .

exe: api cssl csml atomic-exe
atomic-exe: $(objects)
	$(LINKER) $(objects) -o atomic $(LFLAGS) $(LIBS)

lib: api cssl csml atomic-lib
atomic-lib: $(objects)
	$(ARCHIVER) libatomic.a $(objects)

f2py:
	fpp -P -DF2PY atomic_api.f90 > py_api.f90

pyjasmine: f2py lib
	f2py -c --fcompiler=intelem py_api.f90 -m pyjasmine -L. -latomic -lMM $(LIBS)

.f90.o:
	$(F90) $(FFLAGS) $*.f90

clean:
	rm -f *.mod
	rm -f *.o
	rm -f *.so
	rm -f atomic
	rm -f libatomic.a
	rm -f libMM.a
	rm -f py_api.f90

clean-dat:
	rm -f *.dat
	rm -f *.bin.*
	rm -f *.out

clean-all: clean clean-dat
