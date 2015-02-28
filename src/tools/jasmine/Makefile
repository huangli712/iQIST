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
## 02/28/2015 by li huang
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

api:
	cp ../../api/japi.mod .
	cp ../../api/atomic_api.o .

cssl:
	cp ../../common/libMM.a .

csml:
	cp ../../common/constants.mod .
	cp ../../common/parser.mod .
	cp ../../common/linkedlist.mod .

touch:
	touch atomic_main.f90

exe: api cssl csml atomic-exe
atomic-exe: $(objects)
	$(LINKER) $(objects) -o atomic $(LFLAGS) $(LIBS)

lib: api cssl csml touch atomic-lib
atomic-lib: $(objects)
	$(ARCHIVER) libatomic.a $(objects)

.f90.o:
	$(F90) $(FFLAGS) $*.f90

clean:
	rm -f *.mod
	rm -f *.o
	rm -f atomic
	rm -f libatomic.a
	rm -f libMM.a

clean-dat:
	rm -f *.cix
	rm -f *.dat
	rm -f *.bin.*
	rm -f *.out

clean-all: clean clean-dat
