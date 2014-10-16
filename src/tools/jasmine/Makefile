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
MAPI = libMM.a

objects = $(modc) $(core) $(lev1) $(lev2) $(lev3) $(main) $(MAPI)

default: all

all: exe

cssl:
	cp ../../common/libMM.a .

csml:
	cp ../../common/constants.mod .
	cp ../../common/parser.mod .
	cp ../../common/linkedlist.mod .

exe: cssl csml atomic-exe
atomic-exe: $(objects)
	$(LINKER) $(objects) -o atomic $(LFLAGS) $(LIBS)

lib: cssl csml atomic-lib
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
	rm -f *.dat
	rm -f *.bin.*
	rm -f *.out

clean-all: clean clean-dat
