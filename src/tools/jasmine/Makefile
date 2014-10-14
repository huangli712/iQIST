##
##
## Introduction
## ============
##
## Makefile for iqist/src/tools/jasmine directory
##
## Usage
## =====
##
## 'make', 'make lib', 'make clean', 'make clean-dat' or 'make clean-all'
##
##

.SUFFIXES: .f90

include ../../build/make.sys

default: all

modn = atomic_control.o atomic_context.o
lev1 = atomic_driver.o
lev2 = atomic_basis.o atomic_stream.o
lev3 = atomic_dump.o atomic_print.o
main = atomic_main.o
MAPI = libMM.a

objects = $(modn) $(lev1) $(lev2) $(lev3) $(main) $(MAPI)

default: all

all: exe

cssl:
	cp ../../common/libMM.a .

csml:
	cp ../../common/constants.mod .
	cp ../../common/parser.mod .
	cp ../../common/linkedlist.mod .

exe: cssl csml atomic_exe
atomic_exe: $(objects)
	$(LINKER) $(objects) -o atomic $(LFLAGS) $(LIBS)

.f90.o:
	$(F90) $(FFLAGS) $*.f90

clean:
	rm -f *.mod
	rm -f *.o
	rm -f atomic
	rm -f libMM.a

clean-dat:
	rm -f *.dat
	rm -f *.bin.*
	rm -f *.out
