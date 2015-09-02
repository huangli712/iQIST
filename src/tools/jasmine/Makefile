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

include ../../build/make.sys

modc = atomic_control.o atomic_context.o
core = atomic_driver.o
lev1 = atomic_full.o atomic_sector.o
lev2 = atomic_stream.o atomic_util.o
lev3 = atomic_dump.o atomic_print.o
main = atomic_main.o
open = atomic_open.o

CAPI = atomic_api.o
CLIB = libMM.a

objects = $(modc) $(core) $(lev1) $(lev2) $(lev3)
exe_obj = $(objects) $(main) $(CLIB)
lib_obj = $(objects) $(open) $(CAPI)

default: all

all: exe

api:
	cp ../../capi/japi.mod .
	cp ../../capi/atomic_api.o .

cssl:
	cp ../../common/libMM.a .

csml:
	cp ../../common/constants.mod .
	cp ../../common/parser.mod .
	cp ../../common/linkedlist.mod .

exe: api cssl csml atomic-exe
atomic-exe: $(exe_obj)
	$(LINKER) $(exe_obj) -o atomic $(LFLAGS) $(LIBS)

lib: api cssl csml atomic-lib
atomic-lib: $(lib_obj)
	$(ARCHIVER) libatomic.a $(lib_obj)

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
