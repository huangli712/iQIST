##
##
## Introduction
## ============
##
## Makefile for iqist/src/app/azalea directory.
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
## 01/17/2017 by li huang (created)
## 02/01/2017 by li huang (last modified)
##
##

.SUFFIXES: .f90

include ../../../build/make.sys

modc = df_control.o df_context.o
dual = df_diag.o
core = df_core.o
lev1 = df_solver.o df_susc.o
lev2 = df_stream.o df_util.o
lev3 = df_dump.o df_print.o
main = df_main.o

CLIB = libMM.a

objects = $(modc) $(dual) $(core) $(lev1) $(lev2) $(lev3)
exe_obj = $(objects) $(main) $(CLIB)
lib_obj = $(objects)

default: all

all: exe

cssl:
	cp ../../base/libMM.a .

csml:
	cp ../../base/constants.mod .
	cp ../../base/mmpi.mod .
	cp ../../base/parser.mod .
	cp ../../base/linkedlist.mod .
	cp ../../base/version.mod .

exe: cssl csml df-exe
df-exe: $(exe_obj)
	$(LINKER) $(exe_obj) -o df $(LFLAGS) $(LIBS)

lib: cssl csml df-lib
df-lib: $(lib_obj)
	$(ARCHIVER) libdf.a $(lib_obj)

.f90.o:
	$(F90) $(FFLAGS) $*.f90

clean:
	rm -f *.mod
	rm -f *.o
	rm -f *.so
	rm -f df
	rm -f libdf.a
	rm -f libMM.a

clean-dat:
	rm -f *.dat
	rm -f *.out

clean-all: clean clean-dat
