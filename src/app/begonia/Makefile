##
##
## Introduction
## ============
##
## Makefile for iqist/src/app/begonia directory.
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
## 01/17/2017 by li huang (last modified)
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

objects = $(modc) $(dmft) $(core) $(lev1) $(lev2) $(lev3)
exe_obj = $(objects) $(main) $(CLIB)
lib_obj = $(objects) $(open) $(CAPI)

default: all

all: exe

api:
	cp ../../capi/capi.mod .
	cp ../../capi/ctqmc_api.o .

cssl:
	cp ../../base/libMM.a .

csml:
	cp ../../base/constants.mod .
	cp ../../base/mmpi.mod .
	cp ../../base/spring.mod .
	cp ../../base/stack.mod .
	cp ../../base/parser.mod .
	cp ../../base/linkedlist.mod .

exe: api cssl csml ctqmc-exe
ctqmc-exe: $(exe_obj)
	$(LINKER) $(exe_obj) -o ctqmc $(LFLAGS) $(LIBS)

lib: api cssl csml ctqmc-lib
ctqmc-lib: $(lib_obj)
	$(ARCHIVER) libctqmc.a $(lib_obj)

pylib: lib
	cpp -E -DPYAPI ctqmc_open.f90 > ctqmc_open_py.f90
	$(F2PY) -c --f90exec=$(F2PYC) -m pyiqist -L. -lctqmc -lMM $(F2PYL) ctqmc_open_py.f90
	rm -fr *.dSYM ctqmc_open_py.f90

.f90.o:
	$(F90) $(FFLAGS) $*.f90

clean:
	rm -f *.mod
	rm -f *.o
	rm -f *.so
	rm -f ctqmc
	rm -f libctqmc.a
	rm -f libMM.a

clean-dat:
	rm -f *.dat
	rm -f *.out

clean-all: clean clean-dat
