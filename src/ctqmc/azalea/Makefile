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

modc = ctqmc_control.o ctqmc_context.o
dmft = ctqmc_dmft.o
core = ctqmc_solver.o
lev1 = ctqmc_flavor.o ctqmc_update.o
lev2 = ctqmc_record.o ctqmc_status.o ctqmc_stream.o ctqmc_util.o
lev3 = ctqmc_dump.o ctqmc_print.o
main = ctqmc_main.o
open = ctqmc_open.o
CAPI = ctqmc_api.o
CLIB = libMM.a

objects = $(modc) $(dmft) $(core) $(lev1) $(lev2) $(lev3)
exe_obj = $(objects) $(main) $(CAPI) $(CLIB)
lib_obj = $(objects) $(open) $(CAPI) $(CLIB)

default: all

all: exe

api:
	cp ../../api/capi.mod .
	cp ../../api/ctqmc_api.o .

cssl:
	cp ../../common/libMM.a .

csml:
	cp ../../common/constants.mod .
	cp ../../common/mmpi.mod .
	cp ../../common/spring.mod .
	cp ../../common/stack.mod .
	cp ../../common/parser.mod .
	cp ../../common/linkedlist.mod .

exe: api cssl csml ctqmc-exe
ctqmc-exe: $(exe_obj)
	$(LINKER) $(exe_obj) -o ctqmc $(LFLAGS) $(LIBS)

lib: api cssl csml ctqmc-lib
ctqmc-lib: $(lib_obj)
	$(ARCHIVER) libctqmc.a $(lib_obj)

.f90.o:
	$(F90) $(FFLAGS) $*.f90

clean:
	rm -f *.mod
	rm -f *.o
	rm -f ctqmc
	rm -f libctqmc.a
	rm -f libMM.a

clean-dat:
	rm -f *.dat
	rm -f *.out

clean-all: clean clean-dat
