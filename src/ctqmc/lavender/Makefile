##
##
## Introduction
## ============
##
## Makefile for iqist/src/ctqmc/lavender directory.
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
## 11/06/2014 by li huang
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
MAPI = ctqmc_api.o libMM.a

objects = $(modc) $(dmft) $(core) $(lev1) $(lev2) $(lev3) $(main) $(MAPI)

default: all

all: exe

api:
	cp ../api/api.mod .
	cp ../api/ctqmc_api.o .

cssl:
	cp ../../common/libMM.a .

csml:
	cp ../../common/constants.mod .
	cp ../../common/mmpi.mod .
	cp ../../common/sparse.mod .
	cp ../../common/spring.mod .
	cp ../../common/stack.mod .
	cp ../../common/parser.mod .
	cp ../../common/linkedlist.mod .

exe: api cssl csml ctqmc-exe
ctqmc-exe: $(objects)
	$(LINKER) $(objects) -o ctqmc $(LFLAGS) $(LIBS)

lib: api cssl csml ctqmc-lib
ctqmc-lib: $(objects)
	$(ARCHIVER) libctqmc.a $(objects)

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
	rm -f *.bin.*
	rm -f *.out

clean-all: clean clean-dat
