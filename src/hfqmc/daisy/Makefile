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

modn = hfqmc_constants.o hfqmc_control.o hfqmc_context.o
dmft = hfqmc_dmft.o
core = hfqmc_solver.o
lev1 = hfqmc_update.o
lev2 = hfqmc_record.o hfqmc_stream.o
lev3 = hfqmc_fourier.o hfqmc_util.o
dump = hfqmc_dump.o hfqmc_print.o
main = hfqmc_main.o

objects = $(modm) $(modn) $(dmft) $(core) $(lev1) $(lev2) $(lev3) $(dump) $(main)

default: all
all: hfqmc

hfqmc: $(objects)
	$(LINKER) $(objects) -o hfqmc $(LFLAGS) $(LIBS)

.f90.o:
	$(F90) $(FFLAGS) $*.f90

clean:
	rm -f *.mod
	rm -f *.o
	rm -f hfqmc

clean-dat:
	rm -f *.dat
	rm -f *.bin.*
	rm -f *.out

clean-all: clean clean-dat
