##
##
## Introduction
## ============
##
## Makefile for iqist/src/atomic directory.
##
## Usage
## =====
##
## 'make', 'make clean', 'make clean-dat' or 'make clean-all'
##
## Author
## ======
##
## This building system is designed, created, and maintained by
##
## Li Huang // email: huangli@caep.cn
##
## History
## =======
##
## 02/28/2015 by li huang (created)
## 01/02/2024 by li huang (last modified)
##
##

.SUFFIXES: .f90

include ../../build/make.inc

modc = atomic_control.o atomic_context.o
core = atomic_driver.o
lev1 = atomic_full.o atomic_sector.o
lev2 = atomic_stream.o atomic_util.o
lev3 = atomic_dump.o atomic_print.o
main = atomic_main.o
mlib = libflink.a

objects = $(modc) $(core) $(lev1) $(lev2) $(lev3) $(main) $(mlib)

default: all

all: exe

flink: flink_lib flink_mod

flink_lib:
	cp $(FLINK)/libflink.a .

flink_mod:
	cp $(FLINK)/constants.mod .
	cp $(FLINK)/parser.mod .
	cp $(FLINK)/linkedlist.mod .

exe: flink $(objects)
	$(LINKER) $(objects) -o atomic $(LFLAGS) $(LIBS)

.f90.o:
	$(F90) $(FFLAGS) $*.f90

clean:
	rm -f *.mod
	rm -f *.o
	rm -f atomic
	rm -f libflink.a

clean-dat:
	rm -f *.dat
	rm -f *.out
	rm -f *.cix
	rm -f solver.umat.in

clean-all: clean clean-dat
