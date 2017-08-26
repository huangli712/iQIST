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
## 'make', 'make clean', 'make clean-dat' or 'make clean-all'
##
## Author
## ======
##
## This building system is designed, created, and maintained by
##
## Li Huang // email: lihuang.dmft@gmail.com
##
## History
## =======
##
## 01/17/2017 by li huang (created)
## 01/10/2018 by li huang (last modified)
##
##

.SUFFIXES: .f90

include ../../../build/make.sys

modc = dt_control.o dt_context.o
lev1 = dt_run.o dt_df.o
lev2 = dt_stream.o dt_util.o
lev3 = dt_dump.o dt_print.o
main = dt_main.o
mlib = libMM.a

objects = $(modc) $(lev1) $(lev2) $(lev3) $(main) $(mlib)

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

exe: cssl csml $(objects)
	$(LINKER) $(objects) -o dmft $(LFLAGS) $(LIBS) -L/opt/local/lib/ -lfftw3

.f90.o:
	$(F90) -I/opt/local/include $(FFLAGS) $*.f90

clean:
	rm -f *.mod
	rm -f *.o
	rm -f dmft
	rm -f libMM.a

clean-dat:
	rm -f *.dat
	rm -f *.out

clean-all: clean clean-dat
