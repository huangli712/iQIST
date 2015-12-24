##
##
## Introduction
## ============
##
## Makefile for iqist/src/tools/hibiscus/stoch directory.
##
## Usage
## =====
##
## 'make' or 'make clean'
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
## 09/30/2014 by li huang (created)
## 08/17/2015 by li huang (last modified)
##
##

.SUFFIXES: .f90

include ../../../../build/make.sys

modc = sac_control.o sac_context.o
core = sac_stream.o sac_toolbox.o sac_update.o sac_record.o
dump = sac_dump.o sac_print.o
main = sac_main.o

CLIB = libMM.a

objects = $(modc) $(core) $(dump) $(main) $(CLIB)

default: all

all: exe

cssl:
	cp ../../../base/libMM.a .

csml:
	cp ../../../base/constants.mod .
	cp ../../../base/mmpi.mod .
	cp ../../../base/spring.mod .
	cp ../../../base/parser.mod .
	cp ../../../base/linkedlist.mod .

exe: cssl csml sac-exe
sac-exe: $(objects)
	$(LINKER) $(objects) -o sac $(LFLAGS) $(LIBS)

.f90.o:
	$(F90) $(FFLAGS) $*.f90

clean:
	rm -f *.mod
	rm -f *.o
	rm -f sac
	rm -f libMM.a

clean-dat:
	rm -f *.dat
	rm -f *.out

clean-all: clean clean-dat
