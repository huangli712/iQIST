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
## 09/30/2014 by li huang
##
##

.SUFFIXES: .f90

include ../../../build/make.sys

modc = sac_control.o sac_context.o
core = sac_stream.o sac_toolbox.o sac_update.o sac_record.o
dump = sac_dump.o sac_print.o
main = sac_main.o
MAPI = libMM.a

objects = $(modc) $(core) $(dump) $(main) $(MAPI)

default: all

all: exe

cssl:
	cp ../../../common/libMM.a .

csml:
	cp ../../../common/constants.mod .
	cp ../../../common/mmpi.mod .
	cp ../../../common/spring.mod .
	cp ../../../common/parser.mod .
	cp ../../../common/linkedlist.mod .

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
	rm -f *.bin
	rm -f *.out

clean-all: clean clean-dat
