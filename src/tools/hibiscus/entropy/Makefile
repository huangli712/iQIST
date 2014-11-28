##
##
## Introduction
## ============
##
## Makefile for iqist/src/tools/hibiscus/entropy directory.
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
## Li Huang // email: huangli712@gmail.com
##
## History
## =======
##
## 09/30/2014 by li huang
##
##

.SUFFIXES: .f90

include ../../../build/make.sys

modc = entropy_control.o entropy_context.o
core = entropy_stream.o entropy_toolbox.o entropy_update.o
dump = entropy_dump.o entropy_print.o
main = entropy_main.o
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

exe: cssl csml entropy-exe
entropy-exe: $(objects)
	$(LINKER) $(objects) -o entropy $(LFLAGS) $(LIBS)

.f90.o:
	$(F90) $(FFLAGS) $*.f90

clean:
	rm -f *.mod
	rm -f *.o
	rm -f entropy
	rm -f libMM.a

clean-dat:
	rm -f *.dat
	rm -f *.bin
	rm -f *.out

clean-all: clean clean-dat
