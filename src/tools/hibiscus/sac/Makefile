##
##
## Introduction
## ============
##
## Makefile for iqist/src/tools/hibiscus/sac directory.
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

modc = sai_control.o sai_context.o
core = sai_stream.o sai_toolbox.o sai_update.o sai_record.o
dump = sai_dump.o sai_print.o
main = sai_main.o
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

exe: cssl csml sai-exe
sai-exe: $(objects)
	$(LINKER) $(objects) -o sai $(LFLAGS) $(LIBS)

.f90.o:
	$(F90) $(FFLAGS) $*.f90

clean:
	rm -f *.mod
	rm -f *.o
	rm -f sai
	rm -f libMM.a

clean-dat:
	rm -f *.dat
	rm -f *.bin
	rm -f *.out

clean-all: clean clean-dat
