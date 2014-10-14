##
##
## Introduction
## ============
##
## Makefile for iqist/src/tools/hibiscus/entropy2 directory.
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

modc = maxent_control.o maxent_context.o
util = maxent_utils.o 
lev1 = maxent_stream.o maxent_model.o
lev2 = maxent_makeaw.o
lev3 = maxent_driver.o
dump = maxent_dump.o maxent_print.o
main = maxent_main.o
MAPI = libMM.a

objects = $(modc) $(util) $(lev1) $(lev2) $(lev3) $(dump) $(main) $(MAPI)

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

exe: cssl csml maxent-exe
maxent-exe: $(objects)
	$(LINKER) $(objects) -o maxent $(LFLAGS) $(LIBS)

.f90.o:
	$(F90) $(FFLAGS) $*.f90

clean:
	rm -f *.mod
	rm -f *.o
	rm -f maxent
	rm -f libMM.a

clean-dat:
	rm -f *.dat
	rm -f *.bin
	rm -f *.out

clean-all: clean clean-dat
