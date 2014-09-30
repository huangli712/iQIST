##
##
## Introduction
## ============
##
## Makefile for iqist/src/tools/hibiscus/entropy1 directory.
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

include ./make.sys

default: all

modm = mod_mpi.o mod_spring.o
modn = entropy_constants.o entropy_control.o entropy_context.o
core = entropy_stream.o entropy_toolbox.o entropy_update.o
dump = entropy_dump.o entropy_print.o
main = entropy_main.o

objects = $(modm) $(modn) $(core) $(dump) $(main)

all: entropy

entropy: $(objects)
	$(LINKER) $(objects) -o entropy $(LFLAGS) $(LIBS)

.f90.o:
	$(F90) $(FFLAGS) $*.f90

clean:
	rm -f *.mod
	rm -f *.o
	rm -f entropy

clean-dat:
	rm -f *.dat
	rm -f *.bin
	rm -f *.out

clean-all: clean clean-dat
