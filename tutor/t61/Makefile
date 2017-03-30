##
##
## Introduction
## ============
##
## Makefile for iqist/tutor/t61 directory.
##
## Usage
## =====
##
## 'make', 'make clean', 'make clean-dat' or 'make clean-all'
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
## 01/24/2015 by li huang (created)
## 03/30/2017 by li huang (last modified)
##
##

.SUFFIXES: .f90

include ../../build/make.sys

main = template.o

CLIB = libctqmc.a libMM.a

objects = $(main) $(CLIB)

default: all

all: exe

api:
	cp ../../src/capi/capi.mod .
	cp ../../src/ctqmc/narcissus/libctqmc.a .

cssl:
	cp ../../src/base/libMM.a .

csml:
	cp ../../src/base/constants.mod .
	cp ../../src/base/mmpi.mod .

exe: api cssl csml main-exe
main-exe: $(objects)
	$(LINKER) $(objects) -o main $(LFLAGS) $(LIBS)

.f90.o:
	$(F90) $(FFLAGS) $*.f90

clean:
	rm -f *.mod
	rm -f *.o
	rm -f main
	rm -f libctqmc.a
	rm -f libMM.a

clean-dat:
	rm -f *.dat
	rm -f *.out

clean-all: clean clean-dat
