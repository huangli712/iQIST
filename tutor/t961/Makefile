##
##
## Introduction
## ============
##
## Makefile for iqist/tutor/t961 directory.
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
## Li Huang // email: huangli712@gmail.com
##
## History
## =======
##
## 01/24/2015 by li huang
##
##

.SUFFIXES: .f90

include ../../src/build/make.sys

main = template.o
MAPI = libctqmc.a libMM.a

objects = $(main) $(MAPI)

default: all

all: exe

api:
	cp ../../src/api/api.mod .
	cp ../../src/ctqmc/azalea/libctqmc.a .

cssl:
	cp ../../src/common/libMM.a .

csml:
	cp ../../src/common/constants.mod .
	cp ../../src/common/mmpi.mod .

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
	rm -f *.bin.*
	rm -f *.out

clean-all: clean clean-dat
