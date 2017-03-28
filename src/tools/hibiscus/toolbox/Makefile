##
##
## Introduction
## ============
##
## Makefile for iqist/src/tools/hibiscus/toolbox directory.
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
## 10/02/2014 by li huang (created)
## 03/29/2017 by li huang (last modified)
##
##

.SUFFIXES: .f90

include ../../../../build/make.sys

CLIB = libMM.a

objects = $(CLIB)

default: all

all: exe

cssl:
	cp ../../../base/libMM.a .

csml:
	cp ../../../base/constants.mod .

exe: cssl csml makechi makedos makescr

makechi: cssl csml makechi.o
	$(LINKER) makechi.o $(objects) -o mchi $(LFLAGS) $(LIBS)

makedos: cssl csml makedos.o
	$(LINKER) makedos.o $(objects) -o mdos $(LFLAGS) $(LIBS)

makescr: cssl csml makescr.o
	$(LINKER) makescr.o $(objects) -o mscr $(LFLAGS) $(LIBS)

.f90.o:
	$(F90) $(FFLAGS) $*.f90

clean:
	rm -f *.mod
	rm -f *.o
	rm -f mchi mdos mscr
	rm -f libMM.a

clean-dat:
	rm -f dos.*
	rm -f scr.*
	rm -f *.dat

clean-all: clean clean-dat
