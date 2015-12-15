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
## 08/17/2015 by li huang (last modified)
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

exe: cssl csml makedos makescr maketau makeups makekra makesig makestd makechi

makedos: cssl csml makedos.o
	$(LINKER) makedos.o $(objects) -o mdos $(LFLAGS) $(LIBS)

makescr: cssl csml makescr.o
	$(LINKER) makescr.o $(objects) -o mscr $(LFLAGS) $(LIBS)

maketau: cssl csml maketau.o
	$(LINKER) maketau.o $(objects) -o mtau $(LFLAGS) $(LIBS)

makeups: cssl csml makeups.o
	$(LINKER) makeups.o $(objects) -o mups $(LFLAGS) $(LIBS)

makekra: cssl csml makekra.o
	$(LINKER) makekra.o $(objects) -o mkra $(LFLAGS) $(LIBS)

makesig: cssl csml makesig.o
	$(LINKER) makesig.o $(objects) -o msig $(LFLAGS) $(LIBS)

makestd: cssl csml makestd.o
	$(LINKER) makestd.o $(objects) -o mstd $(LFLAGS) $(LIBS)

makechi: cssl csml makechi.o
	$(LINKER) makechi.o $(objects) -o mchi $(LFLAGS) $(LIBS)

.f90.o:
	$(F90) $(FFLAGS) $*.f90

clean:
	rm -f *.mod
	rm -f *.o
	rm -f mdos mscr mtau mups mkra msig mstd mchi
	rm -f libMM.a

clean-dat:
	rm -f dos.*
	rm -f scr.*
	rm -f tau.*
	rm -f ups.*
	rm -f kra.*
	rm -f sig.*
	rm -f std.*
	rm -f *.dat

clean-all: clean clean-dat
