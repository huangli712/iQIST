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
## Li Huang // email: huangli712@gmail.com
##
## History
## =======
##
## 10/02/2014 by li huang
##
##

.SUFFIXES: .f90

include ../../../build/make.sys

MAPI = libMM.a

objects = $(MAPI)

default: all

all: exe

cssl:
	cp ../../../common/libMM.a .

csml:
	cp ../../../common/constants.mod .

exe: cssl csml makedos makescr maketau makeups makekra makesig makestd makechi

makedos: makedos.o
	$(LINKER) makedos.o $(objects) -o mdos $(LFLAGS) $(LIBS)

makescr: makescr.o
	$(LINKER) makescr.o $(objects) -o mscr $(LFLAGS) $(LIBS)

maketau: maketau.o
	$(LINKER) maketau.o $(objects) -o mtau $(LFLAGS) $(LIBS)

makeups: makeups.o
	$(LINKER) makeups.o $(objects) -o mups $(LFLAGS) $(LIBS)

makekra: makekra.o
	$(LINKER) makekra.o $(objects) -o mkra $(LFLAGS) $(LIBS)

makesig: makesig.o
	$(LINKER) makesig.o $(objects) -o msig $(LFLAGS) $(LIBS)

makestd: makestd.o
	$(LINKER) makestd.o $(objects) -o mstd $(LFLAGS) $(LIBS)

makechi: makechi.o
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

clean-all: clean clean-dat
