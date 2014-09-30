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
## 09/30/2014 by li huang
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
	$(LINKER) $(objects) makedos.o -o mdos.x $(LFLAGS) $(LIBS)

makescr: makescr.o
	$(LINKER) $(objects) makescr.o -o mscr.x $(LFLAGS) $(LIBS)

maketau: maketau.o
	$(LINKER) $(objects) maketau.o -o mtau.x $(LFLAGS) $(LIBS)

makeups: makeups.o
	$(LINKER) $(objects) makeups.o -o mups.x $(LFLAGS) $(LIBS)

makekra: makekra.o
	$(LINKER) $(objects) makekra.o -o mkra.x $(LFLAGS) $(LIBS)

makesig: makesig.o
	$(LINKER) $(objects) makesig.o -o msig.x $(LFLAGS) $(LIBS)

makestd: makestd.o
	$(LINKER) $(objects) makestd.o -o mstd.x $(LFLAGS) $(LIBS)

makechi: makechi.o
	$(LINKER) $(objects) makechi.o -o mchi.x $(LFLAGS) $(LIBS)

.f90.o:
	$(F90) $(FFLAGS) $*.f90

clean:
	rm -f *.mod
	rm -f *.o
	rm -f *.x
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
