##
##
## Introduction
## ============
##
## Makefile for iqist/src/tools/jasmine directory
##
## Usage
## =====
##
## 'make', 'make lib', 'make clean', 'make clean-dat' or 'make clean-all'
##
##

.SUFFIXES: .f90

include ../../build/make.sys



default: all

modn = atomic_control.o atomic_sector.o atomic_global.o
lev1 = atomic_tran.o atomic_print.o atomic_config.o atomic_check.o atomic_gaunt.o atomic_read.o atomic_write.o atomic_basis.o
lev2 = atomic_mksector.o atomic_mkspmat.o atomic_natural.o 
lev3 = atomic_hmat.o atomic_fmat.o atomic_occu.o atomic_driver.o
main = atomic_main.o
MAPI = libMM.a

objects = $(modn) $(lev1) $(lev2) $(lev3) $(main) $(MAPI)

default: all

all: exe

cssl:
	cp ../../common/libMM.a .

csml:
	cp ../../common/constants.mod .
	cp ../../common/parser.mod .
	cp ../../common/linkedlist.mod .

exe: cssl csml atomic_exe
atomic_exe: $(objects)
	$(LINKER) $(objects) -o atomic $(LFLAGS) $(LIBS)

.f90.o:
	$(F90) $(FFLAGS) $*.f90

clean:
	rm -f *.mod
	rm -f *.o
	rm -f atomic
	rm -f libMM.a

clean-dat:
	rm -f *.dat
	rm -f *.bin.*
	rm -f *.out
