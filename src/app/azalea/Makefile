##
##
## Introduction
## ============
##
## Makefile for iqist/src/app/azalea directory.
##
## Usage
## =====
##
## 'make', 'make clean', 'make clean-dat' or 'make clean-all'
##
## Author
## ======
##
## This building system is designed, created, and maintained by
##
## Li Huang // email: lihuang.dmft@gmail.com
##
## History
## =======
##
## 01/17/2017 by li huang (created)
## 08/29/2017 by li huang (last modified)
##
##

.SUFFIXES: .f90

include ../../../build/make.sys

modc = dmft_control.o dmft_context.o
impl = impl_dmft.o impl_df.o
lev1 = dmft_stream.o dmft_util.o
lev2 = dmft_dump.o dmft_print.o
main = dmft_main.o
mlib = libMM.a

objects = $(modc) $(impl) $(lev1) $(lev2) $(main) $(mlib)

default: all

all: exe

cssl:
	cp ../../base/libMM.a .

csml:
	cp ../../base/constants.mod .
	cp ../../base/mmpi.mod .
	cp ../../base/parser.mod .
	cp ../../base/linkedlist.mod .
	cp ../../base/version.mod .

exe: cssl csml df-exe
df-exe: $(exe_obj)
	$(LINKER) $(exe_obj) -o df $(LFLAGS) $(LIBS) -L/opt/local/lib/ -lfftw3

lib: cssl csml df-lib
df-lib: $(lib_obj)
	$(ARCHIVER) libdf.a $(lib_obj)

.f90.o:
	$(F90) -I/opt/local/include $(FFLAGS) $*.f90

clean:
	rm -f *.mod
	rm -f *.o
	rm -f *.so
	rm -f df
	rm -f libdf.a
	rm -f libMM.a

clean-dat:
	rm -f *.dat
	rm -f *.out

clean-all: clean clean-dat
