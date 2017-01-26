##
##
## Introduction
## ============
##
## Makefile for iqist/src/base directory.
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
## 09/18/2014 by li huang (created)
## 01/26/2017 by li huang (last modified)
##
##

.SUFFIXES: .f90

include ../../build/make.sys

mod = m_constants.o m_linkedlist.o m_mpi.o m_parser.o m_sparse.o m_spring.o m_stack.o m_version.o
sub = s_error.o s_fft.o s_integrator.o s_matrix.o s_spline.o s_util.o s_vector.o
objects = $(mod) $(sub)

default: all

all: build-mod build-sub build-lib

build-mod: $(mod)
build-sub: $(sub)
build-lib: $(objects)
	$(ARCHIVER) libMM.a $(objects)

.f90.o:
	$(F90) $(FFLAGS) $*.f90

clean:
	rm -f *.mod
	rm -f *.o
	rm -f libMM.a

clean-dat:
	rm -f *.dat
	rm -f *.out

clean-all: clean clean-dat
