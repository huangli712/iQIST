##
##
## Introduction
## ============
##
## Makefile for iqist/src/common directory.
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
## 09/18/2014 by li huang
##
##

.SUFFIXES: .f90

include ../build/make.sys

mod = m_constants.o m_linkedlist.o m_mpi.o m_parser.o m_skynet.o m_sparse.o m_spring.o m_stack.o
sub = s_error.o s_fft.o s_integrator.o s_matrix.o s_spline.o s_util.o s_vector.o
objects = $(mod) $(sub)

default: all

all: build_mod build_sub build_lib

build_mod: $(mod)
build_sub: $(sub)
build_lib: $(objects)
	$(ARCHIVER) libMM.a $(objects)

.f90.o:
	$(F90) $(FFLAGS) $*.f90

clean:
	rm -f *.mod
	rm -f *.o
	rm -f libMM.a

clean-dat:
	rm -f *.dat
	rm -f *.bin.*
	rm -f *.out

clean-all: clean clean-dat
