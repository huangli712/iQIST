.SUFFIXES: .f90

include ../build/make.sys

default: all

mod = mod_const.o mod_mpi.o mod_parser.o mod_skynet.o mod_sparse.o mod_spring.o mod_stack.o
src =
objects = $(mod) $(src)

all: build_mod build_src build_lib

build_mod: $(mod)
build_src: $(src)
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
