.SUFFIXES: .f90

include ../build/make.sys

default: all

mod = m_constants.o m_linkedlist.o m_mpi.o m_parser.o m_skynet.o m_sparse.o m_spring.o m_stack.o
src = s_error.o s_fft.o s_math.o s_matrix.o s_spline.o s_util.o
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
