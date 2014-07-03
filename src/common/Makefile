.SUFFIXES: .f90

include ../build/make.sys

default: all

mod = mod_const.o mod_mpi.o mod_parser.o mod_skynet.o mod_sparse.o mod_spring.o mod_stack.o

all: ctqmc_exe

ctqmc_exe: $(objects)
	$(LINKER) $(objects) -o ctqmc $(LFLAGS) $(LIBS)

ctqmc_lib: $(objects)
	$(ARCHIVER) libctqmc.a $(objects)

.f90.o:
	$(F90) $(FFLAGS) $*.f90

clean:
	rm -f *.mod
	rm -f *.o
	rm -f ctqmc
	rm -f libctqmc.a

clean-dat:
	rm -f *.dat
	rm -f *.bin.*
	rm -f *.out

clean-all: clean clean-dat
