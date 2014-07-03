.SUFFIXES: .f90

include ./make.sys

default: all

modn = mod_constants.o mod_control.o mod_subspace.o mod_mkgood.o mod_fmat.o mod_global.o
lev1 = atomic_util.o atomic_mksubs.o
lev2 = atomic_stream.o atomic_basis.o atomic_hmat.o
lev3 = atomic_mkfmat.o
dump = atomic_dumper.o atomic_print.o atomic_write.o
main = atomic_main.o

objects = $(modn) $(lev1) $(lev2) $(lev3) $(dump) $(main)

all: atomic

atomic: $(objects)
	$(LINKER) $(objects) -o atomic $(LFLAGS) $(LIBS) 

.f90.o:
	$(F90) $(FFLAGS) $*.f90

clean:
	rm -f *.mod
	rm -f *.o
	rm -f atomic
	rm -rf *genmod*
clean-dat:
	rm -f *.dat
	rm -f *.bin
	rm -f *.out

clean-all: clean clean-dat
