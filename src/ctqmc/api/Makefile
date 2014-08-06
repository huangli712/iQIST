
.SUFFIXES: .f90

include ../../build/make.sys

mod = ctqmc_api.o

default: all

all: build_mod

build_mod : $(mod)

.f90.o:
	$(F90) $(FFLAGS) $*.f90

clean:
	rm -f *.mod
	rm -f *.o

clean-dat:
	rm -f *.dat
	rm -f *.bin.*
	rm -f *.out

clean-all: clean clean-dat
