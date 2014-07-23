.SUFFIXES: .f90

include ./make.sys

default: all

modn = mod_constants.o mod_control.o mod_sector.o mod_global.o
lev1 = atomic_util.o atomic_tran.o atomic_print.o atomic_config.o atomic_check.o atomic_gaunt.o atomic_read.o atomic_write.o atomic_basis.o
lev2 = atomic_mksector.o atomic_mkspmat.o atomic_natural.o 
lev3 = atomic_hmat.o atomic_fmat.o atomic_occu.o atomic_driver.o
main = atomic_main.o

objects = $(modn) $(lev1) $(lev2) $(lev3) $(main)

all: atomic

atomic: $(objects)
	$(LINKER) $(objects) -o atomic $(LFLAGS) $(LIBS) 

.f90.o:
	$(F90) $(FFLAGS) $*.f90

clean:
	rm -f *.mod
	rm -f *.o
	rm -f atomic
	rm -f *genmod*
	rm -f *.i90
clean-dat:
	rm -f *.dat
	rm -f *.bin
	rm -f *.out

clean-all: clean clean-dat
