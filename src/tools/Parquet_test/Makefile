.SUFFIXES: .f90

include ./make.sys

default: all

modn = analyze_constants.o analyze_control.o analyze_context.o
dump = analyze_print.o
readin = readin.o
calculate = analyze_calculate.o
utility = analyze_util.o
lavender = analyze_lavender.o

objects = $(modn) $(dump) $(readin) $(calculate) $(utility) $(lavender)

all: analyze_lavender

analyze_lavender: $(objects)
	$(LINKER) $(objects) -o analyze_lavender $(LFLAGS) $(LIBS)

.f90.o:
	$(F90) $(FFLAGS) $*.f90

clean:
	rm -f *.mod
	rm -f *.o
	rm -f *~
	rm -f analyze_lavender

clean-dat:
	rm -f *.dat
	rm -f *.bin.*
	rm -f *.out

clean-all: clean clean-dat
