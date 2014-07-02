.SUFFIXES: .f90

include ./make.sys

default: makedos makescr maketau makeups makekra makesig makestd makechi

makedos: mconst.o makedos.o
	$(LINKER) mconst.o makedos.o -o mdos.x $(LFLAGS) $(LIBS)

makescr: mconst.o makescr.o
	$(LINKER) mconst.o makescr.o -o mscr.x $(LFLAGS) $(LIBS)

maketau: mconst.o maketau.o
	$(LINKER) mconst.o maketau.o -o mtau.x $(LFLAGS) $(LIBS)

makeups: mconst.o makeups.o
	$(LINKER) mconst.o makeups.o -o mups.x $(LFLAGS) $(LIBS)

makekra: mconst.o makekra.o
	$(LINKER) mconst.o makekra.o -o mkra.x $(LFLAGS) $(LIBS)

makesig: mconst.o makesig.o
	$(LINKER) mconst.o makesig.o -o msig.x $(LFLAGS) $(LIBS)

makestd: mconst.o makestd.o
	$(LINKER) mconst.o makestd.o -o mstd.x $(LFLAGS) $(LIBS)

makechi: mconst.o makechi.o
	$(LINKER) mconst.o makechi.o -o mchi.x $(LFLAGS) $(LIBS)

.f90.o:
	$(F90) $(FFLAGS) $*.f90

clean:
	rm -f *.mod
	rm -f *.o
	rm -f *.x

clean-dat:
	rm -f dos.*
	rm -f scr.*
	rm -f tau.*
	rm -f ups.*
	rm -f kra.*
	rm -f sig.*
	rm -f std.*

clean-all: clean clean-dat
