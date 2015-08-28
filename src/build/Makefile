##
##
## Introduction
## ============
##
## This is the top Makefile for the iQIST software package. You can use
## it to build all of the executable programs and libraries. Please do
## not modify it by yourself even you are very familiar with iQIST.
##
## Usage
## =====
##
## 1. get help
## -----------
##
## 'make help' or 'make help-more'
##
## 2. quick build
## --------------
##
## 'make all'
##
## 3. build impurity solver
## ------------------------
##
## 'make solver'
##
## 4. build auxiliary tools
## ------------------------
##
## 'make tool'
##
## 5. build documents
## ------------------
##
## 'make other'
##
## 6. build specified code
## -----------------------
##
## 'make xxx'
##
## here xxx is the code name, please input 'make help-more' to obtain more
## information about the possible values of xxx
##
## 7. clean the directories
## ------------------------
##
## 'make clean'
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
## 02/28/2015 by li huang (created)
## 08/17/2015 by li huang (last modified)
##
##

# targets for help
help:
	@echo '  Interacting Quantum Impurity Solver Toolkit Building System'
	@echo
	@echo
	@echo '  targets:'
	@echo '     help              print short usage information                  '
	@echo '     help-more         print full usage information                   '
	@echo
	@echo '     all               build all components, apis, and documents, etc '
	@echo '     lib               build all components as libs                   '
	@echo '     clean             clean all directories                          '
	@echo
	@echo '     solver            build all ctqmc and hfqmc components           '
	@echo '     solver-lib        build all ctqmc and hfqmc components as libs   '
	@echo '     clean-solver      clean objects for impurity solvers             '
	@echo
	@echo '     tool              build auxiliary tool components                '
	@echo '     tool-lib          build auxiliary tool components as libs        '
	@echo '     clean-tool        clean objects for auxiliary tools              '
	@echo
	@echo '     other             build documents                                '
	@echo '     other-lib         build documents                                '
	@echo '     clean-other       clean doc directories                          '

help-more:
	@echo '  Interacting Quantum Impurity Solver Toolkit Building System'
	@echo
	@echo
	@echo '  targets:'
	@echo '     help              print short usage information        '
	@echo '     help-more         print full usage information         '
	@echo
	@echo '     azalea            build azalea code                    '
	@echo '     azalea-lib        build azalea library for Fortran     '
	@echo '     azalea-pylib      build azalea library for Python      '
	@echo '     clean-azalea      clean azalea directory               '
	@echo
	@echo '     gardenia          build gardenia code                  '
	@echo '     gardenia-lib      build gardenia library for Fortran   '
	@echo '     gardenia-pylib    build gardenia library for Python    '
	@echo '     clean-gardenia    clean gardenia directory             '
	@echo
	@echo '     narcissus         build narcissus code                 '
	@echo '     narcissus-lib     build narcissus library for Fortran  '
	@echo '     narcissus-pylib   build narcissus library for Python   '
	@echo '     clean-narcissus   clean narcissus directory            '
	@echo
	@echo '     begonia           build begonia code                   '
	@echo '     begonia-lib       build begonia library for Fortran    '
	@echo '     begonia-pylib     build begonia library for Python     '
	@echo '     clean-begonia     clean begonia directory              '
	@echo
	@echo '     lavender          build lavender code                  '
	@echo '     lavender-lib      build lavender library for Fortran   '
	@echo '     lavender-pylib    build lavender library for Python    '
	@echo '     clean-lavender    clean lavender directory             '
	@echo
	@echo '     pansy             build pansy code                     '
	@echo '     pansy-lib         build pansy library for Fortran      '
	@echo '     pansy-pylib       build pansy library for Python       '
	@echo '     clean-pansy       clean pansy directory                '
	@echo
	@echo '     manjushaka        build manjushaka code                '
	@echo '     manjushaka-lib    build manjushaka library for Fortran '
	@echo '     manjushaka-pylib  build manjushaka library for Python  '
	@echo '     clean-manjushaka  clean manjushaka directory           '
	@echo
	@echo '     daisy             build daisy code                     '
	@echo '     daisy-lib         build daisy library for Fortran      '
	@echo '     daisy-pylib       build daisy library for Python       '
	@echo '     clean-daisy       clean daisy directory                '
	@echo
	@echo '     jasmine           build jasmine code                   '
	@echo '     jasmine-lib       build jasmine library for Fortran    '
	@echo '     jasmine-pylib     build jasmine library for Python     '
	@echo '     clean-jasmine     clean jasmine directory              '
	@echo
	@echo '     hibiscus          build hibiscus code                  '
	@echo '     hibiscus-lib      build hibiscus code                  '
	@echo '     clean-hibiscus    clean hibiscus directory             '
	@echo
	@echo '     capi              generate capi library                '
	@echo '     capi-lib          generate capi library                '
	@echo '     clean-capi        clean capi directory                 '
	@echo
	@echo '     common            build csml modules and cssl library  '
	@echo '     common-lib        build csml modules and cssl library  '
	@echo '     clean-common      clean common directory               '
	@echo
	@echo '     doc               generate document                    '
	@echo '     doc-lib           generate document                    '
	@echo '     clean-doc         clean doc directory                  '

# all-in-one building target
all: capi common solver tool other

# all-in-one building target (lib mode)
lib: capi-lib common-lib solver-lib tool-lib other-lib

# all-in-one cleaning target
clean: clean-capi clean-common clean-solver clean-tool clean-other

# targets for solver
solver: capi common \
        azalea gardenia narcissus \
        begonia lavender \
        pansy manjushaka \
        daisy
solver-lib: capi-lib common-lib \
	    azalea-lib gardenia-lib narcissus-lib \
	    begonia-lib lavender-lib \
	    pansy-lib manjushaka-lib \
	    daisy-lib
clean-solver: clean-capi clean-common \
              clean-azalea clean-gardenia clean-narcissus \
              clean-begonia clean-lavender \
              clean-pansy clean-manjushaka \
              clean-daisy

# targets for tool
tool: capi common jasmine hibiscus
tool-lib: capi-lib common-lib jasmine-lib hibiscus-lib
clean-tool: clean-capi clean-common clean-jasmine clean-hibiscus

# targets for other
other: doc
other-lib: doc-lib
clean-other: clean-doc

# targets for ctqmc codes: segment version (azalea, gardenia, narcissus)
azalea:
	cd ../ctqmc/azalea/; pwd; make

azalea-lib:
	cd ../ctqmc/azalea/; pwd; make lib

azalea-pylib:
	cd ../ctqmc/azalea/; pwd; make pylib

clean-azalea:
	cd ../ctqmc/azalea/; pwd; make clean

gardenia:
	cd ../ctqmc/gardenia/; pwd; make

gardenia-lib:
	cd ../ctqmc/gardenia/; pwd; make lib

gardenia-pylib:
	cd ../ctqmc/gardenia/; pwd; make pylib

clean-gardenia:
	cd ../ctqmc/gardenia/; pwd; make clean

narcissus:
	cd ../ctqmc/narcissus/; pwd; make

narcissus-lib:
	cd ../ctqmc/narcissus/; pwd; make lib

narcissus-pylib:
	cd ../ctqmc/narcissus/; pwd; make pylib

clean-narcissus:
	cd ../ctqmc/narcissus/; pwd; make clean

# targets for ctqmc codes: general version (begonia, lavender)
begonia:
	cd ../ctqmc/begonia/; pwd; make

begonia-lib:
	cd ../ctqmc/begonia/; pwd; make lib

begonia-pylib:
	cd ../ctqmc/begonia/; pwd; make pylib

clean-begonia:
	cd ../ctqmc/begonia/; pwd; make clean

lavender:
	cd ../ctqmc/lavender/; pwd; make

lavender-lib:
	cd ../ctqmc/lavender/; pwd; make lib

lavender-pylib:
	cd ../ctqmc/lavender/; pwd; make pylib

clean-lavender:
	cd ../ctqmc/lavender/; pwd; make clean

# targets for ctqmc codes: general version (pansy, manjushaka)
pansy:
	cd ../ctqmc/pansy/; pwd; make

pansy-lib:
	cd ../ctqmc/pansy/; pwd; make lib

pansy-pylib:
	cd ../ctqmc/pansy/; pwd; make pylib

clean-pansy:
	cd ../ctqmc/pansy/; pwd; make clean

manjushaka:
	cd ../ctqmc/manjushaka/; pwd; make

manjushaka-lib:
	cd ../ctqmc/manjushaka/; pwd; make lib

manjushaka-pylib:
	cd ../ctqmc/manjushaka/; pwd; make pylib

clean-manjushaka:
	cd ../ctqmc/manjushaka/; pwd; make clean

# targets for hfqmc codes: daisy
daisy:
	cd ../hfqmc/daisy/; pwd; make

daisy-lib:
	cd ../hfqmc/daisy/; pwd; make lib

daisy-pylib:
	cd ../hfqmc/daisy/; pwd; make pylib

clean-daisy:
	cd ../hfqmc/daisy/; pwd; make clean

# targets for atomic codes: jasmine
jasmine:
	cd ../tools/jasmine/; pwd; make

jasmine-lib:
	cd ../tools/jasmine/; pwd; make lib

jasmine-pylib:
	cd ../tools/jasmine/; pwd; make pylib

clean-jasmine:
	cd ../tools/jasmine/; pwd; make clean

# targets for toolkit codes: hibiscus
hibiscus:
	cd ../tools/hibiscus/; pwd; make

hibiscus-lib:
	cd ../tools/hibiscus/; pwd; make

clean-hibiscus:
	cd ../tools/hibiscus/; pwd; make clean

# targets for capi: capi
capi:
	cd ../capi/; pwd; make

capi-lib:
	cd ../capi/; pwd; make

clean-capi:
	cd ../capi/; pwd; make clean

# targets for common: common
common:
	cd ../common/; pwd; make

common-lib:
	cd ../common/; pwd; make

clean-common:
	cd ../common/; pwd; make clean

# targets for doc: doc
doc:
	cd ../../doc/manual/; pwd; make full

clean-doc:
	cd ../../doc/manual/; pwd; make clean
