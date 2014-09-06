##
##
## Introduction
## ============
##
## This is the top Makefile for the iqist software package. You can use
## it to build all of the executable programs and libraries. Please do
## not modify it by yourself.
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
## 5. build documents and api
## --------------------------------
##
## 'make other'
##
## 6. build specified code
## -----------------------
##
## 'make xxx'
##
## 7. clean the directories
## ------------------------
##
## 'make clean'
##
## here xxx is the code name, please input 'make help-more' to obtain more
## information about the possible values of xxx
##
## Author
## ======
##
## This building system is designed, created, implemented, and maintained
## by Li Huang (email: huangli712@gmail.com).
##
## Department of Physics, University of Fribourg, Fribourg CH-1700, Switzerland
##
## History
## =======
##
## 08/19/2014 by li huang
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
	@echo '     clean             clean all directories                          '
	@echo
	@echo '     solver            build all ctqmc and hfqmc impurity solvers     '
	@echo '     clean-solver      clean objects for impurity solvers             '
	@echo
	@echo '     tool              build auxiliary tools                          '
	@echo '     clean-tool        clean objects for auxiliary tools              '
	@echo
	@echo '     other             build api and documents                        '
	@echo '     clean-other       clean api and doc directories                  '

help-more:
	@echo '  Interacting Quantum Impurity Solver Toolkit Building System'
	@echo
	@echo
	@echo '  targets:'
	@echo '     help              print short usage information        '
	@echo '     help-more         print full usage information         '
	@echo
	@echo '     common            build common modules and util library'
	@echo '     clean-common      clean common directory               '
	@echo
	@echo '     azalea            build azalea code                    '
	@echo '     azalea-lib        build azalea library                 '
	@echo '     clean-azalea      clean azalea directory               '
	@echo
	@echo '     gardenia          build gardenia code                  '
	@echo '     gardenia-lib      build gardenia library               '
	@echo '     clean-gardenia    clean gardenia directory             '
	@echo
	@echo '     narcissus         build narcissus code                 '
	@echo '     narcissus-lib     build narcissus library              '
	@echo '     clean-narcissus   clean narcissus directory            '
	@echo
	@echo '     begonia           build begonia code                   '
	@echo '     begonia-lib       build begonia library                '
	@echo '     clean-begonia     clean begonia directory              '
	@echo
	@echo '     lavender          build lavender code                  '
	@echo '     lavender-lib      build lavender library               '
	@echo '     clean-lavender    clean lavender directory             '
	@echo
	@echo '     pansy             build pansy code                     '
	@echo '     pansy-lib         build pansy library                  '
	@echo '     clean-pansy       clean pansy directory                '
	@echo
	@echo '     manjushaka        build manjushaka code                '
	@echo '     manjushaka-lib    build manjushaka library             '
	@echo '     clean-manjushaka  clean manjushaka directory           '
	@echo
	@echo '     daisy             build daisy code                     '
	@echo '     clean-daisy       clean daisy directory                '
	@echo
	@echo '     jasmine           build jasmine code                   '
	@echo '     clean-jasmine     clean jasmine directory              '
	@echo
	@echo '     hibiscus          build hibiscus code                  '
	@echo '     clean-hibiscus    clean hibiscus directory             '
	@echo
	@echo '     api               generate api library                 '
	@echo '     clean-api         clean api directory                  '
	@echo
	@echo '     doc               generate document                    '
	@echo '     clean-doc         clean doc directory                  '

# all-in-one target
all: build-all
build-all: common solver tool other

clean: clean-all
clean-all: clean-common clean-solver clean-tool clean-other

# targets for common
common: build-common
build-common:
	cd ../common/; pwd; make

clean-common:
	cd ../common/; pwd; make clean

# targets for solver
solver: common build-solver
build-solver: azalea gardenia narcissus begonia lavender pansy manjushaka daisy
clean-solver: clean-azalea clean-gardenia clean-narcissus clean-begonia clean-lavender clean-pansy clean-manjushaka clean-daisy

# targets for tool
tool: common build-tool
build-tool: jasmine hibiscus
clean-tool: clean-jasmine clean-hibiscus

# targets for other
other: common build-other
build-other: api doc
clean-other: clean-api clean-doc

# tagets for ctqmc codes: segment version (azalea, gardenia, narcissus)
azalea: build-azalea
build-azalea:
	cd ../ctqmc/azalea/; pwd; make

azalea-lib: build-azalea-lib
build-azalea-lib:
	cd ../ctqmc/azalea/; pwd; make lib

clean-azalea:
	cd ../ctqmc/azalea/; pwd; make clean

gardenia: build-gardenia
build-gardenia:
	cd ../ctqmc/gardenia/; pwd; make

gardenia-lib: build-gardenia-lib
build-gardenia-lib:
	cd ../ctqmc/gardenia/; pwd; make lib

clean-gardenia:
	cd ../ctqmc/gardenia/; pwd; make clean

narcissus: build-narcissus
build-narcissus:
	cd ../ctqmc/narcissus/; pwd; make

narcissus-lib: build-narcissus-lib
build-narcissus-lib:
	cd ../ctqmc/narcissus/; pwd; make lib

clean-narcissus:
	cd ../ctqmc/narcissus/; pwd; make clean

# tagets for ctqmc codes: general version (begonia, lavender)
begonia: build-begonia
build-begonia:
	cd ../ctqmc/begonia/; pwd; make

begonia-lib: build-begonia-lib
build-begonia-lib:
	cd ../ctqmc/begonia/; pwd; make lib

clean-begonia:
	cd ../ctqmc/begonia/; pwd; make clean

lavender: build-lavender
build-lavender:
	cd ../ctqmc/lavender/; pwd; make

lavender-lib: build-lavender-lib
build-lavender-lib:
	cd ../ctqmc/lavender/; pwd; make lib

clean-lavender:
	cd ../ctqmc/lavender/; pwd; make clean

# tagets for ctqmc codes: general version (pansy, manjushaka)
pansy: build-pansy
build-pansy:
	cd ../ctqmc/pansy/; pwd; make

pansy-lib: build-pansy-lib
build-pansy-lib:
	cd ../ctqmc/pansy/; pwd; make lib

clean-pansy:
	cd ../ctqmc/pansy/; pwd; make clean

manjushaka: build-manjushaka
build-manjushaka:
	cd ../ctqmc/manjushaka/; pwd; make

manjushaka-lib: build-manjushaka-lib
build-manjushaka-lib:
	cd ../ctqmc/manjushaka/; pwd; make lib

clean-manjushaka:
	cd ../ctqmc/manjushaka/; pwd; make clean

# tagets for hfqmc codes: daisy
daisy: build-daisy
build-daisy:
	cd ../hfqmc/daisy/; pwd; make

clean-daisy:
	cd ../hfqmc/daisy/; pwd; make clean

# tagets for atomic codes: jasmine
jasmine: build-jasmine
build-jasmine:
	cd ../tools/jasmine/; pwd; make

clean-jasmine:
	cd ../tools/jasmine/; pwd; make clean

# tagets for toolkit codes: hibiscus
hibiscus: build-hibiscus
build-hibiscus:
	cd ../tools/hibiscus/; pwd; make

clean-hibiscus:
	cd ../tools/hibiscus/; pwd; make clean

# tagets for api: api
api: build-api
build-api:
	cd ../ctqmc/api/; pwd; make

clean-api:
	cd ../ctqmc/api/; pwd; make clean

# tagets for doc: doc
doc: build-doc
build-doc:
	cd ../../doc/guide/; pwd; make

clean-doc:
	cd ../../doc/guide/; pwd; make clean
