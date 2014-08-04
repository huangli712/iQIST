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
## 5. build documents and libraries
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
## 04/08/2014 by li huang
##
##

# targets for help
help:
	@echo '  Interacting Quantum Impurity Solver Toolkit Building System'
	@echo
	@echo
	@echo '  targets:'
	@echo '     help              print short usage information                      '
	@echo '     help-more         print full usage information                       '
	@echo '     all               build all components, libraries, and documents, etc'
	@echo '     clean             clean all directories                              '
	@echo '     solver            build all ctqmc and hfqmc impurity solvers         '
	@echo '     clean-solver      clean objects for impurity solvers                 '
	@echo '     tool              build auxiliary tools                              '
	@echo '     clean-tool        clean objects for auxiliary tools                  '
	@echo '     other             build api and documents                            '
	@echo '     clean-other       clean api and doc directories                      '

help-more:
	@echo '  Interacting Quantum Impurity Solver Toolkit Building System'
	@echo
	@echo
	@echo '  targets:'
	@echo '     help              print short usage information        '
	@echo '     help-more         print full usage information         '
	@echo '     common            build common modules and util library'
	@echo '     clean-common      clean common directory               '
	@echo '     azalea            build azalea code                    '
	@echo '     clean-azalea      clean azalea directory               '
	@echo '     gardenia          build gardenia code                  '
	@echo '     clean-gardenia    clean gardenia directory             '
	@echo '     narcissus         build narcissus code                 '
	@echo '     clean-narcissus   clean narcissus directory            '
	@echo '     begonia           build begonia code                   '
	@echo '     clean-begonia     clean begonia directory              '
	@echo '     lavender          build lavender code                  '
	@echo '     clean-lavender    clean lavender directory             '
	@echo '     pansy             build pansy code                     '
	@echo '     clean-pansy       clean pansy directory                '
	@echo '     manjushaka        build manjushaka code                '
	@echo '     clean-manjushaka  clean manjushaka directory           '
	@echo '     jasmine           build jasmine code                   '
	@echo '     clean-jasmine     clean jasmine directory              '
	@echo '     hibiscus          build hibiscus code                  '
	@echo '     clean-hibiscus    clean hibiscus directory             '
	@echo '     api               generate api library                 '
	@echo '     clean-api         clean api directory                  '
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
solver: build-solver
build-solver: azalea gardenia narcissus begonia lavender pansy manjushaka daisy
clean-solver: clean-azalea clean-gardenia clean-narcissus clean-begonia clean-lavender clean-pansy clean-manjushaka clean-daisy

# targets for tool
tool: build-tool
build-tool: jasmine hibiscus
clean-tool: clean-jasmine clean-hibiscus

# targets for other
other: build-other
build-other: api doc
clean-other: clean-api clean-doc

# tagets for ctqmc codes: segment version (azalea, gardenia, narcissus)
azalea: build-azalea
build-azalea:
	cd ../ctqmc/azalea/; pwd; make

clean-azalea:
	cd ../ctqmc/azalea/; pwd; make clean

gardenia: build-gardenia
build-gardenia:
	cd ../ctqmc/gardenia/; pwd; make

clean-gardenia:
	cd ../ctqmc/gardenia/; pwd; make clean

narcissus: build-narcissus
build-narcissus:
	cd ../ctqmc/narcissus/; pwd; make

clean-narcissus:
	cd ../ctqmc/narcissus/; pwd; make clean

# tagets for ctqmc codes: general version (begonia, lavender)
begonia: build-begonia
build-begonia:
	cd ../ctqmc/begonia/; pwd; make

clean-begonia:
	cd ../ctqmc/begonia/; pwd; make clean

lavender: build-lavender
build-lavender:
	cd ../ctqmc/lavender/; pwd; make

clean-lavender:
	cd ../ctqmc/lavender/; pwd; make clean

# tagets for ctqmc codes: general version (pansy, manjushaka)
pansy: build-pansy
build-pansy:
	cd ../ctqmc/pansy/; pwd; make

clean-pansy:
	cd ../ctqmc/pansy/; pwd; make clean

manjushaka: build-manjushaka
build-manjushaka:
	cd ../ctqmc/manjushaka/; pwd; make

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
	cd ../../doc/; pwd; make

clean-doc:
	cd ../../doc/; pwd; make clean
