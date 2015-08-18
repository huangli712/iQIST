##
##
## Introduction
## ============
##
## This is the top Makefile for the iqist software package. You can use
## it to build all of the executable programs and libraries. Please do
## not modify it by yourself even you are very familiar with iqist.
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
	@echo '     daisy-lib         build daisy library                  '
	@echo '     clean-daisy       clean daisy directory                '
	@echo
	@echo '     jasmine           build jasmine code                   '
	@echo '     jasmine-lib       build jasmine library                '
	@echo '     clean-jasmine     clean jasmine directory              '
	@echo
	@echo '     hibiscus          build hibiscus code                  '
	@echo '     clean-hibiscus    clean hibiscus directory             '
	@echo
	@echo '     common            build common modules and util library'
	@echo '     clean-common      clean common directory               '
	@echo
	@echo '     api               generate api library                 '
	@echo '     clean-api         clean api directory                  '
	@echo
	@echo '     doc               generate document                    '
	@echo '     clean-doc         clean doc directory                  '

# all-in-one building target
all: common api solver tool other

# all-in-one building target (lib mode)
lib: common api solver-lib tool-lib other

# all-in-one cleaning target
clean: clean-common clean-api clean-solver clean-tool clean-other

# targets for solver
solver: common api \
        azalea gardenia narcissus \
        begonia lavender \
        pansy manjushaka \
        daisy
solver-lib: common api \
	    azalea-lib gardenia-lib narcissus-lib \
	    begonia-lib lavender-lib \
	    pansy-lib manjushaka-lib \
	    daisy-lib
clean-solver: clean-azalea clean-gardenia clean-narcissus \
              clean-begonia clean-lavender \
              clean-pansy clean-manjushaka \
              clean-daisy

# targets for tool
tool: common api jasmine hibiscus
tool-lib: common api jasmine-lib hibiscus
clean-tool: clean-jasmine clean-hibiscus

# targets for other
other: doc
clean-other: clean-doc

# targets for ctqmc codes: segment version (azalea, gardenia, narcissus)
azalea:
	cd ../ctqmc/azalea/; pwd; make

azalea-lib:
	cd ../ctqmc/azalea/; pwd; make lib

clean-azalea:
	cd ../ctqmc/azalea/; pwd; make clean

gardenia:
	cd ../ctqmc/gardenia/; pwd; make

gardenia-lib:
	cd ../ctqmc/gardenia/; pwd; make lib

clean-gardenia:
	cd ../ctqmc/gardenia/; pwd; make clean

narcissus:
	cd ../ctqmc/narcissus/; pwd; make

narcissus-lib:
	cd ../ctqmc/narcissus/; pwd; make lib

clean-narcissus:
	cd ../ctqmc/narcissus/; pwd; make clean

# targets for ctqmc codes: general version (begonia, lavender)
begonia:
	cd ../ctqmc/begonia/; pwd; make

begonia-lib:
	cd ../ctqmc/begonia/; pwd; make lib

clean-begonia:
	cd ../ctqmc/begonia/; pwd; make clean

lavender:
	cd ../ctqmc/lavender/; pwd; make

lavender-lib:
	cd ../ctqmc/lavender/; pwd; make lib

clean-lavender:
	cd ../ctqmc/lavender/; pwd; make clean

# targets for ctqmc codes: general version (pansy, manjushaka)
pansy:
	cd ../ctqmc/pansy/; pwd; make

pansy-lib:
	cd ../ctqmc/pansy/; pwd; make lib

clean-pansy:
	cd ../ctqmc/pansy/; pwd; make clean

manjushaka:
	cd ../ctqmc/manjushaka/; pwd; make

manjushaka-lib:
	cd ../ctqmc/manjushaka/; pwd; make lib

clean-manjushaka:
	cd ../ctqmc/manjushaka/; pwd; make clean

# targets for hfqmc codes: daisy
daisy:
	cd ../hfqmc/daisy/; pwd; make

daisy-lib:
	cd ../hfqmc/daisy/; pwd; make lib

clean-daisy:
	cd ../hfqmc/daisy/; pwd; make clean

# targets for atomic codes: jasmine
jasmine:
	cd ../tools/jasmine/; pwd; make

jasmine-lib:
	cd ../tools/jasmine/; pwd; make lib

clean-jasmine:
	cd ../tools/jasmine/; pwd; make clean

# targets for toolkit codes: hibiscus
hibiscus:
	cd ../tools/hibiscus/; pwd; make

clean-hibiscus:
	cd ../tools/hibiscus/; pwd; make clean

# targets for common: common
common:
	cd ../common/; pwd; make

clean-common:
	cd ../common/; pwd; make clean

# targets for api: api
api:
	cd ../api/; pwd; make

clean-api:
	cd ../api/; pwd; make clean

# targets for doc: doc
doc:
	cd ../../doc/manual/; pwd; make full

clean-doc:
	cd ../../doc/manual/; pwd; make clean
