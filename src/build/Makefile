##
##
##

help:
	@echo 'hello world!'

all: common solver tool other
clean: clean-common clean-solver clean-tool clean-other

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
