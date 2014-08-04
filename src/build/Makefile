##
##
##

help:
	echo 'hello world!'

all: common solver tool other
clean: clean-common clean-solver clean-tool clean-other

common: build-common
build-common:
	cd ../common/; pwd; make

clean-common:
	cd ../common/; pwd; make clean

solver: build-solver
build-solver: azalea gardenia narcissus begonia lavender pansy manjushaka daisy
clean-solver: clean-azalea clean-gardenia clean-narcissus clean-begonia clean-lavender clean-pansy clean-manjushaka clean-daisy

tool: build-tool
build-tool: jasmine hibiscus
clean-tool: clean-jasmine clean-hibiscus

other: build-other
build-other: api doc
clean-other: clean-api clean-doc

azalea: build-azalea
build-azalea:
	cd ../ctqmc/azalea/; pwd; make

clean-azalea:
	cd ../ctqmc/azalea/; pwd; make clean

gardenia: build-gardenia
build-gardenia:
	cd ../ctqmc/azalea/; pwd; make
clean-gardenia:

narcissus: build-narcissus
build-narcissus:
	cd ../ctqmc/azalea/; pwd; make
clean-narcissus:

begonia: build-begonia
build-begonia:
	cd ../ctqmc/azalea/; pwd; make
clean-begonia:

lavender: build-lavender
build-lavender:
	cd ../ctqmc/azalea/; pwd; make
clean-lavender:

pansy: build-pansy
build-pansy:
	cd ../ctqmc/azalea/; pwd; make
clean-pansy:

manjushaka: build-manjushaka
build-manjushaka:
	cd ../ctqmc/azalea/; pwd; make
clean-manjushaka:

daisy: build-daisy
build-daisy:
	cd ../ctqmc/azalea/; pwd; make
clean-daisy:

jasmine: build-jasmine
build-jasmine:
	cd ../ctqmc/azalea/; pwd; make
clean-jasmine:

hibiscus: build-hibiscus
build-hibiscus:
	cd ../ctqmc/azalea/; pwd; make
clean-hibiscus:

api: build-api
build-api:
	cd ../ctqmc/azalea/; pwd; make
clean-api:

doc: build-doc
build-doc:
	cd ../ctqmc/azalea/; pwd; make
clean-doc:
