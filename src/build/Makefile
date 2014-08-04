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
	cd ../ctqmc/gardenia/; pwd; make

clean-gardenia:
	cd ../ctqmc/gardenia/; pwd; make clean

narcissus: build-narcissus
build-narcissus:
	cd ../ctqmc/narcissus/; pwd; make

clean-narcissus:
	cd ../ctqmc/narcissus/; pwd; make clean

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

daisy: build-daisy
build-daisy:
	cd ../ctqmc/daisy/; pwd; make

clean-daisy:
	cd ../ctqmc/daisy/; pwd; make clean

jasmine: build-jasmine
build-jasmine:
	cd ../ctqmc/jasmine/; pwd; make

clean-jasmine:
	cd ../ctqmc/jasmine/; pwd; make clean

hibiscus: build-hibiscus
build-hibiscus:
	cd ../ctqmc/hibiscus/; pwd; make

clean-hibiscus:
	cd ../ctqmc/hibiscus/; pwd; make clean

api: build-api
build-api:
	cd ../ctqmc/api/; pwd; make

clean-api:
	cd ../ctqmc/api/; pwd; make clean

doc: build-doc
build-doc:
	cd ../ctqmc/doc/; pwd; make

clean-doc:
	cd ../ctqmc/doc/; pwd; make clean
