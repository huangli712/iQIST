##
##
## Introduction
## ============
##
## Makefile for iqist/src/tools/hibiscus directory.
##
## Usage
## =====
##
## 'make' or 'make clean'
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
## 11/10/2014 by li huang (created)
## 08/17/2015 by li huang (last modified)
##
##

default: all

all: build-entropy build-stoch build-script build-swing build-toolbox
clean: clean-entropy clean-stoch clean-script clean-swing clean-toolbox

# targets for entropy
build-entropy:
	cd ./entropy; pwd; make

clean-entropy:
	cd ./entropy; pwd; make clean

# targets for stoch
build-stoch:
	cd ./stoch; pwd; make

clean-stoch:
	cd ./stoch; pwd; make clean

# targets for script
build-script:
	cd ./script; pwd; make

clean-script:
	cd ./script; pwd; make clean

# targets for swing
build-swing:
	cd ./swing; pwd; make

clean-swing:
	cd ./swing; pwd; make clean

# targets for toolbox
build-toolbox:
	cd ./toolbox; pwd; make

clean-toolbox:
	cd ./toolbox; pwd; make clean
