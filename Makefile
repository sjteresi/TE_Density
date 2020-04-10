# scripts for development
# __file__ Makefile
# __author__ Michael Teresi

ROOT_DIR:=$(shell dirname $(realpath $(firstword $(MAKEFILE_LIST))))
DEV_DATA := $(ROOT_DIR)/../TE_Data
DEV_CACHE := $(DEV_DATA)/cache
DEV_GENES := $(DEV_DATA)/camarosa_gtf_data.gtf
DEV_TES := $(DEV_DATA)/camarosa_gff_data.gff
DEV_GENOME := "Camarosa"

.PHONY: dev help clean test flake8 lint

dev: | clean        ## execute with default testing arguments
	@echo executing dev main
	mkdir -p $(DEV_CACHE)
	-python3 $(ROOT_DIR)/transposon/density.py $(DEV_GENES) $(DEV_TES) $(DEV_GENOME) -vv -s $(DEV_CACHE)

help:               ## Show this help.
	fgrep -h "##" $(MAKEFILE_LIST) | fgrep -v fgrep | sed -e 's/\\$$//' | sed -e 's/##//'

clean:              ## remove temporary output files
	@echo cleaning temporary files
	rm -f $(DEV_CACHE)/*.h5

test:               ## run the tests
	@echo running all tests
	pytest $(ROOT_DIR)

flake8:             ## run style guide
	flake8 $(ROOT_DIR)

lint:               ## run linter
	pylint $(ROOT_DIR)/transposon
