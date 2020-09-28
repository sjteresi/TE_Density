# scripts for development
# __file__ Makefile
# __author__ Michael Teresi

ROOT_DIR:=$(shell dirname $(realpath $(firstword $(MAKEFILE_LIST))))
DEV_DATA := $(ROOT_DIR)/../TE_Data
DEV_CACHE := $(DEV_DATA)/cache
DEV_GENES := $(DEV_DATA)/Camarosa_Genes.gtf
DEV_TES := $(DEV_DATA)/Camarosa_EDTA_TEs.gff
DEV_GENOME := "Camarosa"

.PHONY: dev help clean test flake8 lint

dev: | tags         ## execute with default testing arguments
	mkdir -p $(DEV_CACHE)
	$(ROOT_DIR)/process_genome.py $(DEV_GENES) $(DEV_TES) $(DEV_GENOME) -vv -s $(DEV_CACHE)

help:               ## Show this help.
	fgrep -h "##" $(MAKEFILE_LIST) | fgrep -v fgrep | sed -e 's/\\$$//' | sed -e 's/##//'

clean:              ## remove temporary output files
	@echo cleaning temporary files
	rm -f $(DEV_CACHE)/*.h5

test:               ## run the tests
	pytest $(ROOT_DIR)

flake8:             ## run style guide
	flake8 $(ROOT_DIR)

lint:               ## run linter
	pylint $(ROOT_DIR)/transposon

tags:               ## run ctags
	ctags \
		$(ROOT_DIR)/*.py \
		$(ROOT_DIR)/transposon/*.py

