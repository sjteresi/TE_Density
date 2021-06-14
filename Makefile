# scripts for development
# __file__ Makefile
# __author__ Michael Teresi

ROOT_DIR:=$(shell dirname $(realpath $(firstword $(MAKEFILE_LIST))))
DEV_DATA := $(realpath $(ROOT_DIR)/../TE_Data)
DEV_CACHE := $(DEV_DATA)/tmp                                # MAGIC default
DEV_CACHE_OVERLAP := $(addsuffix /overlap,$(DEV_CACHE))     # MAGIC default
DEV_GENES := $(DEV_DATA)/Camarosa_Genes.gtf
DEV_TES := $(DEV_DATA)/Camarosa_EDTA_TEs.gff
DEV_PROCESSED_GENES := $(DEV_DATA)/filtered_input_data/Cleaned_Camarosa_Genes.tsv
DEV_PROCESSED_TES := $(DEV_DATA)/filtered_input_data/Cleaned_Camarosa_EDTA_TEs.tsv
DEV_PROD_CONF := $(ROOT_DIR)/config/production_run_config.ini
DEV_GENOME := "Camarosa"

.PHONY: dev help clean test flake8 lint

dev: | tags         ## execute with default testing arguments
	mkdir -p $(DEV_CACHE)
	$(ROOT_DIR)/process_genome.py $(DEV_GENES) $(DEV_TES) $(DEV_GENOME) -vv


# Use this command to make the processed/cleaned input files
filtered_data:
	python $(ROOT_DIR)/examples/Strawberry/src/import_strawberry_gene_anno.py $(DEV_GENES)
	python $(ROOT_DIR)/examples/Strawberry/src/import_strawberry_EDTA.py $(DEV_TES)

production: | tags         ## execute with default production arguments
	mkdir -p $(DEV_CACHE)
	$(ROOT_DIR)/process_genome.py $(DEV_PROCESSED_GENES) $(DEV_PROCESSED_TES) $(DEV_GENOME) -c $(DEV_PROD_CONF) -n 16 -vv

help:               ## Show this help.
	fgrep -h "##" $(MAKEFILE_LIST) | fgrep -v fgrep | sed -e 's/\\$$//' | sed -e 's/##//'

clean_overlap:      ## remove overlap temp files
	@echo cleaning temporary files
	cd $(DEV_CACHE_OVERLAP) && rm *.h5

test:               ## run the tests
	mkdir -p $(ROOT_DIR)/tests/test_h5_cache_loc
	mkdir -p $(ROOT_DIR)/tests/output_data
	pytest $(ROOT_DIR)

flake8:             ## run style guide
	flake8 $(ROOT_DIR)

lint:               ## run linter
	pylint $(ROOT_DIR)/transposon

tags:               ## run ctags
	ctags \
		$(ROOT_DIR)/*.py \
		$(ROOT_DIR)/transposon/*.py

