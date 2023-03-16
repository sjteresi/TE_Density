# scripts for development
# __file__ Makefile
# __author__ Michael Teresi

ROOT_DIR:=$(shell dirname $(realpath $(firstword $(MAKEFILE_LIST))))

SYS_TEST_DIR := tests/system_test_input_data
SYS_TEST_GENES := $(ROOT_DIR)/$(SYS_TEST_DIR)/Cleaned_TAIR10_GFF3_genes_main_chromosomes.tsv
SYS_TEST_TES := $(ROOT_DIR)/$(SYS_TEST_DIR)/Cleaned_TAIR10_chr_main_chromosomes.fas.mod.EDTA.TEanno.tsv


.PHONY: help
help:               ## Show this help
	@grep -E '^[a-z_A-Z0-9^.(]+:.*?## .*$$' $(MAKEFILE_LIST) | awk 'BEGIN {FS = ":.*?## "}; {printf "\033[36m%-30s\033[0m %s\n", $$1, $$2}'


.PHONY: system_test
system_test:        ## run system test on sample data
	mkdir -p ./tmp
	python $(ROOT_DIR)/process_genome.py $(SYS_TEST_GENES) $(SYS_TEST_TES) Test -o ./tmp


.PHONY: system_clean
system_clean:       ## clean the system test
	rm -rf ./tmp


.PHONY: test
test:               ## run the tests
	mkdir -p $(ROOT_DIR)/tests/test_h5_cache_loc
	mkdir -p $(ROOT_DIR)/tests/output_data
	pytest $(ROOT_DIR)


.PHONY: flake8
flake8:             ## run style guide
	flake8 $(ROOT_DIR)


.PHONY: lint
lint:               ## run linter
	pylint $(ROOT_DIR)/transposon


.PHONY: tags
tags:               ## run ctags
	ctags \
		$(ROOT_DIR)/*.py \
		$(ROOT_DIR)/transposon/*.py
