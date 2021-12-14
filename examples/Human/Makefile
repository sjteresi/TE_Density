# scripts for running human TE Density examples
# __file__ Makefile
# __author__ Scott Teresi

ROOT_DIR:=$(shell dirname $(realpath $(firstword $(MAKEFILE_LIST))))
DEV_GFF_READ_EXECUTABLE := /home/scott/Documents/Uni/Research/gffread
#DEV_GFF_READ_EXECUTABLE := /mnt/research/edgerpat_lab/Scotty/gffread
DEV_DATA := $(realpath $(ROOT_DIR)/../../../TE_Density_Example_Data/Human)
#DEV_GENES := $(DEV_DATA)/Genes/Chr7_13_Human_Genes.tsv
DEV_GENES := $(DEV_DATA)/Genes/gencode.v37.annotation.gff3
#DEV_TES := $(DEV_DATA)/TEs/Chr7_13_Human_TEs.tsv
DEV_TES := $(DEV_DATA)/TEs/Human_TEs.tsv
DEV_PROD_CONF := $(ROOT_DIR)/../../config/production_run_config.ini

DEV_FILTERED := $(realpath $(ROOT_DIR)/../../../TE_Data/filtered_input_data)
DEV_HDF5 := $(realpath $(ROOT_DIR)/../../../TE_Data/finalized_data/10KB)
DEV_RESULTS:= $(realpath $(ROOT_DIR)/results)

filter_genes:
	@echo Filtering human genes into appropriate format for TE Density
	python $(ROOT_DIR)/src/import_human_gene_anno.py $(DEV_GENES)

filter_TEs:
	python $(ROOT_DIR)/src/import_human_te_anno.py $(DEV_TES)

examine_TE_levels:
	@echo
	python $(ROOT_DIR)/src/retrieve_info_of_genes.py $(DEV_HDF5)/Human/Human_Chrom_Subset_chr7.h5 $(DEV_FILTERED)/Cleaned_Chr7_13_Human_Genes.tsv 7 Human_Chr7 $(ROOT_DIR)/list_of_chr7_genes.txt
	@echo
	python $(ROOT_DIR)/src/retrieve_info_of_genes.py $(DEV_HDF5)/Human/Human_Chrom_Subset_chr13.h5 $(DEV_FILTERED)/Cleaned_Chr7_13_Human_Genes.tsv 13 Human_Chr13 $(ROOT_DIR)/list_of_chr13_genes.txt
help:               ## Show this help.
	fgrep -h "##" $(MAKEFILE_LIST) | fgrep -v fgrep | sed -e 's/\\$$//' | sed -e 's/##//'
