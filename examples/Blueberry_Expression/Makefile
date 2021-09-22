# scripts for running blueberry TE Density examples
# __file__ Makefile
# __author__ Scott Teresi

ROOT_DIR:=$(shell dirname $(realpath $(firstword $(MAKEFILE_LIST))))
DEV_DATA := $(realpath $(ROOT_DIR)/../../../TE_Density_Example_Data/Blueberry)
DEV_GENES := $(DEV_DATA)/Genes/Blueberry_Genes.gff
DEV_TES := $(DEV_DATA)/TEs/Blueberry_EDTA_TEs.gff
DEV_GENE_EXPRESSION := $(DEV_DATA)/Genes/Blueberry_TPM_All.tsv
DEV_DENSITY_FILES := $(DEV_DATA)/../../TE_Data/finalized_data/10KB/Draper
DEV_GENOME := "BlueberryExample"

filter_genes:
	@echo Filtering blueberry genes into appropriate format for TE Density
	@echo 
	python $(ROOT_DIR)/src/import_blueberry_gene_anno.py $(DEV_GENES)

filter_TEs:
	@echo Filtering blueberry TEs into appropriate format for TE Density
	@echo 
	python $(ROOT_DIR)/src/import_blueberry_EDTA.py $(DEV_TES)

calculate_TE_Density:
	@echo Running TE Density for blueberry
	sbatch $(ROOT_DIR)/src/TE_Density_Blueberry.sb

generate_expression_graphs:
	@echo Generating TE density vs. gene expression graphs
	@echo 
	mkdir -p $(ROOT_DIR)/results/graphs
	python $(ROOT_DIR)/src/compare_expression.py $(DEV_GENE_EXPRESSION) $(DEV_GENES) $(DEV_DENSITY_FILES) -o $(ROOT_DIR)/results/graphs
