# scripts for running rice synteny TE Density examples
# __file__ Makefile
# __author__ Scott Teresi

ROOT_DIR:=$(shell dirname $(realpath $(firstword $(MAKEFILE_LIST))))
DEV_GFF_READ_EXECUTABLE := /home/scott/Documents/Uni/Research/gffread
#DEV_GFF_READ_EXECUTABLE := /mnt/research/edgerpat_lab/Scotty/gffread
DEV_DATA := $(realpath $(ROOT_DIR)/../../../TE_Density_Example_Data/Arabidopsis)
DEV_ARAB_FASTA_DATA := $(DEV_DATA)/Sequences/TAIR10_chr_main_chromosomes.fas
DEV_ARAB_GENES := $(DEV_DATA)/Genes/TAIR10_GFF3_genes_main_chromosomes.gff
DEV_ARAB_EXP := $(DEV_DATA)/Genes/At-Expression_AtBrCv_4Wisecaver.csv
DEV_ARAB_TEs := $(DEV_DATA)/TEs/TAIR10_chr_main_chromosomes.fas.mod.EDTA.TEanno.gff3
DEV_FILTERED := $(realpath $(ROOT_DIR)/../../../TE_Data/filtered_input_data)
DEV_HDF5 := $(realpath $(ROOT_DIR)/../../../TE_Data/finalized_data/10KB/Arabidopsis)
DEV_RESULTS:= $(realpath $(ROOT_DIR)/results)

.PHONY: create_CDS fix_fasta_names fix_CDS_names

create_CDS:
	@echo
	@echo Creating CDS from GFF and fasta file for Arabidopsis
	$(DEV_GFF_READ_EXECUTABLE)/gffread -x $(DEV_DATA)/Sequences/Arabidopsis_CDS.fasta -g $(DEV_ARAB_FASTA_DATA) $(DEV_ARAB_GENES)
	@echo

run_EDTA_HPCC:
	@echo Running EDTA for Arabidopsis
	sbatch $(ROOT_DIR)/src/Annotate_EDTA_Arabidopsis_thaliana.sb

filter_genes:
	@echo Filtering Arabidopsis genes into appropriate format for TE Density
	python $(ROOT_DIR)/src/import_Arabidopsis_gene_anno.py $(DEV_ARAB_GENES)

filter_TEs:
	@echo Filtering Arabidopsis TEs into appropriate format for TE Density
	python $(ROOT_DIR)/src/import_Arabidopsis_EDTA.py $(DEV_ARAB_TEs)

run_TE_Density_HPCC:
	@echo Running TE Density for Glaberrima
	@echo sbatch file contains paths inside
	sbatch $(ROOT_DIR)/src/TE_Density_Arabidopsis.sb

generate_dotplots:
	@echo Generating TE density dotplot for Arabidopsis
	mkdir -p $(DEV_RESULTS)/graphs
	python $(ROOT_DIR)/src/generate_dotplots.py $(DEV_HDF5)/Arabidopsis_Chr1.h5 $(DEV_FILTERED)/Cleaned_TAIR10_GFF3_genes_main_chromosomes.tsv 1 -o $(DEV_RESULTS)/graphs

compare_density_upstream_downstream:
	@echo Comparing upstream and downstream density arrays with chi-squared
	mkdir -p $(DEV_RESULTS)/chi_squared
	python $(ROOT_DIR)/src/chi_squared_density_comparison.py $(DEV_FILTERED)/Cleaned_TAIR10_GFF3_genes_main_chromosomes.tsv $(DEV_HDF5) -o $(DEV_RESULTS)/chi_squared

compare_centromeric:
	@echo Comparing centromeric/pericentromeric and regular density arrays
	mkdir -p $(DEV_RESULTS)/centromeric
	python $(ROOT_DIR)/src/compare_centromeric_densities.py $(DEV_FILTERED)/Cleaned_TAIR10_GFF3_genes_main_chromosomes.tsv $(DEV_HDF5) $(DEV_ARAB_EXP) -o $(DEV_RESULTS)/centromeric
