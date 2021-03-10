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

create_CDS:
	@echo
	@echo Creating CDS from GFF and fasta file for Human 
	$(DEV_GFF_READ_EXECUTABLE)/gffread -x $(DEV_DATA)/Human_CDS.fasta -g $(DEV_DATA)/Fasta/hg38_chr1_chr7.fasta $(DEV_DATA)/Genes/Human_Gene_Annotation_Chr1_Chr7.gff3
	@echo

fix_fasta_names:
	@echo
	@echo Fixing the fasta names for Human so that they are not too long for EDTA
	python $(ROOT_DIR)/src/fix_fasta_names.py $(DEV_DATA)/Fasta/hg38_chr1_chr7.fasta Human
	@echo

fix_CDS_names:
	@echo
	@echo Fixing the CDS fasta names for Human so that they are not too long for EDTA
	python $(ROOT_DIR)/src/fix_cds_names.py $(DEV_DATA)/Human_CDS.fasta Human
	@echo

run_EDTA_HPCC:
	@echo Running EDTA for Human
	sbatch $(ROOT_DIR)/src/Annotate_EDTA_Human.sb


################################
# NOTE all of the above is used for generating an EDTA annotation for Human
# I have tried to not use that for my analyses, and instead rely on the community
# annotations, however if I continue to have trouble with formatting the community
# annotations I will utilize the EDTA version
################################


filter_genes:
	@echo Filtering human genes into appropriate format for TE Density
	python $(ROOT_DIR)/src/import_human_gene_anno.py $(DEV_GENES)

filter_TEs:
	python $(ROOT_DIR)/src/import_human_te_anno.py $(DEV_TES)


help:               ## Show this help.
	fgrep -h "##" $(MAKEFILE_LIST) | fgrep -v fgrep | sed -e 's/\\$$//' | sed -e 's/##//'
