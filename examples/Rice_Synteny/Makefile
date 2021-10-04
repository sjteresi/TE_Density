# scripts for running rice synteny TE Density examples
# __file__ Makefile
# __author__ Scott Teresi

ROOT_DIR:=$(shell dirname $(realpath $(firstword $(MAKEFILE_LIST))))
DEV_GFF_READ_EXECUTABLE := /home/scott/Documents/Uni/Research/gffread
#DEV_GFF_READ_EXECUTABLE := /mnt/research/edgerpat_lab/Scotty/gffread  # use if on HPCC
DEV_DATA := $(realpath $(ROOT_DIR)/../../../TE_Density_Example_Data/Rice)
DEV_GLABERRIMA_DATA := $(DEV_DATA)/Oryza_Glaberrima
DEV_GLABERRIMA_FASTA_DATA := $(DEV_GLABERRIMA_DATA)/Sequences
DEV_GLABERRIMA_GENES := $(DEV_GLABERRIMA_DATA)/Genes/Oryza_glaberrima.Oryza_glaberrima_V1.50.gff3
DEV_GLABERRIMA_TEs := $(DEV_GLABERRIMA_DATA)/TEs/Oryza_Glaberrima_NewNames.fasta.mod.EDTA.TEanno.gff3
DEV_SATIVA_DATA := $(DEV_DATA)/Oryza_Sativa
DEV_SATIVA_FASTA_DATA := $(DEV_SATIVA_DATA)/Sequences
DEV_SATIVA_GENES := $(DEV_SATIVA_DATA)/Genes/Oryza_sativa.IRGSP-1.0.50.gff3
DEV_SATIVA_TEs := $(DEV_SATIVA_DATA)/TEs/Oryza_Sativa_NewNames.fasta.mod.EDTA.TEanno.gff3
DEV_FILTERED := $(realpath $(ROOT_DIR)/../../../TE_Data/filtered_input_data)
DEV_HDF5 := $(realpath $(ROOT_DIR)/../../../TE_Data/finalized_data/10KB)
DEV_RESULTS:= $(realpath $(ROOT_DIR)/results)

.PHONY: create_CDS fix_fasta_names fix_CDS_names

create_CDS:
	@echo
	@echo Creating CDS from GFF and fasta file for glaberrima
	$(DEV_GFF_READ_EXECUTABLE)/gffread -x $(DEV_GLABERRIMA_FASTA_DATA)/Oryza_glaberrima_CDS.fasta -g $(DEV_GLABERRIMA_FASTA_DATA)/Oryza_glaberrima.Oryza_glaberrima_V1.dna.toplevel.fa $(DEV_GLABERRIMA_GENES)
	@echo Creating CDS from GFF and fasta file for sativa
	$(DEV_GFF_READ_EXECUTABLE)/gffread -x $(DEV_SATIVA_FASTA_DATA)/Oryza_sativa_CDS.fasta -g $(DEV_SATIVA_FASTA_DATA)/Oryza_sativa.IRGSP-1.0.dna.toplevel.fa $(DEV_SATIVA_GENES)
	@echo


fix_fasta_names:
	@echo
	@echo Fixing the fasta names for glaberrima so that they are not too long for EDTA
	python $(ROOT_DIR)/src/fix_fasta_names.py $(DEV_GLABERRIMA_FASTA_DATA)/Oryza_glaberrima.Oryza_glaberrima_V1.dna.toplevel.fa Oryza_Glaberrima
	@echo Fixing the fasta names for sativa so that they are not too long for EDTA
	python $(ROOT_DIR)/src/fix_fasta_names.py $(DEV_SATIVA_FASTA_DATA)/Oryza_sativa.IRGSP-1.0.dna.toplevel.fa Oryza_Sativa
	@echo



fix_CDS_names:
	@echo
	@echo Fixing the CDS fasta names for glaberrima so that they are not too long for EDTA
	python $(ROOT_DIR)/src/fix_cds_names.py $(DEV_GLABERRIMA_FASTA_DATA)/Oryza_glaberrima_CDS.fasta Oryza_Glaberrima
	@echo Fixing the CDS fasta names for sativa so that they are not too long for EDTA
	python $(ROOT_DIR)/src/fix_cds_names.py $(DEV_SATIVA_FASTA_DATA)/Oryza_sativa_CDS.fasta Oryza_Sativa
	@echo

run_EDTA_HPCC:
	@echo Running EDTA for glaberrima
	sbatch $(ROOT_DIR)/src/Annotate_EDTA_Rice_Glaberrima.sb
	@echo Running EDTA for sativa
	sbatch $(ROOT_DIR)/src/Annotate_EDTA_Rice_Sativa.sb


filter_genes:
	@echo Filtering glaberrima genes into appropriate format for TE Density
	python $(ROOT_DIR)/src/import_rice_gene_anno.py $(DEV_GLABERRIMA_GENES)
	@echo Filtering sativa genes into appropriate format for TE Density
	python $(ROOT_DIR)/src/import_rice_gene_anno.py $(DEV_SATIVA_GENES)

filter_TEs:
	@echo Filtering blueberry TEs into appropriate format for TE Density
	python $(ROOT_DIR)/src/import_rice_EDTA.py $(DEV_GLABERRIMA_TEs)
	python $(ROOT_DIR)/src/import_rice_EDTA.py $(DEV_SATIVA_TEs)

run_TE_Density_HPCC:
	@echo Running TE Density for Glaberrima
	sbatch $(ROOT_DIR)/src/TE_Density_Glaberrima.sb
	sbatch $(ROOT_DIR)/src/TE_Density_Sativa.sb

filter_syntelogs:
	@echo Filtering syntelog file from SynMap to a cleaner version for downstream analysis
	mkdir -p $(ROOT_DIR)/results
	python $(ROOT_DIR)/src/generate_pairs.py $(DEV_DATA)/SynMap_Results/Glaberrima_VS_Sativa_SynMap.tsv -o $(ROOT_DIR)/results


compare_syntelog_TE_differences:
	@echo  Generate graphs of syntelog TE differences
	@echo This is for chromosome 1
	mkdir -p $(ROOT_DIR)/results/graphs
	python $(ROOT_DIR)/src/compare_density.py $(DEV_RESULTS)/set_syntelogs.tsv $(DEV_HDF5)/Sativa/Sativa_1.h5 $(DEV_HDF5)/Glaberrima/Glaberrima_1.h5 $(DEV_FILTERED)/Cleaned_Oryza_sativa.IRGSP-1.0.50.tsv $(DEV_FILTERED)/Cleaned_Oryza_glaberrima.Oryza_glaberrima_V1.50.tsv -o $(ROOT_DIR)/results

# TODO Rename inputs for Rice
identify_interesing_genes:
	@echo Reporting genes with interesting levels of TE density
	mkdir -p $(DEV_RESULTS)/tables
	python $(ROOT_DIR)/src/find_abnormal_genes.py $(DEV_FILTERED)/Cleaned_Oryza_sativa.IRGSP-1.0.50.tsv $(DEV_HDF5)/Sativa -o $(DEV_RESULTS)/tables
