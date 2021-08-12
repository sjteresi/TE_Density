# Examples:
The purpose of this repository is to demonstrate the usage of TE Density on multiple different genomes as well as show interesting application and data visualization opportunities using the TE Density's output data sets.
The makefile in each genome-specific directory contains the sequence of commands we used to perform each analysis.


## Example 1 - Genome-Wide Trends of TE Presence Relative to Genes (Arabidopsis):
The TE Density toolkit may be used to investigate the relationship of genes and TE presence genome-wide.
Here, we create a dotplot of the average TE density values of all genes by TE type, upstream, intragenically, and downstream for chromosome 1 of the Arabidopsis genome.

### Step 1 - Create a CDS FASTA file for use in creating a TE annotation:
Please reference the makefile command `create_CDS`, here we use the GFFRead tool to generate a CDS FASTA file which is useful in creating a TE annotation with EDTA.

### Step 2 - Generate a TE Annotation with EDTA:
Please reference the makefile command `run_EDTA_HPCC`, here we run EDTA on MSU's high-performance computing cluster (HPCC).

### Step 3 - Preprocess each annotation file prior to running TE Density:
Please reference the makefile commands `filter_genes` and `filter_TEs`, here we use the Python package Pandas to perform reformat the annotation files. This part will likely need to be custom tailored to the user's own annotation files.

### Step 4 - Run TE Density:
Please reference the makefile command `run_TE_Density_HPCC`, here we call `process_genome.py` for the Arabidopsis genome but do so in an SBATCH script because we are using the resources of the HPCC. It is better to err on the side of more RAM if the genome is taking too long to compute, it may have stalled out due to insufficient RAM.

### Step 5 - Begin analysis of TE Density data:
#### Generate dotplots of average TE Density values for all genes as the window changes:
Please reference the makefile command `generate_dotplots` and its python file `generate_dotplots.py`. Here we initialize the DensityData class to aid in accessing the data in the output HDF5 files.

