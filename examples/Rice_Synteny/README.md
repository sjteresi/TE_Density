# TE Density Comparisons of Syntelogs in Rice:
This README contains information pertaining to how we compared TE Density levels of syntelogs using two closely related rice genomes. Please refer to the `Makefile` for our explicit commands we used from start to finish in creating and using the data relevant to this example.

# Downloading of Rice FASTA and Gene Annotation Files:

# Generating a TE Annotation for Rice using EDTA:
First, we need to generate a TE annotation for the two rice genomes. We do this by creating a CDS FASTA for each genome and then fixing (shortening) the names of the sequences in both the regular FASTA and CDS FASTA to run EDTA without any warnings (it does not like long sequence ID names).
We use three main code blocks in the `Makefile` to get the necessary files prior to running EDTA:
1. `create_CDS` utilizes `gffread` to create a CDS FASTA file using the regular FASTA file and a gene annotation.
2. `fix_fasta_names` alters the sequence ID names in the FASTA file so that they aren't too long for EDTA.
3. `fix_CDS_names` alters the sequence ID names in the FASTA file so that they aren't too long for EDTA.
4. `run_EDTA_HPCC` executes the commands required to run EDTA on the HPCC (the computing cluster of MSU). Please refer to the `Annotate_EDTA_Rice_Glaberrima.sb` and `Annotate_EDTA_Rice_Sativa.sb` files to see the commands and resources we used in generating the EDTA annotations.

# Perform the Pre-processing Steps of TE Density, Filter Both Annotation Files:
Here please refer to `filter_genes` and `filter_TEs` to view the commands we invoked to pre-process our data. The `import_rice_gene_anno.py` and `import_rice_EDTA.py` are the primary scripts the user will need to edit for their own purposes. 
1. Rice gene annotation...



# Running SynMap:
This section describes the methods to run [SynMap](https://genomevolution.org/CoGe/SynMap.pl) on CoGe. I ran SynMap with mostly [default options](https://genomevolution.org/wiki/index.php/SynMap), I did change one option: under *Merge Syntenic Blocks* I set it to `Quota Align Merge`. Here is the [link](https://genomevolution.org/r/1how2) for *Glaberrima vs Sativa*.
