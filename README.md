# Summary:
This software pipeline calculates TE density.
TE density is a metric that we define as the number of TE-occupied base-pairs in a given window (a base-pair range), divided by that window value, relative to each gene in the genome.
TE density values are calculated for every gene in the genome for the combination of TE *(TE superfamily || TE order) x (upstream || intragenic || downstream) wrt* a window length.
The pipeline requires two main inputs, a gene annotation and a TE annotation.
The pipeline outputs an [HDF5](https://portal.hdfgroup.org/display/HDF5/HDF5) file of TE density values for each pseudomolecule in the genome.
The pipeline can be invoked through `process_genome.py`.
We provide a class structure in a Python script (`transposon/density_data.py`) to access the sub-arrays of the results and aid the user in analyzing their data.

# Purpose:
Calculate TE density for every gene in a genome, along sliding windows, upstream, intragenically, and downstream of genes,  and correctly allocate TE density values to a gene along TE Order and Superfamily identities.
Our main goal in writing this tool was to provide a clear and concise pipeline/metric to aid in the study of TE presence relative to genes.

# Usage:
This code requires two input files, both a gene and a TE annotation.
We divide the pipeline into three stages described below, preprocessing, processing, and postprocessing.
Preprocessing requires the user to reformat their annotation files for the general pipeline, processing requires the user the user to executing the `process_genome.py` script with their chosen input arguments, and postprocessing aids the user in accessing the subarrays of the resulting TE density datafiles.

## Datasets for Examples:
The datasets can be accessed on Dryad at the following [link](https://datadryad.org/stash/share/mFjpHlP53Y-BUP4nKI0LQGVmLkXevNatnz8MLlK36zw)

## Installation Requirements:
Please install *pip* so that you may easily install Python packages.
To install all required packages for this tool run `pip install -r requirements.txt`.
We suggest that you use a virtual environment for this tool in order to better keep track of the packages being installed.
Please refer to this [guide](https://packaging.python.org/guides/installing-using-pip-and-virtual-environments/) for directions on how to use a virtual environment in conjunction with *pip*.

## System Test:
Once you have the Python packages installed and the correct version of Python
(3.10), please try the system test with `make system_test` from within the root
directory. This runs TE Density on a slightly modified *Arabidopsis thaliana*
genome. I ran the TE Density on the Arabidopsis system test files (`tests/system_test_input_data/`) in under 10 minutes on a Dell Latitude 5490 Laptop running Ubuntu 20.04.5 with an Intel i7-8650U and 32 GB of RAM.


## Preprocessing:
The user is responsible for preprocessing both of their annotation files into the format described below.
To aid in this endeavor, we provide the scripts we used to reformat our input data.
The user is encouraged to reference these scripts during their own preprocessing stage, but their explicit use is not required as long the user is able to reformat each annotation into their respective formats shown below.
Since annotation files can be slightly different, depending on the software used to create them, we use two custom scripts to preprocess/import the annotations into the format required for the pipeline.
We also use an optional helper script to reclassify TE identities.
We provide examples of preprocessing the annotation files for several different genomes in the `examples` directory; i.e: `examples/Arabidopsis/src/import_Arabidopsis_EDTA.py` and `examples/Arabidopsis/src/import_Arabidopsis_gene_anno.py`.
An example helper script may be found at `examples/Arabidopsis/src/replace_names_Arabidopsis.py`.

### Preprocessing of the Gene Annotation:
Importing the gene annotation is easy, the main difficulty arises when trying to acquire the correct gene name from the last relevant column of the annotation.
We found that the gene name frequently requires using substring or regex methods to extract it correctly.
The import gene script should output a tab-separated file with a filename prefix of "Cleaned\_" with format:

| Gene\_Name | Chromosome | Feature | Start | Stop | Strand | Length |
| :--: | :--: | :--: | :--: | :--: | :--: | :--: |
| Gene\_XYZ | Chromosome\_1 | gene | 500.0 | 600.0 | + | 101.0 |
| Gene\_ABC | Chromosome\_1 | gene | 400.0 | 700.0 | + | 301.0 |

The length is calculated by (Stop - Start + 1) because the gene is an inclusive range.

### Preprocessing of the TE Annotation:
Importing the TE annotation is another matter because the user may want to recategorize some groupings of TEs to simplify analysis later on, or because some text processing is required to acquire the TE Order and Superfamily information.
For the examples used in the publication, we chose to use [EDTA](https://github.com/oushujun/EDTA) to generate our TE annotations due to its ease-of-use, accuracy, and the fact that it outputs a TE annotation ideal for our pipeline.
We also provide our scripts that we used to generate our TE annotations using EDTA in our example directories.
The import TE script should output a tab-separated file with a filename prefix of "Cleaned\_" with format:
| Chromosome | Start | Stop | Strand | Order | SuperFamily | Length |
| :--: | :--: | :--: | :--: | :--: | :--: | :--: |
| Chromosome\_1 | 500.0 | 600.0 | + | LTR | Copia | 101.0 |

### (Optional) Redefining the TE Groupings:
Sometimes the user may want to rename groupings of TEs in their annotation. For example in `examples/Arabidosis/src/replace_names_Arabidopsis.py` we reclassified TEs to better fit the naming scheme proposed in [Wicker et al. (2007)](https://www.nature.com/articles/nrg2165).
We also renamed the superfamily values of TEs whose order value is known.
For example, if a TE's order is LTR, but its superfamily is unknown, we renamed the superfamily to "Unknown\_LTR\_Superfam".
We also did this process for TIR elements, renaming the unknown superfamily to "Unknown\_TIR\_Superfam" so that we could better distinguish between unknown LTR and TIR elements.
The replace names stage is part of the general preprocessing workflow, and is invoked during the import TEs step.
If the user does not wish to change any names or reclassify any TEs they may just leave the `te_annot_renamer` function blank, save for the return statement, which is required for the import script.

## General Pipeline (Processing):
Once the two annotation files are filtered and saved to disk, they can then be used as primary inputs to the general pipeline.
For an example on how to invoke the TE Density pipepline please refer to `examples/Arabidopsis/src/TE_Density_Arabidopsis.sb`.
It shouldn't take more than one or a few hours to run `process_genome.py` once the revised annotation is created, so if you find your command stalling out and taking a long time, please consider requesting more RAM for each processor.

### Sample Command-Line Usage:
`python process_genome.py $DATA_DIR/Cleaned_TAIR10_GFF3_genes_main_chromosomes.tsv $DATA_DIR/Cleaned_TAIR10_chr_main_chromosomes.fas.mod.EDTA.TEanno.tsv $GENOME -c $ROOT_DIR/config/production_run_config.ini -n 5 -o $OUTPUT_DIR`


### Configuration File:
A configuration file is required by the density algorithm to set the beginning window size, window delta, and max window size variables. Two configuration files are located in the `config` folder.
`test_run_config.ini` contains a setup for quickly testing the code.
 `production_run_config.ini` contains the paradigm that we use for our general implementation of the code.
 When running `process_genome.py` the code defaults to `test_run_config.ini`, however with the `-c` option the path to an optional configuration file may be supplied.
More on the configuration below in the **Optional Arguments** section.

### Command-line Arguments:
#### Required Arguments:
1. **Gene input file**: The absolute path to the cleaned gene annotation file, created during the preprocessing stage.
2. **TE input file**: The absolute path to the cleaned TE annotation file, created during the preprocessing stage
3. **Genome ID**: String representing the genome ID, used for file naming later.

#### Optional Arguments:
* `--num_threads`: The number of processors to use when running TE Density. For maximum speed, use 1:1 ratio of processors to pseudomolecules. TE density is only calculated between genes and TEs belonging to the same pseudomolecule (chromosomes for chromosome-scale assemblies).
* `--config_file`: The path to a configuration file containing the sequence of base-pair windows that will be used during the calculation of TE density. Please refer to `config/production_run_config.ini` for an example. Values must be integers. The user may provide their own configuration file as long as the format is preserved, only the integers may be changed. Defaults to `test_run_config.ini` which defines the range [400, 800] with a step value of 200.
* `--reset_h5`: Defaults to `False`. When `True`, resets the cache of wrapped GeneData and TransposonData input files. If the `--reset_h5` flag is provided as an argument to `process_genome.py`, then the boolean `True` is supplied. Useful if the user has previously initiated the general pipeline (the wrapped GeneData and TransposonData files are cached to disk at the beginning), but has since altered the preprocessed input files. This must be done otherwise the pipeline will continue to use and recognize the wrapped dataframes that were previously created, making the updated preprocessed input files irrelevant.
* `--revise_anno`: Defaults to `False`. When `True`, forces recreation of the revised TE annotation data files. If the `--revise_anno` flag is provided as an argument, then the boolean `True` is supplied. Useful if the user has previously initiated the general pipeline, and the pipeline was successful in creating a revised annotation, but the user has since altered the preprocessed input files and wants the revised annotation to reflect that. This must be done otherwise the pipeline will continue to use the cached version of the revised TE annotation, making the updated preprocessed input files irrelevant.
* `--output_dir`: The path to the directory in which all output will be placed. The output directory defaults to `../TE_Data/` relative to the source code directory.
* `--single_process`: If supplied, the flag will be `True` and the code will run without multiprocessing.

## Postprocessing:
We provide a Python script containing a data class useful in accessing the HDF5 output files.
This python script, `transposon/density_data.py`, contains the DensityData class structure which provides getter methods to access the sub-arrays of the output data (HDF5).
We also provide multiple examples of how to access the HDF5 data within the data analysis example directories.
Please refer to `examples/general_read_density_data.py` for a barebones implementation of how to initialize DensityData for a single output file.
### Example Script:
`python examples/general_read_density_data.py CLEANED_GENE_ANNOTATION.tsv DENSITY_DATA_FOLDER "Arabidopsis_(.*?).h5"`.
The `DENSITY_DATA_FOLDER` argument must **only** contain density data results (not GeneData or TEData)


## Dependencies

Tested on `Python 3.10` and `Ubuntu 22.04`.

SEE ./requirements/requirements.txt


### Troubleshooting

HDF5 and pytables requires the HDF5 runtime and development files.
e.g. on Ubuntu 22.04:
```bash
# apt install python3.10-dev
# apt install python3-distutils
# apt install libhdf5-serial-dev
```

You may require additional BLAS dependencies to support NUMPY / HDF5 etc.
```bash
# apt install libblas3 liblapack3 liblapack-dev libblas-dev
# apt install libatlas-base-dev
```

## Performance:
I ran the TE Density on the Arabidopsis system test files (`tests/system_test_input_data/`) in under 10 minutes on a Dell Latitude 5490 Laptop running Ubuntu 20.04.5 with an Intel i7-8650U and 32 GB of RAM.
